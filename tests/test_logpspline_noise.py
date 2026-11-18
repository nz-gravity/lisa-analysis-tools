"""Full complex XYZ spectra must reach the existing likelihood unchanged."""

import copy
import pickle
import unittest
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import numpy as np

from lisatools.analysiscontainer import AnalysisContainer
from lisatools.domains import FDSettings, FDSignal, WDMSettings
from lisatools.logpsplines import LogPSplineNoise
from lisatools.sensitivity import (
    CompositeSensitivityBackend,
    CompositeSensitivityMatrix,
)


def spectrum(frequency, scale=1e-40):
    """Unequal complex correlations with distinguishable frequency bins."""
    root = np.broadcast_to(np.eye(3, dtype=complex), (len(frequency), 3, 3)).copy()
    root[:, 1, 0] = 0.25 + 0.2j
    root[:, 2, 0] = -0.15 + 0.1j
    root[:, 2, 1] = 0.3 - 0.12j
    root[:, 2, 2] = 1.5
    covariance = root @ root.conj().swapaxes(-1, -2)
    return scale * (1 + 20 * frequency[:, None, None]) * covariance


class LogPSplineNoiseTest(unittest.TestCase):
    def setUp(self):
        self.settings = FDSettings(33, 1e-4, min_freq=1e-4, force_backend="cpu")
        self.frequency = np.asarray(self.settings.f_arr)
        self.matrix = spectrum(self.frequency)
        self.noise = LogPSplineNoise(self.frequency, self.matrix)

    def test_full_complex_matrix_inverse_and_determinant(self):
        backend = CompositeSensitivityBackend(
            self.settings, extra_components=[self.noise], force_backend="cpu"
        )
        matrix = backend("logpspline", psd_params=None)
        np.testing.assert_allclose(matrix.sens_mat, self.matrix.transpose(1, 2, 0), atol=0)
        np.testing.assert_allclose(
            matrix.invC.transpose(2, 0, 1), np.linalg.inv(self.matrix), rtol=2e-14, atol=0
        )
        np.testing.assert_allclose(matrix.detC.real, np.linalg.det(self.matrix), rtol=2e-13, atol=0)
        self.assertEqual(len(matrix.components), 1)  # No second instrument model.
        self.assertGreater(np.max(np.abs(matrix.sens_mat.imag)), 0)

    def test_analysis_container_matches_independent_complex_gaussian(self):
        rng = np.random.default_rng(29)
        # LDC continuous-transform coefficients, not unit-Whittle coefficients.
        residual = 1e-19 * (
            rng.normal(size=(len(self.frequency), 3))
            + 1j * rng.normal(size=(len(self.frequency), 3))
        )
        matrix = CompositeSensitivityMatrix(self.settings, [self.noise])
        container = AnalysisContainer(FDSignal(residual.T, self.settings), matrix)
        quadratic = np.einsum(
            "fi,fi->", residual.conj(), np.linalg.solve(self.matrix, residual[..., None])[..., 0]
        ).real
        expected = -2 * self.settings.df * quadratic - np.linalg.slogdet(self.matrix)[1].sum()
        self.assertAlmostEqual(float(container.likelihood()), expected, places=9)
        self.assertAlmostEqual(
            float(container.likelihood(noise_only=True)),
            -np.linalg.slogdet(self.matrix)[1].sum(),
            places=9,
        )
        diagonal = np.zeros_like(self.matrix)
        diagonal[:, np.arange(3), np.arange(3)] = self.matrix[:, np.arange(3), np.arange(3)]
        q_diagonal = np.einsum(
            "fi,fi->", residual.conj(), np.linalg.solve(diagonal, residual[..., None])[..., 0]
        ).real
        diagonal_loglike = -2 * self.settings.df * q_diagonal - np.linalg.slogdet(diagonal)[1].sum()
        self.assertGreater(abs(expected - diagonal_loglike), 1.0)

    def test_active_subband_and_roundoff_match_without_interpolation(self):
        subset = FDSettings(33, 1e-4, min_freq=8e-4, max_freq=2e-3, force_backend="cpu")
        indices = np.rint(np.asarray(subset.f_arr) / subset.df).astype(int) - 1
        np.testing.assert_array_equal(
            self.noise.covariance(subset), self.matrix[indices].transpose(1, 2, 0)
        )
        rounded = LogPSplineNoise(self.frequency * (1 + 2e-16), self.matrix)
        np.testing.assert_array_equal(
            rounded.covariance(self.settings), self.noise.covariance(self.settings)
        )

    def test_rejects_missing_bins_dc_and_other_domains(self):
        missing = LogPSplineNoise(self.frequency[::2], self.matrix[::2])
        with self.assertRaisesRegex(ValueError, "every active FD frequency"):
            missing.covariance(self.settings)
        outside = FDSettings(34, 1e-4, min_freq=1e-4, force_backend="cpu")
        with self.assertRaisesRegex(ValueError, "extrapolation"):
            self.noise.covariance(outside)
        with self.assertRaisesRegex(ValueError, "exclude DC"):
            self.noise.covariance(FDSettings(33, 1e-4, force_backend="cpu"))
        with self.assertRaisesRegex(TypeError, "FDSettings"):
            self.noise.covariance(WDMSettings(8, 8, 5.0, force_backend="cpu"))

    def test_rejects_invalid_matrices_at_physical_psd_scales(self):
        invalid = {}
        invalid["shape"] = self.matrix[:, :2, :2]
        invalid["finite"] = self.matrix.copy()
        invalid["finite"][0, 0, 0] = np.nan
        invalid["Hermitian"] = self.matrix.copy()
        invalid["Hermitian"][0, 0, 1] += 1e-41j
        invalid["positive definite"] = self.matrix.copy()
        invalid["positive definite"][0, 0, 0] = -1e-40
        for message, value in invalid.items():
            with self.subTest(message=message), self.assertRaisesRegex(ValueError, message):
                LogPSplineNoise(self.frequency, value)
        for frequency in (self.frequency[::-1], self.frequency * np.nan, self.frequency * 0):
            with self.subTest(frequency=frequency), self.assertRaises(ValueError):
                LogPSplineNoise(frequency, self.matrix)
        with self.assertRaisesRegex(ValueError, "channel order"):
            LogPSplineNoise(self.frequency, self.matrix, channels=("A", "E", "T"))

    def test_copies_inputs_and_survives_settings_serialization(self):
        self.matrix[:] = 0
        self.assertGreater(np.max(np.abs(self.noise.covariance(self.settings))), 0)
        for restored in (copy.deepcopy(self.noise), pickle.loads(pickle.dumps(self.noise))):
            np.testing.assert_array_equal(
                restored.covariance(self.settings), self.noise.covariance(self.settings)
            )

    def test_result_mean_draw_and_preview(self):
        import xarray as xr

        values = np.stack((self.matrix, 2 * self.matrix))[None]
        data = xr.DataArray(
            values,
            dims=("chain", "draw", "frequency", "channel", "channel_aux"),
            coords={
                "chain": [0],
                "draw": [0, 1],
                "frequency": self.frequency,
                "channel": [0, 1, 2],
                "channel_aux": [0, 1, 2],
            },
        )
        result = SimpleNamespace(
            spectrum=data, posterior=xr.Dataset(coords={"chain": [0], "draw": [0, 1]})
        )
        mean = LogPSplineNoise.from_result(result)
        np.testing.assert_allclose(mean.spectral_density, 1.5 * self.matrix, rtol=1e-14, atol=0)
        draw = LogPSplineNoise.from_result(result, chain=0, draw=1)
        np.testing.assert_allclose(draw.spectral_density, 2 * self.matrix, rtol=1e-14, atol=0)
        with self.assertRaisesRegex(ValueError, "both chain and draw"):
            LogPSplineNoise.from_result(result, draw=0)
        result.spectrum = data.isel(draw=[1])
        with self.assertRaisesRegex(ValueError, "preview"):
            LogPSplineNoise.from_result(result)
        result.spectrum = data.assign_coords(channel=["A", "E", "T"])
        with self.assertRaisesRegex(ValueError, "XYZ"):
            LogPSplineNoise.from_result(result)


class EreborLogPSplineConfigurationTest(unittest.TestCase):
    def test_real_erebor_walkers_match_unwindowed_fd_likelihood(self):
        from eryn.state import BranchSupplemental
        from mpi4py import MPI

        from lisatools.globalfit.run import GlobalFit
        from lisatools.globalfit.stock import erebor
        from lisatools.globalfit.stock.erebor.injections import UserDataProcessor

        n, dt = 64, 5.0
        settings = FDSettings(
            n // 2 + 1, 1 / (n * dt), min_freq=0.025, max_freq=0.08, force_backend="cpu"
        )
        frequency = np.asarray(settings.f_arr)
        covariance = spectrum(frequency)
        noise = LogPSplineNoise(frequency, covariance)
        data = 1e-21 * np.random.default_rng(7).normal(size=(3, n))
        with TemporaryDirectory(prefix="erebor-logpspline-") as output:
            fit = erebor.blank(
                nwalkers=2,
                ntemps=2,
                use_gpu=False,
                gpus=None,
                dt=dt,
                nf=8,
                nt=8,
                Tobs=n * dt,
                min_freq=0.025,
                max_freq=0.08,
                # The reference here describes rectangular-window data. The
                # inherited wavelet taper must also be disabled explicitly.
                window_taper_wavelets=0,
                window_taper_duration=0.0,
                file_store_dir=output + "/",
                data_processor_class=UserDataProcessor,
                processor_init_kwargs={"data": data, "dt": dt},
            )
            erebor.configure_logpspline_noise(fit, noise, settings=settings)
            fit.build()
            runtime = GlobalFit(fit, MPI.COMM_SELF)
            priors = {}
            for name in fit.branch_names:
                priors.update(fit.source_info[name].priors)
            state = runtime.load_info(priors)
            nt, nw = runtime.ntemps, runtime.nwalkers
            state.supplemental = BranchSupplemental(
                {"walker_inds": np.tile(np.arange(nw), (nt, 1))},
                base_shape=(nt, nw),
                copy=True,
            )
            containers = runtime.setup_acs(state)
            residual = (dt * np.fft.rfft(data, axis=-1))[:, settings.active_slice].T
            quadratic = np.einsum(
                "fi,fi->",
                residual.conj(),
                np.linalg.solve(covariance, residual[..., None])[..., 0],
            ).real
            expected = -2 * settings.df * quadratic - np.linalg.slogdet(covariance)[1].sum()
            self.assertFalse(containers[0].likelihood_source_only)
            np.testing.assert_allclose(
                containers[0].sens_mat.sens_mat, covariance.transpose(1, 2, 0), atol=0
            )
            np.testing.assert_allclose(
                np.asarray(containers.likelihood()), expected, rtol=1e-12, atol=1e-9
            )

    def test_fixed_noise_configuration_and_guards(self):
        from lisatools.globalfit.stock.erebor.noise import configure_logpspline_noise

        settings = FDSettings(33, 1e-4, min_freq=1e-4, force_backend="cpu")
        frequency = np.asarray(settings.f_arr)
        noise = LogPSplineNoise(frequency, spectrum(frequency))
        general = SimpleNamespace(sensitivity_init_kwargs=None)
        fit = SimpleNamespace(general=general, built=False, branches={})
        self.assertIs(configure_logpspline_noise(fit, noise, settings=settings), fit)
        self.assertIs(general.domain_settings, settings)
        self.assertIs(general.sensitivity_backend_class, CompositeSensitivityBackend)
        self.assertEqual(general.fixed_psd_kwargs, {"psd_params": None, "galfor_params": None})
        self.assertIsNone(general.fixed_psd_params)
        self.assertFalse(general.likelihood_source_only)
        backend = general.sensitivity_backend_class(
            settings, force_backend="cpu", **general.sensitivity_init_kwargs
        )
        np.testing.assert_allclose(
            backend("erebor", **general.fixed_psd_kwargs).sens_mat,
            spectrum(frequency).transpose(1, 2, 0),
            atol=0,
        )
        fit.branches = {"psd": object()}
        with self.assertRaisesRegex(ValueError, "sampled noise"):
            configure_logpspline_noise(fit, noise, settings=settings)
        fit.branches = {}
        fit.built = True
        with self.assertRaisesRegex(ValueError, "before building"):
            configure_logpspline_noise(fit, noise, settings=settings)


if __name__ == "__main__":
    unittest.main()
