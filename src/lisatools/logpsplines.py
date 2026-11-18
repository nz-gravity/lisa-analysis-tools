"""Import reconstructed LogPSplinePSD XYZ spectra into the noise likelihood.

The fitting library is optional: this boundary consumes physical, one-sided
spectral matrices, not spline coefficients, precision matrices, or whitened
corrections. The first integration supports stationary FD covariances only.
"""

from __future__ import annotations

import numpy as np

from .domains import FDSettings
from .sensitivity import NoiseComponent
from .utils.utility import asnumpy


class LogPSplineNoise(NoiseComponent):
    """A fixed full XYZ spectral matrix reconstructed by LogPSplinePSD.

    Args:
        frequency: Strictly increasing positive physical frequencies in Hz.
        spectral_density: Complex ``(Nf, 3, 3)`` one-sided spectral density,
            with ``S[i,j] = E[z_i conj(z_j)]``. Use the same TDI generation,
            observables, units, and window convention as the residual data.
        channels: Input channel order; must be ``("X", "Y", "Z")``.

    Small Hermitian roundoff (relative tolerance ``1e-10`` per frequency)
    is symmetrized. Non-finite, non-Hermitian, or non-positive-definite
    matrices are rejected, including at physical PSD scales. No diagonal
    loading, frequency interpolation, extrapolation, or FFT rescaling is
    applied. Every requested FD frequency must already occur in the table.

    This is a conditional, fixed-noise likelihood. A posterior mean or a
    chosen draw does not marginalize noise uncertainty in Erebor.
    """

    name = "logpspline_xyz"

    def __init__(self, frequency, spectral_density, *, channels=("X", "Y", "Z")):
        if tuple(channels) != ("X", "Y", "Z"):
            raise ValueError("LogPSplineNoise requires channel order (X, Y, Z).")
        frequency = np.array(frequency, dtype=float, copy=True)
        covariance = np.array(spectral_density, dtype=np.complex128, copy=True)
        if (
            frequency.ndim != 1
            or frequency.size == 0
            or not np.isfinite(frequency).all()
            or np.any(frequency <= 0)
            or np.any(np.diff(frequency) <= 0)
        ):
            raise ValueError("frequency must be increasing, finite, positive Hz bins.")
        if covariance.shape != (frequency.size, 3, 3):
            raise ValueError("spectral_density must have shape (Nf, 3, 3).")
        if not np.isfinite(covariance).all():
            raise ValueError("spectral_density must contain only finite values.")
        adjoint = covariance.conj().swapaxes(-1, -2)
        scale = np.max(np.abs(covariance), axis=(-1, -2))
        error = np.max(np.abs(covariance - adjoint), axis=(-1, -2))
        if np.any(error > 1e-10 * scale):
            raise ValueError("spectral_density must be Hermitian at every frequency.")
        covariance = covariance * 0.5 + adjoint * 0.5
        try:
            np.linalg.cholesky(covariance)
        except np.linalg.LinAlgError as exc:
            raise ValueError(
                "spectral_density must be positive definite at every frequency."
            ) from exc
        frequency.flags.writeable = False
        covariance.flags.writeable = False
        self.frequency = frequency
        self.spectral_density = covariance

    @classmethod
    def from_result(cls, result, *, chain=None, draw=None, channels=("X", "Y", "Z")):
        """Use a stationary ``PSDResult`` mean or one reconstructed draw.

        With no indices, average complete reconstructed matrices over every
        chain/draw; incomplete spectrum previews are refused. With both
        ``chain`` and ``draw``, select zero-based positional indices instead.
        The result must use physical Hz and original data units. For a fit
        to reference-whitened coefficients or log-frequency coordinates,
        reconstruct ``L0 R L0^H`` on physical bins and use the constructor.
        Time-dependent results require an explicit physical reconstruction.

        No import of LogPSplinePSD is needed to consume its result object.
        """
        spectrum = result.spectrum
        dims = ("chain", "draw", "frequency", "channel", "channel_aux")
        if set(spectrum.dims) != set(dims):
            raise ValueError("Expected a stationary matrix PSDResult spectrum.")
        for axis in ("channel", "channel_aux"):
            labels = tuple(spectrum.coords[axis].values)
            if labels not in ((0, 1, 2), ("X", "Y", "Z")):
                raise ValueError("PSDResult channels must identify XYZ in that order.")
        if (chain is None) != (draw is None):
            raise ValueError("Specify both chain and draw, or neither for the mean.")
        if chain is None:
            for axis in ("chain", "draw"):
                posterior = result.posterior
                if spectrum.sizes[axis] != posterior.sizes[axis] or not np.array_equal(
                    spectrum.coords[axis], posterior.coords[axis]
                ):
                    raise ValueError(
                        "Posterior mean requires all reconstructed draws, not a preview."
                    )
            spectrum = spectrum.mean(dim=("chain", "draw"), skipna=False)
        else:
            spectrum = spectrum.isel(chain=chain, draw=draw)
        values = spectrum.transpose("frequency", "channel", "channel_aux").values
        return cls(spectrum.coords["frequency"].values, values, channels=channels)

    def covariance(self, settings):
        """Return ``(3, 3, N_active)`` on a matching positive-frequency FD grid.

        WDM/STFT folding and cross-frequency gap covariance need separate
        validation and are deliberately unsupported by this initial adapter.
        """
        if not isinstance(settings, FDSettings):
            raise TypeError("LogPSplineNoise currently requires FDSettings.")
        target = asnumpy(settings.f_arr)
        if np.any(target <= 0):
            raise ValueError("Select a positive-frequency active FD band (exclude DC).")
        # Match roundoff from equivalent Hz-grid calculations without ever
        # substituting a neighbouring frequency or evaluating between bins.
        indices = np.searchsorted(self.frequency, target)
        right = np.minimum(indices, len(self.frequency) - 1)
        left = np.maximum(indices - 1, 0)
        use_left = np.abs(self.frequency[left] - target) < np.abs(self.frequency[right] - target)
        indices = np.where(use_left, left, right)
        if not np.all(np.isclose(self.frequency[indices], target, rtol=1e-12, atol=0)):
            raise ValueError(
                "Reconstruct LogPSplinePSD on every active FD frequency; "
                "interpolation and extrapolation are not supported."
            )
        return settings.xp.asarray(
            np.ascontiguousarray(np.moveaxis(self.spectral_density[indices], 0, -1))
        )
