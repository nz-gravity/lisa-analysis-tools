Sensitivity Information
~~~~~~~~~~~~~~~~~~~~~~~~~~

The analytic sensitivity code is heavily based on an original code by Stas Babak, Antoine Petiteau for the LDC team.

The module has four layers: single-channel PSD classes (:class:`~lisatools.sensitivity.Sensitivity`
subclasses), sensitivity (covariance) matrices, the composite noise model of the global fit, and
the native XYZ backend. The 9-month and 12-month global-fit runs use the composite model
(:class:`~lisatools.sensitivity.CompositeSensitivityBackend`) on XYZ TDI 2.0 channels with the
unequal-arm instrument noise (:class:`~lisatools.sensitivity.UnequalArmInstrumentNoise`,
``UNEQUAL_ARM=1``) on the WDM domain (``WDM_PSD_METHOD=layer_calibrated``). The analytic
equal-arm curves and matrices below are the quick-study building blocks; the global-fit noise
model is the unequal-arm composite one.

References for noise models:
  * 'Proposal': LISA Consortium Proposal for L3 mission: LISA_L3_20170120 (https://atrium.in2p3.fr/13414ec1-c9ac-44b4-bace-7004468f684c)
  * 'SciRDv1': Science Requirement Document: ESA-L3-EST-SCI-RS-001 14/05/2018 (https://atrium.in2p3.fr/f5a78d3e-9e19-47a5-aa11-51c81d370f5f)
  * 'MRDv1': Mission Requirement Document: ESA-L3-EST-MIS-RS-001 08/12/2017
  * 'sangria': Model used for LDC2A data generation.

.. autofunction:: lisatools.sensitivity.get_sensitivity

.. autofunction:: lisatools.sensitivity.get_stock_sensitivity_options

.. autofunction:: lisatools.sensitivity.get_stock_sensitivity_matrix_options

.. autofunction:: lisatools.sensitivity.get_stock_sensitivity_from_str

.. autofunction:: lisatools.sensitivity.check_sensitivity

Sensitivity Matrix
---------------------
The sensitivity matrix is designed to house sensitivity information that will enter into the Likelihood function. 
It is a matrix because the 3 channels in LISA will be correlated. Therefore, in XYZ channels, this is a 3x3 matrix that includes correlated cross-terms between the channels. 
When working with AET channels (uncorrelated in idealized noise situations), then the sensitivity matrix will be an array of length 2 for AE or 3 for AET. 
The :func:`lisatools.diagnostic.inner_product` will adjust its computation based on the shape of the sensitivity matrix input. 
The user does not have to do anything special for this change to work. It happens under the hood.

.. autoclass:: lisatools.sensitivity.SensitivityMatrixBase
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.SensitivityMatrix
    :members:
    :show-inheritance:

Stock Sensitivity Matrices
******************************

.. autoclass:: lisatools.sensitivity.XYZ2SensitivityMatrix
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.XYZ1SensitivityMatrix
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.AET2SensitivityMatrix
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.AET1SensitivityMatrix
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.AE2SensitivityMatrix
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.AE1SensitivityMatrix
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.LISASensSensitivityMatrix
    :members:
    :show-inheritance:


Sensitivity Base Class
--------------------------

.. autoclass:: lisatools.sensitivity.Sensitivity
    :members:
    :show-inheritance:

Stock Sensitivity Models
*************************

.. autoclass:: lisatools.sensitivity.LISASens
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.CornishLISASens
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.FlatPSDFunction
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.X1TDISens
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.Y1TDISens
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.Z1TDISens
    :members:
    :show-inheritance:  

.. autoclass:: lisatools.sensitivity.XY1TDISens
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.YZ1TDISens
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.ZX1TDISens
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.A1TDISens
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.E1TDISens
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.T1TDISens
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.X2TDISens
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.Y2TDISens
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.Z2TDISens
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.XY2TDISens
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.YZ2TDISens
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.ZX2TDISens
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.A2TDISens
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.E2TDISens
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.T2TDISens
    :members:
    :show-inheritance:

Unequal-arm TDI 2.0 elements
*****************************

The XYZ TDI 2.0 covariance elements for six independent link light travel times (the
production instrument model, through :class:`~lisatools.sensitivity.UnequalArmInstrumentNoise`).

.. autoclass:: lisatools.sensitivity.UnequalArmXX2TDISens
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.UnequalArmYY2TDISens
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.UnequalArmZZ2TDISens
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.UnequalArmXY2TDISens
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.UnequalArmXZ2TDISens
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.UnequalArmYZ2TDISens
    :members:
    :show-inheritance:

Composite Noise Model (global fit)
------------------------------------

A :class:`~lisatools.sensitivity.CompositeSensitivityMatrix` sums noise components: the
instrument (:class:`~lisatools.sensitivity.UnequalArmInstrumentNoise` with its
:class:`~lisatools.sensitivity.LinkDelayTable` in production), the galactic foreground
(:class:`~lisatools.sensitivity.GalacticForeground`, optionally time-modulated by
:class:`~lisatools.sensitivity.GalForTimeModulation`), an optional
:class:`~lisatools.sensitivity.SGWB` and the tabulated
:class:`~lisatools.sensitivity.MojitoNoiseEstimates`.
:class:`~lisatools.sensitivity.CompositeSensitivityBackend` builds one such matrix per walker
from the sampled ``psd`` and ``galfor`` parameters.

.. autoclass:: lisatools.sensitivity.SensitivityBackendBase
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.CompositeSensitivityBackend
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.CompositeSensitivityMatrix
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.MojitoNoiseSensitivityMatrix
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.NoiseComponent
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.SeparableComponent
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.InstrumentNoise
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.UnequalArmInstrumentNoise
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.LinkDelayTable
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.GalacticForeground
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.GalForTimeModulation
    :members:
    :show-inheritance:

.. autodata:: lisatools.sensitivity.GALFOR_MODULATION_TABLE

.. autofunction:: lisatools.sensitivity.packaged_galfor_modulation_path

.. autofunction:: lisatools.sensitivity.modulation_from_elements

.. autoclass:: lisatools.sensitivity.SGWB
    :members:
    :show-inheritance:

.. autoclass:: lisatools.sensitivity.MojitoNoiseEstimates
    :members:
    :show-inheritance:

.. autofunction:: lisatools.sensitivity.estimate_noise_params_from_file

.. autofunction:: lisatools.sensitivity.tdi_generation_from_channel

LogPSplinePSD XYZ covariance (experimental)
******************************************

:class:`~lisatools.logpsplines.LogPSplineNoise` imports a reconstructed,
one-sided complex XYZ spectral density through the existing
:class:`~lisatools.sensitivity.CompositeSensitivityBackend`. This initial
integration supports a **fixed stationary FD covariance**. It preserves
the complex off-diagonal elements and uses the existing matrix inverse,
quadratic form and log-determinant. LogPSplinePSD/JAX is not required to
consume an already reconstructed matrix.

The input is ``S[f, channel, channel_aux]`` in XYZ order, with
``S_ij = E[z_i conj(z_j)]``, original data units and physical Hz bins.
The adapter changes the layout to ``(3, 3, N_active)``. It checks
Hermiticity and positive definiteness and refuses missing frequencies,
DC, interpolation and extrapolation. Reconstruct on every active native
FD bin; a set of plotting frequencies or adaptive bin centres is insufficient.
Match the residual's TDI generation, channel definition, units, preprocessing
and observation window before importing a covariance.

For example, given physical ``frequency`` and ``S_xyz`` arrays:

.. code-block:: python

    from lisatools.domains import FDSettings
    from lisatools.logpsplines import LogPSplineNoise
    from lisatools.sensitivity import CompositeSensitivityBackend

    # N is the total rFFT bin count; the active band excludes DC/Nyquist.
    fd = FDSettings(N, df, min_freq=fmin, max_freq=fmax, force_backend="cpu")
    noise = LogPSplineNoise(frequency, S_xyz)
    backend = CompositeSensitivityBackend(
        fd, extra_components=[noise], force_backend="cpu"
    )
    sensitivity = backend("logpspline", psd_params=None)

``psd_params=None`` is essential: the imported covariance supplies the noise;
passing the stock instrument levels would add another instrument component.
An ordinary stationary LogPSplinePSD ``PSDResult`` in physical units and Hz
can instead be imported with ``LogPSplineNoise.from_result(result)``. This
averages complete reconstructed matrices, preserving positive definiteness;
it rejects incomplete spectrum previews. Use ``chain=0, draw=0`` for an
explicit reconstructed draw. An elementwise complex median is not a
positive-definite matrix summary.

Reference-whitened case-study fits
++++++++++++++++++++++++++++++++

The `short XYZ gaps notebook
<https://nz-gravity.github.io/logpsplines_lisa_case_studies/notebooks/04_short_xyz_gaps.html>`_
fits a dimensionless matrix correction ``R`` after whitening native XYZ
coefficients. Its log-frequency spline coordinates are not physical Hz.
Use its ``matrix_mean`` reconstruction, including restored conditional theta
fields, to evaluate ``R`` on the native grid. Then undo the reference
whitening before constructing the adapter:

.. code-block:: python

    # Existing notebook objects: selected, components, frequency, reference, bins.
    # R is evaluated at native bins, not at the adaptive aggregate centres.
    coordinate = np.log(frequency / frequency[0]) / np.log(frequency[-1] / frequency[0])
    R = matrix_mean(selected["gapped"].posterior, components, coordinate)
    L0 = np.linalg.cholesky(reference["gapped"][bins])
    S_window = L0 @ R @ L0.conj().swapaxes(-1, -2)
    noise = LogPSplineNoise(frequency, S_window)

Here ``L0`` must describe the same observed window as the residual. For the
notebook's ``z = rfft(window * xyz) * sqrt(2*dt/n)``, the corresponding LAT
coefficients are ``x = dt * rfft(window * xyz)`` and ``z = sqrt(2*df) * x``.
The existing FD log likelihood is therefore
``-2*df*sum(x^H S_window^-1 x) - sum(log(det(S_window)))``, up to common
constants. No additional factor of two or frequency-bin width belongs in
``S_window``.

This retains a Whittle approximation: a same-frequency 3x3 matrix does not
represent gap-induced correlations between frequencies, pseudo-covariance,
or the notebook's exact selected-feature Gaussian benchmark. An underlying
time-dependent physical spectrum also needs an explicit observation-window
model; silently averaging epochs is not supported.

Using a fixed covariance in Erebor
+++++++++++++++++++++++++++++++++

For an unbuilt fixed-noise fit with user-supplied XYZ time data:

.. code-block:: python

    from lisatools.globalfit.stock import erebor
    from lisatools.globalfit.stock.erebor.injections import UserDataProcessor

    # xyz_td has shape (3, n); fd describes its natural rFFT grid.
    fit = erebor.blank(
        dt=dt,
        Tobs=1 / df,
        window_taper_wavelets=0,
        window_taper_duration=0.0,
        data_processor_class=UserDataProcessor,
        processor_init_kwargs={"data": xyz_td, "dt": dt},
    )
    erebor.configure_logpspline_noise(fit, noise, settings=fd)
    # Add source branches/moves, then build/run using the usual Erebor API.

The helper retains both likelihood terms and disables the stock instrument
component. It rejects already-built fits, existing extra noise components,
and sampled ``psd``, ``galfor`` or ``sgwb`` branches. Supply the complete
desired fixed noise covariance. It leaves data conditioning to the caller;
use the same window for the data and imported covariance. The example
disables Erebor's additional taper: for a window-aware covariance, supply
the already windowed/gapped time stream. Otherwise, reconstruct the covariance
for the actual Erebor analysis window, including its configured taper duration.

This first draft conditions on a fitted matrix or one selected draw. It does
not propagate the LogPSplinePSD posterior or implement joint spline updates
in ``PSDMove``. WDM/STFT projection and time-dependent covariances need
separate normalization and likelihood tests before production integration.

.. autoclass:: lisatools.logpsplines.LogPSplineNoise
    :members:
    :show-inheritance:

.. autofunction:: lisatools.globalfit.stock.erebor.noise.configure_logpspline_noise

Native XYZ Backend
----------------------

An alternative ``sensitivity_backend`` computed by the C++/CUDA ``SensitivityMatrixWrap``
kernel; the production runs use :class:`~lisatools.sensitivity.CompositeSensitivityBackend`.

.. autoclass:: lisatools.sensitivity.XYZSensitivityBackend
    :members:
    :show-inheritance:

Synthetic Noise Generation
-----------------------------

Time-domain noise and foreground realizations, used by the stock fits' synthetic data streams
(:mod:`lisatools.globalfit.preprocessing`, :mod:`lisatools.globalfit.stock.erebor.injections`).
The instrument-noise draw uses the equal-arm :class:`~lisatools.sensitivity.InstrumentNoise`
covariance.

.. autofunction:: lisatools.sensitivity.generate_correlated_instrument_noise_td

.. autofunction:: lisatools.sensitivity.generate_foreground_td

.. autofunction:: lisatools.sensitivity.annual_amplitude_envelope

.. autofunction:: lisatools.sensitivity.annual_modulation_matrix
