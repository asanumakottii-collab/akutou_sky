"""Scalar delta-M transport with an explicit single-scattering correction.

The backend's native SO forcing uses the untruncated single-scatter phase.
For narrow fog peaks this aliases into a coarse incoming quadrature. Instead
solve the consistently truncated problem, then replace its observer single
scatter with the backend's full-phase, scaled-optical-depth single scatter
(TMS, Nakajima & Tanaka 1988). No colour smoothing or intensity clipping.
"""
import numpy as np
import sasktran2 as sk


def truncated_coefficients(scaled_coefficients, order):
    """Remove the forward delta consistently from normalized scalar moments.

    Input is the backend's TMS-scaled beta_l/(1-f), including the (2l+1)
    convention. The resulting zeroth moment must be exactly one.
    """
    coefficients = np.asarray(scaled_coefficients)
    if order < 1 or coefficients.shape[0] < order:
        raise ValueError('Insufficient phase coefficients for delta-M order')
    degrees = (2*np.arange(order)+1).reshape((order,) + (1,)*(coefficients.ndim-1))
    result = coefficients[:order] - degrees*(coefficients[0]-1.)
    result[0] = 1.
    return result


class DeltaMTransport:
    def __init__(self, full_config, geometry, viewing, order):
        self.order = order
        self.geometry = geometry
        self.config = sk.Config()
        for name in ('num_threads', 'num_stokes', 'num_sza',
                     'num_successive_orders_incoming', 'num_successive_orders_outgoing',
                     'successive_orders_reduced_horizon_quadrature',
                     'successive_orders_altitude_grid_m',
                     'successive_orders_horizontal_angle_grid_radians',
                     'num_successive_orders_iterations', 'successive_orders_relative_tolerance',
                     'successive_orders_absolute_tolerance', 'solar_refraction',
                     'los_refraction', 'multiple_scatter_refraction', 'wavelength_batch_size'):
            setattr(self.config, name, getattr(full_config, name))
        self.config.num_streams = order
        self.config.num_singlescatter_moments = order
        self.config.delta_m_scaling = False  # full atmosphere is already scaled
        self.config.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrders
        self.multiple = sk.Engine(self.config, geometry, viewing)
        single_config = sk.Config()
        for name in ('num_threads', 'num_stokes', 'solar_refraction',
                     'los_refraction', 'multiple_scatter_refraction', 'wavelength_batch_size'):
            setattr(single_config, name, getattr(self.config, name))
        single_config.num_streams = order
        single_config.num_singlescatter_moments = order
        single_config.delta_m_scaling = False
        single_config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
        self.single = sk.Engine(single_config, geometry, viewing)

    def correction(self, full_atmosphere):
        """Return truncated total minus truncated observer single scatter.

        Call after the full-phase single-scatter engine has populated/scaled
        full_atmosphere. Both raw solvers use identical scaled media.
        """
        raw = sk.Atmosphere(self.geometry, self.config,
                            wavelengths_nm=full_atmosphere.wavelengths_nm,
                            calculate_derivatives=False)
        source = full_atmosphere.storage
        raw.storage.total_extinction[:] = source.total_extinction
        raw.storage.ssa[:] = source.ssa
        raw.storage.solar_irradiance[:] = source.solar_irradiance
        raw.surface.albedo[:] = full_atmosphere.surface.albedo
        # Native scaling stores beta_l/(1-f); beta_0=1/(1-f).
        # Subtract (2l+1) f/(1-f), leaving beta'_0=1 exactly.
        raw.storage.leg_coeff[:] = truncated_coefficients(source.leg_coeff, self.order)
        multiple = self.multiple.calculate_radiance(raw).radiance
        single = self.single.calculate_radiance(raw).radiance
        return multiple - single
