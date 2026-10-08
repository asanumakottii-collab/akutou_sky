"""Idealized aerosol water uptake; separate from condensed fog droplets.

Petters & Kreidenweis (2007), doi:10.5194/acp-7-1961-2007.
Kelvin curvature and deliquescence hysteresis are neglected. A single uniform
growth factor shifts the lognormal distribution without changing particle count.
"""
from .config import ModelConfig

WATER_INDEX = complex(1.333, -1e-9)  # visible constant-index approximation

def aerosol_water_uptake(config: ModelConfig):
    """Return growth factor, wet median radius (nm), and volume-mixed index.

    aod550 remains the DRY column extinction when this effect is active. The
    transport module converts it using wet/dry Mie cross sections at 550 nm.
    """
    config.validate()
    growth_cubed = 1 + config.hygroscopicity_kappa * config.relative_humidity / (1 - config.relative_humidity)
    factor = growth_cubed ** (1 / 3)
    dry_index = complex(config.refractive_index_real, -config.refractive_index_imag)
    wet_index = (dry_index + WATER_INDEX * (growth_cubed - 1)) / growth_cubed
    return factor, config.median_radius_nm * factor, wet_index
