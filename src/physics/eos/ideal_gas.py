"""Calorically perfect ideal-gas thermodynamics."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class IdealGasModel:
    """Constant-gamma ideal-gas model used by the current solver."""

    gamma: float
    gas_constant: float
    name: str = "ideal_gas"


def pressure_from_internal_energy(rho, e, gas_model: IdealGasModel):
    """p = (gamma - 1) * rho * e"""
    return (gas_model.gamma - 1.0) * rho * e


def sound_speed(rho, p, gas_model: IdealGasModel):
    """a = sqrt(gamma * p / rho)"""
    return np.sqrt(gas_model.gamma * p / rho)


def temperature_from_rho_p(rho, p, gas_model: IdealGasModel):
    """T = p / (rho * R)"""
    return p / (rho * gas_model.gas_constant)


def total_energy_density(rho, u, p, gas_model: IdealGasModel, v=0.0):
    """rhoE = p / (gamma - 1) + 0.5 * rho * (u^2 + v^2)"""
    return p / (gas_model.gamma - 1.0) + 0.5 * rho * (u**2 + v**2)


def pressure_from_conserved(rho, rhou, rhov, rhoE, gas_model: IdealGasModel):
    """Recover pressure from the conservative state."""
    kinetic = 0.5 * (rhou**2 + rhov**2) / rho
    return (gas_model.gamma - 1.0) * (rhoE - kinetic)
