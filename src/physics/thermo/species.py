"""Minimal species bookkeeping for future reacting-flow support."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class SpeciesSet:
    """Collection of species metadata and a default inert composition."""

    names: tuple[str, ...]
    molecular_weights: np.ndarray
    default_mass_fractions: np.ndarray

    @property
    def nspecies(self) -> int:
        return len(self.names)


def create_species_set(names, molecular_weights, default_mass_fractions):
    """Construct and validate species metadata."""
    names_tuple = tuple(names)
    molecular_weights_array = np.asarray(molecular_weights, dtype=np.float64)
    default_mass_fractions_array = np.asarray(default_mass_fractions, dtype=np.float64)

    if molecular_weights_array.ndim != 1:
        raise ValueError("Species molecular weights must be a 1D array.")
    if default_mass_fractions_array.ndim != 1:
        raise ValueError("Default mass fractions must be a 1D array.")
    if len(names_tuple) == 0:
        raise ValueError("At least one species is required.")
    if len(names_tuple) != molecular_weights_array.size:
        raise ValueError("Species names and molecular weights must have matching length.")
    if len(names_tuple) != default_mass_fractions_array.size:
        raise ValueError("Species names and default mass fractions must have matching length.")
    if np.any(molecular_weights_array <= 0.0):
        raise ValueError("Species molecular weights must be positive.")

    mass_fraction_sum = np.sum(default_mass_fractions_array)
    if not np.isclose(mass_fraction_sum, 1.0):
        raise ValueError("Default mass fractions must sum to 1.0.")
    if np.any(default_mass_fractions_array < 0.0):
        raise ValueError("Default mass fractions must be non-negative.")

    return SpeciesSet(
        names=names_tuple,
        molecular_weights=molecular_weights_array,
        default_mass_fractions=default_mass_fractions_array,
    )


def make_inert_air_species():
    """Single-species placeholder for today's non-reacting solver."""
    return create_species_set(
        names=("air",),
        molecular_weights=(28.97e-3,),
        default_mass_fractions=(1.0,),
    )


def initialize_mass_fractions(ncells: int, species_set: SpeciesSet):
    """Allocate per-cell species mass fractions using the default composition."""
    mass_fractions = np.zeros((ncells, species_set.nspecies), dtype=np.float64)
    mass_fractions[:] = species_set.default_mass_fractions
    return mass_fractions
