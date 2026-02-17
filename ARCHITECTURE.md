# ARCHITECTURE.md — System Design & Invariants

## Core Philosophy
This project follows a **Layered Domain Architecture**. To maintain legibility and prevent "spaghetti physics," dependencies must only flow **downward**. A lower layer must never import or depend on a layer above it.

## The Layers
1. **App/CLI (`app/`)**: Thin wrappers that handle user input via `argparse`. This layer translates human intent into simulation parameters.
2. **Solver (`src/solver/`)**: The execution engine. It manages the time-stepping loop, checks for convergence, and coordinates between the Physics and Mesh layers.
3. **Physics (`src/physics/`)**: The numerical kernels. These are "pure functions" that implement the Navier-Stokes equations or other fluid laws using NumPy.
4. **Mesh (`src/core/`)**: The foundation. This layer defines the grid, cell connectivity, and geometric properties like areas and volumes.

## Taste Invariants (Rules for the Agent)
These rules are enforced to ensure the codebase remains "agent-legible":

* **Array Orientation**: All spatial arrays representing the physical domain must follow the `(Nx, Ny)` shape convention to match Cartesian mapping.
* **Pure Physics**: Physics kernels must be side-effect free. They take data, compute fluxes/sources, and return data. They do not perform I/O or printing.
* **No Global State**: All simulation data (velocity, pressure, time) must be encapsulated in a `State` object. We do not use global variables.
* **Boundary Handling**: Boundary conditions must be handled by a dedicated provider or wrapper. We avoid "if/else" logic for boundaries inside the main numerical loops to keep the math clean.
* **Mechanical Enforcement**: Any architectural violation (e.g., a circular import or a Physics kernel trying to call the Solver) should be caught by a linter or a structural test.
