# GIFPY — Liquid Simulation Formula Reference

This document is the mathematical source of truth for the initial **3D weakly-compressible Smoothed Particle Hydrodynamics (WCSPH)** solver. It is a small reference implementation, not yet a production-quality free-surface renderer.

## 1. State and units

For particle i: position x_i in metres, velocity v_i in m/s, mass m_i in kg, density rho_i in kg/m^3, pressure p_i in Pa, acceleration a_i in m/s^2, smoothing length h in metres, and timestep dt in seconds.

Let r_ij = x_i - x_j and r = norm(r_ij). Kernels are zero for r >= h. All kernels below use 3D normalization.

## 2. 3D kernels

### 2.1 Poly6 density kernel

W_poly6(r,h) = 315/(64 pi h^9) * (h^2-r^2)^3, for 0 <= r < h; otherwise 0.

Its gradient is grad_i W_poly6 = -945/(32 pi h^9) * (h^2-r^2)^2 * r_ij. The implementation uses the kernel for density; this gradient is documented for completeness.

### 2.2 Spiky pressure-gradient kernel

grad_i W_spiky = -45/(pi h^6) * ((h-r)^2/r) * r_ij, for 0 < r < h. At r = 0, the implementation returns zero to avoid division by zero.

### 2.3 Viscosity Laplacian

laplacian W_visc = 45/(pi h^6) * (h-r), for 0 <= r < h; otherwise 0.

## 3. Density estimation

rho_i = sum_j m_j W_poly6(norm(x_i-x_j), h).

The sum includes the particle itself. The solver uses a 3D uniform spatial hash to find candidate neighbors.

For cell width h, the integer cell coordinate is
c(x) = (floor(x/h), floor(y/h), floor(z/h)).

Each particle is inserted into the cell keyed by c(x_i). To query neighbors of particle i, inspect its own cell and the 26 immediately adjacent cells (offsets delta_x, delta_y, delta_z in {-1,0,1}). Because every cell has width h, any particle with Euclidean distance less than h must be in one of those 27 cells. The kernel is still evaluated with the exact distance, so cell membership only filters candidates; it does not approximate the SPH equations.

Expected work is O(N + sum_i k_i) for grid construction and queries, where k_i is the number of candidate particles in the 27 cells around particle i. This is typically much cheaper than O(N^2) for spatially distributed particles, but can still approach O(N^2) when many particles occupy the same few cells. Grid coordinates are checked for integer-range overflow.

## 4. Weakly-compressible equation of state

p_i = c_0^2 (rho_i - rho_0).

The current baseline clamps pressure to non-negative values: p_i = max(0, c_0^2 (rho_i-rho_0)). This is a stabilization/modeling choice, not a universal SPH rule.

## 5. Pressure acceleration

a_i^p = -sum_(j != i) m_j (p_i/rho_i^2 + p_j/rho_j^2) grad_i W_spiky(r_ij,h).

This symmetric pair form supports linear-momentum conservation when pair evaluations are consistent. Density denominators use a positive numerical floor.

## 6. Viscous acceleration

a_i^nu = mu sum_(j != i) (m_j/rho_j) (v_j-v_i) laplacian W_visc(r_ij,h).

Here mu is the coefficient used by this simplified discretization; its physical interpretation depends on the SPH formulation. Do not treat it as interchangeable with dynamic viscosity without deriving the units and density factors.

## 7. Momentum, gravity, and integration

dv_i/dt = a_i^p + a_i^nu + g; dx_i/dt = v_i.

Semi-implicit Euler updates velocity first, then position:
v_i^(n+1) = v_i^n + dt*a_i^n
x_i^(n+1) = x_i^n + dt*v_i^(n+1)

Gravity is a 3D acceleration vector. The default convention is positive Y upward, so default gravity is (0,-9.81,0).

## 8. Domain boundaries

After integration, clamp each coordinate to the corresponding axis-aligned domain bounds. If a particle crosses a wall, reverse the outward velocity component using v_n = -e*v_n, where restitution e is in [0,1]. This is a simple baseline; robust liquid-wall interaction needs validated boundary particles or another boundary formulation.

## 9. Safeguards and validation

- Reject non-finite parameters and state, non-positive smoothing length, rest density, timestep, and particle mass.
- Avoid evaluating the spiky gradient at zero separation.
- Protect density denominators from zero.
- Check finite positions, velocities, densities, pressures, and accelerations after each step.
- Test kernel support, self-density, gravity, 3D interactions, boundaries, and timestep sensitivity.
- Track density deviation abs(rho_i-rho_0)/rho_0 and total momentum in later validation.
- The uniform spatial hash is built from current positions once per solver step and reused for both density and force evaluation; positions are only updated after force accumulation. The standalone `compute_density_and_pressure()` call builds its own grid. This reuse reduces grid-construction overhead without changing the kernel support test or SPH equations.

## 10. Deferred work

Surface tension, incompressible pressure projection, adaptive timesteps, advanced boundary conditions, GPU compute, and accelerated neighbor search remain future work. Add the equations, assumptions, and tests to this file before implementing each method.

## 11. Repeatable scaling and neighbor-search benchmark

The `gifpy_benchmark` executable measures wall-clock time for complete solver steps at several particle counts. Initial positions come from a deterministic jittered 3D lattice, so repeated runs use the same starting distribution. The timer includes spatial-grid rebuilding, density/pressure evaluation, force accumulation, integration, and boundary handling.

For each particle count, it also compares the production spatial-hash density/pressure pass against a local brute-force density reference using the exact same initial particle positions, masses, smoothing length, and Poly6 kernel. It reports both timings, the reference/grid timing ratio, and the maximum relative density error. The benchmark fails if that error exceeds `1e-12`; the brute-force path is only a validation/performance reference and is not used by the production solver.

Run it from a configured CMake build directory:

```sh
./gifpy_benchmark
./gifpy_benchmark --counts 64,256,1024,4096 --steps 20
```

The output includes full-step timing CSV rows and density-comparison rows. Compare timings on the same machine/build configuration; they are measurements rather than performance guarantees, and CI only runs a small smoke configuration. The CLI caps particle counts at 8,192 because the brute-force density reference has O(N^2) work. For larger simulations, use a separate production-solver benchmark without the reference comparison.
