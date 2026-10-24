# GIFPY — Liquid Simulation Formula Reference

This document is the mathematical source of truth for the first solver: a **2D weakly-compressible Smoothed Particle Hydrodynamics (WCSPH)** reference implementation. The first implementation is intentionally small and testable; it is not yet a production-quality free-surface renderer.

## 1. State and units

For particle \(i\):
- position \(\mathbf{x}_i\): m
- velocity \(\mathbf{v}_i\): m s\(^{-1}\)
- mass \(m_i\): kg
- density \(\rho_i\): kg m\(^{-2}\) in this 2D model
- pressure \(p_i\): kg s\(^{-2}\) (2D pressure units)
- acceleration \(\mathbf{a}_i\): m s\(^{-2}\)
- smoothing length \(h\): m
- timestep \(\Delta t\): s

The 2D solver uses consistent 2D kernel normalization. Parameters such as reference density and viscosity are model inputs; timestep and kernel support are numerical choices. Do not change physical constants merely to hide numerical instability.

## 2. Kernel functions

Let \(\mathbf{r}_{ij}=\mathbf{x}_i-\mathbf{x}_j\), \(r=\|\mathbf{r}_{ij}\|\), and \(q=r/h\). Kernels are zero for \(r\ge h\).

### 2.1 Poly6 density kernel (2D)

\[
W_{poly6}(r,h)=\frac{4}{\pi h^8}(h^2-r^2)^3,\quad 0\le r<h
\]

Its gradient:
\[
\nabla_i W_{poly6}=-\frac{24}{\pi h^8}(h^2-r^2)^2\mathbf{r}_{ij}
\]

### 2.2 Spiky pressure-gradient kernel (2D)

\[
\nabla_i W_{spiky}= -\frac{30}{\pi h^5}\frac{(h-r)^2}{r}\mathbf{r}_{ij},\quad 0<r<h
\]

The gradient is defined as zero at \(r=0\) to avoid division by zero.

### 2.3 Viscosity Laplacian (2D)

\[
\nabla^2 W_{visc}=\frac{40}{\pi h^5}(h-r),\quad 0\le r<h
\]

These kernels have distinct roles: poly6 estimates density, spiky gradient computes pressure force, and the viscosity Laplacian computes viscous diffusion.

## 3. Density estimation

Include the particle itself in the neighbor sum:
\[
\rho_i=\sum_j m_j W_{poly6}(\|\mathbf{x}_i-\mathbf{x}_j\|,h)
\]

This sum is over all particles within support radius \(h\), including \(i=j\). The initial reference implementation uses a simple all-pairs neighbor loop for clarity and small test cases; replace it with a uniform grid or spatial hash before scaling to large particle counts.

## 4. Equation of state

The initial linear weakly-compressible equation of state is:
\[
p_i=c_0^2(\rho_i-\rho_0)
\]
where \(\rho_0\) is rest density and \(c_0\) is artificial speed of sound. For stability in free-surface cases, pressure may be clamped to non-negative values in this first baseline:
\[
p_i=\max(0,c_0^2(\rho_i-\rho_0))
\]
This clamp is a modeling/stabilization choice, not a universal SPH law. The C++ implementation documents it explicitly.

## 5. Pressure acceleration

A symmetric pressure-force discretization is:
\[
\mathbf{a}^{p}_i=-\sum_{j\ne i}m_j\left(\frac{p_i}{\rho_i^2}+\frac{p_j}{\rho_j^2}\right)\nabla_i W_{spiky}(r_{ij},h)
\]
The symmetric pair form helps conserve linear momentum when pairs and gradients are evaluated consistently. Density denominators are protected by a small positive numerical floor.

## 6. Viscous acceleration

\[
\mathbf{a}^{\nu}_i=\mu\sum_{j\ne i}\frac{m_j}{\rho_j}(\mathbf{v}_j-\mathbf{v}_i)\nabla^2 W_{visc}(r_{ij},h)
\]
Here \(\mu\) is the model's kinematic-viscosity-like coefficient used by this discretization. Its precise interpretation depends on the selected SPH formulation; do not assume it is interchangeable with dynamic viscosity without deriving the density factors and units for a chosen formulation.

## 7. Momentum and gravity

\[
\frac{d\mathbf{v}_i}{dt}=\mathbf{a}^{p}_i+\mathbf{a}^{\nu}_i+\mathbf{g}
\]
\[
\frac{d\mathbf{x}_i}{dt}=\mathbf{v}_i
\]

Gravity is an acceleration vector in m s\(^{-2}\). In a screen-like coordinate system where positive y points downward, use \(\mathbf{g}=(0,g)\); otherwise use the coordinate convention consistently across simulation and rendering.

## 8. Time integration

The baseline uses semi-implicit (symplectic) Euler:
\[
\mathbf{v}^{n+1}_i=\mathbf{v}^{n}_i+\Delta t\,\mathbf{a}^{n}_i
\]
\[
\mathbf{x}^{n+1}_i=\mathbf{x}^{n}_i+\Delta t\,\mathbf{v}^{n+1}_i
\]
It is generally more stable for mechanical systems than explicit position Euler, but does not guarantee stability for arbitrary \(\Delta t\). Compute forces from the current state, update velocity, then position.

## 9. Domain boundaries

The initial demo uses axis-aligned box boundaries. After position integration, clamp each particle to the domain and reflect only inward-violating velocity:
\[
v_n\leftarrow -e v_n
\]
where \(e\in[0,1]\) is coefficient of restitution. This simple boundary treatment is for baseline tests only; robust liquid-wall interaction needs boundary particles or a validated boundary formulation.

## 10. Numerical safeguards and validation

- Reject non-finite input parameters and non-positive \(h\), \(\rho_0\), and \(\Delta t\).
- Never evaluate the spiky gradient at \(r=0\).
- Protect density denominators from zero.
- Require finite positions, velocities, densities, pressures, and accelerations after each step.
- Test kernel support (zero at/outside \(h\)), self-density contribution, finite forces, gravity-only acceleration, and bounded-domain behavior.
- Track timestep sensitivity, density deviation \(\epsilon_\rho=|\rho_i-\rho_0|/\rho_0\), and total momentum in later validation.
- The first implementation uses all-pairs neighbor search, \(O(N^2)\); it is for a small reference solver only. A uniform grid/spatial hash is a required optimization before large scenes.

## 11. Explicitly deferred

Surface tension, incompressible pressure projection, adaptive timesteps, advanced boundary conditions, 3D kernels, GPU compute, and high-performance neighbor structures are not part of this baseline. Each must have its equations, assumptions, and tests documented before implementation.
