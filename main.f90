program gfoil36
  use cufft
  implicit none



## Grid preprocessing (run once, offline)
- Load Plot3D `grid.xyz` — read NXI×NETA node coordinates x(j,i), y(j,i)
- Compute forward metrics: $x_\xi, x_\eta, y_\xi, y_\eta$ via 2nd-order central differences
- Compute Jacobian: $J = x_\xi y_\eta - x_\eta y_\xi$
- Compute inverse metrics: $\xi_x, \xi_y, \eta_x, \eta_y$
- Compute metric tensor: $g^{11}, g^{22}, g^{12}$
- Compute physical spacings: $d\xi_s(j,i)$, $d\eta_s(j,i)$, $\Delta z$ uniform
- Validate: J > 0, skewness < 0.15, y+ ~ 1–2, wake cut gap < 1e-6
- Estimate dt from min cell size

## Initialisation (once per run)
- Allocate 3D arrays (NZ, NETA, NXI) for u, v, w, p and temporaries
- Set freestream IC: u=cos(AoA), v=sin(AoA), w=0, p=0 at all (k,j,i)
- Add small random perturbation over full 3D domain
- Zero RK accumulators ru, rv, rw
- Compute Poisson eigenvalues: DST-I in ξ (NXI modes), DFT in z (NZ modes)
- Create cuFFT plans (batched 1D transforms)
- All arrays already on GPU via managed memory

## Time loop (repeat NSTEPS times)

### RK3 substeps (inner loop, 3 times per timestep)

- **Momentum RHS** — OpenACC kernel over (k,j,i) interior
  - Contravariant velocities: $U = \xi_x u + \xi_y v$, $V = \eta_x u + \eta_y v$, $W = w$
  - Convection u,v,w: 3rd-order upwind in ξ and η (using U,V), central in z (using W)
  - Diffusion u,v,w: FD2 central with $g^{11}$, $g^{22}$, $g^{12}$ cross terms, plus $\Delta z^2$ in z
  - Pressure gradient in physical space: $\partial p/\partial x = \xi_x \partial p/\partial\xi + \eta_x \partial p/\partial\eta$, $\partial p/\partial z$ central
  - RK accumulator update: ru = rka\*ru + RHS\_u (and same for v, w)
  - Intermediate velocity: u\* = u + dt\*rkb\*ru (same for v\*, w\*)

- **Tripping** (first NTRIP steps, first substep only)
  - Broadband spanwise-varying forcing near x = TRIP\_X
  - Applied to u\* and w\* only

- **Boundary conditions on u\*, v\*, w\***
  - Wall j=1: no-slip, u\*=v\*=w\*=0
  - Far field j=NETA: freestream u\*=cos(AoA), v\*=sin(AoA), w\*=0
  - Inlet i=1: freestream
  - Outlet i=NXI: zero-gradient (copy from i=NXI-1)
  - Spanwise k=1,NZ: periodic
  - Wake cut i=1/i=NXI: match values across cut

- **Pressure Poisson solve** (final RK substep only)
  - Compute Rhie-Chow face velocities: $\hat{U}_{i+1/2}$, $\hat{V}_{j+1/2}$, $\hat{W}_{k+1/2}$ — interpolate cell velocities to faces with pressure gradient correction to suppress checkerboard
  - Compute divergence of u\* using Rhie-Chow face velocities: $\nabla \cdot \mathbf{u}^* / \Delta t$
  - Forward transform RHS: DST-I in ξ, DFT in z (cuFFT)
  - Spectral division by $(eig\_\xi(i) + eig\_z(k))$ — approximate, uniform eigenvalues
  - Iterative metric correction (NPCORR=3 iterations):
    - Compute cross-term residual using current φ: $-2g^{12} \partial^2\phi/\partial\xi\partial\eta$
    - Add residual to RHS
    - Re-solve with forward/inverse FFT
  - Inverse transform: IDFT in z, IDST-I in ξ (cuFFT)
  - Normalise by FFT scaling factors

- **Velocity correction**
  - $u = u^* - \Delta t (\xi_x \partial\phi/\partial\xi + \eta_x \partial\phi/\partial\eta)$
  - $v = v^* - \Delta t (\xi_y \partial\phi/\partial\xi + \eta_y \partial\phi/\partial\eta)$
  - $w = w^* - \Delta t \partial\phi/\partial z$
  - Pressure update: $p = p + \phi$
  - Re-apply BCs on corrected u, v, w

- **Advance time**: t = t + dt (after all 3 substeps)

## Force monitoring (every NPRINT steps)
- Download wall pressure from GPU: p(k,1,i)
- Z-average wall pressure: $\bar{p}(i) = \frac{1}{N_z}\sum_k p(k,1,i)$
- Integrate CL, CD via wall-normal pressure: outward normal from $(\eta_x, \eta_y)$ at j=1
- Count CL zero-crossings for Strouhal number estimate
- Print step, t, CL, CD to stdout

## Finalisation
- Destroy cuFFT plans
- Deallocate all arrays
- Print Strouhal number estimate




end program gfoil3d
