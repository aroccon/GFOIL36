
!## Grid preprocessing (run once, offline)
!- Load Plot3D `grid.xyz` — read NXI×NETA node coordinates x(j,i), y(j,i)
!- Compute forward metrics: $x_\xi, x_\eta, y_\xi, y_\eta$ via 2nd-order central differences
!- Compute Jacobian: $J = x_\xi y_\eta - x_\eta y_\xi$
!- Compute inverse metrics: $\xi_x, \xi_y, \eta_x, \eta_y$
!- Compute metric tensor: $g^{11}, g^{22}, g^{12}$
!- Compute physical spacings: $d\xi_s(j,i)$, $d\eta_s(j,i)$, $\Delta z$ uniform
!- Validate: J > 0, skewness < 0.15, y+ ~ 1–2, wake cut gap < 1e-6
!- Estimate dt from min cell size

!## Initialisation (once per run)
!!- Allocate 3D arrays (NZ, NETA, NXI) for u, v, w, p and temporaries
!- Set freestream IC: u=cos(AoA), v=sin(AoA), w=0, p=0 at all (k,j,i)
!- Add small random perturbation over full 3D domain
!- Zero RK accumulators ru, rv, rw
!- Compute Poisson eigenvalues: DST-I in ξ (NXI modes), DFT in z (NZ modes)
!- Create cuFFT plans (batched 1D transforms)
!- All arrays already on GPU via managed memory

!## Time loop (repeat NSTEPS times)

!### RK3 substeps (inner loop, 3 times per timestep)

!- **Momentum RHS** — OpenACC kernel over (k,j,i) interior
!  - Contravariant velocities: $U = \xi_x u + \xi_y v$, $V = \eta_x u + \eta_y v$, $W = w$
!  - Convection u,v,w: 3rd-order upwind in ξ and η (using U,V), central in z (using W)
!  - Diffusion u,v,w: FD2 central with $g^{11}$, $g^{22}$, $g^{12}$ cross terms, plus $\Delta z^2$ in z
!  - Pressure gradient in physical space: $\partial p/\partial x = \xi_x \partial p/\partial\xi + \eta_x \partial p/\partial\eta$, $\partial p/\partial z$ central
!  - RK accumulator update: ru = rka\*ru + RHS\_u (and same for v, w)
!  - Intermediate velocity: u\* = u + dt\*rkb\*ru (same for v\*, w\*)

!- **Tripping** (first NTRIP steps, first substep only)
!  - Broadband spanwise-varying forcing near x = TRIP\_X
!  - Applied to u\* and w\* only

!- **Boundary conditions on u\*, v\*, w\***
!  - Wall j=1: no-slip, u\*=v\*=w\*=0
!  - Far field j=NETA: freestream u\*=cos(AoA), v\*=sin(AoA), w\*=0
!  - Inlet i=1: freestream
!  - Outlet i=NXI: zero-gradient (copy from i=NXI-1)
!  - Spanwise k=1,NZ: periodic
!  - Wake cut i=1/i=NXI: match values across cut

!- **Pressure Poisson solve** (final RK substep only)
!  - Compute Rhie-Chow face velocities: $\hat{U}_{i+1/2}$, $\hat{V}_{j+1/2}$, $\hat{W}_{k+1/2}$ — interpolate cell velocities to faces with pressure gradient correction to suppress checkerboard
!  - Compute divergence of u\* using Rhie-Chow face velocities: $\nabla \cdot \mathbf{u}^* / \Delta t$
!  - Forward transform RHS: DST-I in ξ, DFT in z (cuFFT)
!  - Spectral division by $(eig\_\xi(i) + eig\_z(k))$ — approximate, uniform eigenvalues
!  - Iterative metric correction (NPCORR=3 iterations):
!    - Compute cross-term residual using current φ: $-2g^{12} \partial^2\phi/\partial\xi\partial\eta$
!    - Add residual to RHS
!    - Re-solve with forward/inverse FFT
!  - Inverse transform: IDFT in z, IDST-I in ξ (cuFFT)
!  - Normalise by FFT scaling factors

!- **Velocity correction**
!  - $u = u^* - \Delta t (\xi_x \partial\phi/\partial\xi + \eta_x \partial\phi/\partial\eta)$
!  - $v = v^* - \Delta t (\xi_y \partial\phi/\partial\xi + \eta_y \partial\phi/\partial\eta)$
!  - $w = w^* - \Delta t \partial\phi/\partial z$
!  - Pressure update: $p = p + \phi$
!  - Re-apply BCs on corrected u, v, w

!- **Advance time**: t = t + dt (after all 3 substeps)!

!## Force monitoring (every NPRINT steps)
!- Download wall pressure from GPU: p(k,1,i)
!- Z-average wall pressure: $\bar{p}(i) = \frac{1}{N_z}\sum_k p(k,1,i)$
!- Integrate CL, CD via wall-normal pressure: outward normal from $(\eta_x, \eta_y)$ at j=1
!- Count CL zero-crossings for Strouhal number estimate
!- Print step, t, CL, CD to stdout

!## Finalisation
!- Destroy cuFFT plans
!- Deallocate all arrays
!- Print Strouhal number estimate


program gfoil36
  implicit none
  ! parameters
  integer, parameter :: nx = 512,       ! xi  (streamwise / C-grid wrap)
  integer, parameter :: ny = 256        ! eta (wall-normal, j=1 wall, j=NJ far)
  integer, parameter :: nz = 128        ! zeta (spanwise, periodic)
  real(8), parameter :: re   = 50000d0  ! Reynolds number
  real(8), parameter :: aoa  = 5d0      ! angle of attack (degrees)
  real(8), parameter :: lz   = 0.2d0    ! spanwise length
  integer, parameter :: nsteps = 100    ! steps to run 
  integer, parameter :: nprint = 10     ! print interval
  real(8), parameter :: PI      = acos(-1d0)
  real(8), parameter :: nu      = 1d0 / RE
  real(8), parameter :: dz      = LZ / dble(NK)
  real(8), parameter :: dzi     = 1d0/dz
  real(8), parameter :: ddzi    = 1d0/dz/fz
  real(8), parameter :: AOA_RAD = AOA * PI / 180d0
  real(8), parameter :: u0      = dcos(AOA_RAD)
  real(8), parameter :: v0      = dsin(AOA_RAD)
  real(8), parameter :: rka(3) = [ 0d0,        -17d0/60d0, -5d0/12d0 ]
  real(8), parameter :: rkb(3) = [ 8d0/15d0,    5d0/12d0,   3d0/4d0  ]
  double precision, allocatable :: xg(:,:), yg(:,:)
  double precision, allocatable :: xi_x(:,:), xi_y(:,:)               ! d(xi)/dx,  d(xi)/dy
  double precision, allocatable :: et_x(:,:), et_y(:,:)               ! d(eta)/dx, d(eta)/dy
  double precision, allocatable :: jac(:,:)                           ! Jacobian J = x_xi*y_et - x_et*y_xi
  double precision, allocatable :: g11(:,:), g22(:,:), g12(:,:)       ! xi_x^2  + xi_y^2 ---  et_x^2  + et_y^2 -- xi_x*et_x + xi_y*et_y
  double precision, allocatable :: dxi(:,:),det(:,:)                  ! physical xi  cell spacing and physical eta cell spacing
  double precision, allocatable :: u(:,:,:), v(:,:,:), w(:,:,:)
  double precision, allocatable :: p(:,:,:)
  double precision, allocatable :: ru(:,:,:), rv(:,:,:), rw(:,:,:)
  double precision, allocatable :: us(:,:,:), vs(:,:,:), ws(:,:,:)
  double precision, allocatable :: phi(:,:,:), div(:,:,:)
  integer :: i, j, k, rk, step, t, dt, stage
    ! Metric temporaries inside loops
  real(8) :: x_xi, x_et, y_xi, y_et, jac_loc
  real(8) :: dxip, dxim, detp, detm     ! neighbour spacings
  ! RHS components
  double precision :: rhs_u, rhs_v, rhs_w, conv_u, conv_v, conv_w
  double precision :: diff_u, diff_v, diff_w, cross_u, cross_v, cross_w
  ! Contravariant velocities
  real(8) :: U_con, V_con               ! U = xi_x*u + xi_y*v, V = et_x*u + et_y*v
  ! Rhie-Chow face velocities
  real(8) :: Uf_ip, Uf_im               ! contravariant U at i+1/2 and i-1/2
  real(8) :: Vf_jp, Vf_jm               ! contravariant V at j+1/2 and j-1/2
  real(8) :: Wf_kp, Wf_km               ! w at k+1/2 and k-1/2
  ! Pressure gradient components
  double precision :: dp_dx, dp_dy, dp_dz, dphi_xi, dphi_et, dphi_dz, dphi_dx, dphi_dy, dmin
 
  ! allocate arrays
  allocate( xg(nx,ny), yg(nx,ny) )
  allocate( xi_x(nx,ny), xi_y(nx,ny) )
  allocate( et_x(nx,ny), et_y(nx,ny) )
  allocate( jac(nx,ny) )
  allocate( g11(nx,ny), g22(nx,ny), g12(nx,ny) )
  allocate( dxi(nx,ny), det(nx,ny) )
  allocate( u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz) )
  allocate( p(nx,ny,nz) )
  allocate( ru(nx,ny,nz), rv(nx,ny,nz), rw(nx,ny,nz) )
  allocate( us(nx,ny,nz), vs(nx,ny,nz), ws(nx,ny,nz) )
  allocate( phi(nx,ny,nz), div(nx,ny,nz) )

  ! 1. PLACEHOLDER GRID  (uniform Cartesian — replace with Construct2D reader)
  ! This makes the code testable on a channel flow immediately, at AoA=0 on a Cartesian grid: xi_x=1/dx, et_y=1/dy, g12=0 exactly.
  block
    real(8) :: dx_c, dy_c
    integer :: ii, jj
    dx_c = 1d0 / dble(nx)
    dy_c = 1d0 / dble(ny)
    do ii = 1, nx
      do jj = 1, ny
        xg(ii,jj) = (dble(ii) - 0.5d0) * dx_c
        yg(ii,jj) = (dble(jj) - 0.5d0) * dy_c
      end do
    end do
  end block

  ! 2. COMPUTE METRICS FROM x(i,j), y(i,j)
  do i = 1, nx
    do j = 1, ny
      ! d/dxi — central, one-sided at i=1 and i=NI
      if (i == 1) then
        x_xi = xg(2,j)  - xg(1,j)
        y_xi = yg(2,j)  - yg(1,j)
      else if (i == nx) then
        x_xi = xg(nx,j) - xg(nx-1,j)
        y_xi = yg(nx,j) - yg(nx-1,j)
      else
        x_xi = (xg(i+1,j) - xg(i-1,j)) * 0.5d0
        y_xi = (yg(i+1,j) - yg(i-1,j)) * 0.5d0
      end if
      ! d/deta — central, one-sided at j=1 and j=NJ
      if (j == 1) then
        x_et = xg(i,2)  - xg(i,1)
        y_et = yg(i,2)  - yg(i,1)
      else if (j == NJ) then
        x_et = xg(i,NJ) - xg(i,NJ-1)
        y_et = yg(i,NJ) - yg(i,NJ-1)
      else
        x_et = (xg(i,j+1) - xg(i,j-1)) * 0.5d0
        y_et = (yg(i,j+1) - yg(i,j-1)) * 0.5d0
      end if
      ! Jacobian
      jac_loc   = x_xi * y_et - x_et * y_xi
      jac(i,j)  = jac_loc
      ! Inverse metrics
      xi_x(i,j) =  y_et   / jac_loc
      xi_y(i,j) = -x_et   / jac_loc
      et_x(i,j) = -y_xi   / jac_loc
      et_y(i,j) =  x_xi   / jac_loc
      ! Contravariant metric tensor
      g11(i,j)  = xi_x(i,j)**2 + xi_y(i,j)**2
      g22(i,j)  = et_x(i,j)**2 + et_y(i,j)**2
      g12(i,j)  = xi_x(i,j)*et_x(i,j) + xi_y(i,j)*et_y(i,j)
      ! Physical cell spacings (arc length)
      dxi(i,j)  = sqrt(x_xi**2 + y_xi**2)
      det(i,j)  = sqrt(x_et**2 + y_et**2)
    end do
  end do


  ! 4. initialize flow field (freestream + tiny random perturbation)
  call random_seed()
  block
    real(8) :: noise
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          !call random_number(noise)
          u(i,j,k) = U0 !+ (noise - 0.5d0) * 1d-3
          !call random_number(noise)
          v(i,j,k) = V0 !+ (noise - 0.5d0) * 1d-3
          !call random_number(noise)
          w(i,j,k) =   !   (noise - 0.5d0) * 1d-3
          p(i,j,k) = 0d0
          ru(i,j,k) = 0d0
          rv(i,j,k) = 0d0
          rw(i,j,k) = 0d0
        end do
      end do
    end do
  end block


  ! start temporal loop
  do step = 1, nsteps
    ! start RK3 substeps, projection 
    ! compute rhs = - advection + nu * diffusion
    ! advection: 3rd-order upwind in xi and eta, 2nd-order central in z
    ! diffusion:  2nd-order central with metric tensor
    do stage = 1, 3
      do k = 1, nz
        do j = 2, nj-1
          do i = 2, nx-1
            ! Local spacings
            dxip = 0.5d0 * (dxi(i,j) + dxi(i+1,j))
            dxim = 0.5d0 * (dxi(i,j) + dxi(i-1,j))
            detp = 0.5d0 * (det(i,j) + det(i,j+1))
            detm = 0.5d0 * (det(i,j) + det(i,j-1))

            ! Contravariant velocities at cell centre
            U_con = xi_x(i,j)*u(i,j,k) + xi_y(i,j)*v(i,j,k)
            V_con = et_x(i,j)*u(i,j,k) + et_y(i,j)*v(i,j,k)

            !--- advection of u ---
            ! xi direction: 3rd-order upwind
            if (U_con >= 0d0) then
              conv_u=U_con*(2d0*u(i+1,j,k)-6d0*u(i,j,k)+3d0*u(i-1,j,k)+u(i-2,j,k))/(6d0*dxi(i,j))
            else
              conv_u=U_con*(-u(i+2,j,k)+3d0*u(i+1,j,k)+6d0*u(i,j,k)-2d0*u(i-1,j,k))/(6d0*dxi(i,j))*(-1d0)
            end if
            ! eta direction: 3rd-order upwind
            if (V_con >= 0d0) then
              conv_u=conv_u+V_con*( 2d0*u(i,j+1,k)-6d0*u(i,j,k)+3d0*u(i,j-1,k)+u(i,j-2,k))/(6d0*det(i,j))
            else
              conv_u=conv_u+V_con*(-u(i,j+2,k)+3d0*u(i,j+1,k)+6d0*u(i,j,k)-2d0*u(i,j-1,k))/(6d0*det(i,j))*(-1d0)
            end if
            ! z direction: 2nd-order central (periodic)
            conv_u =conv_u+w(i,j,k)*(u(i,j,kp1(k))-u(i,j,km1(k)))*0.5d0*dzi

            !--- Convection of v (same pattern) ---
            if (U_con >= 0d0) then
              conv_v=U_con*( 2d0*v(i+1,j,k)-6d0*v(i,j,k)+3d0*v(i-1,j,k)+v(i-2,j,k))/(6d0*dxi(i,j))
            else
              conv_v=U_con*(-v(i+2,j,k)+3d0*v(i+1,j,k)+6d0*v(i,j,k)-2d0*v(i-1,j,k))/(6d0*dxi(i,j))*(-1d0)
            end if
            if (V_con >= 0d0) then
              conv_v=conv_v+V_con*( 2d0*v(i,j+1,k)-6d0*v(i,j,k)+3d0*v(i,j-1,k)+v(i,j-2,k))/(6d0*det(i,j))
            else
              conv_v=conv_v+V_con*(-v(i,j+2,k)+3d0*v(i,j+1,k)+6d0*v(i,j,k)-2d0*v(i,j-1,k))/(6d0*det(i,j))*(-1d0)
            end if
            conv_v=conv_v+w(i,j,k)*(v(i,j,kp1(k))-v(i,j,km1(k)))*0.5d0*dzi

            !--- Convection of w ---
            if (U_con >= 0d0) then
              conv_w=U_con*(2d0*w(i+1,j,k)-6d0*w(i,j,k)+3d0*w(i-1,j,k)+w(i-2,j,k))/(6d0*dxi(i,j))
            else
              conv_w=U_con*(-w(i+2,j,k)+3d0*w(i+1,j,k)+6d0*w(i,j,k)-2d0*w(i-1,j,k))/(6d0*dxi(i,j))*(-1d0)
            end if
            if (V_con >= 0d0) then
              conv_w=conv_w+V_con*(2d0*w(i,j+1,k)-6d0*w(i,j,k)+3d0*w(i,j-1,k)+w(i,j-2,k))/(6d0*det(i,j))
            else
              conv_w=conv_w+V_con*(-w(i,j+2,k)+3d0*w(i,j+1,k)+6d0*w(i,j,k)-2d0*w(i,j-1,k))/(6d0*det(i,j))*(-1d0)
            end if
            conv_w=conv_w+w(i,j,k)*(w(i,j,kp1(k))-w(i,j,km1(k)))*0.5d0*dzi

            !--- Diffusion of u ---
            diff_u =          g11(i,j)*(u(i+1,j,k) - 2d0*u(i,j,k) + u(i-1,j,k))/dxi(i,j)**2                                    ! xi term:  g11 * d2u/dxi2
            diff_u = diff_u + g22(i,j)*(u(i,j+1,k) - 2d0*u(i,j,k) + u(i,j-1,k))/det(i,j)**2                                    ! eta term: g22 * d2u/deta2
            diff_u = diff_u + (u(i,j,kp1(k)) - 2d0*u(i,j,k) + u(i,j,km1(k)))*ddzi                                              ! z term: d2u/dz2 (Cartesian)
            cross_u = 2d0*g12(i,j)*(u(i+1,j+1,k) - u(i+1,j-1,k) - u(i-1,j+1,k) + u(i-1,j-1,k))/(4d0 * dxi(i,j) * det(i,j))     ! cross metric term: 2*g12 * d2u/dxi/deta
            diff_u = diff_u + cross_u
            !--- Diffusion of v ---
            diff_v =          g11(i,j)*(v(i+1,j,k) - 2d0*v(i,j,k) + v(i-1,j,k))/dxi(i,j)**2 
            diff_v = diff_v + g22(i,j)*(v(i,j+1,k) - 2d0*v(i,j,k) + v(i,j-1,k))/det(i,j)**2 
            diff_v = diff_v + (v(i,j,kp1(k)) - 2d0*v(i,j,k) + v(i,j,km1(k)))*ddzi 
            cross_v = 2d0*g12(i,j)*(v(i+1,j+1,k) - v(i+1,j-1,k) - v(i-1,j+1,k) + v(i-1,j-1,k))/(4d0 * dxi(i,j) * det(i,j))
            diff_v = diff_v + cross_v
            !--- Diffusion of w ---
            diff_w =          g11(i,j)*(w(i+1,j,k) - 2d0*w(i,j,k) + w(i-1,j,k))/dxi(i,j)**2 
            diff_w = diff_w + g22(i,j)*(w(i,j+1,k) - 2d0*w(i,j,k) + w(i,j-1,k))/det(i,j)**2 
            diff_w = diff_w + (w(i,j,kp1(k)) - 2d0*w(i,j,k) + w(i,j,km1(k)))*ddzi 
            cross_w = 2d0*g12(i,j)*(w(i+1,j+1,k) - w(i+1,j-1,k) - w(i-1,j+1,k) + w(i-1,j-1,k))/(4d0 * dxi(i,j) * det(i,j))
            diff_w = diff_w + cross_w
            !--- RHS assembly ---
            rhs_u = -conv_u + nu*diff_u
            rhs_v = -conv_v + nu*diff_v
            rhs_w = -conv_w + nu*diff_w
            !--- RK3 accumulator update ---
            ru(i,j,k) = rka(stage)*ru(i,j,k) + rhs_u
            rv(i,j,k) = rka(stage)*rv(i,j,k) + rhs_v
            rw(i,j,k) = rka(stage)*rw(i,j,k) + rhs_w
            ! compute prijectd velocity u*, v*, w* for Poisson RHS
            us(i,j,k) = u(i,j,k) + dt*rkb(stage)*ru(i,j,k)
            vs(i,j,k) = v(i,j,k) + dt*rkb(stage)*rv(i,j,k)
            ws(i,j,k) = w(i,j,k) + dt*rkb(stage)*rw(i,j,k)
          end do
        end do
      end do
      !$acc end parallel loop

      ! 5b. Apply BCs on u*, v*, w*
      call apply_bcs(us, vs, ws, nx, ny, nz, U0, V0)

    end do  ! rk loop

    !-----------------------------------------------------------------------
    ! 5c. PRESSURE POISSON (once per time step, after all RK substeps)
    !     Solve: laplacian(phi) = (1/dt) * div(u*)
    !     where div is computed with Rhie-Chow face velocities
    !-----------------------------------------------------------------------

    !--- Rhie-Chow divergence ---
    call rhie_chow_div(us, vs, ws, p, phi, div, &
                       xi_x, xi_y, et_x, et_y, jac, g11, g22, g12, &
                       dxi, det, nx, ny, nz, DZI, dt)

    !--- Poisson solve (placeholder) ---
    call poisson_solve(phi, div, nx, ny, nz, DZ)

    !--- Velocity correction ---
    do k = 1, nz
      do j = 2, ny-1
        do i = 2, nx-1
          ! Pressure correction gradient in computational space
          dphi_xi = (phi(i+1,j,k) - phi(i-1,j,k)) / (2d0 * dxi(i,j))
          dphi_et = (phi(i,j+1,k) - phi(i,j-1,k)) / (2d0 * det(i,j))
          dphi_dz = (phi(i,j,kp1(k)) - phi(i,j,km1(k))) * 0.5d0 * DZI
          ! Transform to physical space
          dphi_dx = xi_x(i,j)*dphi_xi + et_x(i,j)*dphi_et
          dphi_dy = xi_y(i,j)*dphi_xi + et_y(i,j)*dphi_et
          ! Correct velocity
          u(i,j,k) = us(i,j,k) - dt * dphi_dx
          v(i,j,k) = vs(i,j,k) - dt * dphi_dy
          w(i,j,k) = ws(i,j,k) - dt * dphi_dz
          ! Update pressure
          p(i,j,k) = p(i,j,k) + phi(i,j,k)

        end do
      end do
    end do
    !$acc end parallel loop

    ! Re-apply BCs on corrected velocity
    call apply_bcs(u, v, w, nx, ny, nz, U0, V0)

    t = t + dt
  end do  ! time loop

  write(*,'(A)') '  Done.'






contains
  !==========================================================================
  ! Periodic index helpers (inline, no overhead)
  !==========================================================================
  pure integer function kp1(k)
    integer, intent(in) :: k
    kp1 = mod(k, NK) + 1
  end function

  pure integer function km1(k)
    integer, intent(in) :: k
    km1 = mod(k + NK - 2, NK) + 1
  end function

  !==========================================================================
  ! BOUNDARY CONDITIONS
  !   Wall    j=1:    no-slip   u=v=w=0
  !   Farfield j=NJ:  freestream u=U0, v=V0, w=0
  !   Inlet   i=1:    freestream
  !   Outlet  i=NI:   zero-gradient (copy from i=NI-1)
  !   Spanwise:       periodic (handled by kp1/km1)
  !   Wake cut i=1/NI: match values (C-grid topology)
  !==========================================================================
  subroutine apply_bcs(uu, vv, ww, ni, nj, nk, u0, v0)
    integer,  intent(in)    :: ni, nj, nk
    real(8),  intent(in)    :: u0, v0
    real(8),  intent(inout) :: uu(ni,nj,nk), vv(ni,nj,nk), ww(ni,nj,nk)
    integer :: i, j, k

    ! Wall: no-slip
    !$acc parallel loop collapse(2) present(uu,vv,ww)
    do k = 1, nk
      do i = 1, ni
        uu(i,1,k) = 0d0;  vv(i,1,k) = 0d0;  ww(i,1,k) = 0d0
      end do
    end do
    !$acc end parallel loop

    ! Far field: freestream
    !$acc parallel loop collapse(2) present(uu,vv,ww)
    do k = 1, nk
      do i = 1, ni
        uu(i,nj,k) = u0;  vv(i,nj,k) = v0;  ww(i,nj,k) = 0d0
      end do
    end do
    !$acc end parallel loop

    ! Inlet: freestream
    !$acc parallel loop collapse(2) present(uu,vv,ww)
    do k = 1, nk
      do j = 1, nj
        uu(1,j,k) = u0;  vv(1,j,k) = v0;  ww(1,j,k) = 0d0
      end do
    end do
    !$acc end parallel loop

    ! Outlet: zero-gradient
    !$acc parallel loop collapse(2) present(uu,vv,ww)
    do k = 1, nk
      do j = 1, nj
        uu(ni,j,k) = uu(ni-1,j,k)
        vv(ni,j,k) = vv(ni-1,j,k)
        ww(ni,j,k) = ww(ni-1,j,k)
      end do
    end do
    !$acc end parallel loop

    ! Wake cut: match i=1 and i=NI (C-grid — same physical line)
    !$acc parallel loop collapse(2) present(uu,vv,ww)
    do k = 1, nk
      do j = 1, nj
        uu(1,j,k)  = 0.5d0 * (uu(1,j,k)  + uu(ni,j,k))
        vv(1,j,k)  = 0.5d0 * (vv(1,j,k)  + vv(ni,j,k))
        ww(1,j,k)  = 0.5d0 * (ww(1,j,k)  + ww(ni,j,k))
        uu(ni,j,k) = uu(1,j,k)
        vv(ni,j,k) = vv(1,j,k)
        ww(ni,j,k) = ww(1,j,k)
      end do
    end do
    !$acc end parallel loop

  end subroutine apply_bcs

  !==========================================================================
  ! RHIE-CHOW DIVERGENCE
  !   Computes div = (1/dt) * divergence(u*) using face velocities that
  !   include a pressure gradient correction to suppress the checkerboard.
  !
  !   Face velocity (Rhie-Chow, xi direction at i+1/2):
  !     U_face = 0.5*(U_i + U_{i+1})
  !             - 0.5*dt*[ (dp/dxi)_{i+1} - (dp/dxi)_i ] / dxi_face
  !             + 0.5*dt*[ (dp/dx)_{i+1} + (dp/dx)_i ] interpolated
  !   The correction damps the odd-even pressure decoupling.
  !==========================================================================
  subroutine rhie_chow_div(uu, vv, ww, pp, phi, dv, &
                           xi_x, xi_y, et_x, et_y, jc, g11, g22, g12, &
                           dxi, det, ni, nj, nk, dzi, dt)
    integer,  intent(in)  :: ni, nj, nk
    real(8),  intent(in)  :: dzi, dt
    real(8),  intent(in)  :: uu(ni,nj,nk), vv(ni,nj,nk), ww(ni,nj,nk)
    real(8),  intent(in)  :: pp(ni,nj,nk)
    real(8),  intent(out) :: phi(ni,nj,nk), dv(ni,nj,nk)
    real(8),  intent(in)  :: xi_x(ni,nj), xi_y(ni,nj)
    real(8),  intent(in)  :: et_x(ni,nj), et_y(ni,nj)
    real(8),  intent(in)  :: jc(ni,nj), g11(ni,nj), g22(ni,nj), g12(ni,nj)
    real(8),  intent(in)  :: dxi(ni,nj), det(ni,nj)
    integer  :: i, j, k
    real(8)  :: Uf_ip, Uf_im, Vf_jp, Vf_jm, Wf_kp, Wf_km
    real(8)  :: U_i, U_ip1, V_j, V_jp1
    real(8)  :: dpdxi_i, dpdxi_ip1, dpdet_j, dpdet_jp1
    real(8)  :: dUdxi, dVdet, dWdz

    phi = 0d0   ! initialise output

    !$acc parallel loop collapse(3) independent &
    !$acc&   private(Uf_ip,Uf_im,Vf_jp,Vf_jm,Wf_kp,Wf_km, &
    !$acc&           U_i,U_ip1,V_j,V_jp1, &
    !$acc&           dpdxi_i,dpdxi_ip1,dpdet_j,dpdet_jp1, &
    !$acc&           dUdxi,dVdet,dWdz)
    do k = 1, nk
      do j = 2, nj-1
        do i = 2, ni-1

          ! Contravariant velocity at cell centres
          U_i   = xi_x(i,j)  *uu(i,j,k)   + xi_y(i,j)  *vv(i,j,k)
          U_ip1 = xi_x(i+1,j)*uu(i+1,j,k) + xi_y(i+1,j)*vv(i+1,j,k)
          V_j   = et_x(i,j)  *uu(i,j,k)   + et_y(i,j)  *vv(i,j,k)
          V_jp1 = et_x(i,j+1)*uu(i,j+1,k) + et_y(i,j+1)*vv(i,j+1,k)

          ! Pressure gradient in xi at cells i and i+1 (for RC correction)
          dpdxi_i   = (pp(i+1,j,k) - pp(i-1,j,k)) / (2d0*dxi(i,j))
          dpdxi_ip1 = (pp(i+2,j,k) - pp(i,  j,k)) / (2d0*dxi(i+1,j))
          dpdet_j   = (pp(i,j+1,k) - pp(i,j-1,k)) / (2d0*det(i,j))
          dpdet_jp1 = (pp(i,j+2,k) - pp(i,j,  k)) / (2d0*det(i,j+1))

          ! Rhie-Chow face velocity at i+1/2 (xi direction)
          Uf_ip = 0.5d0*(U_i + U_ip1) &
                - 0.5d0*dt*(dpdxi_ip1 - dpdxi_i) / dxi(i,j)

          ! Rhie-Chow face velocity at i-1/2
          Uf_im = 0.5d0*(xi_x(i-1,j)*uu(i-1,j,k)+xi_y(i-1,j)*vv(i-1,j,k) + U_i) &
                - 0.5d0*dt*(dpdxi_i - (pp(i,j,k)-pp(i-2,j,k))/(2d0*dxi(i-1,j))) &
                / dxi(i,j)

          ! Rhie-Chow face velocity at j+1/2 (eta direction)
          Vf_jp = 0.5d0*(V_j + V_jp1) &
                - 0.5d0*dt*(dpdet_jp1 - dpdet_j) / det(i,j)

          ! Rhie-Chow face velocity at j-1/2
          Vf_jm = 0.5d0*(et_x(i,j-1)*uu(i,j-1,k)+et_y(i,j-1)*vv(i,j-1,k) + V_j) &
                - 0.5d0*dt*(dpdet_j - (pp(i,j,k)-pp(i,j-2,k))/(2d0*det(i,j-1))) &
                / det(i,j)

          ! Face velocities in z (no RC correction needed — Cartesian)
          Wf_kp = 0.5d0 * (ww(i,j,k) + ww(i,j,kp1(k)))
          Wf_km = 0.5d0 * (ww(i,j,k) + ww(i,j,km1(k)))

          ! Divergence: (1/J) * [d(J*U)/dxi + d(J*V)/deta] + dW/dz
          dUdxi = (jc(i+1,j)*Uf_ip - jc(i-1,j)*Uf_im) &
                  / (2d0 * dxi(i,j) * jc(i,j))
          dVdet = (jc(i,j+1)*Vf_jp - jc(i,j-1)*Vf_jm) &
                  / (2d0 * det(i,j) * jc(i,j))
          dWdz  = (Wf_kp - Wf_km) * dzi

          dv(i,j,k) = (dUdxi + dVdet + dWdz) / dt

        end do
      end do
    end do
    !$acc end parallel loop

  end subroutine rhie_chow_div

  !==========================================================================
  ! POISSON SOLVE  —  PLACEHOLDER
  !   Solves: laplacian(phi) = rhs
  !   Replace with cuFFT-based solver (DST-I in xi, DFT in z)
  !   For now: returns phi = 0 (no pressure correction)
  !   This lets you verify momentum + BCs before adding Poisson.
  !==========================================================================
  subroutine poisson_solve(phi, rhs, ni, nj, nk, dz)
    integer, intent(in)    :: ni, nj, nk
    real(8), intent(in)    :: dz
    real(8), intent(inout) :: phi(ni,nj,nk)
    real(8), intent(in)    :: rhs(ni,nj,nk)

    ! TODO: replace with cuFFT DST-I / DFT solver
    ! For debugging: zero phi means no pressure correction (pure advection test)
    phi = 0d0

    ! Uncomment to use a simple iterative solver for CPU debugging:
    ! call gauss_seidel(phi, rhs, ni, nj, nk, dz)

  end subroutine poisson_solve

end program gfoil36