!=============================================================================
! GFOIL36 — Incompressible NS on a C-grid
! Collocated grid + Rhie-Chow, fractional step, RK3

program gfoil36
  implicit none
  integer, parameter :: nx = 512
  integer, parameter :: ny = 256
  integer, parameter :: nz = 128
  real(8), parameter :: re   = 50000d0
  real(8), parameter :: aoa  = 5d0
  real(8), parameter :: lz   = 0.2d0
  integer, parameter :: nsteps = 100
  integer, parameter :: nprint = 10
  real(8), parameter :: PI      = acos(-1d0)
  real(8), parameter :: nu      = 1d0 / re
  real(8), parameter :: dz_c    = lz / dble(nz)
  real(8), parameter :: dzi     = 1d0/dz_c
  real(8), parameter :: ddzi    = 1d0/dz_c**2
  real(8), parameter :: aoa_rad = aoa * PI / 180d0
  real(8), parameter :: u0      = cos(aoa_rad)
  real(8), parameter :: v0      = sin(aoa_rad)
  real(8), parameter :: rka(3) = [ 0d0,        -17d0/60d0, -5d0/12d0 ]
  real(8), parameter :: rkb(3) = [ 8d0/15d0,    5d0/12d0,   3d0/4d0  ]
  real(8), allocatable :: xg(:,:), yg(:,:)
  real(8), allocatable :: xi_x(:,:), xi_y(:,:)
  real(8), allocatable :: et_x(:,:), et_y(:,:)
  real(8), allocatable :: jac(:,:)
  real(8), allocatable :: g11(:,:), g22(:,:), g12(:,:)
  real(8), allocatable :: dxi(:,:), det(:,:)
  real(8), allocatable :: u(:,:,:), v(:,:,:), w(:,:,:)
  real(8), allocatable :: p(:,:,:)
  real(8), allocatable :: ru(:,:,:), rv(:,:,:), rw(:,:,:)
  real(8), allocatable :: us(:,:,:), vs(:,:,:), ws(:,:,:)
  real(8), allocatable :: phi(:,:,:), div(:,:,:)
  integer :: i, j, k, step, stage
  real(8) :: t, dt, dmin
  real(8) :: x_xi, x_et, y_xi, y_et, jac_loc
  real(8) :: rhs_u, rhs_v, rhs_w
  real(8) :: conv_u, conv_v, conv_w
  real(8) :: diff_xi_u, diff_xi_v, diff_xi_w
  real(8) :: diff_z_u,  diff_z_v,  diff_z_w
  real(8) :: cross_u, cross_v, cross_w
  real(8) :: U_con, V_con
  real(8) :: dphi_xi, dphi_et, dphi_dz, dphi_dx, dphi_dy
  real(8) :: noise

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

  do i = 1, nx
    do j = 1, ny
      if (i == 1) then
        x_xi = xg(2,j) - xg(1,j);     y_xi = yg(2,j) - yg(1,j)
      else if (i == nx) then
        x_xi = xg(nx,j) - xg(nx-1,j); y_xi = yg(nx,j) - yg(nx-1,j)
      else
        x_xi = (xg(i+1,j) - xg(i-1,j)) * 0.5d0
        y_xi = (yg(i+1,j) - yg(i-1,j)) * 0.5d0
      end if
      if (j == 1) then
        x_et = xg(i,2) - xg(i,1);     y_et = yg(i,2) - yg(i,1)
      else if (j == ny) then
        x_et = xg(i,ny) - xg(i,ny-1); y_et = yg(i,ny) - yg(i,ny-1)
      else
        x_et = (xg(i,j+1) - xg(i,j-1)) * 0.5d0
        y_et = (yg(i,j+1) - yg(i,j-1)) * 0.5d0
      end if
      jac_loc   = x_xi*y_et - x_et*y_xi
      jac(i,j)  = jac_loc
      xi_x(i,j) =  y_et/jac_loc;  xi_y(i,j) = -x_et/jac_loc
      et_x(i,j) = -y_xi/jac_loc;  et_y(i,j) =  x_xi/jac_loc
      g11(i,j)  = xi_x(i,j)**2 + xi_y(i,j)**2
      g22(i,j)  = et_x(i,j)**2 + et_y(i,j)**2
      g12(i,j)  = xi_x(i,j)*et_x(i,j) + xi_y(i,j)*et_y(i,j)
      dxi(i,j)  = sqrt(x_xi**2 + y_xi**2)
      det(i,j)  = sqrt(x_et**2 + y_et**2)
    end do
  end do

  dmin = minval(dxi(:,:))
  dt   = 0.4d0 * dmin / 2d0

  write(*,'(A)') repeat('=',52)
  write(*,'(A,I0,A,I0,A,I0)') '  Grid  : ', nx,' x ',ny,' x ',nz
  write(*,'(A,F8.1)')          '  Re    : ', re
  write(*,'(A,ES10.3)')        '  dt    : ', dt
  write(*,'(A,ES10.3)')        '  dy1   : ', minval(det(:,1))
  write(*,'(A,ES10.3)')        '  dt_diff_expl: ', re*minval(det(:,1))**2/2d0
  write(*,'(A)') repeat('=',52)

  call random_seed()
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        call random_number(noise)
        u(i,j,k) = u0 + (noise-0.5d0)*1d-3
        call random_number(noise)
        v(i,j,k) = v0 + (noise-0.5d0)*1d-3
        call random_number(noise)
        w(i,j,k) =      (noise-0.5d0)*1d-3
        p(i,j,k)  = 0d0
        ru(i,j,k) = 0d0; rv(i,j,k) = 0d0; rw(i,j,k) = 0d0
      end do
    end do
  end do

  t = 0d0
  do step = 1, nsteps
    do stage = 1, 3
      do k = 1, nz
        do j = 2, ny-1    ! full domain interior
          do i = 2, nx-1  ! diffusion always valid (needs ±1 only)

            U_con = xi_x(i,j)*u(i,j,k) + xi_y(i,j)*v(i,j,k)
            V_con = et_x(i,j)*u(i,j,k) + et_y(i,j)*v(i,j,k)

            !--- Convection of u ---
            ! xi: 3rd-order upwind in interior, 2nd-order central at i=2 and i=nx-1
            if (i >= 3 .and. i <= nx-2) then
              if (U_con >= 0d0) then
                conv_u=U_con*(2d0*u(i+1,j,k)-6d0*u(i,j,k)+3d0*u(i-1,j,k)+u(i-2,j,k))/(6d0*dxi(i,j))
              else
                conv_u=U_con*(-u(i+2,j,k)+3d0*u(i+1,j,k)+6d0*u(i,j,k)-2d0*u(i-1,j,k))/(6d0*dxi(i,j))*(-1d0)
              end if
            else
              conv_u = U_con*(u(i+1,j,k)-u(i-1,j,k))/(2d0*dxi(i,j))
            end if
            ! eta: 3rd-order upwind in interior, 2nd-order central at j=2 and j=ny-1
            if (j >= 3 .and. j <= ny-2) then
              if (V_con >= 0d0) then
                conv_u=conv_u+V_con*(2d0*u(i,j+1,k)-6d0*u(i,j,k)+3d0*u(i,j-1,k)+u(i,j-2,k))/(6d0*det(i,j))
              else
                conv_u=conv_u+V_con*(-u(i,j+2,k)+3d0*u(i,j+1,k)+6d0*u(i,j,k)-2d0*u(i,j-1,k))/(6d0*det(i,j))*(-1d0)
              end if
            else
              conv_u = conv_u + V_con*(u(i,j+1,k)-u(i,j-1,k))/(2d0*det(i,j))
            end if
            ! z: always 2nd-order central (periodic, no boundary issue)
            conv_u = conv_u + w(i,j,k)*(u(i,j,kp1(k))-u(i,j,km1(k)))*0.5d0*dzi

            !--- Convection of v ---
            if (i >= 3 .and. i <= nx-2) then
              if (U_con >= 0d0) then
                conv_v=U_con*(2d0*v(i+1,j,k)-6d0*v(i,j,k)+3d0*v(i-1,j,k)+v(i-2,j,k))/(6d0*dxi(i,j))
              else
                conv_v=U_con*(-v(i+2,j,k)+3d0*v(i+1,j,k)+6d0*v(i,j,k)-2d0*v(i-1,j,k))/(6d0*dxi(i,j))*(-1d0)
              end if
            else
              conv_v = U_con*(v(i+1,j,k)-v(i-1,j,k))/(2d0*dxi(i,j))
            end if
            if (j >= 3 .and. j <= ny-2) then
              if (V_con >= 0d0) then
                conv_v=conv_v+V_con*(2d0*v(i,j+1,k)-6d0*v(i,j,k)+3d0*v(i,j-1,k)+v(i,j-2,k))/(6d0*det(i,j))
              else
                conv_v=conv_v+V_con*(-v(i,j+2,k)+3d0*v(i,j+1,k)+6d0*v(i,j,k)-2d0*v(i,j-1,k))/(6d0*det(i,j))*(-1d0)
              end if
            else
              conv_v = conv_v + V_con*(v(i,j+1,k)-v(i,j-1,k))/(2d0*det(i,j))
            end if
            conv_v = conv_v + w(i,j,k)*(v(i,j,kp1(k))-v(i,j,km1(k)))*0.5d0*dzi

            !--- Convection of w ---
            if (i >= 3 .and. i <= nx-2) then
              if (U_con >= 0d0) then
                conv_w=U_con*(2d0*w(i+1,j,k)-6d0*w(i,j,k)+3d0*w(i-1,j,k)+w(i-2,j,k))/(6d0*dxi(i,j))
              else
                conv_w=U_con*(-w(i+2,j,k)+3d0*w(i+1,j,k)+6d0*w(i,j,k)-2d0*w(i-1,j,k))/(6d0*dxi(i,j))*(-1d0)
              end if
            else
              conv_w = U_con*(w(i+1,j,k)-w(i-1,j,k))/(2d0*dxi(i,j))
            end if
            if (j >= 3 .and. j <= ny-2) then
              if (V_con >= 0d0) then
                conv_w=conv_w+V_con*(2d0*w(i,j+1,k)-6d0*w(i,j,k)+3d0*w(i,j-1,k)+w(i,j-2,k))/(6d0*det(i,j))
              else
                conv_w=conv_w+V_con*(-w(i,j+2,k)+3d0*w(i,j+1,k)+6d0*w(i,j,k)-2d0*w(i,j-1,k))/(6d0*det(i,j))*(-1d0)
              end if
            else
              conv_w = conv_w + V_con*(w(i,j+1,k)-w(i,j-1,k))/(2d0*det(i,j))
            end if
            conv_w = conv_w + w(i,j,k)*(w(i,j,kp1(k))-w(i,j,km1(k)))*0.5d0*dzi
            !--- Explicit diffusion: xi + z + g12 cross (eta handled implicitly)
            diff_xi_u = g11(i,j)*(u(i+1,j,k)-2d0*u(i,j,k)+u(i-1,j,k))/dxi(i,j)**2
            diff_xi_v = g11(i,j)*(v(i+1,j,k)-2d0*v(i,j,k)+v(i-1,j,k))/dxi(i,j)**2
            diff_xi_w = g11(i,j)*(w(i+1,j,k)-2d0*w(i,j,k)+w(i-1,j,k))/dxi(i,j)**2
            diff_z_u  = (u(i,j,kp1(k))-2d0*u(i,j,k)+u(i,j,km1(k)))*ddzi
            diff_z_v  = (v(i,j,kp1(k))-2d0*v(i,j,k)+v(i,j,km1(k)))*ddzi
            diff_z_w  = (w(i,j,kp1(k))-2d0*w(i,j,k)+w(i,j,km1(k)))*ddzi
            ! g12 cross term: needs ±1 in both i and j — valid at i=2..nx-1, j=2..ny-1
            cross_u=2d0*g12(i,j)*(u(i+1,j+1,k)-u(i+1,j-1,k)-u(i-1,j+1,k)+u(i-1,j-1,k))/(4d0*dxi(i,j)*det(i,j))
            cross_v=2d0*g12(i,j)*(v(i+1,j+1,k)-v(i+1,j-1,k)-v(i-1,j+1,k)+v(i-1,j-1,k))/(4d0*dxi(i,j)*det(i,j))
            cross_w=2d0*g12(i,j)*(w(i+1,j+1,k)-w(i+1,j-1,k)-w(i-1,j+1,k)+w(i-1,j-1,k))/(4d0*dxi(i,j)*det(i,j))

            rhs_u = -conv_u + nu*(diff_xi_u + diff_z_u + cross_u)
            rhs_v = -conv_v + nu*(diff_xi_v + diff_z_v + cross_v)
            rhs_w = -conv_w + nu*(diff_xi_w + diff_z_w + cross_w)

            ru(i,j,k) = rka(stage)*ru(i,j,k) + rhs_u
            rv(i,j,k) = rka(stage)*rv(i,j,k) + rhs_v
            rw(i,j,k) = rka(stage)*rw(i,j,k) + rhs_w
            us(i,j,k) = u(i,j,k) + dt*rkb(stage)*ru(i,j,k)
            vs(i,j,k) = v(i,j,k) + dt*rkb(stage)*rv(i,j,k)
            ws(i,j,k) = w(i,j,k) + dt*rkb(stage)*rw(i,j,k)
          end do
        end do
      end do
      !$acc end parallel loop
      call implicit_eta(us, u, g22, det, nx, ny, nz, nu, dt, rkb(stage), u0)
      call implicit_eta(vs, v, g22, det, nx, ny, nz, nu, dt, rkb(stage), v0)
      call implicit_eta(ws, w, g22, det, nx, ny, nz, nu, dt, rkb(stage), 0d0)
      call apply_bcs(us, vs, ws, nx, ny, nz, u0, v0)
    end do

    call rhie_chow_div(us, vs, ws, p, phi, div, &
                       xi_x, xi_y, et_x, et_y, jac, g11, g22, g12, &
                       dxi, det, nx, ny, nz, dzi, dt)
    call poisson_solve(phi, div, nx, ny, nz, dz_c)

    do k = 1, nz
      do j = 2, ny-1
        do i = 2, nx-1
          dphi_xi = (phi(i+1,j,k)-phi(i-1,j,k))/(2d0*dxi(i,j))
          dphi_et = (phi(i,j+1,k)-phi(i,j-1,k))/(2d0*det(i,j))
          dphi_dz = (phi(i,j,kp1(k))-phi(i,j,km1(k)))*0.5d0*dzi
          dphi_dx = xi_x(i,j)*dphi_xi + et_x(i,j)*dphi_et
          dphi_dy = xi_y(i,j)*dphi_xi + et_y(i,j)*dphi_et
          u(i,j,k) = us(i,j,k) - dt*dphi_dx
          v(i,j,k) = vs(i,j,k) - dt*dphi_dy
          w(i,j,k) = ws(i,j,k) - dt*dphi_dz
          p(i,j,k) = p(i,j,k)  + phi(i,j,k)
        end do
      end do
    end do
    call apply_bcs(u, v, w, nx, ny, nz, u0, v0)
    t = t + dt
    if (mod(step,nprint)==0) write(*,'(A,I6,A,F9.4,A,ES10.3)') &
      '  step=',step,'  t=',t,'  |u|max=',maxval(abs(u))
  end do
  write(*,'(A)') '  Done.'
  deallocate(xg,yg,xi_x,xi_y,et_x,et_y,jac,g11,g22,g12,dxi,det)
  deallocate(u,v,w,p,ru,rv,rw,us,vs,ws,phi,div)

contains

  pure integer function kp1(k)
    integer, intent(in) :: k
    kp1 = mod(k, nz) + 1
  end function

  pure integer function km1(k)
    integer, intent(in) :: k
    km1 = mod(k + nz - 2, nz) + 1
  end function

  subroutine implicit_eta(us_out, u_old, g22, det, ni, nj, nk, nu, dt, rkb, far_val)
    integer, intent(in)    :: ni, nj, nk
    real(8), intent(in)    :: nu, dt, rkb, far_val
    real(8), intent(inout) :: us_out(ni,nj,nk)
    real(8), intent(in)    :: u_old(ni,nj,nk)
    real(8), intent(in)    :: g22(ni,nj), det(ni,nj)
    integer :: i, j, k
    real(8) :: alpha, fac
    real(8) :: lo(nj), diag(nj), up(nj), rr(nj), sol(nj)

    do k = 1, nk
      do i = 1, ni
        do j = 2, nj-1
          alpha   = nu * dt * rkb * g22(i,j) / (2d0 * det(i,j)**2)
          lo(j)   = -alpha
          diag(j) =  1d0 + 2d0*alpha
          up(j)   = -alpha
          rr(j)   = us_out(i,j,k) &
                  + alpha*(u_old(i,j+1,k) - 2d0*u_old(i,j,k) + u_old(i,j-1,k))
        end do
        ! Wall BC j=1: no-slip
        diag(1) = 1d0;  up(1) = 0d0;  lo(1) = 0d0;  rr(1) = 0d0
        ! Far-field BC j=nj: freestream
        diag(nj) = 1d0;  lo(nj) = 0d0;  up(nj) = 0d0;  rr(nj) = far_val
        ! Thomas forward sweep
        do j = 2, nj
          fac     = lo(j) / diag(j-1)
          diag(j) = diag(j) - fac * up(j-1)
          rr(j)   = rr(j)   - fac * rr(j-1)
        end do
        ! Back substitution
        sol(nj) = rr(nj) / diag(nj)
        do j = nj-1, 1, -1
          sol(j) = (rr(j) - up(j)*sol(j+1)) / diag(j)
        end do
        do j = 1, nj
          us_out(i,j,k) = sol(j)
        end do
      end do
    end do
  end subroutine implicit_eta

  subroutine apply_bcs(uu, vv, ww, ni, nj, nk, u0, v0)
    integer, intent(in)    :: ni, nj, nk
    real(8), intent(in)    :: u0, v0
    real(8), intent(inout) :: uu(ni,nj,nk), vv(ni,nj,nk), ww(ni,nj,nk)
    integer :: i, j, k
    do k = 1, nk
      do i = 1, ni
        uu(i,1,k)  = 0d0;  vv(i,1,k)  = 0d0;  ww(i,1,k)  = 0d0
        uu(i,nj,k) = u0;   vv(i,nj,k) = v0;   ww(i,nj,k) = 0d0
      end do
    end do
    do k = 1, nk
      do j = 1, nj
        uu(1,j,k) = u0;  vv(1,j,k) = v0;  ww(1,j,k) = 0d0
        uu(ni,j,k) = uu(ni-1,j,k)
        vv(ni,j,k) = vv(ni-1,j,k)
        ww(ni,j,k) = ww(ni-1,j,k)
      end do
    end do
    do k = 1, nk
      do j = 1, nj
        uu(1,j,k)  = 0.5d0*(uu(1,j,k) + uu(ni,j,k))
        vv(1,j,k)  = 0.5d0*(vv(1,j,k) + vv(ni,j,k))
        ww(1,j,k)  = 0.5d0*(ww(1,j,k) + ww(ni,j,k))
        uu(ni,j,k) = uu(1,j,k)
        vv(ni,j,k) = vv(1,j,k)
        ww(ni,j,k) = ww(1,j,k)
      end do
    end do
  end subroutine apply_bcs

  subroutine rhie_chow_div(uu, vv, ww, pp, phi, dv, &
                           xi_x, xi_y, et_x, et_y, jc, g11, g22, g12, &
                           dxi, det, ni, nj, nk, dzi, dt)
    integer, intent(in)  :: ni, nj, nk
    real(8), intent(in)  :: dzi, dt
    real(8), intent(in)  :: uu(ni,nj,nk), vv(ni,nj,nk), ww(ni,nj,nk)
    real(8), intent(in)  :: pp(ni,nj,nk)
    real(8), intent(out) :: phi(ni,nj,nk), dv(ni,nj,nk)
    real(8), intent(in)  :: xi_x(ni,nj), xi_y(ni,nj)
    real(8), intent(in)  :: et_x(ni,nj), et_y(ni,nj)
    real(8), intent(in)  :: jc(ni,nj), g11(ni,nj), g22(ni,nj), g12(ni,nj)
    real(8), intent(in)  :: dxi(ni,nj), det(ni,nj)
    integer  :: i, j, k
    real(8)  :: Uf_ip, Uf_im, Vf_jp, Vf_jm, Wf_kp, Wf_km
    real(8)  :: U_i, U_ip1, V_j, V_jp1
    real(8)  :: dpdxi_i, dpdxi_ip1, dpdet_j, dpdet_jp1
    real(8)  :: dUdxi, dVdet, dWdz
    phi = 0d0
    do k = 1, nk
      do j = 3, nj-2
        do i = 3, ni-2
          U_i   = xi_x(i,j)  *uu(i,j,k)   + xi_y(i,j)  *vv(i,j,k)
          U_ip1 = xi_x(i+1,j)*uu(i+1,j,k) + xi_y(i+1,j)*vv(i+1,j,k)
          V_j   = et_x(i,j)  *uu(i,j,k)   + et_y(i,j)  *vv(i,j,k)
          V_jp1 = et_x(i,j+1)*uu(i,j+1,k) + et_y(i,j+1)*vv(i,j+1,k)
          dpdxi_i   = (pp(i+1,j,k)-pp(i-1,j,k))/(2d0*dxi(i,j))
          dpdxi_ip1 = (pp(i+2,j,k)-pp(i,  j,k))/(2d0*dxi(i+1,j))
          dpdet_j   = (pp(i,j+1,k)-pp(i,j-1,k))/(2d0*det(i,j))
          dpdet_jp1 = (pp(i,j+2,k)-pp(i,j,  k))/(2d0*det(i,j+1))
          Uf_ip = 0.5d0*(U_i+U_ip1) - 0.5d0*dt*(dpdxi_ip1-dpdxi_i)/dxi(i,j)
          Uf_im = 0.5d0*(xi_x(i-1,j)*uu(i-1,j,k)+xi_y(i-1,j)*vv(i-1,j,k)+U_i) &
                - 0.5d0*dt*(dpdxi_i-(pp(i,j,k)-pp(i-2,j,k))/(2d0*dxi(i-1,j)))/dxi(i,j)
          Vf_jp = 0.5d0*(V_j+V_jp1) - 0.5d0*dt*(dpdet_jp1-dpdet_j)/det(i,j)
          Vf_jm = 0.5d0*(et_x(i,j-1)*uu(i,j-1,k)+et_y(i,j-1)*vv(i,j-1,k)+V_j) &
                - 0.5d0*dt*(dpdet_j-(pp(i,j,k)-pp(i,j-2,k))/(2d0*det(i,j-1)))/det(i,j)
          Wf_kp = 0.5d0*(ww(i,j,k)+ww(i,j,kp1(k)))
          Wf_km = 0.5d0*(ww(i,j,k)+ww(i,j,km1(k)))
          dUdxi = (jc(i+1,j)*Uf_ip-jc(i-1,j)*Uf_im)/(2d0*dxi(i,j)*jc(i,j))
          dVdet = (jc(i,j+1)*Vf_jp-jc(i,j-1)*Vf_jm)/(2d0*det(i,j)*jc(i,j))
          dWdz  = (Wf_kp-Wf_km)*dzi
          dv(i,j,k) = (dUdxi+dVdet+dWdz)/dt
        end do
      end do
    end do
  end subroutine rhie_chow_div

  subroutine poisson_solve(phi, rhs, ni, nj, nk, dz)
    integer, intent(in)    :: ni, nj, nk
    real(8), intent(in)    :: dz
    real(8), intent(inout) :: phi(ni,nj,nk)
    real(8), intent(in)    :: rhs(ni,nj,nk)
    phi = 0d0
  end subroutine poisson_solve

end program gfoil36
