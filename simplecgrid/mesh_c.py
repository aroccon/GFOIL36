import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mclr


# ---------------------------------------------------------------------------
# 1. NACA 4-digit airfoil generator
# ---------------------------------------------------------------------------
def naca4(code="0012", n_points=80):
    """
    Generate NACA 4-digit airfoil coordinates.

    Points go from trailing edge -> upper surface -> leading edge ->
    lower surface -> trailing edge (standard CCW ordering for a C-grid).

    Parameters
    ----------
    code : str     4-digit NACA designation, e.g. "0012".
    n_points : int number of points per surface (upper / lower).

    Returns
    -------
    x, y : 1D arrays of airfoil coordinates.
    """
    m = int(code[0]) / 100.0       # max camber
    p = int(code[1]) / 10.0        # position of max camber
    t = int(code[2:]) / 100.0      # thickness

    # half-cosine spacing -> clustering at LE only, coarser at TE
    beta = np.linspace(0.0, np.pi, n_points)
    x = 1.0 - np.cos(beta / 2.0)

    # thickness distribution (sharp TE coefficient = -0.1036)
    yt = 5.0 * t * (0.2969 * np.sqrt(x)
                    - 0.1260 * x
                    - 0.3516 * x**2
                    + 0.2843 * x**3
                    - 0.1036 * x**4)

    # mean camber line
    yc = np.where(x < p,
                  m / p**2 * (2.0 * p * x - x**2) if p > 0 else np.zeros_like(x),
                  m / (1.0 - p)**2 * ((1.0 - 2.0 * p) + 2.0 * p * x - x**2)
                  if p > 0 else np.zeros_like(x))
    dyc = np.where(x < p,
                   2.0 * m / p**2 * (p - x) if p > 0 else np.zeros_like(x),
                   2.0 * m / (1.0 - p)**2 * (p - x) if p > 0 else np.zeros_like(x))
    theta = np.arctan(dyc)

    xu = x - yt * np.sin(theta);  yu = yc + yt * np.cos(theta)
    xl = x + yt * np.sin(theta);  yl = yc - yt * np.cos(theta)

    # TE -> upper -> LE -> lower -> TE
    X = np.concatenate([xu[::-1], xl[1:]])
    Y = np.concatenate([yu[::-1], yl[1:]])
    return X, Y


# ---------------------------------------------------------------------------
# 2. Normal-direction stretching
# ---------------------------------------------------------------------------
def stretch(n, y_stretch=0.5, y_offset=2.95):
    """Exponential clustering near the wall (eta=0)."""
    iy = np.linspace(0.0, 1.0, n)
    return (np.exp((y_stretch * (y_offset + iy))**2) - np.exp((y_stretch * y_offset)**2)) \
           / (np.exp((y_stretch * (y_offset + 1))**2) - np.exp((y_stretch * y_offset)**2))


# ---------------------------------------------------------------------------
# 3. C-grid construction
# ---------------------------------------------------------------------------
def _resample_curve(x, y, n):
    """Resample a polyline at n points uniformly spaced in arclength."""
    seg = np.hypot(np.diff(x), np.diff(y))
    s = np.concatenate([[0.0], np.cumsum(seg)])
    s_new = np.linspace(0.0, s[-1], n)
    return np.interp(s_new, s, x), np.interp(s_new, s, y)


def build_cgrid(px, py, R=1.25, x_wake=2.0, n_wake=20, n_normal=15):
    """
    Build a C-grid around the airfoil.

    Inner boundary traversal:
        downstream (top of wake cut) -> TE -> upper surface -> LE ->
        lower surface -> TE -> downstream (bottom of wake cut)
    Outer boundary mirrors this: top straight line going upstream ->
    front semicircle -> bottom straight line going downstream.
    """
    le_idx = int(np.argmin(px))
    xe = px[0]                                       # trailing edge x

    # --- inner boundary -----------------------------------------------------
    up_wake_x = np.linspace(x_wake, xe, n_wake)
    up_wake_y = np.zeros_like(up_wake_x)

    # airfoil (skip first point, it's at the TE already in up_wake)
    af_x = px[1:]
    af_y = py[1:]

    # lower wake (skip first point, already at the TE)
    lo_wake_x = np.linspace(xe, x_wake, n_wake)[1:]
    lo_wake_y = np.zeros_like(lo_wake_x)

    inner_x = np.concatenate([up_wake_x, af_x, lo_wake_x])
    inner_y = np.concatenate([up_wake_y, af_y, lo_wake_y])

    # --- outer boundary (built densely, then resampled to match) -----------
    cx = px[le_idx]

    # top straight: (x_wake, +R) -> (cx, +R)
    top_x = np.linspace(x_wake, cx, 60)
    top_y = np.full_like(top_x, R)

    # front semicircle: (cx, +R) -> (cx-R, 0) -> (cx, -R)
    phi = np.linspace(0.5 * np.pi, 1.5 * np.pi, 80)
    arc_x = R * np.cos(phi) + cx
    arc_y = R * np.sin(phi)

    # bottom straight: (cx, -R) -> (x_wake, -R)
    bot_x = np.linspace(cx, x_wake, 60)
    bot_y = np.full_like(bot_x, -R)

    outer_x = np.concatenate([top_x, arc_x[1:], bot_x[1:]])
    outer_y = np.concatenate([top_y, arc_y[1:], bot_y[1:]])

    # resample outer to match inner point count
    outer_x, outer_y = _resample_curve(outer_x, outer_y, len(inner_x))

    assert len(inner_x) == len(outer_x)

    # --- transfinite interpolation in the normal direction ------------------
    mx = len(inner_x)
    eta = stretch(n_normal)

    xk = np.zeros((mx, n_normal))
    yk = np.zeros((mx, n_normal))
    for i in range(mx):
        xk[i, :] = inner_x[i] + eta * (outer_x[i] - inner_x[i])
        yk[i, :] = inner_y[i] + eta * (outer_y[i] - inner_y[i])

    return xk, yk, (inner_x, inner_y), (outer_x, outer_y)


# ---------------------------------------------------------------------------
# 4. Plot3D (P3D) writer — 2D, single-block, formatted ASCII
# ---------------------------------------------------------------------------
def write_plot3d(filename, xk, yk):
    """
    Write a 2D single-block mesh in formatted Plot3D ASCII.

    File layout:
        nx  ny
        x[0,0] x[1,0] ... x[nx-1,0]  x[0,1] ...  (nx*ny values, i fastest)
        y[0,0] y[1,0] ... y[nx-1,0]  y[0,1] ...  (nx*ny values, i fastest)

    One value per line, %18.8E, matching the convention of the reference file.

    NOTE on orientation: the i-direction is reversed relative to how
    build_cgrid() constructs the curve. build_cgrid() walks the C clockwise
    (upper-wake-outlet -> TE -> LE -> TE -> lower-wake-outlet), which gives a
    *negative* Jacobian x_xi*y_et - x_et*y_xi when eta points outward from the
    wall. Reversing i here produces the standard right-handed orientation
    (positive Jacobian) that solvers expect.
    """
    xk = xk[::-1, :]
    yk = yk[::-1, :]
    nx, ny = xk.shape
    with open(filename, "w") as f:
        f.write(f"{nx:12d}{ny:12d}\n")
        for v in xk.flatten(order="F"):
            f.write(f"{v:18.8E}\n")
        for v in yk.flatten(order="F"):
            f.write(f"{v:18.8E}\n")


# ---------------------------------------------------------------------------
# 5. Run + plot
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    # airfoil
    px, py = naca4("0012", n_points=30)

    # C-grid
    xk, yk, inner, outer = build_cgrid(
        px, py,
        R=1.25,        # far-field radius
        x_wake=2.0,    # downstream extent
        n_wake=20,     # points along each wake cut
        n_normal=15,   # points in wall-normal direction
    )

    # plot — draw actual grid lines, not pcolormesh cells
    fig, ax = plt.subplots(figsize=(14, 8))

    # j-lines: one curve per normal-direction index (wrap around the C)
    for j in range(xk.shape[1]):
        ax.plot(xk[:, j], yk[:, j], "k-", lw=0.4)
    # i-lines: one curve per streamwise index (radial spokes)
    for i in range(xk.shape[0]):
        ax.plot(xk[i, :], yk[i, :], "k-", lw=0.4)

    # highlight the airfoil itself
    ax.fill(px, py, color="crimson", zorder=10, lw=0)
    ax.plot(px, py, "k-", lw=1.2, zorder=11)

    ax.set_aspect("equal")
    ax.set_title("C-grid around NACA 0012", fontsize=14)
    ax.set_xlabel("x / c")
    ax.set_ylabel("y / c")
    ax.grid(False)
    plt.tight_layout()

    out = "cgrid_naca0012.png"
    plt.savefig(out, dpi=140, bbox_inches="tight")
    print(f"Saved: {out}")
    print(f"Mesh shape: {xk.shape}")

    # Plot3D export
    p3d_out = "cgrid_naca0012.p3d"
    write_plot3d(p3d_out, xk, yk)
    print(f"Saved: {p3d_out}")
    print()
    print("=" * 52)
    print(f"  For main.f90, set at the top of the program:")
    print(f"    integer, parameter :: nx = {xk.shape[0]}")
    print(f"    integer, parameter :: ny = {xk.shape[1]}")
    print("=" * 52)

    plt.show()