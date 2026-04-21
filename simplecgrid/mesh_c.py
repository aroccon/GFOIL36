"""
C-mesh generation around a NACA 0012 airfoil, with sharp trailing edge
and coincident wake cut (solver-ready).

Based on: https://alpynepyano.github.io/healthyNumerics/posts/cfd-03-grids-for-airfoils.html

The airfoil is generated analytically with the sharp-TE variant of the NACA
00xx equation (last coefficient -0.1036 so y(1)=0 exactly).

Exports to PLOT3D in construct2D's 2D formatted single-block layout.

Index convention (same as construct2D):
  j = 0        -> wall + wake cut (inner C boundary)
  j = nj-1     -> far-field outer boundary
  i = 0        -> wake outflow, upper side
  i = ni-1     -> wake outflow, lower side
  i ~ ni/2     -> leading-edge nose
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mclr

# ---------------------------------------------------------------
# STEP 0 — Helpers
# ---------------------------------------------------------------
def scale01(z):
    return (z - np.min(z)) / (np.max(z) - np.min(z))

def fStretch(my):
    yStretch = 0.5
    yOffset = 2.95
    iy = np.linspace(0, 1, my)
    return scale01(np.exp((yStretch * (yOffset + iy))**2))

def coord(a, b, xi):
    return a + xi * (b - a)


# ---------------------------------------------------------------
# STEP 1 — NACA 0012 airfoil (analytic, sharp TE)
# ---------------------------------------------------------------
# Symmetric NACA 00xx profile:
#   y = 5t * [ 0.2969 sqrt(x) - 0.1260 x - 0.3516 x^2
#              + 0.2843 x^3 - 0.1036 x^4 ]
# The last coefficient is -0.1036 (not the original -0.1015) so that y(1) = 0
# exactly, giving a genuinely sharp trailing edge.

def naca00xx(t, n_surface=20, le_clustering=0.0, te_clustering=0.0):
    """
    Generate NACA 00xx coordinates with sharp trailing edge.
    Returns points ordered TE (upper) -> LE -> TE (lower).

    t              : thickness ratio (e.g. 0.12 for NACA 0012)
    n_surface      : number of points per surface
    le_clustering  : 0 = uniform, 1 = full cosine (dense at LE)
    te_clustering  : 0 = uniform, 1 = full cosine (dense at TE)

    With le_clustering=1.0, te_clustering=0.3 you get tight LE spacing
    (needed to resolve the nose curvature) and modest TE density.
    """
    # Build x-distribution by blending cosine (dense at both ends)
    # with uniform (same spacing everywhere). Separate weights for each end.
    u = np.linspace(0.0, 1.0, n_surface)

    # Fully clustered at LE: x = 0.5*(1 - cos(pi*u))  -> dense near u=0 and u=1
    # We want to weight the "dense at u=0" part with le_clustering
    # and the "dense at u=1" part with te_clustering.
    # Trick: split the cosine into its LE-bias and TE-bias halves.
    x_uniform = u
    x_cosine_le = 1.0 - np.cos(0.5 * np.pi * u)       # dense near u=0 only
    x_cosine_te = np.sin(0.5 * np.pi * u)             # dense near u=1 only

    # Blend: normalize so x(0)=0 and x(1)=1 in all cases.
    # Each of the three base distributions already satisfies that.
    # Weight them: more le_clustering -> more x_cosine_le,
    #              more te_clustering -> more x_cosine_te,
    #              remainder goes to uniform.
    w_le = le_clustering
    w_te = te_clustering
    w_un = max(0.0, 1.0 - w_le - w_te)
    # If the user requests both > 0.5 we just renormalize
    total = w_le + w_te + w_un
    x = (w_le * x_cosine_le + w_te * x_cosine_te + w_un * x_uniform) / total

    y = 5.0 * t * (
        0.2969 * np.sqrt(x)
        - 0.1260 * x
        - 0.3516 * x**2
        + 0.2843 * x**3
        - 0.1036 * x**4
    )

    # Walk: TE -> upper -> LE  then  LE -> lower -> TE.
    xu = x[::-1];   yu =  y[::-1]
    xl = x[1:];     yl = -y[1:]

    px = np.concatenate([xu, xl])
    py = np.concatenate([yu, yl])
    return px, py


px, py = naca00xx(t=0.12, n_surface=30,
                  le_clustering=1.0,   # tight LE spacing
                  te_clustering=0.3)   # moderate TE density
nn = px.size
print(f"NACA 0012, sharp TE: {nn} surface points")
print(f"  First point (TE upper): ({px[0]:.6f}, {py[0]:.6f})")
print(f"  Middle point (LE)     : ({px[nn//2]:.6f}, {py[nn//2]:.6f})")
print(f"  Last  point (TE lower): ({px[-1]:.6f}, {py[-1]:.6f})")


# ---------------------------------------------------------------
# STEP 2 — Build the C-grid boundaries (with coincident wake cut)
# ---------------------------------------------------------------
ix1 = nn // 4
ix2 = nn - ix1

R = 1.25          # far-field radius
xback = 2.0       # downstream extent of wake / outflow position
nw = 20           # number of wake points on each side

# --- INNER BOUNDARY (j = 0): wake-top + airfoil + wake-bottom ------
# Upper wake runs from outflow back to the TE, at y = 0.
# Then we walk the airfoil (upper surface -> nose -> lower surface).
# Then the lower wake runs from TE back out to outflow, at y = 0.
# BOTH wake lines share y = 0 so the cut is coincident (1-to-1).

q1x = np.linspace(xback, 1.0, nw);  q1y = np.zeros(nw)     # upper wake
q2x = px[0:ix1];                    q2y = py[0:ix1]        # upper rear (incl. TE)
q3x = px[ix1:ix2];                  q3y = py[ix1:ix2]      # nose
q4x = px[ix2:nn];                   q4y = py[ix2:nn]       # lower rear (incl. closing TE)
q5x = np.linspace(1.0, xback, nw);  q5y = np.zeros(nw)     # lower wake

# Drop duplicate TE points at the wake/airfoil junctions:
#   q1 ends at (1,0), q2 starts at (1,0)  -> drop q2[0]
#   q4 ends at (1,0), q5 starts at (1,0)  -> drop q4[-1]
qtx = np.concatenate([q1x, q2x[1:], q3x, q4x[:-1], q5x])
qty = np.concatenate([q1y, q2y[1:], q3y, q4y[:-1], q5y])

# --- OUTER BOUNDARY (j = nj-1): far-field C-shape -----------------
# --- OUTER BOUNDARY (j = nj-1): far-field C-shape -----------------
# The outer boundary must be partitioned with EXACTLY the same segment
# sizes as the inner boundary, otherwise radial lines connecting inner
# to outer cross each other near the TE. Map:
#   inner q1 (upper wake)  -> outer wake-top strip:   x from R+1 down to 1, y=R
#   inner q2[1:] (upper)   -> outer upper arc from (1, R) around to nose top
#   inner q3 (nose)        -> outer left semicircle
#   inner q4[:-1] (lower)  -> outer lower arc from nose bottom to (1, -R)
#   inner q5 (lower wake)  -> outer wake-bottom strip: x from 1 up to R+1, y=-R
# Here we use a simple construction: far-field is a rectangle to the right
# of the semicircle (centered at the nose cx).

n_wake_top = q1x.size               # nw
n_upper    = q2x[1:].size           # ix1 - 1
n_nose     = q3x.size               # ix2 - ix1
n_lower    = q4x[:-1].size          # nn - ix2 - 1
n_wake_bot = q5x.size               # nw

cx = q3x[0]                         # x of the nose (airfoil leading edge)

# Q1a: outer above the upper wake — horizontal at y=R, from (xback, R) down to (1, R)
Q1ax = np.linspace(xback, 1.0, n_wake_top)
Q1ay = R * np.ones(n_wake_top)

# Q1b: outer above the upper airfoil surface — horizontal at y=R,
# from just left of x=1 down to the "corner" where the rectangle meets the nose semicircle
Q1bx = np.linspace(1.0, cx, n_upper + 1)[1:]    # skip the duplicate x=1 point
Q1by = R * np.ones(n_upper)

# Q2: nose semicircle, from top corner (cx, R) around through (cx-R, 0) to bottom corner (cx, -R)
phi = np.linspace(0, 1, n_nose)
phi = phi * np.pi + 0.5 * np.pi     # from pi/2 to 3pi/2
Q2x = R * np.cos(phi) + cx
Q2y = R * np.sin(phi)

# Q3a: outer below the lower airfoil surface — horizontal at y=-R, from corner out to (1, -R)
Q3ax = np.linspace(cx, 1.0, n_lower + 1)[:-1]   # skip duplicate x=1 point
Q3ay = -R * np.ones(n_lower)

# Q3b: outer below the lower wake — horizontal at y=-R, from (1, -R) out to (xback, -R)
Q3bx = np.linspace(1.0, xback, n_wake_bot)
Q3by = -R * np.ones(n_wake_bot)

Qtx = np.concatenate([Q1ax, Q1bx, Q2x, Q3ax, Q3bx])
Qty = np.concatenate([Q1ay, Q1by, Q2y, Q3ay, Q3by])


# ---------------------------------------------------------------
# STEP 3 — Fill the interior by interpolating between boundaries
# ---------------------------------------------------------------
mx = qtx.size         # ni — points around the C
my = 15               # nj — points radially outward
gsi = fStretch(my)    # wall-clustering

xk = np.zeros((mx, my))
yk = np.zeros((mx, my))
for t in range(mx):
    xk[t, :] = coord(qtx[t], Qtx[t], gsi)
    yk[t, :] = coord(qty[t], Qty[t], gsi)

ni, nj = xk.shape
print(f"Mesh dimensions: ni={ni}, nj={nj}  (cells: {(ni-1)*(nj-1)})")

# --- sanity check: wake cut at j=0 should be coincident on upper & lower ---
print(f"Wake cut check (j=0): "
      f"upper start ({xk[0,0]:.4f}, {yk[0,0]:.4f})  "
      f"lower end   ({xk[-1,0]:.4f}, {yk[-1,0]:.4f})")


# ---------------------------------------------------------------
# STEP 4 — Export to PLOT3D (construct2D-compatible format)
# ---------------------------------------------------------------
def write_plot3d_2d(filename, x, y):
    """
    2D formatted single-block PLOT3D (construct2D style):
      - Header:  "<ni><nj>" (two integers, width-12 fields)
      - Body:    all x values (one per line, %18.8E), then all y values
      - Fortran column-major ordering: i varies fastest
    """
    ni, nj = x.shape
    with open(filename, 'w') as f:
        f.write(f"{ni:>12d}{nj:>12d}\n")
        for j in range(nj):
            for i in range(ni):
                f.write(f"{x[i, j]:18.8E}\n")
        for j in range(nj):
            for i in range(ni):
                f.write(f"{y[i, j]:18.8E}\n")
    print(f"Wrote {filename}  ({ni} x {nj} = {ni*nj} points)")


write_plot3d_2d('cmesh.p3d', xk, yk)


# ---------------------------------------------------------------
# STEP 5 — Plot
# ---------------------------------------------------------------
fig = plt.figure(figsize=(18, 9))
myCmap = mclr.ListedColormap(['white', 'white'])
ax = fig.add_subplot(111)
ax.pcolormesh(xk, yk, np.zeros_like(xk),
              edgecolors='k', linewidths=0.5, cmap=myCmap, shading='auto')
#ax.fill(px, py, 'r', zorder=10)
ax.set_aspect('equal', 'datalim')
ax.set_title(f'C-mesh around NACA 0012, sharp TE   ({ni} x {nj})')
plt.tight_layout()
plt.savefig('cmesh.png', dpi=130, bbox_inches='tight')
print("Saved cmesh.png")
plt.show()