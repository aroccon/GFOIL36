"""
Convert a 2D C-grid (NI x NJ) in PLOT3D-style single-block ASCII format
into a 3D vtkStructuredGrid (.vts) by extruding uniformly in Z.

Input file layout:
    line 1: "NI NJ"
    then NI*NJ X values (j-major: all i for j=0, then j=1, ...)
    then NI*NJ Y values

Usage: python p3d_to_vts.py input.dat output.vts [nz] [zmin] [zmax]
"""
import sys
import numpy as np
import vtk
from vtk.util import numpy_support as nps


def read_plot3d_2d(path):
    with open(path) as f:
        toks = f.read().split()
    ni, nj = int(toks[0]), int(toks[1])
    vals = np.array(toks[2:], dtype=float)
    n = ni * nj
    assert vals.size >= 2 * n, f"expected {2*n} values, got {vals.size}"
    X = vals[0:n].reshape(nj, ni)       # row-major: j outer, i inner
    Y = vals[n:2*n].reshape(nj, ni)
    return ni, nj, X, Y


def build_structured_grid(X2d, Y2d, z_levels):
    """X2d,Y2d shape (nj, ni); returns vtkStructuredGrid with dims (ni, nj, nk)."""
    nj, ni = X2d.shape
    nk = len(z_levels)

    # VTK expects points in k-major, j-middle, i-fastest order
    pts = np.empty((nk * nj * ni, 3), dtype=np.float64)
    idx = 0
    for k in range(nk):
        z = z_levels[k]
        for j in range(nj):
            for i in range(ni):
                pts[idx, 0] = X2d[j, i]
                pts[idx, 1] = Y2d[j, i]
                pts[idx, 2] = z
                idx += 1

    vtk_points = vtk.vtkPoints()
    vtk_points.SetData(nps.numpy_to_vtk(pts, deep=True))

    sg = vtk.vtkStructuredGrid()
    sg.SetDimensions(ni, nj, nk)   # (i, j, k) order
    sg.SetPoints(vtk_points)
    return sg


def main():
    in_path  = sys.argv[1] if len(sys.argv) > 1 else "cmesh.p3d"
    out_path = sys.argv[2] if len(sys.argv) > 2 else "grid.vts"
    nz       = int(sys.argv[3])   if len(sys.argv) > 3 else 20
    zmin     = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
    zmax     = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0

    ni, nj, X, Y = read_plot3d_2d(in_path)
    print(f"read 2D grid: NI={ni}, NJ={nj}")
    print(f"  X range: [{X.min():.4f}, {X.max():.4f}]")
    print(f"  Y range: [{Y.min():.4f}, {Y.max():.4f}]")

    z_levels = np.linspace(zmin, zmax, nz)
    print(f"extruding to NK={nz} uniform levels in [{zmin}, {zmax}]")

    sg = build_structured_grid(X, Y, z_levels)

    w = vtk.vtkXMLStructuredGridWriter()
    w.SetFileName(out_path)
    w.SetInputData(sg)
    w.SetDataModeToBinary()
    w.Write()
    print(f"wrote {out_path}  (dims {ni} x {nj} x {nz})")


if __name__ == "__main__":
    main()