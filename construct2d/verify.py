import numpy as np

with open('naca0012_sharp.p3d', 'rb') as f:
    def read_record():
        n = np.frombuffer(f.read(4), dtype=np.int32)[0]
        data = f.read(n)
        f.read(4)  # trailing marker
        return data

    nblocks = np.frombuffer(read_record(), dtype=np.int32)[0]
    dims    = np.frombuffer(read_record(), dtype=np.int32)
    ni, nj  = dims[0], dims[1]

    coords  = np.frombuffer(read_record(), dtype=np.float64)
    x = coords[:ni*nj].reshape(nj, ni)   # j outer, i inner
    y = coords[ni*nj:2*ni*nj].reshape(nj, ni)

import matplotlib.pyplot as plt
plt.plot(x[0,:], y[0,:], 'k-', lw=1.5)   # wall (j=0)
plt.axis('equal')
plt.show()