import scipy
import numpy as np

mesh = np.array(
[
    [
        [True, True, True, ],
        [True, False, True, ],
        [True, True, True, ],
    ],
]
)

print(mesh.shape)
mesh = mesh[0,:,:]


mesh[scipy.ndimage.binary_fill_holes(mesh.astype(int))] = True

mesh = np.array([mesh])

print(mesh)