from mayavi import mlab
from scipy.ndimage import convolve
import numpy as np


def extract_surface(arr):
    """
    Extract the surface of a 3D array representing a voxelized mesh.

    This function identifies the surface voxels of a 3D array. A surface voxel is defined as a filled
    voxel (value True) that is adjacent to at least one empty voxel (value False). The function pads the
    input array to handle edge cases, identifies the surface voxels, and then removes the padding before
    returning the surface mask.

    Parameters:
    arr (numpy.ndarray): A 3D boolean array representing a voxelized mesh.

    Returns:
    numpy.ndarray: A 3D boolean array where surface voxels are marked as True.
    """
    # Pad the array by one voxel on all sides to handle edge cases during surface detection
    padded_array = np.pad(arr, pad_width=1, mode="constant", constant_values=0)

    # Initialize a mask for surface voxels
    surface_mask = np.zeros_like(padded_array, dtype=bool)

    # Iterate over each dimension (x, y, z) to identify surface voxels
    for dim in range(3):
        # Create slices for the current dimension
        slice_front = [slice(None)] * 3
        slice_back = [slice(None)] * 3
        slice_front[dim] = slice(1, None)
        slice_back[dim] = slice(None, -1)

        # Identify surface voxels by checking adjacent voxels in the current dimension
        surface_mask[tuple(slice_front)] |= (
            padded_array[tuple(slice_front)] & ~padded_array[tuple(slice_back)]
        )
        surface_mask[tuple(slice_back)] |= (
            padded_array[tuple(slice_back)] & ~padded_array[tuple(slice_front)]
        )

    # Remove padding and return the mask with only surface voxels marked
    return surface_mask[1:-1, 1:-1, 1:-1]


def vis(mesh):
    """
    Visualize the surface of a voxelized mesh.

    This function visualizes only the surface of the provided voxelized mesh. It identifies the
    surface voxels (those adjacent to at least one empty voxel) and plots them using Mayavi's mlab.

    Parameters:
    mesh (numpy.ndarray): A 3D boolean array representing the voxelized mesh.
    """
    # Adjust the scale factor for visualization; can be modified for denser meshes
    scale_factor = 1

    # Extract the coordinates of surface voxels in the mesh
    surface = extract_surface(mesh)
    x, y, z = np.where(surface)

    # Visualize the surface voxels using cube glyphs
    mlab.points3d(
        x,
        y,
        z,
        mode="cube",
        color=(0, 0, 1),  # Blue color for the cubes
        scale_mode="none",
        scale_factor=scale_factor,
    )

    # Configure the viewing angle for better perception
    mlab.view(azimuth=45, elevation=45)

    # Display the 3D visualization
    mlab.show()


def array_to_textfile(array, filename):
    """
    Write the contents of a 3D NumPy array to a text file.

    This function takes a 3D NumPy array and writes its contents to a specified file.
    Each line in the file contains the indices (x, y, z) and the value at that position
    in the array, with the value being converted from boolean to integer (1 for True, 0 for False).

    Parameters:
    array (numpy.ndarray): A 3D NumPy array whose contents are to be written to the file.
    filename (str): The name of the file where the array data will be written.

    Raises:
    ValueError: If the input is not a 3D NumPy array.
    """
    # Check that the input is a 3D NumPy array
    if not isinstance(array, np.ndarray) or len(array.shape) != 3:
        raise ValueError("The input must be a 3D NumPy array.")

    # Open the file for writing
    with open(filename, "w") as file:
        # Iterate over each element in the array
        for x in range(array.shape[0]):
            for y in range(array.shape[1]):
                for z in range(array.shape[2]):
                    # Convert the boolean value to an integer
                    value = int(array[x, y, z])
                    # Write the indices and value to the file
                    file.write(f"{x} {y} {z} {value}\n")


if __name__ == "__main__":
    pass
