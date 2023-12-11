from mayavi import mlab
import numpy as np
import os
import fnmatch


def delete_stl_files(root_dir, exclude_files=None):
    """
    Delete .stl files from a specified directory and its subdirectories,
    with an option to exclude certain files from deletion.

    This function recursively traverses through the given directory and
    all of its subdirectories to find and delete files with the .stl extension.
    It provides an option to exclude specific .stl files from deletion based on their filenames.

    Parameters:
    root_dir (str): The root directory path where the search and deletion of .stl files begins.
    exclude_files (list of str, optional): A list of .stl filenames that should be excluded from deletion.
                                           Defaults to None, meaning no files are excluded.

    Returns:
    None: This function does not return anything. It directly deletes files from the filesystem.
    """

    # Check if exclude_files is None and, if so, initialize it as an empty list
    if exclude_files is None:
        exclude_files = []

    # Walk through all directories and subdirectories starting from root_dir
    for dirpath, dirnames, filenames in os.walk(root_dir):
        # Filter for .stl files in the current directory
        for filename in fnmatch.filter(filenames, "*.stl"):
            # Proceed only if the file is not in the exclude list
            if filename not in exclude_files:
                file_path = os.path.join(dirpath, filename)
                print("Deleted: ", filename)
                os.remove(file_path)  # Delete the file


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


def vis(mesh, dim=3):
    """
    Visualize the surface of a voxelized mesh.

    This function visualizes only the surface of the provided voxelized mesh. It identifies the
    surface voxels (those adjacent to at least one empty voxel) and plots them using Mayavi's mlab.

    Parameters:
    mesh (numpy.ndarray): A 3D boolean array representing the voxelized mesh.
    dim (int, optional): Dimension of the input mesh.
                         Default 3 = input mesh is a 3D mesh.
    """
    # Adjust the scale factor for visualization; can be modified for denser meshes
    scale_factor = 1

    # Extract the coordinates of surface voxels in the mesh
    if dim == 3:
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

    else:
        # Convert bool values to 1 and 0 and replace 0 values with NaN so that they are not visualized
        mesh = mesh.astype(int)
        mesh = np.where(mesh == 0, np.nan, mesh)

        # Create a mesh grid (x, y coordinates)
        x, y = np.mgrid[0 : mesh.shape[0], 0 : mesh.shape[1]]

        # Visualize using Mayavi
        mlab.surf(x, y, mesh)

        # Configure the viewing angle for better perception
        mlab.view(azimuth=0, elevation=0)

        # Display the 2D visualization
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
