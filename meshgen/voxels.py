import trimesh as trm
import numpy as np
from scipy import ndimage
import utilities
import mesher

from tqdm import tqdm  # Import tqdm
from concurrent.futures import ProcessPoolExecutor


def load_mesh(path):
    """
    Load a mesh from a specified file path using trimesh.

    This function takes a file path as input, and uses the trimesh library to load the mesh.
    It supports various mesh formats that trimesh can handle.

    Parameters:
    path (str): The file path of the mesh to be loaded.

    Returns:
    trimesh.Trimesh: A trimesh object representing the loaded mesh.
    """
    # Load the mesh from the given file path using trimesh
    mesh = trm.load(path)
    return mesh


def get_leading_direction(shape):
    """
    Identify the leading direction of a 3D shape based on its dimensions.

    This function calculates which of the three dimensions (x, y, z) of a shape
    is the largest and returns the index of that dimension. It uses numpy's argmax
    function to find the index of the maximum value in the shape dimensions.

    Parameters:
    shape (tuple or list): A 3-element iterable representing the dimensions of the shape.
                           Typically, this would be in the format (x, y, z).

    Returns:
    int: The index of the leading dimension (0 for x, 1 for y, 2 for z).
    """
    # Find the index of the maximum value in shape dimensions
    direction_idx = np.argmax([shape[0], shape[1], shape[2]])

    return direction_idx


def fill_extended_mesh(original_mesh, new_mesh):
    """
    Embed an original mesh into a new, larger mesh.

    This function takes an original mesh and a new mesh (which is larger in size).
    It calculates the leading direction (dimension) of the original mesh and embeds it
    into the center of the new mesh along the remaining dimensions. The function assumes
    that the meshes are 3D and represented in a NumPy array-like format.

    Parameters:
    original_mesh (numpy.ndarray): The original mesh to be embedded.
    new_mesh (numpy.ndarray): The new, larger mesh where the original mesh will be embedded.

    Returns:
    numpy.ndarray: The new mesh with the original mesh embedded in it.
    """
    # Determine the shapes of the original and new meshes
    original_shape = original_mesh.shape
    new_shape = new_mesh.shape

    # Identify the leading direction (dimension) of the original mesh
    direction_idx = get_leading_direction(original_mesh.shape)
    remaining_indices = [i for i in range(3) if i != direction_idx]

    # Fill the original mesh inside its surface
    original_mesh = fill_mesh_inside_surface(original_mesh)

    # TODO: this is just to make things work
    if (new_shape[0] * new_shape[1] * new_shape[2]) < (
        original_shape[0] * original_shape[1] * original_shape[2]
    ):
        return original_mesh

    # Embed the original mesh into the new mesh based on the leading direction
    if remaining_indices == [0, 1]:
        # Calculate starting points for embedding
        start_0 = (new_shape[0] - original_shape[0]) // 2
        start_1 = (new_shape[1] - original_shape[1]) // 2

        # Embed the original mesh into the new mesh
        new_mesh[
            start_0 : start_0 + original_shape[0],
            start_1 : start_1 + original_shape[1],
            :,
        ] = original_mesh

    elif remaining_indices == [0, 2]:
        start_0 = (new_shape[0] - original_shape[0]) // 2
        start_1 = (new_shape[2] - original_shape[2]) // 2

        new_mesh[
            start_0 : start_0 + original_shape[0],
            :,
            start_1 : start_1 + original_shape[2],
        ] = original_mesh

    else:
        start_0 = (new_shape[1] - original_shape[1]) // 2
        start_1 = (new_shape[2] - original_shape[2]) // 2

        new_mesh[
            :,
            start_0 : start_0 + original_shape[1],
            start_1 : start_1 + original_shape[2],
        ] = original_mesh

    return new_mesh


def get_lbm_shape(original_shape):
    """
    Adjust the shape of a 3D object to make it suitable for Lattice Boltzmann Method (LBM) simulations.

    This function checks if the leading dimension of the original shape is a multiple of 128.
    If not, it raises a ValueError. Then, it adjusts the remaining dimensions to be the smallest
    multiple of 32 that is greater than or equal to the maximum of these dimensions, keeping the
    leading dimension unchanged.

    Parameters:
    original_shape (tuple or list): A 3-element iterable representing the original shape dimensions (x, y, z).

    Returns:
    tuple: The new shape, adjusted for LBM simulations.

    Raises:
    ValueError: If the leading dimension of the original shape is not a multiple of 128.
    """
    # Identify the leading direction (dimension) of the shape
    direction_idx = get_leading_direction(original_shape)

    # Check if the leading dimension is a multiple of 128, raise an error if not
    if (original_shape[direction_idx] % 128) != 0:
        error_msg = f"Shape {original_shape} is not suitable for LBM simulation."
        raise ValueError(error_msg)

    # Find the maximum dimension size that is not the leading direction
    # and calculate the smallest multiple of 32 that is greater than or equal to that size
    max_remaining_size = max(
        original_shape[i] for i in range(len(original_shape)) if i != direction_idx
    )
    rounded_size = ((max_remaining_size + 31) // 32) * 32

    # Construct the new shape, adjusting non-leading dimensions
    new_shape = [
        original_shape[i] if i == direction_idx else rounded_size
        for i in range(len(original_shape))
    ]
    new_shape = tuple(new_shape)

    return new_shape


def complete_mesh(original_mesh):
    """
    Adjust an original mesh to make it suitable for Lattice Boltzmann Method (LBM) simulations.

    This function first determines a new shape for the mesh that is compatible with LBM simulations
    using the 'get_lbm_shape' function. It then creates an empty mesh of this new shape and embeds
    the original mesh into it using the 'fill_extended_mesh' function.

    Parameters:
    original_mesh (numpy.ndarray): The original mesh to be adjusted.

    Returns:
    numpy.ndarray: The new mesh, adjusted in shape and filled, suitable for LBM simulations.
    """
    # Get the shape of the original mesh
    original_shape = original_mesh.shape

    # Determine the new shape suitable for LBM simulations
    new_shape = get_lbm_shape(original_shape)

    # Create an empty mesh with the new shape
    empty_mesh = np.zeros(new_shape, dtype=bool)

    # Fill the original mesh into the new, empty mesh
    new_mesh = fill_extended_mesh(original_mesh, empty_mesh)

    return new_mesh


def fill_mesh_inside_surface(mesh):
    """
    Fill the internal voids of a 3D mesh.

    This function takes a 3D mesh, represented as a NumPy array, and fills its internal voids.
    It uses the 'binary_fill_holes' function from SciPy's ndimage module to identify and fill
    these voids. The input mesh is expected to be a boolean array, where True indicates the
    presence of the mesh and False indicates empty space.

    Parameters:
    mesh (numpy.ndarray): A 3D boolean array representing the mesh.

    Returns:
    numpy.ndarray: The modified mesh with internal voids filled.
    """
    # Fill the internal voids of the mesh using SciPy's binary_fill_holes function
    mesh[ndimage.binary_fill_holes(mesh)] = True

    return mesh


def is_leading_dir(direction, leading_direction):
    """
    Check if a given direction is the leading direction.

    Parameters:
    direction (int): The direction to be checked (0 for x, 1 for y, 2 for z).
    leading_direction (int): The leading direction against which to check.

    Returns:
    bool: True if 'direction' is the same as 'leading_direction', otherwise False.
    """
    # Compare the given direction with the leading direction
    return direction == leading_direction


def calculate_margin(original_bounds, new_bounds, leading_direction):
    """
    Calculate the margin between the original and new bounds of a mesh.

    This function calculates the margin for each dimension (x, y, z). For the leading direction,
    the margin is set to [0, 0]. For the other two directions, it calculates the absolute difference
    between the original and new bounds.

    Parameters:
    original_bounds (array-like): The original bounds of the mesh.
    new_bounds (array-like): The new bounds of the mesh after processing.
    leading_direction (int): The index of the leading direction (0 for x, 1 for y, 2 for z).

    Returns:
    numpy.ndarray: A 2D array representing the margins for each direction.
    """
    margin = []
    # Calculate margin for each dimension
    for j in range(3):
        if is_leading_dir(j, leading_direction):
            margin.append([0, 0])
        else:
            margin.append(
                [
                    abs(new_bounds[0][j] - original_bounds[0][j]),
                    abs(new_bounds[1][j] - original_bounds[1][j]),
                ]
            )

    return np.array(margin)


def calculate_segment_direction_and_increment(leading_direction, voxel_size):
    """
    Calculate the direction vector and increment for segmenting a mesh.

    This function determines the direction vector and increment value for segmenting a mesh along
    its leading direction. The direction vector is a unit vector along the leading direction, and
    the increment is calculated based on the voxel size and the leading direction.

    Parameters:
    leading_direction (int): The index of the leading direction (0 for x, 1 for y, 2 for z).
    voxel_size (int): The size of a single voxel.

    Returns:
    tuple: A tuple containing the direction vector and the increment value.
    """
    # Create a direction vector based on the leading direction
    segment_direction = np.array(
        [int(is_leading_dir(i, leading_direction)) for i in range(3)]
    )

    # Calculate increment based on voxel size and leading direction
    segment_increment = np.array(
        [voxel_size * int(is_leading_dir(i, leading_direction)) for i in range(3)]
    )

    return segment_direction, segment_increment


def slice_mesh(mesh, direction, plane_origin, increment):
    """
    Slice a mesh in a specified direction using two parallel planes.

    This function slices a 3D mesh using two parallel planes. The first plane is defined by the
    'plane_origin' point and the 'direction' vector. The second plane is parallel to the first and
    offset by the 'increment' value. The section of the mesh that lies between these two planes is returned.

    Parameters:
    mesh (trimesh.Trimesh): The mesh to be sliced.
    direction (array-like): A 3-element array representing the normal vector of the slicing plane.
    plane_origin (array-like): A 3-element array representing a point on the first slicing plane.
    increment (float): The distance between the two parallel slicing planes.

    Returns:
    trimesh.Trimesh: The portion of the mesh that lies between the two slicing planes.
    """
    # Slice the mesh with the first plane in the given direction
    mesh_one_direction = trm.intersections.slice_mesh_plane(
        mesh, direction, plane_origin
    )

    # Slice the resulting mesh with a second plane, in the opposite direction and offset by increment
    mesh_second_direction = trm.intersections.slice_mesh_plane(
        mesh_one_direction, direction * (-1), plane_origin + increment
    )

    return mesh_second_direction


def split_mesh(mesh, voxel_size, n_segments):
    """
    Split a mesh into several segments along its leading direction.

    This function divides a given mesh into a specified number of segments, ensuring that each
    segment has a uniform size along the mesh's leading direction. The leading direction is
    determined based on the mesh bounds. The function calculates the direction and increment
    for segmenting, then slices the mesh accordingly, and computes margins for each segment.

    Parameters:
    mesh (trimesh.Trimesh): The mesh to be split.
    voxel_size (int): The size of a single voxel, used in calculating margins.
    n_segments (int): The number of segments to divide the mesh into.

    Returns:
    tuple: A tuple containing two lists -
           1. submeshes: a list of trimesh.Trimesh objects, each representing a segment.
           2. margins: a list of integers, each representing the margin of a corresponding segment.
    """
    # Determine the bounds of the mesh and identify the leading direction
    bounds = mesh.bounds
    leading_direction = get_leading_direction(tuple(bounds[1] - bounds[0]))

    # Calculate the direction and increment for segmenting the mesh
    segment_direction, segment_increment = calculate_segment_direction_and_increment(
        leading_direction, voxel_size
    )

    submeshes = []  # List to store the segmented parts of the mesh
    margins = []  # List to store the margin of each segment

    # Iterate over the number of segments to split the mesh
    for i in range(n_segments):
        # Calculate the origin for slicing the mesh
        plane_origin = bounds[0] + i * segment_increment

        # Slice the mesh to create a segment
        sliced_segment = slice_mesh(
            mesh, segment_direction, plane_origin, segment_increment
        )
        sliced_bounds = sliced_segment.bounds

        # Calculate the margin for the current segment
        margins.append(
            calculate_margin(bounds, sliced_bounds, leading_direction) // voxel_size
        )
        submeshes.append(sliced_segment)

    return submeshes, margins


def complete_segment(segment, margins):
    """
    Complete a segment of a mesh by adding margins in each direction.

    This function takes a segment of a mesh and adds specified margins around it in all three
    dimensions (x, y, z). The margins are added by creating a larger array and embedding the
    original segment in it, centered with the given margins.

    Parameters:
    segment (numpy.ndarray): A 3D array representing a segment of the mesh.
    margins (numpy.ndarray): A 2D array where each row represents the margins [start, end]
                             for each dimension (x, y, z).

    Returns:
    numpy.ndarray: The segment with added margins, represented as a larger 3D array.
    """
    a1 = margins[0][0]  # offsets in the x direction
    a2 = margins[0][1]

    b1 = margins[1][0]  # offsets in the y direction
    b2 = margins[1][1]

    c1 = margins[2][0]  # offsets in the z direction
    c2 = margins[2][1]

    # Calculate the new shape including the margins
    new_shape = (
        segment.shape[0] + a1 + a2,
        segment.shape[1] + b1 + b2,
        segment.shape[2] + c1 + c2,
    )

    # Create a new array with the new shape, initialized as empty (False)
    new_array = np.zeros(new_shape, dtype=bool)

    # Embed the original segment into the center of the new array, surrounded by the margins
    new_array[
        a1 : a1 + segment.shape[0],
        b1 : b1 + segment.shape[1],
        c1 : c1 + segment.shape[2],
    ] = segment

    return new_array


def fill_slice(mesh_slice, leading_direction):
    """
    Fill a mesh slice and expand it along the leading direction.

    This function takes a slice of a mesh (represented as a 2D array), fills its internal voids,
    and then expands this slice into a 3D array by adding an extra dimension along the specified
    leading direction. This step is typically used to prepare the slice for processes that require
    3D input or to maintain consistency in the mesh's dimensionality across operations.

    Parameters:
    mesh_slice (numpy.ndarray): A 2D array representing a slice of the mesh.
    leading_direction (int): The axis along which to expand the slice (0 for x, 1 for y, 2 for z).

    Returns:
    numpy.ndarray: A 3D array representing the filled and expanded slice of the mesh.
    """
    # Fill the internal voids of the slice
    mesh_slice = fill_mesh_inside_surface(mesh_slice)

    # Expand the 2D slice into a 3D array along the specified leading direction
    mesh_slice = np.expand_dims(mesh_slice, axis=leading_direction)

    return mesh_slice


def calculate_voxel_size(mesh, res):
    """
    Calculate the voxel size for a mesh based on the specified resolution.

    This function calculates the voxel size needed to represent a mesh at a given resolution. It
    finds the maximum dimension of the mesh and divides it by a factor derived from the resolution.
    The function is used to determine the appropriate voxel size for voxelizing the mesh.

    Parameters:
    mesh (trimesh.Trimesh): The mesh for which the voxel size is being calculated.
    res (float): The desired resolution, used as a factor in the calculation.

    Returns:
    float: The calculated voxel size.
    """
    # Calculate the maximum size of the mesh in any dimension
    maxima = [np.max(mesh.bounds[:, i]) - np.min(mesh.bounds[:, i]) for i in range(3)]
    h = np.max(maxima)

    # Calculate the voxel size based on the maximum size and desired resolution
    voxel_size = h / (128 * res)

    return voxel_size


def voxelize_without_splitting(mesh, voxel_size):
    """
    Convert a mesh into a voxelized representation without splitting.

    This function takes a mesh and converts it into a voxelized form using the specified voxel size.
    The voxelization process divides the mesh into cubic cells of the given size. It does not split
    the mesh into separate segments but rather represents the entire mesh in a voxelized format.

    Parameters:
    mesh (trimesh.Trimesh): The mesh to be voxelized.
    voxel_size (float): The edge length of each voxel in the voxelized representation.

    Returns:
    numpy.ndarray: A 3D array representing the voxelized mesh, where each cell is a voxel.
    """
    # Convert the mesh into a voxelized form with the specified voxel size
    voxelized_object = mesh.voxelized(voxel_size).matrix

    return voxelized_object


def process_submesh(submsh, margin, voxel_size, leading_direction):
    """
    Process a single submesh for voxelization.

    This function handles the voxelization of a single mesh segment. It voxelizes the segment,
    completes it by adding necessary margins, and adjusts the slicing based on the leading direction.
    Finally, it fills the slice to ensure a solid segment.

    Parameters:
    submsh (trimesh.Trimesh): The submesh segment to be voxelized.
    margin (numpy.ndarray): The margin to be added around the voxelized segment.
    voxel_size (float): The size of each voxel.
    leading_direction (int): The leading direction for slicing (0 for x, 1 for y, 2 for z).

    Returns:
    numpy.ndarray: A 3D array representing the processed and filled submesh segment.
    """
    # Voxelizing the submesh without splitting
    voxelized_segment = voxelize_without_splitting(submsh, voxel_size)

    # Completing the segment by adding necessary margins
    comp = complete_segment(voxelized_segment, margin.astype(int))

    # Adjusting slicing based on the leading direction
    if leading_direction == 0:
        comp = comp[0, 1:-1, 1:-1]  # Slicing for the x-direction
    elif leading_direction == 1:
        comp = comp[1:-1, 0, 1:-1]  # Slicing for the y-direction
    else:  # leading_direction == 2
        comp = comp[1:-1, 1:-1, 0]  # Slicing for the z-direction

    # Filling the slice for further processing
    return fill_slice(comp, leading_direction)


def voxelize_with_splitting(mesh, voxel_size, split, num_processes=1, dim=3):
    """
    Voxelizes a mesh with splitting along the leading direction.

    This function voxelizes a mesh by first splitting it into segments along its leading direction.
    Each segment is voxelized independently, and then they are concatenated back together. The function
    also handles the edge layers separately to ensure a complete voxelization.

    Parameters:
    mesh (trimesh.Trimesh): The mesh to be voxelized.
    voxel_size (float): The size of each voxel.
    split (int): The number of segments to split the mesh into.
    num_processes(int, optional): The number of subprocesses the main process should be divided into.
                                  Default 1 = no subprocesses.
    dim (int, optional): Desired dimension of the output mesh.
                         Default 3 = output mesh is a 3D mesh.

    Returns:
    numpy.ndarray: A 3D array representing the voxelized mesh.
    """
    # Determine the bounds of the mesh and identify the leading direction
    bounds = mesh.bounds
    leading_direction = get_leading_direction(tuple(bounds[1] - bounds[0]))

    if dim == 2 and leading_direction == 2:
        error_msg = f"2D mesh not supported (not located in XY plane)."
        raise ValueError(error_msg)

    # Split the mesh into segments along the leading direction
    submeshes, margins = split_mesh(mesh, voxel_size, n_segments=split)

    # Parallel processing of submeshes
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        # Prepare arguments for each submesh processing
        futures = [
            executor.submit(
                process_submesh, submsh, margin, voxel_size, leading_direction
            )
            for submsh, margin in zip(submeshes, margins)
        ]

        # Collecting results with progress display
        comps = [
            future.result()
            for future in tqdm(futures, total=len(submeshes), desc="Voxelizing")
        ]

    # Concatenate the processed segments along the leading direction
    ary = np.concatenate(comps, axis=leading_direction)

    # Adjust the concatenated array based on the leading direction
    if leading_direction == 0:
        ary = ary[1:-1, :-1, :]
    elif leading_direction == 1:
        ary = ary[:-1, 1:-1, :]
    else:
        ary = ary[:-1, :, 1:-1]

    # Adjust edge layers based on the leading direction
    if leading_direction == 0:
        first_layer = fill_slice(ary[0, :, :], 0)
        last_layer = fill_slice(ary[-1, :, :], 0)
        ary = np.concatenate([first_layer, ary, last_layer], axis=0)
    elif leading_direction == 1:
        first_layer = fill_slice(ary[:, 0, :], 1)
        last_layer = fill_slice(ary[:, -1, :], 1)
        ary = np.concatenate([first_layer, ary, last_layer], axis=1)
    else:  # leading_direction == 2
        first_layer = fill_slice(ary[:, :, 0], 2)
        last_layer = fill_slice(ary[:, :, -1], 2)
        ary = np.concatenate([first_layer, ary, last_layer], axis=2)

    # If output mesh is a 2D mesh, decrease the array in dimension
    if dim == 2:
        ary = np.squeeze(ary, axis=2)

    return ary


def voxelize_mesh(name, res=1, split=None, num_processes=1, dim=3, **kwargs):
    """
    Load a mesh from a file, voxelize it with or without splitting, and adjust it for LBM simulations.

    This function loads a mesh from the specified path and calculates an appropriate voxel size based
    on the mesh's bounds and a resolution factor. It then voxelizes the mesh, either as a whole or by
    splitting it into segments, depending on the 'split' parameter. The voxelized mesh is adjusted
    using the 'complete_mesh' function to make it suitable for LBM simulations.

    Parameters:
    path (str): The file path of the mesh to be loaded and voxelized.
    res (int, optional): A resolution factor for voxelization. Default is 1.
    split (int, optional): The number of segments to split the mesh into before voxelization.
                           If None, the mesh is voxelized without splitting. Default is None.
    num_processes(int, optional): The number of subprocesses the main process should be divided into.
                                  Default 1 = no subprocesses.
    dim (int, optional): Desired dimension of the output mesh.
                         Default 3 = output mesh is a 3D mesh.

    Returns:
    numpy.ndarray: The voxelized and adjusted mesh, suitable for LBM simulations.
    """

    modified_kwargs = kwargs
    modified_kwargs["resolution"] = res

    stl_file_path = mesher.gmsh_surface(name, **modified_kwargs)

    # Load the mesh from the specified path
    mesh = load_mesh(stl_file_path)

    # Calculate the voxel size based on the mesh's bounds and resolution factor
    voxel_size = calculate_voxel_size(mesh, res)

    if split is None:
        # Voxelize the mesh without splitting
        output = voxelize_without_splitting(mesh, voxel_size)
    else:
        # Calculate the bounding box dimensions of the mesh
        bounds = mesh.bounds
        x, y, z = bounds[0]
        dx, dy, dz = bounds[1] - bounds[0]

        # Generate a box-shaped mesh based on the bounding box dimensions and voxel size
        name = mesher.box_stl(name, x, y, z, dx, dy, dz, voxel_size)
        boxed_mesh = trm.load(name)
        # Voxelize the mesh with splitting along the leading direction
        output = voxelize_with_splitting(
            boxed_mesh, voxel_size, split, num_processes=num_processes, dim=dim
        )

    utilities.delete_stl_files("../meshgen/")

    # If the mesh is meant for 2D simulations, no further operations have to be done
    if dim == 2:
        lbm_mesh = output
        return lbm_mesh

    # Adjust the voxelized mesh for 3D LBM simulations
    lbm_mesh = complete_mesh(output)

    return lbm_mesh


if __name__ == "__main__":
    import time

    start_time = time.time()

    # name = "tcpc_classic"
    name = "2dtcpc"
    msh = voxelize_mesh(
        name, res=2, split=2 * 128, num_processes=6, dim=2, offset=0.3, h=0.01
    )
    print(msh.shape)

    print(time.time() - start_time)

    # utilities.array_to_textfile(~msh, 'output')
    utilities.vis(msh, dim=2)
#     TODO: optimize output array - delete coordinates
#     TODO: try benchmark for easy tcpc
