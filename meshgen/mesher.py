import gmsh
import shutil


def modify_geo_file(input_file_path, output_file_path, **kwargs):
    """
    Modify a GEO file by replacing placeholders with specified values.

    This function reads a GEO file, replaces specified placeholders with their corresponding values,
    and writes the modified content to a new file. Placeholders in the GEO file are expected to be
    in the format 'DEFINE_VARIABLE', where 'VARIABLE' is the name of the variable to be replaced.

    Parameters:
    input_file_path (str): The path to the original GEO file.
    output_file_path (str): The path where the modified GEO file will be saved.
    **kwargs: Arbitrary keyword arguments where keys are the variable names (as they appear in the
              'DEFINE_VARIABLE' placeholders) and values are the replacement values.

    Example:
    modify_geo_file('original.geo', 'modified.geo', angle=45, length=10)
    """

    # Open and read the content of the original GEO file
    with open(input_file_path, "r") as file:
        file_data = file.read()

    # Iterate through each keyword argument to replace placeholders
    for variable, value in kwargs.items():
        # Construct the placeholder string
        placeholder_variable = "DEFINE_" + variable.upper()
        replace_string = str(value)

        # Replace the placeholder with the actual value
        file_data = file_data.replace(placeholder_variable, replace_string)

    # Write the modified data to the new GEO file
    with open(output_file_path, "w") as file:
        file.write(file_data)


def box_stl(stl_file, x, y, z, dx, dy, dz, voxel_size):
    """
    Generate a boxed STL file from a template and specified parameters.

    This function creates a boxed version of an STL file by modifying a template GEO file with
    specific dimensions and then using Gmsh to generate the mesh and export it as an STL file.
    The template GEO file should contain placeholders for the STL file name and box dimensions.

    Parameters:
    stl_file (str): The name of the STL file (without the '.stl' extension).
    x, y, z (float): The coordinates of the bottom-left-front corner of the box.
    dx, dy, dz (float): The dimensions of the box in the x, y, and z directions, respectively.
    voxel_size (float): The size of the voxels for the mesh.

    Returns:
    str: The path to the generated boxed STL file.
    """
    # Paths for the input template and output GEO file
    input_file_path = f"../meshgen/geo_templates/boxed_stl_template.geo"
    output_file_path = (
        f"../meshgen/geo_templates/temp/boxed_{stl_file}_template_filled.geo"
    )

    # Modify the GEO template with the provided dimensions and voxel size
    modify_geo_file(
        input_file_path,
        output_file_path,
        stl_file=stl_file + ".stl",
        x=x,
        y=y,
        z=z,
        dx=dx,
        dy=dy,
        dz=dz,
        voxel_size=voxel_size,
    )

    # Copy the original STL file to the working directory
    shutil.copyfile(
        f"../meshgen/stl_models/{stl_file}.stl",
        f"../meshgen/geo_templates/temp/{stl_file}.stl",
    )

    # Initialize Gmsh for mesh generation
    gmsh.initialize()

    # Open and process the modified GEO file
    gmsh.open(output_file_path)

    # Generate a 2D mesh from the GEO file
    gmsh.model.mesh.generate(2)

    # Export the generated mesh as an STL file
    stl_file_path = f"../meshgen/stl_models/temp/boxed_{stl_file}.stl"
    gmsh.write(stl_file_path)

    # Finalize Gmsh to free resources
    gmsh.finalize()

    return stl_file_path


def gmsh_surface(name_geo, **kwargs):
    """
    Generate a surface mesh using Gmsh based on a template GEO file and custom parameters.

    This function modifies a template GEO file with provided parameters, then uses Gmsh to generate
    a 2D mesh and exports it as an STL file. It's designed to work with GEO files containing placeholders
    that are replaced by the specified keyword arguments.

    Parameters:
    name_geo (str): The base name of the GEO file (without the '.geo' extension).
    **kwargs: Arbitrary keyword arguments where keys are the variable names (as they appear in the
              placeholders of the GEO file) and values are the replacement values.

    Returns:
    str: The path to the generated STL file.
    """
    # Initialize Gmsh for mesh generation
    gmsh.initialize()

    # Paths for the input template and output GEO file
    geo_file_path = f"../meshgen/geo_templates/{name_geo}_template.geo"
    output_geo_file_path = (
        f"../meshgen/geo_templates/temp/{name_geo}_template_filled.geo"
    )

    # Modify the GEO template with the provided parameters
    modify_geo_file(geo_file_path, output_geo_file_path, **kwargs)

    # Open and process the modified GEO file
    gmsh.open(output_geo_file_path)

    # Generate a 2D mesh from the GEO file
    gmsh.model.mesh.generate(2)

    # Export the generated mesh as an STL file
    stl_file_path = f"../meshgen/stl_models/{name_geo}.stl"
    gmsh.write(stl_file_path)

    # Finalize Gmsh to free resources
    gmsh.finalize()

    return stl_file_path


if __name__ == "__main__":
    pass
