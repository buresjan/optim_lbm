# Project Overview

This project is composed of three integral parts, each pertaining to a different aspect of computational modeling and optimization. Below is an outline of the project's structure and instructions on how to navigate and use the components.

## Directory Structure

The project is organized into three main directories, each housing a distinct segment of the codebase:

- `nomad`
- `meshgen`
- `lbm`

### NOMAD - Nonlinear Optimization with the MADS Algorithm

- **Location**: `nomad/`
- **Description**: This directory contains the NOMAD C++ package which implements the Mesh Adaptive Direct Search (MADS) optimization method. NOMAD is designed to handle blackbox optimization tasks where the objective and constraint functions are provided as external computer codes.

### Mesh Generation for LBM

- **Location**: `meshgen/`
- **Description**: Here, you will find custom Python code for mesh generation tailored for Lattice Boltzmann Methods (LBM) simulations. The code utilizes `gmsh` and `trimesh` libraries to create and manipulate mesh data effectively.

### Lattice Boltzmann Code

- **Location**: `lbm/`
- **Description**: This directory contains the Lattice Boltzmann code developed by the Department of Mathematics at the Faculty of Nuclear Sciences and Physical Engineering (FNSPE) of the Czech Technical University (CTU) in Prague. The Lattice Boltzmann Method (LBM) is a class of computational fluid dynamics (CFD) algorithms for fluid simulation.

## Environment Setup

To ensure consistent execution of code across different environments, an `environment.yml` file is provided at the root of the project directory. This file specifies all the necessary dependencies required to run the code.

To create the conda environment, follow these steps:

```bash
$ conda env create -f environment.yml
```

This command will create a new conda environment that includes all the libraries and versions specified in the `environment.yml` file.

## Getting Started

After setting up the environment, you can navigate to each directory to work with individual components of the project. Each directory may have its own set of instructions or README files for more detailed guidance on how to use the tools and code provided.
