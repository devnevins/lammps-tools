# lammps-tools
As part of my graduate work I created some tools for the molecular dynamics simulation program [LAMMPS](https://www.lammps.org/). These utilities greatly increased the utility of using LAMMPS in the geologic domain. There are utilities for creating intial configurations using the skew-start methodology along with analysis code for computing coordination numbers, their associated cutoffs, 

# Citing this code
If you use or modify any of this code please cite this code with the following reference. This is a reference to my Dissertation which contains the details of what the code did and what algorithms were used.

```
Nevins, D. I. (2009). Understanding silicate geoliquids at high temperatures and pressures through molecular dynamics simulations. University of California, Santa Barbara.
```

# Why am I doing this?
I've been asked multiple times for this code so I decided to put it here so people can get to it. Enjoy!

# Note
- This code is OLD (circa 2009) and is for Python 2.x. You'll have to port to Python 3. If you do, please issue a pull request!
- It was written by a graduate student. You have been warned!
- This is not written to any specific recommended approach to Python (e.g. no proper directory structure). See previous statement.
- This code is specific to a custom version of LAMMPS that was called LAMMPSG. Things may have changed significantly.
- Utilizes (probably an old version of) Numarray for numerical computation and manipulation of arrays.
- Defines cutoff distances for computing interactions between different particle types.

# Included code

## Coordination Statistics (cn.py)

`cn.py` computes coordination statistics for LAMMPSG. Key features include initialization and configuration using `ConfigParser`, computation of coordination statistics based on distances and specified cutoffs, handling of periodic boundary conditions, and efficient numerical operations using Numarray. The script processes position files, supports skipping initial runs, and outputs coordination statistics for each species pair. It also includes features for handling different species, error checking, and optimizations for performance. The output includes coordination numbers and additional statistics for specific scenarios (e.g., T around O, O around T). The script is tailored for LAMMPSG simulations with dependencies on external configuration files.

## Cutoff (cutoff.py)

`cutoff.py` computes cutoff values for the coordination number program (`cn.py`). It reads an RDF file, extracts relevant columns, analyzes the data, and generates cutoff values for coordination numbers between different material pairs. The computed cutoff values are then outputted in a format compatible with the coordination number program. The script is designed to be executed with the command: `./cutoff.py materialfilename`. 

### Note
- There's a missing import statement for the `re` module that you may need to add at the beginning of the script.

## Diffusion (diffusion.py)

`diffusion.py` calculates the diffusion coefficient from mean squared displacement (MSD) data. It includes a Material class for handling information about materials, a conversion function for mass units, and a linear regression function. In the main block, the script parses command-line arguments, reads MSD data, performs linear regression, converts the slope to the diffusion coefficient, and writes the results to a CSV file.

## Radial Distribution Function (rdf.py)

`rdf.py` calculates the Radial Distribution Function (RDF) in molecular dynamics simulations. It initializes parameters, defines classes for particles and materials, includes a function to convert mass units, and computes RDF based on particle positions over multiple time steps. The main script reads configuration and position files, normalizes the RDF, and outputs the results to a file. The script is well-commented and handles cases where parameters are missing. Overall, it is a tool for analyzing molecular dynamics simulation data, focusing on the RDF to describe particle density variations with distance.

## Skew-start (skewstart.py)

`skewstart.py` generates initial conditions for molecular dynamics simulations using the skew-start method proposed by Keith Refson. The code takes configuration parameters from a specified file, including details about materials, box size, temperature, and output format. It then employs the skew-start method to generate positions for particles, assigns velocities based on a specified temperature, and removes system-wide momentum. The final step involves outputting the generated initial conditions in a chosen format (IODYN, LAMMPS, or MOLDY). The code includes classes for handling different output formats, abstract functions for position and velocity generation, and a main section that orchestrates the entire process. Overall, the code serves the purpose of preparing the initial state for molecular dynamics simulations with flexibility in output format selection.

## Viscosity Calculation (viscosity.py)

`viscosity.py` performs viscosity calculations using numerical integration methods. It includes classes for window spacing and sliding windows for data processing. The main execution reads configuration parameters, processes viscosity calculations using rectangular, trapezoidal, and Simpson's integration, and outputs estimated viscosity values. The script is structured with clear sections for imports, classes, functions, and the main execution. It requires command line arguments for a configuration file and input data file. The final output includes detailed results and a table with timestep, autocorrelation function (ACF) values, and integrated values for each method.

















