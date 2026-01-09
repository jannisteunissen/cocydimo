cocydimo: Conducting Cylinder Discharge Model
=======

This code implements a reduced discharge model for streamer discharges, in which the channels are represented by a collection of conducting cylindrical segments. The conductivity and electric potential together evolve on a numerical mesh, as described in this [paper](https://doi.org/10.1016/j.cpc.2025.109733).

Getting the code and compiling
==

To compile the code, the following packages are required:

* Recent versions of `gcc` and `gfortran`. Other compilers can probably be used as well by modifying the build process.
* `f2py` and `meson` and python3 development files (e.g., `python3-dev' on Debian/Ubuntu) to create a Python module

The code (including the [afivo library](https://github.com/MD-CWI/afivo)) can be downloaded with:

    git clone --recurse-submodules https://github.com/jannisteunissen/cocydimo.git

Compiling for the first time can take some time:

    cd cocydimo
    make

Running
==

The model can be run using

    ./2d_model # 2D axisymmetric, no grid refinement
    ./3d_model # 3D with adaptive mesh refinement and branching

To see a full list of options, use:

    ./2d_model.py -h
    ./3d_model.py -h

Example
==

![alt text](https://github.com/jannisteunissen/cocydimo/blob/main/media/ICPIG_2025_example.png?raw=true)
