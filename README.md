## ParVoro++: A Scalable Parallel Algorithm for Constructing 3D Voronoi Tessellations Based on KD-Tree Decomposition

## Licensing

ParVoro++ is released as open source software under a BSD style [license](./LICENSE).

## Dependencies

CMake-3.10+, MPICH-3.0, HDF5-1.8.8, zlib-1.2.8, Silo-4.10.2

## Installation

You will need a C++ compiler to compile the code and the CMake makefile generator tool to 
generate a Unix Makefile out of the project file. 

To build the executable, create an empty build directory first. In this directory, run

    edit CMakeLists.txt
    
    cmake [path_to_source] 

where [path_to_source] should point to the src/ directory in this source distribution.
After the CMake tool has successfully created the makefile, run

    make 

    make install

to build the executable.

or 

    make install DESTDIR=/PATH_TO_INSTALL 

## Execution

cd path/to/ParVoro++/install/bin/

```
Edit KD_AUTO.sh or KD_GHOST.sh

```
./KD_AUTO.sh

You can visualize the Voronoi tessellation with VisIt or ParaView.

## Acknowledgments

This work was supported by National Natural Science Fund of China  under grant No.61403036.

## Questions / contact info

Guoqing Wu
Institute of Applied Physics and Computational Mathematics, Beijing, China
wu_guoqing@iapcm.ac.cn




