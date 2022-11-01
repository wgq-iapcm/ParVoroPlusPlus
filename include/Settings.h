///////////////////////////////////////////////////////////////////////////////
//
// ParVoro++: A scalable parallel algorithm for constructing 3D Voronoi tessellations based on KD-tree decomposition
//
// Guoqing Wu
// High Performance Computing Center, Institute of Applied Physics and Computational Mathematics, Beijing
// wu_guoqing@iapcm.ac.cn
//
///////////////////////////////////////////////////////////////////////////////
#pragma once

#define PROGRAM_VERSION_STRING "1.0.0"

#define EPSILON 0.00001

#define TOLERANCE 1e-10

#define isPowerOf2(n) (n>0? (n&(n-1))==0 : 0)

#define MAX_TIMES 10

