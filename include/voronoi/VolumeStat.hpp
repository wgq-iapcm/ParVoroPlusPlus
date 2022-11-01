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

#include "StandardIncludes.h"
#include <voro++/voro++.hh>

double 
volume_stat(diy::Master& master)
{
    size_t nblocks = master.size();
    double total = 0;

    for(int i = 0; i < nblocks; i++)
    {
	auto b = (ParticleBlock*)master.block(i);

	voro::voronoicell c;
        voro::c_loop_order clo(*(b->con), b->po);
	if(clo.start()) do if((b->con->compute_cell(c,clo))) {
		total += c.volume();
	}while (clo.inc());
    }

    return total;
}



