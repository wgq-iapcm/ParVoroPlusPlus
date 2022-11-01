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

void 
local_cells_ghost(ParticleBlock* b,
                  const diy::Master::ProxyWithLink& cp,
                  bool wrap,
                  int nx,
                  int ny,
                  int nz,
                  bool test)
{
    int nx_,ny_,nz_;
    voro::pre_container pcon(b->global_bounds.min[0], b->global_bounds.max[0], b->global_bounds.min[1], b->global_bounds.max[1], b->global_bounds.min[2], b->global_bounds.max[2], wrap, wrap, wrap);
    for(int i = 0; i <  b->num_orig_particles; i++)
        pcon.put(i,  b->particles[i][0],  b->particles[i][1],  b->particles[i][2]);
    pcon.guess_optimal(nx_,ny_,nz_);

    if(nx != 0 && ny != 0 && nz != 0)
    {
        nx_ = nx;
        ny_ = ny;
        nz_ = nz;
    }
    
    b->con = new voro::container(b->global_bounds.min[0], b->global_bounds.max[0], b->global_bounds.min[1], b->global_bounds.max[1], b->global_bounds.min[2], b->global_bounds.max[2], nx_, ny_, nz_, wrap, wrap, wrap, 8);

    // add walls here as you wish
    {
        //voro::wall_sphere* sphere = new voro::wall_sphere(0.5, 0.5, 0.5, 0.3);
        //b->con->add_wall(*sphere);
    }

    pcon.setup(b->po,*(b->con));

    for(int j = b->num_orig_particles; j < b->num_particles; j++)
        b->con->put(j, b->particles[j][0], b->particles[j][1], b->particles[j][2]);

    if(test)
    {
        voro::voronoicell c;
        voro::c_loop_order clo(*b->con, b->po);
        if(clo.start()) do if(b->con->compute_cell(c, clo));
        while(clo.inc());
    }
}

void Voronoi_ghost(diy::Master& master,
                   bool wrap,
                   int nx,
                   int ny,
                   int nz,
                   bool test)
{
    master.foreach([&](ParticleBlock* b, const diy::Master::ProxyWithLink & cp)
                  {local_cells_ghost(b, cp, wrap, nx, ny, nz, test);});
}


