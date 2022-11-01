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

#include "particle.hpp"
#include <diy/reduce.hpp>

// get total number of particles
size_t
get_total_num(diy::mpi::communicator& world, 
              diy::Master& master)
{
    size_t nparticles = 0;
    for(unsigned i = 0; i < master.size(); i++)
    	nparticles += ((ParticleBlock*)master.block(i))->particles.size();

    size_t tot_particles;
    diy::mpi::all_reduce(world, nparticles, tot_particles, std::plus<size_t>());

    return tot_particles;
}

void
wrap_pt(diy::Point<float,3>& rp, 
        diy::Direction wrap_dir, 
        diy::ContinuousBounds& domain)
{
    rp[0] -= wrap_dir[0] * (domain.max[0] - domain.min[0]);
    rp[1] -= wrap_dir[1] * (domain.max[1] - domain.min[1]);
    rp[2] -= wrap_dir[2] * (domain.max[2] - domain.min[2]);
}

double distance(double* u, double* v)
{
    double n = 0;
    for (int i = 0; i < 3; ++i)
    {
        n += (u[i] - v[i]) * (u[i] - v[i]);
    }
    return sqrt(n);
}



