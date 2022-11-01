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
#include <diy/types.hpp>

size_t  
gen_particles_uniform(const diy::ContinuousBounds& domain, // global data bounds
                      size_t n,
                      std::vector<diy::Point<float,3> >&  particles)
{
    particles.resize(n);

    std::random_device rd;
    std::mt19937 gen(rd());
    for(size_t i = 0; i < n; ++i)
        for(unsigned j = 0; j < 3; ++j)
        {
            std::uniform_real_distribution<float> dis(domain.min[j], domain.max[j]);
            particles[i][j] = dis(gen);
        }

    return n;
}


