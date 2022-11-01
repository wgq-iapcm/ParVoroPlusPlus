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

#include <diy/types.hpp>
#include <diy/log.hpp>
#include <diy/point.hpp>
#include <voro++/voro++.hh> 
#include "StandardIncludes.h"
#include "import/random/genRandomParticles.hpp"
#include "import/hdf5/readHDF5.hpp"

struct RemoteParticleInfo
{
    diy::BlockID  bid;     
    size_t        index;   
};

struct RemoteParticle
{
    diy::Point<float,3>  particle;
    RemoteParticleInfo   info;
};

struct OutCell
{
    int ijk;
    int q;
    int pid;
};

struct ParticleBlock
{
    typedef  diy::Point<float,3>  SimpleParticle;

    ParticleBlock()  {}
    ParticleBlock(const diy::ContinuousBounds& bounds_): bounds(bounds_)  {}

    static void* create_block()
    {
        return new ParticleBlock;
    }

    static void destroy_block(void* b_)
    {
        auto b = static_cast<ParticleBlock*>(b_);
        if(b->sphere != nullptr)
            delete b->sphere;
        if(b->cylinder != nullptr)
            delete b->cylinder;
        if(b->cone != nullptr)        
            delete b->cone;
        if(b->plane != nullptr)
            delete b->plane;
        if(b->con != nullptr)
            delete b->con;
        delete b;
    }

    // block data
    diy::BlockID                       bid;              // BlockID of the block
    int                                lid;              // the local id of the local block
    diy::ContinuousBounds              core{3};          // local block extents without ghost
    diy::ContinuousBounds              bounds{3};        // local block extents with ghost
    diy::ContinuousBounds              global_bounds{3}; // global data extents
    diy::ContinuousBounds              box{3};           // box in current round of particle redistribution
    // particle data
    std::vector<SimpleParticle>        particles;          // all particles, original plus those received from neighbors
    std::vector<RemoteParticleInfo>    infos;              // information of remote particles
    size_t                             num_orig_particles; // number of original particles in this block before any neighbor exchange
    size_t                             num_particles;      // current number of particles in this block after any neighbor exchange; 
						             // original particles appear first followed by received particles
    // voronoi container
    voro::container*                   con; 
    voro::particle_order               po;  
    std::list<OutCell>                 outcell;						           
    size_t                             enqueued;
    // voronoi walls
    voro::wall_sphere*                 sphere;
    voro::wall_cylinder*               cylinder;
    voro::wall_cone*                   cone;
    voro::wall_plane*                  plane;    
};

struct AddBlock
{
    AddBlock(diy::Master& master_): master(master_) {}
    ParticleBlock* operator()(int gid,                                 
                              const diy::ContinuousBounds& core,       
                              const diy::ContinuousBounds& bounds,    
                              const diy::ContinuousBounds& domain,     
                              const diy::RegularContinuousLink& link) 
    const
    {
    	ParticleBlock* b = new ParticleBlock(core);
    	diy::RegularContinuousLink*    l = new diy::RegularContinuousLink(link);
    	diy::Master&                   m = const_cast<diy::Master&>(master);

    	m.add(gid, b, l);

    	int rank = master.communicator().rank();

    	b->bid.proc           = rank;
    	b->bid.gid            = gid;
    	b->lid                = master.lid(gid);
    	b->core               = core;
    	b->bounds             = bounds;
    	b->global_bounds      = domain;
    	b->box                = domain;
    	b->num_orig_particles = 0;
    	b->num_particles      = 0;
    	b->con                = nullptr;
    	b->sphere             = nullptr;
    	b->cylinder           = nullptr;
    	b->cone               = nullptr;
    	b->plane              = nullptr;

    	return b;
    }

    diy::Master&  master;
};

struct AddAndGenerate: public AddBlock
{
    AddAndGenerate(diy::Master& master_): AddBlock(master_) {}

    void operator()(int gid,
                    const diy::ContinuousBounds& core,
                    const diy::ContinuousBounds& bounds,
                    const diy::ContinuousBounds& domain,
                    const diy::RegularContinuousLink& link) const
    {
    	ParticleBlock* b = AddBlock::operator()(gid, core, bounds, domain, link);
    	if(b->lid == 0)
    	{
      	    size_t sizes[3], n;
      	    sizes[0] = (size_t)(domain.max[0] - domain.min[0]);
      	    sizes[1] = (size_t)(domain.max[1] - domain.min[1]);
      	    sizes[2] = (size_t)(domain.max[2] - domain.min[2]);

      	    if(AddBlock::master.communicator().rank() < AddBlock::master.communicator().size()-1) 
		n = (sizes[0] * sizes[1] * sizes[2]) / AddBlock::master.communicator().size();
      	    else
		n = (sizes[0] * sizes[1] * sizes[2]) - (sizes[0] * sizes[1] * sizes[2]) / AddBlock::master.communicator().size() * (AddBlock::master.communicator().size()-1);

            b->num_particles = gen_particles_uniform(domain, n, b->particles); // generate uniform random particles in global domain;
            b->num_orig_particles = b->num_particles;
    	}
    }
};

struct AddAndReadHDF5: public AddBlock
{
    AddAndReadHDF5(diy::Master&                     master_,
                   const std::string                infile_,
                   int	                             nblocks_,
                   const std::vector<std::string>&  coordinates_):
                   AddBlock(master_),
                   infile(infile_),
                   nblocks(nblocks_),
                   coordinates(coordinates_) {}

    void  operator()(int gid,
                     const diy::ContinuousBounds& core,
                     const diy::ContinuousBounds& bounds,
                     const diy::ContinuousBounds& domain,
                     const diy::RegularContinuousLink& link) const
    {
	ParticleBlock* b = AddBlock::operator()(gid, core, bounds, domain, link);
    	int rank, size;
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    	MPI_Comm_size(MPI_COMM_WORLD, &size);

    	std::vector<std::vector<float> > particles;
    	if(b->lid == 0)
      	    read_particles_32f(MPI_COMM_WORLD,
                              infile,
                              rank,
                              size,
                              particles,
                              coordinates);

    	b->num_particles      = particles.size();
    	b->num_orig_particles = particles.size();
    	b->particles.resize(particles.size());
    	for(size_t i = 0; i < particles.size(); ++i)
      	    for(unsigned j = 0; j < 3; ++j)
        	b->particles[i][j] = particles[i][j];
    }

    int                              nblocks;
    const std::string		      infile;
    const std::vector<std::string>&  coordinates;
};


