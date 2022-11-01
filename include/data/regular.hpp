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

#include "diy/reduce.hpp"
#include "diy/partners/swap.hpp"
#include "diy/pick.hpp"
#include "particle.hpp"
#include "common.hpp"

void 
redistribute(ParticleBlock* b, 
             const diy::ReduceProxy& srp, 
             const diy::RegularSwapPartners& partners)
{
    typedef diy::Point<float,3>  SimpleParticle;

    unsigned round  = srp.round();
    for(unsigned i = 0; i < srp.in_link().size(); ++i)
    {
        int nbr_gid = srp.in_link().target(i).gid;
        if(nbr_gid == srp.gid())
            continue;

        std::vector<SimpleParticle>  in_particles;
        srp.dequeue(nbr_gid, in_particles);
        size_t nparticles = in_particles.size();

        b->particles.resize(b->num_particles + nparticles);
        size_t o = b->num_particles;
        for(size_t j = 0; j < in_particles.size(); ++j)
            b->particles[o++] = in_particles[j];

        b->num_particles += nparticles;
    }
    b->num_orig_particles = b->num_particles;

    if(srp.out_link().size() == 0)  
        return;

    std::vector<std::vector<SimpleParticle> > out_particles(srp.out_link().size());
    unsigned group_size = srp.out_link().size();
    unsigned cur_dim    = partners.dim(round);

    for(size_t i = 0; i < b->particles.size(); ++i)
    {
        int loc = floor((b->particles[i][cur_dim] - b->box.min[cur_dim]) / (b->box.max[cur_dim] - b->box.min[cur_dim]) * group_size);

        if(loc == out_particles.size())
	    loc -= 1;
        if(loc < 0)
	    loc = 0;

        out_particles[loc].push_back(b->particles[i]);
    }
    int pos = -1;
    for(unsigned i = 0; i < group_size; ++i)
    {
        if(srp.out_link().target(i).gid == srp.gid())
        {
            b->particles.resize(out_particles[i].size());
            for(size_t j = 0; j < out_particles[i].size(); j++)
	        b->particles[j] = out_particles[i][j];
            b->num_particles = out_particles[i].size();
            b->num_orig_particles = b->num_particles;
            pos = i;
        }
        else
            srp.enqueue(srp.out_link().target(i), out_particles[i]);
    } 

    float new_min = b->box.min[cur_dim] + (b->box.max[cur_dim] - b->box.min[cur_dim]) / group_size * pos;
    float new_max = b->box.min[cur_dim] + (b->box.max[cur_dim] - b->box.min[cur_dim]) / group_size * (pos + 1);
    b->box.min[cur_dim] = new_min;
    b->box.max[cur_dim] = new_max;
}

void 
regular_redistribute(diy::Master& master, 
                     const diy::Assigner& assigner)
{
    diy::ContinuousBounds domain = master.block<ParticleBlock>(master.loaded_block())->global_bounds;
    diy::RegularDecomposer<diy::ContinuousBounds> decomposer(3, domain, assigner.nblocks());

    int k = 2;
    diy::RegularSwapPartners  partners(decomposer, k, false);
    diy::reduce(master, assigner, partners, &redistribute);
}

void 
regular_send_ghost(ParticleBlock* b, 
                   const diy::Master::ProxyWithLink& cp, 
                   float ghostsize)
{
    diy::RegularContinuousLink*    l = static_cast<diy::RegularContinuousLink*>(cp.link());

    assert(ghostsize < (b->core.max[0] - b->core.min[0])/2.0);
    assert(ghostsize < (b->core.max[1] - b->core.min[1])/2.0);
    assert(ghostsize < (b->core.max[2] - b->core.min[2])/2.0);

    std::vector< std::set<unsigned> > to_send(b->num_orig_particles);
    for(unsigned i = 0; i < l->size(); ++i)
    {
        diy::ContinuousBounds neigh_bounds = l->core(i);
        diy::wrap_bounds(neigh_bounds, l->wrap(i), b->global_bounds);

        for(size_t j = 0; j < b->num_orig_particles; j++)  
        {
            diy::Point<double, 3> p;
            p[0] = b->particles[j][0];
            p[1] = b->particles[j][1];
            p[2] = b->particles[j][2];
            if(diy::distance(neigh_bounds, p) < ghostsize)
                to_send[j].insert(i);
        }
    }

    size_t enqueued = 0;
    RemoteParticle rp;
    for(size_t p = 0; p < b->num_orig_particles; p++)
    {
        for(std::set<unsigned>::iterator it  = to_send[p].begin(); it != to_send[p].end(); it++)
        {
            rp.particle = b->particles[p];
            rp.info.bid = b->bid;
            rp.info.index = p;
            wrap_pt(rp.particle, l->wrap(*it), b->global_bounds);
            cp.enqueue(l->target(*it), rp);
            ++enqueued;
        }
    }
}

void 
recv_particles(ParticleBlock* b, 
               const diy::Master::ProxyWithLink& cp)
{
    std::vector<int> in;
    cp.incoming(in);

    size_t total = 0;
    for(unsigned i = 0; i < in.size(); i++)
    {
        diy::MemoryBuffer& in_queue = cp.incoming(in[i]);
        total += (in_queue.size() - in_queue.position) / sizeof(RemoteParticle);
    }
    b->particles.resize(b->num_particles + total);
    b->infos.resize(b->infos.size() + total);

    for(unsigned i = 0; i < in.size(); ++i)
    {
        diy::MemoryBuffer& incoming = cp.incoming(in[i]);
        size_t incoming_sz  = (incoming.size() - incoming.position) / sizeof(RemoteParticle);
        std::vector<RemoteParticle> pts;
        pts.resize(incoming_sz);
        diy::load(incoming, &pts[0], incoming_sz);

        for(size_t j = 0; j < incoming_sz; j++)
        {
            b->particles[b->num_particles] = pts[j].particle;
            b->infos[b->num_particles - b->num_orig_particles] = pts[j].info;
            b->num_particles++;
        }
    }
}

void 
regular_exchange_ghost(diy::Master& master, 
                       const float ghostsize)
{
    master.foreach([&](ParticleBlock* b, const diy::Master::ProxyWithLink & cp)
                  { regular_send_ghost(b, cp, ghostsize); }); 
    master.exchange();
    master.foreach(&recv_particles);
}



