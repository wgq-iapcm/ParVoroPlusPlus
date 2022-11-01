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

#include <diy/algorithms.hpp>
#include <diy/reduce.hpp>
#include <diy/partners/swap.hpp>
#include <diy/pick.hpp>
#include "particle.hpp"
#include "regular.hpp"
#include "common.hpp"

typedef std::vector<diy::RegularContinuousLink> LinkVector;
typedef std::vector<size_t> LastNeighbors;

void 
kdtree_redistribute(diy::Master& master,
	            const diy::Assigner& assigner,
		    bool wrap)
{
    if(assigner.nblocks() == 1)
        return;

    int bins = 1024;     
    diy::ContinuousBounds domain = master.block<ParticleBlock>(master.loaded_block())->global_bounds;
    diy::kdtree(master, assigner, 3, domain, &ParticleBlock::particles, bins, wrap);

    for(int i = 0; i < master.size(); i++)
    {
        auto b = (ParticleBlock*)master.block(i);
        b->num_orig_particles = b->num_particles = b->particles.size();
        
        auto link = static_cast<diy::RegularContinuousLink*>(master.link(i));
        b->box    = link->bounds();
        b->core   = link->core();
        b->bounds = link->bounds();        
    }
}

size_t 
kdtree_send_ghost(ParticleBlock* b, 
                  const diy::Master::ProxyWithLink& cp,
                  size_t last_neighbor,
                  float ghostsize)
{
    auto l = dynamic_cast<diy::RegularContinuousLink*>(cp.link());

    if(last_neighbor == l->size())
        return 0;
                 
    std::vector<std::set< std::pair<diy::BlockID, diy::Direction> > > to_send;
    to_send.resize(b->num_orig_particles);
    for(size_t i = last_neighbor; i < l->size(); i++)
    {
        diy::ContinuousBounds neigh_bounds = l->bounds(i);
        diy::wrap_bounds(neigh_bounds, l->wrap(i), b->global_bounds);
        for(size_t j = 0; j < b->num_orig_particles; j++)  
        {
            diy::Point<double, 3> p;
            p[0] = b->particles[j][0];
            p[1] = b->particles[j][1];
            p[2] = b->particles[j][2];

            if(diy::distance(neigh_bounds, p) < ghostsize)
                to_send[j].insert(std::make_pair(l->target(i), l->wrap(i)));
        }
    }  

    // enqueue the particles
    size_t enqueued = 0;
    RemoteParticle rp;
    for(size_t p = 0; p < b->num_orig_particles; p++)
    {
	for(auto it = to_send[p].begin(); it != to_send[p].end(); it++)
        {
            rp.particle = b->particles[p];
            rp.info.bid = b->bid;
            rp.info.index = p;
            wrap_pt(rp.particle, it->second, b->global_bounds);
            cp.enqueue(it->first, rp);	
            ++enqueued;
        }
    }

    return enqueued;
}

void 
adaptive_nbrs(ParticleBlock*    b,
              const diy::Master::ProxyWithLink& cp,
              const LinkVector&                 links,
              LastNeighbors&                    neighbors,
              bool                              first,
              float ghostsize)
{
    int lid = cp.master()->lid(cp.gid());
    const diy::RegularContinuousLink& original_link = links[lid];
    size_t& last_neighbor = neighbors[lid];
    auto link  = dynamic_cast<diy::RegularContinuousLink*>(cp.link());

    cp.collectives()->clear();

    LinkVector in_links;
    if(!first)     
    {
        for(size_t i = last_neighbor; i < link->size(); ++i)
        {
            diy::MemoryBuffer& in = cp.incoming(link->target(i).gid);

            auto l = dynamic_cast<diy::RegularContinuousLink*>(diy::LinkFactory::load(in));
            in_links.push_back(*l);
            delete l;
        }
        size_t last_last_neighbor = last_neighbor;
        last_neighbor = link->size(); 

        recv_particles(b, cp);

        std::set< std::pair<diy::BlockID, diy::Direction> > neighbor_blocks;
        for (size_t i = 0; i < link->size(); ++i)
            neighbor_blocks.insert(std::make_pair(link->target(i), link->wrap(i)));
        size_t original_size_unique = link->size_unique();

        for(size_t i = 0; i < in_links.size(); ++i)
        {
            size_t ii = i + last_last_neighbor;
            for(size_t j = 0; j < in_links[i].size(); ++j)
            {
                if(in_links[i].target(j).gid == cp.gid()) continue;        

                diy::Direction wrap = link->wrap(ii);
                for(unsigned k = 0; k < 3; ++k)
                {
                    wrap[k] += in_links[i].wrap(j)[k];
                    if(wrap[k] < -1 || wrap[k] > 1)
                    {
                        fprintf(stderr, "Warning: something is odd with the wrap, "
                                "it exceeds a single wrap-around\n");
                    }
                }

                bool inserted = neighbor_blocks.insert(std::make_pair(in_links[i].target(j),wrap)).second;
                if(!inserted) continue;

                link->add_neighbor(in_links[i].target(j));
                link->add_direction(in_links[i].direction(j));
                link->add_bounds(in_links[i].bounds(j));
                link->add_wrap(wrap);
            }
        }
        cp.master()->add_expected(link->size_unique() - original_size_unique);
    }

    for(size_t i = last_neighbor; i < link->size(); ++i)
    {
        diy::MemoryBuffer& out = cp.outgoing(link->target(i));
        diy::LinkFactory::save(out, &original_link);
    }

    int done = 1;
    if(b->num_orig_particles)
    {        
        size_t num = kdtree_send_ghost(b, cp, last_neighbor, ghostsize);
        done = (num == 0);
    }
    cp.all_reduce(done, std::logical_and<int>());
}

void 
kdtree_exchange(diy::Master& master,
		const diy::Assigner& assigner,
		bool wrap,
		float ghostsize)
{
    if(ghostsize < EPSILON)
        return;

    LinkVector original_links;
    for(size_t i = 0; i < master.size(); ++i)
        original_links.push_back(*dynamic_cast<diy::RegularContinuousLink*>(master.link(i)));

    LastNeighbors last_neighbors(master.size(), 0); 
    bool first = true;
    int done = false;

    while(!done)
    {
        master.foreach([&](ParticleBlock* b, const diy::Master::ProxyWithLink& cp)
                       { adaptive_nbrs(b, cp, original_links, last_neighbors, first, ghostsize); });
        master.exchange();

        first = false;
        done = master.proxy(master.loaded_block()).read<int>();
    }

    for(unsigned i = 0; i < master.size(); ++i)
        master.replace_link(i, new diy::RegularContinuousLink(original_links[i]));
}


