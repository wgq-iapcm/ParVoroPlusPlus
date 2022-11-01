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

#include <voro++/voro++.hh>
#include "StandardIncludes.h"
#include "data/particle.hpp"
#include "data/common.hpp"
#include "utility/timer.h"

typedef std::vector<diy::RegularContinuousLink> LinkVector;
typedef std::vector<size_t> LastNeighbors;

void 
local_cells(ParticleBlock* b,
            const diy::Master::ProxyWithLink& cp,
            bool wrap,
            int nx,
            int ny,
            int nz)
{
    int nx_, ny_, nz_;
    if(nx == 0 || ny == 0 || nz == 0)
    {
        double lamda = pow(b->num_orig_particles / (4.6*(b->global_bounds.max[0]-b->global_bounds.min[0])*(b->global_bounds.max[1]-b->global_bounds.min[1])*(b->global_bounds.max[2]-b->global_bounds.min[2])), 0.3333333); 
        nx_ = ceil(lamda * (b->global_bounds.max[0]-b->global_bounds.min[0]));
        ny_ = ceil(lamda * (b->global_bounds.max[1]-b->global_bounds.min[1]));
        nz_ = ceil(lamda * (b->global_bounds.max[2]-b->global_bounds.min[2]));         
    }
    else
    {
        nx_ = nx;
        ny_ = ny;
        nz_ = nz;
    }
        
    while(nx_ * ny_ * nz_ > 8000000) //TODO: large (nx_,ny_,nz_) may cause std::bad_alloc exception
    {
        nx_ = int(nx_ * 0.9);
        ny_ = int(ny_ * 0.9);
        nz_ = int(nz_ * 0.9);
    }
           
    b->con = new voro::container(b->global_bounds.min[0], b->global_bounds.max[0], b->global_bounds.min[1], b->global_bounds.max[1], b->global_bounds.min[2], b->global_bounds.max[2], nx_, ny_, nz_, wrap, wrap, wrap, 8);

    // add walls here as you wish
    {
        //b->sphere = new voro::wall_sphere(0.5, 0.5, 0.5, 0.3);
        //b->con->add_wall(*b->sphere);
    }
        
    for(int j = 0; j < b->num_orig_particles; j++)
        b->con->put(b->po, j,  b->particles[j][0],  b->particles[j][1],  b->particles[j][2]);  
}

void 
update_cells(ParticleBlock* b)
{
    for(int j = b->con->total_particles(); j < b->num_particles; j++)
        b->con->put(j, b->particles[j][0], b->particles[j][1], b->particles[j][2]);
}

int 
incomplete_cells(ParticleBlock* b,
                 const diy::Master::ProxyWithLink& cp,
                 size_t last_neighbor,
                 bool first)
{
    if(!first && b->enqueued == 0)
        return 0;
        
    auto l = dynamic_cast<diy::RegularContinuousLink*>(cp.link());
    std::vector<std::set<int> > to_send( b->num_orig_particles);
    
    if(first)
    {   
        voro::voronoicell c;
        double x, y, z, r;
        voro::c_loop_order clo(*b->con, b->po);
        if(clo.start()) do if(b->con->compute_cell(c, clo))
        {
            int pid = clo.pid();
            double center[3] = {static_cast<double>(b->particles[pid][0]), static_cast<double>(b->particles[pid][1]), static_cast<double>(b->particles[pid][2])}; 
            double max_rad = sqrt(c.max_radius_squared());        
            int j;
            for(j = 0; j < 3; ++j)
            {
                if(center[j] - b->core.min[j] <= max_rad) break;
                if(b->core.max[j] - center[j] <= max_rad) break;
            }
            if(j == 3)	
                continue;
                            
            for(int i = last_neighbor; i < l->size(); ++i)
            {
                diy::ContinuousBounds neigh_bounds = l->bounds(i);
                diy::wrap_bounds(neigh_bounds, l->wrap(i), b->global_bounds);
               
                if(diy::distance(neigh_bounds, diy::Point<double,3> { center }) <= max_rad)
                    to_send[pid].insert(i);
            }

            OutCell oc;
            oc.ijk = clo.ijk;
            oc.q = clo.q;
            oc.pid = pid;
            b->outcell.push_back(oc);        
        }while(clo.inc());
    }
    else
    {
        for(auto it = b->outcell.begin(); it != b->outcell.end(); it++)
        {
            voro::voronoicell c;
            b->con->compute_cell(c, it->ijk, it->q);
             
            double center[3] = {static_cast<double>(b->particles[it->pid][0]), static_cast<double>(b->particles[it->pid][1]), static_cast<double>(b->particles[it->pid][2])}; 
            double max_rad = sqrt(c.max_radius_squared());
            bool complete = true;
            for(int i = last_neighbor; i < l->size(); ++i)
            {
                diy::ContinuousBounds neigh_bounds = l->bounds(i);
                diy::wrap_bounds(neigh_bounds, l->wrap(i), b->global_bounds);

                if(diy::distance(neigh_bounds, diy::Point<double,3> { center }) <= max_rad)
                {
                    to_send[it->pid].insert(i);
                    complete = false;
                }
            }     
         
            if(complete)
                b->outcell.erase(it--);
        }
    }

    size_t enqueued = 0;
    RemoteParticle rp;
    for(int p = 0; p < b->num_orig_particles; p++)
    {
        for(std::set<int>::iterator it = to_send[p].begin(); it != to_send[p].end(); it++)
        {         
            rp.particle =  b->particles[p];
            rp.info.bid = b->bid;
            rp.info.index = p;
            wrap_pt(rp.particle, l->wrap(*it), b->global_bounds);
            cp.enqueue(l->target(*it), rp);
            ++enqueued;
        }
    }    
    b->enqueued = enqueued;

    return enqueued;
}


void 
voronoi(ParticleBlock*                    b,
        const diy::Master::ProxyWithLink& cp,
        const LinkVector&                 links,
        LastNeighbors&                    neighbors,
        bool                              first,
        bool wrap)
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

    if(b->num_orig_particles)
        update_cells(b);

    for(size_t i = last_neighbor; i < link->size(); ++i)
    {
        diy::MemoryBuffer& out = cp.outgoing(link->target(i));
        diy::LinkFactory::save(out, &original_link);
    }

    int done = 1;
    if(b->num_orig_particles)
    {        
        size_t num = incomplete_cells(b, cp, last_neighbor, first);
        done = (num == 0);
    }
    cp.all_reduce(done, std::logical_and<int>()); 
}

void 
Voronoi_kd(diy::Master& master,                               
           bool wrap,
           int nx,
           int ny,
           int nz)
{ 
    master.foreach([&](ParticleBlock* b, const diy::Master::ProxyWithLink & cp)
                  {local_cells(b, cp, wrap, nx, ny, nz);});
    
    LinkVector original_links;
    for(size_t i = 0; i < master.size(); ++i)
        original_links.push_back(*dynamic_cast<diy::RegularContinuousLink*>(master.link(i)));

    LastNeighbors last_neighbors(master.size(), 0);    
    bool first    = true;
    int done      = false; 
    while(!done)
    {           
        master.foreach([&](ParticleBlock* b, const diy::Master::ProxyWithLink& cp)
                       { voronoi(b, cp, original_links, last_neighbors, first, wrap); });
        master.exchange();
        first = false;
        done = master.proxy(master.loaded_block()).read<int>();             
    }

    for(unsigned i = 0; i < master.size(); ++i)
        master.replace_link(i, new diy::RegularContinuousLink(original_links[i]));       
}



