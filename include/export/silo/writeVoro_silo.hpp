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

#include <pmpio.h>
#include <silo.h>
#include <voro++/voro++.hh>
#include <diy/mpi.hpp>
#include "StandardIncludes.h"
#include "SiloUtils.hpp"

using namespace std;
using namespace voro;

void 
writeVoro(diy::Master& master,
          diy::StaticAssigner& assigner,
          const string& directory = "results_silo", 
          int  timestep = 0,
          double time = 0,
          int numFiles = -1,
          int mpiTag = 0);

static void 
write_block_ucd_lines(DBfile* silofile, 
                      ParticleBlock* b, 			  
                      int timestep, 
                      double time)
{
    DBSetAllowEmptyObjects(1);

    int nproc = 1, rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    char blockName[1024];
    snprintf(blockName, 1024, "block_%d", b->lid);
    DBMkDir(silofile, blockName);
    DBSetDir(silofile, blockName);

    // add timestep/time metadata if needed.
    DBoptlist* optlist = DBMakeOptlist(10);
    double dtime = static_cast<double>(time);
    DBAddOption(optlist, DBOPT_CYCLE, &timestep);
    DBAddOption(optlist, DBOPT_DTIME, &dtime);

    vector<double> vx, vy, vz;
    vector<double> sx,sy,sz; //site coords

    voro::c_loop_order clo(*b->con, b->po);
    voro::voronoicell cell;
    if(clo.start()) do if((b->con->compute_cell(cell, clo)))
    {
	double x,y,z;
	clo.pos(x,y,z);
        sx.push_back(x);
        sy.push_back(y);
        sz.push_back(z);

        int i,j,k;double *ptsp=cell.pts,*pt2;
        for(i=0;i<cell.p;i++,ptsp+=3) 
        {
            for(j=0;j<cell.nu[i];j++) 
            {
                k=cell.ed[i][j];
                if(k<i) 
                {
                    pt2=cell.pts+3*k;
                    vx.push_back(x+*ptsp*0.5);vx.push_back(x+*pt2*0.5);
                    vy.push_back(y+ptsp[1]*0.5);vy.push_back(y+0.5*pt2[1]);
                    vz.push_back(z+ptsp[2]*0.5);vz.push_back(z+0.5*pt2[2]);
                }
            }
	}
    }while(clo.inc());

    double* coords[3];
    coords[0] = &(vx[0]);
    coords[1] = &(vy[0]);
    coords[2] = &(vz[0]);

    int nzones = vx.size()/2;
    vector<int> nodelist(vx.size());
    for(int i = 0; i < vx.size()/2; i++)
    {
        nodelist[2*i] = 2*i;
        nodelist[2*i+1] = 2*i+1;
    }
    int nshapetypes = vx.size()/2;
    vector<int> shapetype(nshapetypes, DB_ZONETYPE_BEAM);
    vector<int> shapesize(nshapetypes, 2);
    vector<int> shapecounts(nshapetypes, 1);

    DBPutZonelist2(silofile, "zonelist", nzones, 3, &nodelist[0], nodelist.size(), 0, 0, 0, &shapetype[0], &shapesize[0], &shapecounts[0], nshapetypes, optlist);
    DBPutUcdmesh(silofile, (char*)"voronoi", 3, 0, coords, vx.size(), nzones, "zonelist", NULL, DB_DOUBLE, optlist); 

    double* sites[3];
    sites[0] = &(sx[0]);
    sites[1] = &(sy[0]);
    sites[2] = &(sz[0]);
    DBPutPointmesh(silofile, "site", 3, sites, sx.size(), DB_DOUBLE, optlist);

    float extent_x[] = {b->bounds.min[0], b->bounds.max[0]};
    float extent_y[] = {b->bounds.min[1], b->bounds.max[1]};
    float extent_z[] = {b->bounds.min[2], b->bounds.max[2]};
    int extent_dims[] = {2, 2, 2};
    float *extent[] = {extent_x, extent_y, extent_z};
    DBPutQuadmesh(silofile, "extent", NULL, extent, extent_dims, 3, DB_FLOAT, DB_COLLINEAR, optlist);

    int ranks[1] = {rank};
    int dims[1] = {1};
    DBPutQuadvar1(silofile, "rank", "extent", ranks, dims, 1, NULL, 0, DB_INT, DB_ZONECENT, optlist); 

    DBSetDir(silofile, "..");
    DBFreeOptlist(optlist);
}

static void 
writeMultiXXXObjects(DBfile* silofile, PMPIO_baton_t* baton, int timestep, double time, vector<vector<int> >& procInGroup, vector<int>& nblockInRank)
{
    int nproc = 1, rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int groupRank = PMPIO_GroupRank(baton, rank);
    int rankInGroup = PMPIO_RankInGroup(baton, rank);
    if(rankInGroup == 0)
    {
        vector<char*> VoronoiNames;
        vector<char*> ExtentNames;
        vector<char*> SiteNames;
        vector<char*> RankNames;

        int count = 0;
        int numChunks = procInGroup[groupRank].size();
        for(int i = 0; i < numChunks; ++i)
        {
            int numBlocks = nblockInRank[procInGroup[groupRank][i]];
            for(int j = 0; j < numBlocks; ++j)
	    {
                count++;
                // Voronoi Mesh
	        char Name1[1024];
	        snprintf(Name1, 1024, "domain_%d/block_%d/voronoi", i,j);
	        VoronoiNames.push_back(strdup(Name1));

                // Local spatial extent
	        char Name2[1024];
	        snprintf(Name2, 1024, "domain_%d/block_%d/extent", i,j);
	        ExtentNames.push_back(strdup(Name2));

                // Voronoi site
	        char Name3[1024];
	        snprintf(Name3, 1024, "domain_%d/block_%d/site", i,j);
	        SiteNames.push_back(strdup(Name3));

                // Block rank
	        char Name4[1024];
	        snprintf(Name4, 1024, "domain_%d/block_%d/rank", i,j);
	        RankNames.push_back(strdup(Name4));
            }
        }
        
        vector<int> VoronoiTypes(count, DB_UCDMESH);
        vector<int> ExtentTypes(count, DB_COLLINEAR);
        vector<int> SiteTypes(count, DB_POINTMESH);
        vector<int> RankTypes(count, DB_QUADVAR);
        // Write the mesh and variable data.
        DBSetDir(silofile, "/");
        DBPutMultimesh(silofile, "voronoi", count, &VoronoiNames[0], &VoronoiTypes[0], NULL);
        DBPutMultimesh(silofile, "extent", count, &ExtentNames[0], &ExtentTypes[0], NULL);
        DBPutMultimesh(silofile, "site", count, &SiteNames[0], &SiteTypes[0], NULL);
        DBPutMultivar(silofile, "rank", count, &RankNames[0], &RankTypes[0], NULL);

        for(int i = 0; i < VoronoiNames.size(); ++i)
        {
            free(VoronoiNames[i]);
            free(ExtentNames[i]);
            free(SiteNames[i]);
            free(RankNames[i]);
        }
    }

    PMPIO_HandOffBaton(baton, (void*)silofile);
    PMPIO_Finish(baton); 
}

static void
write_master(string masterDirName, int timestep, double time, vector<vector<int> >& procInGroup, vector<int>& nblockInRank)
{
    int nproc = 1, rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0)
    {
        char masterFileName[1024];
        snprintf(masterFileName, 1024, "%s/Timestep.%05d/summary.root", masterDirName.c_str(), timestep);

        int driver = DB_HDF5;
        DBfile* silofile = DBCreate(masterFileName, DB_CLOBBER, DB_LOCAL, "Master file", driver);

        vector<char*> VoronoiNames;
        vector<char*> ExtentNames;
        vector<char*> SiteNames;
        vector<char*> RankNames;

        int count = 0;
        int numFiles = procInGroup.size();
        for(int i = 0; i < numFiles; ++i)
        {
            int numChunks = procInGroup[i].size();
            for(int c = 0; c < numChunks; ++c)
            {
                int numBlocks = nblockInRank[procInGroup[i][c]];
                for(int b = 0; b < numBlocks; ++b)
                {
                    count++;
                    // Mesh.
                    char Name1[1024];
                    snprintf(Name1, 1024, "data_T%05d_G%05d.silo:/domain_%d/block_%d/voronoi", timestep,i,c,b);
                    VoronoiNames.push_back(strdup(Name1));

                    // Local spatial extent
	            char Name2[1024];
	            snprintf(Name2, 1024, "data_T%05d_G%05d.silo:/domain_%d/block_%d/extent", timestep, i, c, b);
	            ExtentNames.push_back(strdup(Name2));

                    // Voronoi site
	            char Name3[1024];
	            snprintf(Name3, 1024, "data_T%05d_G%05d.silo:/domain_%d/block_%d/site", timestep, i, c, b);
	            SiteNames.push_back(strdup(Name3));

                    // Block rank
	            char Name4[1024];
	            snprintf(Name4, 1024, "data_T%05d_G%05d.silo:/domain_%d/block_%d/rank", timestep, i, c, b);
	            RankNames.push_back(strdup(Name4));
                }
            }
        }

        vector<int> VoronoiTypes(count, DB_UCDMESH);
        vector<int> ExtentTypes(count, DB_COLLINEAR);
        vector<int> SiteTypes(count, DB_POINTMESH);
        vector<int> RankTypes(count, DB_QUADVAR);

        DBPutMultimesh(silofile, "voronoi", count, &VoronoiNames[0], &VoronoiTypes[0], NULL);
        DBPutMultimesh(silofile, "extent", count, &ExtentNames[0], &ExtentTypes[0], NULL);
        DBPutMultimesh(silofile, "site", count, &SiteNames[0], &SiteTypes[0], NULL);
        DBPutMultivar(silofile, "rank", count, &RankNames[0], &RankTypes[0], NULL);

        DBClose(silofile);
        for(int i = 0; i < VoronoiNames.size(); ++i)
        {
            free(VoronoiNames[i]);
            free(ExtentNames[i]);
            free(SiteNames[i]);
            free(RankNames[i]);
        }
    }
}

void 
writeVoro(diy::Master& master,
          diy::StaticAssigner& assigner,
          const string& directory, 
          int timestep,
          double time,
          int numFiles,
          int mpiTag)
{
    DBfile* silofile;
    PMPIO_baton_t* baton;
    vector<vector<int> > procInGroup;
    vector<int> nblockInRank;

    init_silo(assigner, directory, timestep, time, numFiles, mpiTag, silofile, baton, procInGroup, nblockInRank);

    for(int i = 0; i < master.size(); i++)
    {
        ParticleBlock* b = (ParticleBlock*)master.block(i);
        write_block_ucd_lines(silofile, b, timestep, time);
    }

    writeMultiXXXObjects(silofile, baton, timestep, time, procInGroup, nblockInRank);
    write_master(directory, timestep, time, procInGroup, nblockInRank);
    write_dumps(directory, timestep);
}



