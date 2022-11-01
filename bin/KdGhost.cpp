///////////////////////////////////////////////////////////////////////////////
//
// ParVoro++: A scalable parallel algorithm for constructing 3D Voronoi tessellations based on KD-tree decomposition
//
// Guoqing Wu
// High Performance Computing Center, Institute of Applied Physics and Computational Mathematics, Beijing
// wu_guoqing@iapcm.ac.cn
//
///////////////////////////////////////////////////////////////////////////////

#include <diy/decomposition.hpp>
#include <diy/assigner.hpp>
#include <diy/master.hpp>
#include <diy/mpi.hpp>
#include <diy/reduce-operations.hpp>
#include "StandardIncludes.h"
#include "data/particle.hpp"
#include "data/common.hpp"
#include "data/regular.hpp"
#include "data/kdtree.hpp"
#include "voronoi/VolumeStat.hpp"
#include "voronoi/VoronoiGhost.hpp"
#include "export/silo/writeVoro_silo.hpp"
#include "utility/opts.h"
#include "utility/timer.h"

typedef ParticleBlock               Block;
typedef diy::ContinuousBounds       Bounds;

int main(int argc, char* argv[])
{
  diy::mpi::environment     env(argc, argv);
  diy::mpi::communicator    world;

  using namespace opts;
  Options ops;

  int rank = world.rank();      // MPI usual
  int size = world.size();      // MPI usual
  int tot_blocks = size;        // total number of blocks in the domain
  int num_threads = 1;          // number of threads diy can use
  int mem_blocks = -1;          // number of blocks to keep in memory
  int num_files = -1;           // number of output files
  float cutoffRadius = 0.0;     // distance between particles
  std::string format;           // hdf5 or random particles
  std::string infile;           // input file name
  std::string outfile;          // output file name
  std::vector<std::string>  coordinates; // coordinates to read
  int nx=0, ny=0, nz=0;         // param[in] of voro::container, the number of grid blocks in each of the three coordinate directions  
  bool kdtree, wrap, debug, test, help;

  ops
      >> Option('b', "blocks",       tot_blocks,   "total number of blocks to use")
      >> Option('t', "threads",      num_threads,  "number of threads to use")
      >> Option('n', "num_files",    num_files,    "number of files for parallel i/o")
      >> Option('c', "cutoffRadius", cutoffRadius, "width of ghost layer")
      >> Option('I', "nx",           nx,           "param[in] of voro::container, the number of grid blocks of voro::container in x direction; in default, it is guessed by voro++")
      >> Option('J', "ny",           ny,           "param[in] of voro::container, the number of grid blocks of voro::container in y direction; in default, it is guessed by voro++")
      >> Option('K', "nz",           nz,           "param[in] of voro::container, the number of grid blocks of voro::container in z direction; in default, it is guessed by voro++")  
  ;

  Bounds domain(3);
  domain.min[0] = domain.min[1] = domain.min[2] = 0;
  domain.max[0] = domain.max[1] = domain.max[2] = 20.;
  ops
      >> Option('x', "xlo",  domain.min[0],  "domain min x")
      >> Option('y', "ylo",  domain.min[1],  "domain min y")
      >> Option('z', "zlo",  domain.min[2],  "domain min z")
      >> Option('X', "xhi",  domain.max[0],  "domain max x")
      >> Option('Y', "yhi",  domain.max[1],  "domain max y")
      >> Option('Z', "zhi",  domain.max[2],  "domain max z")
  ;


  ops
      >> Option("kdtree",    kdtree,  "use kdtree decomposition")
      >> Option("wrap",      wrap,    "use periodic boundary")
      >> Option("debug",     debug,   "print debugging info")
      >> Option("test",      test,    "measuring the pure computation time of the Voronoi algorithm")      
      >> Option('h', "help", help,    "show help")
  ;

  std::string varnames;
  if( !ops.parse(argc,argv) || help || !(ops >> PosOption(infile) >> PosOption(outfile) >> PosOption(varnames)) )
  {
      if (world.rank() == 0)
      {
          std::cout << "Usage: " << argv[0] << " [OPTIONS] infile outfile coordinates\n";
          std::cout << "ParVoro++: A Scalable Parallel toolkit for Constructing 3D Voronoi Tessellations.\n";
          std::cout << "Author: Guoqing Wu (wu_guoqing@iapcm.ac.cn)\n";
          std::cout << ops;
      }
      return 0;
  }
  // parse variables
  std::istringstream ss(varnames);
  std::string word;
  while(std::getline(ss, word, ',')) { coordinates.push_back(word); };

  if(kdtree)
  {
    if(wrap && tot_blocks < 64)
    {
       if(rank == 0)
           fprintf(stderr, "Warning: using k-d tree with wrap on and fewer than 64 blocks is likely to fail\n");
       return 1;
    }
  }

  assert(isPowerOf2(tot_blocks) && isPowerOf2(size));
  assert(tot_blocks = size * num_threads);
    
  diy::Master  master(world, num_threads, mem_blocks, &Block::create_block, &Block::destroy_block);
  diy::RoundRobinAssigner  assigner(world.size(), tot_blocks);

  // decompose
  diy::RegularDecomposer<Bounds>::BoolVector        wraps;
  diy::RegularDecomposer<Bounds>::BoolVector        share_face;
  diy::RegularDecomposer<Bounds>::CoordinateVector  ghosts;
  ghosts.assign(3, cutoffRadius); 
  if(wrap)
      wraps.assign(3, true);

  double times[5]; // timing
  timing(times, 0, -1, MPI_COMM_WORLD);
  timing(times, 1, -1, MPI_COMM_WORLD);
  
  if(infile == "!") // generate random particles
  {
      AddAndGenerate genblock(master);
      diy::RegularDecomposer<Bounds> decomposer(3, domain, tot_blocks, share_face, wraps, ghosts);
      decomposer.decompose(world.rank(), assigner, genblock);
  }
  else // read HDF5 file
  {
      AddAndReadHDF5 creat_and_read(master, infile, size, coordinates);
      diy::RegularDecomposer<Bounds> decomposer(3, domain, tot_blocks, share_face, wraps, ghosts);
      decomposer.decompose(world.rank(), assigner, creat_and_read);
  }

  timing(times, -1, 1, MPI_COMM_WORLD);
  timing(times, 2, -1, MPI_COMM_WORLD);

  size_t nblocks = master.size();
  size_t tot_particles = get_total_num(world, master);

  // sort and distribute particles to all blocks
  if(kdtree)
      kdtree_redistribute(master, assigner, wrap);
  else
      regular_redistribute(master, assigner);
  if(rank == 0)
      fprintf(stderr,"Particle redistribution completed.\n");
  
  //neighbor ghost exchange
  if(kdtree)
      kdtree_exchange(master, assigner, wrap, cutoffRadius);
  else
      regular_exchange_ghost(master, cutoffRadius);
    
  world.barrier();    
  if(rank ==0)
      fprintf(stderr,"Ghost particle exchange completed.\n");

  timing(times, -1, 2, MPI_COMM_WORLD);
  timing(times, 3, -1, MPI_COMM_WORLD);

  //voronoi tessellation
  Voronoi_ghost(master, wrap, nx, ny, nz, test);

  world.barrier();
  if(rank == 0)
      fprintf(stderr,"Voronoi tessellation completed.\n");

  timing(times, -1, 3, MPI_COMM_WORLD);
  timing(times, 4, -1, MPI_COMM_WORLD);

  //output
  if(outfile != "!")
  {
      writeVoro(master, assigner, outfile, 0, 0, num_files, 0);
      if(rank == 0)
          fprintf(stderr,"Output completed.\n");
  }

  timing(times, -1, 4, MPI_COMM_WORLD);
  timing(times, -1, 0, MPI_COMM_WORLD);

  float lvol, gvol = 0, dvol = 0;
  if(debug)
  {
      lvol = volume_stat(master);
      world.barrier();
      diy::mpi::all_reduce(master.communicator(), lvol, gvol, std::plus<float>());
      dvol = (domain.max[0]-domain.min[0])*(domain.max[1]-domain.min[1])*(domain.max[2]-domain.min[2]);
  }

  if(rank == 0) 
  {
      fprintf(stderr,"//=================== Input Parameter ==================//\n");
      fprintf(stderr," infile: %s\n outfile: %s\n ranksize: %d\n threads: %d\n tot_blocks: %d\n cutoff: %f\n wrap: %d\n tot_particles: %ld\n",infile.c_str(), outfile.c_str(), world.size(), num_threads, tot_blocks, cutoffRadius, wrap, tot_particles);
      fprintf(stderr,"//================== Analysis Results ==================//\n");
      if(debug)
          fprintf(stderr," Voronoi Volume: %lf, Domain Volume: %lf, ErrorRatio: %lf\n",gvol,dvol,abs(gvol-dvol)/dvol);
      fprintf(stderr, " total time: %.3lf s = %.3lf s input + %.3lf s decomposition + %.3lf s tessellation + %.3lf s output\n", times[0], times[1], times[2], times[3], times[4]);
  }

  return 0;
}



