#!/bin/bash

#PBS -N voronoi
#PBS -j oe
#PBS -l nodes=1:ppn=24
#PBS -l walltime=10:00:00
#PBS -q batch

#----------------------------------------------------------------------------
#
# mpi run script
#
# Wu Guoqing
# Institute of Applied Physics and Computational Mathematics, Beijing
#
# wu_guoqing[at]iapcm.ac.cn
#
#----------------------------------------------------------------------------

ARCH=LINUX
#ARCH=PLUTO

# number of procs
num_procs=4

if [ "$ARCH" == "PLUTO" ]; then
    cd $PBS_O_WORKDIR
fi

# executable
exe=./KdGhost

# input file, "!" indicates generating random particles in domain [xlo,xhi]*[ylo,yhi]*[zlo,zhi]
#infile="!"
infile="./unit-cube.h5"

# output file, "!" indicates no output file
#outfile="!"
outfile="results"

# variables with the delimitation character ',' (no spacing)
vars="x,y,z"

# possible options are --wrap --debug --test --kdtree --blocks <totblocks> --threads <num_threads> --cutoffRadius <cutoffRadius> --files <num_files> --nx <nx> --ny <ny> --nz <nz> --xlo <domain.min[0]> --ylo <domain.min[1]> --zlo <domain.min[2]> --xhi <domain.max[0]> --yhi <domain.max[1]> --zhi <domain.max[2]>
opts="--debug --kdtree --blocks 64 --threads 1 --cutoffRadius 0.05 --xlo 0 --ylo 0 --zlo 0 --xhi 1 --yhi 1 --zhi 1" 

#------
#
# program arguments
#  

args="$opts $infile $outfile $vars"

#------
#
# run commands
#

if [ "$ARCH" = "LINUX" ]; then

    mpiexec -n $num_procs $exe $args
    #mpiexec -n $num_procs valgrind -q $exe $args
    #mpiexec -n $num_procs valgrind -q --tool=memcheck --leak-check=yes $exe $args

fi

if [ "$ARCH" = "PLUTO" ]; then

    NP=`cat $PBS_NODEFILE | wc -l`
    source /public/software/profile.d/compiler_intel-composer_xe_2017.0.098.sh
    source /public/software/profile.d/mpi_intelmpi-2017.sh
    #mpirun -np $NP -machinefile $PBS_NODEFILE $exe $args
    mpiexec -n $num_procs -machinefile $PBS_NODEFILE $exe $args

fi


