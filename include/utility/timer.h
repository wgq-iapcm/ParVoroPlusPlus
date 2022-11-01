
#pragma once

#include "StandardIncludes.h"
#include <mpi.h>

//
// starts / stops timing
// (does a barrier on comm)
//
// times: timing data
// start: index of timer to start (-1 if not used)
// stop: index of timer to stop (-1 if not used)
//
void timing(double* times, int start, int stop, MPI_Comm comm)
{
    if(start < 0 && stop < 0)
    {
        for(int i = 0; i < MAX_TIMES; i++)
            times[i] = 0.0;
    }

    MPI_Barrier(comm);
    if(start >= 0)
        times[start] = MPI_Wtime();
    if(stop >= 0)
        times[stop] = MPI_Wtime() - times[stop];
}


