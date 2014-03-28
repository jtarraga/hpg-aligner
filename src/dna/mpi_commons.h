#ifdef _MPI
#ifndef _MPI_COMMONS_H
#define _MPI_COMMON_H

#include "mpi.h"

//--------------------------------------------------------------------

#define GET_READ_OFFSET     1
#define SET_READ_OFFSET     2
#define UPDATE_READ_OFFSET  3
#define READ_END_OF_FILE    4
#define GET_WRITE_OFFSET    5
#define SET_WRITE_OFFSET    6
#define UPDATE_WRITE_OFFSET 7
#define WRITE_END_OF_FILE   8

//--------------------------------------------------------------------

void mpi_master(int np);


//--------------------------------------------------------------------
//--------------------------------------------------------------------

#endif // _MPI_COMMONS_H
#endif // _MPI
