#include "mpi_main.h"

void mpi_main(options_t *options) {

  int provided;
  MPI_Init_thread(NULL, NULL, MPI_THREAD_SERIALIZED, &provided);
  //  MPI_Init(NULL, NULL);

  int id, np, namelen;
  char name[MPI_MAX_PROCESSOR_NAME];

  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Get_processor_name(name, &namelen);

  printf ("Hello, I'm process %i of %i, and I'm running on node %s\n", id, np, name);
  fflush(stdout);

  if (id == 0) {
    printf("I'm the master (id %i)\n", id);
    mpi_master(options);
  } else {
    printf("I'm a worker (id %i)\n", id);
    mpi_worker(options);
  }
 
  MPI_Finalize ();
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
