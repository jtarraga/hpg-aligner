#include "mp_main.h"

void mp_main(options_t *options) {

  int id = options->id;
  int np = options->np;

  char hostname[1024];
  gethostname(hostname, 1024);

  struct timeval stop, start;
  printf("Starting mapping (process %i (%i) of %i on %s)...\n", id, getpid(), np, hostname);
  gettimeofday(&start, NULL);
  mp_mapper(id, np, options);
  gettimeofday(&stop, NULL);
  printf("...end of mapping (process %i (%i) of %i on %s) in %0.2f s. Done!!\n\n", id, getpid(), np, hostname,
	 (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
