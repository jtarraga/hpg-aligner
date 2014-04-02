#ifdef _MPI

#include "mpi_commons.h"

//--------------------------------------------------------------------


void mpi_master(int np) {
  MPI_Status status;
  MPI_Offset read_offset = 0, write_offset = 0;
  int i, skip, num_readers = np - 1, num_writers = np - 1;

  int read_lock = 0, write_lock = 0;
  int read_lock_owner = 0, write_lock_owner = 0;

  int readers[np], writers[np];
  int first_reader = 0, last_reader = 0;
  int first_writer = 0, last_writer = 0;
  int readers_length = 0, writers_length = 0;
  
  return;

  while (1) {
    // Probe for an incoming message from any process
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    //    printf("master: probe, source = %i, tag = %i\n", status.MPI_SOURCE, status.MPI_TAG);

    if (status.MPI_TAG == GET_READ_OFFSET) {

      // receive the incoming message
      MPI_Recv(&skip, 1, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
      if (read_lock) {
	// ohhh, this reader has to wait for the offset
	if (readers_length == 0) {
	  first_reader = 0;
	  last_reader = 0;
	} else {
	  last_reader++;
	  if (last_reader >= np) last_reader = 0;
	}
	readers_length++;
	readers[last_reader] = status.MPI_SOURCE;
      } else {
	// lock and send the offset
	read_lock = 1;
	read_lock_owner = status.MPI_SOURCE;
	MPI_Send(&read_offset, sizeof(MPI_Offset), MPI_CHAR, read_lock_owner, SET_READ_OFFSET, MPI_COMM_WORLD);
      }

    } else if (status.MPI_TAG == UPDATE_READ_OFFSET) {

      // receive the updated offset
      MPI_Recv(&read_offset, sizeof(MPI_Offset), MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);

      // and then check if there are any reader waiting for the offset
      if (readers_length) {
	read_lock = 1;
	read_lock_owner = readers[first_reader];
	MPI_Send(&read_offset, sizeof(MPI_Offset), MPI_CHAR, read_lock_owner, SET_READ_OFFSET, MPI_COMM_WORLD);
	
	// update waiting readers
	first_reader++;
	if (first_reader >= np) first_reader = 0;
	readers_length--;
      } else {
	read_lock = 0;
      }
      
    } else if (status.MPI_TAG == READ_END_OF_FILE) {

      // the receive incoming message
      MPI_Recv(&skip, 1, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
      read_lock = 0;

      if (readers_length > 0) {
	for (i = 0; i < readers_length; i++) {
	  MPI_Send(&read_offset, sizeof(MPI_Offset), MPI_CHAR, readers[first_reader], 
		   READ_END_OF_FILE, MPI_COMM_WORLD);
	  first_reader++;
	  if (first_reader >= np) first_reader = 0;
	}
	readers_length = 0;
      }
      // and delete
      num_readers--;
      if (num_readers == 0 && num_writers == 0) {
	break;
      }
    } else if (status.MPI_TAG == GET_WRITE_OFFSET) {

      // receive the incoming message
      MPI_Recv(&skip, 1, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
      if (write_lock) {
	// ohhh, this writer has to wait for the offset
	if (writers_length == 0) {
	  first_writer = 0;
	  last_writer = 0;
	} else {
	  last_writer++;
	  if (last_writer >= np) last_writer = 0;
	}
	writers_length++;
	writers[last_writer] = status.MPI_SOURCE;
      } else {
	// lock and send the offset
	write_lock = 1;
	write_lock_owner = status.MPI_SOURCE;
	MPI_Send(&write_offset, sizeof(MPI_Offset), MPI_CHAR, write_lock_owner, 
		 SET_WRITE_OFFSET, MPI_COMM_WORLD);
      }
      
    } else if (status.MPI_TAG == UPDATE_WRITE_OFFSET) {

      // receive the updated offset
      MPI_Recv(&write_offset, sizeof(MPI_Offset), MPI_CHAR, status.MPI_SOURCE, 
	       status.MPI_TAG, MPI_COMM_WORLD, &status);

      // and then check if there are any reader waiting for the offset
      if (writers_length) {
	write_lock = 1;
	write_lock_owner = writers[first_writer];
	MPI_Send(&write_offset, sizeof(MPI_Offset), MPI_CHAR, write_lock_owner, 
		 SET_WRITE_OFFSET, MPI_COMM_WORLD);
	
	// update waiting readers
	first_writer++;
	if (first_writer >= np) first_writer = 0;
	writers_length--;
      } else {
	write_lock = 0;
      }      
      
    } else if (status.MPI_TAG == WRITE_END_OF_FILE) {

      // the receive incoming message
      MPI_Recv(&skip, 1, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);

      // and delete
      num_writers--;
      if (num_readers == 0 && num_writers == 0) {
	break;
      }
    }
  }
}

//--------------------------------------------------------------------

#endif // _MPI
