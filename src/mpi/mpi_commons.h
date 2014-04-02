#ifndef _MPI_COMMONS_H
#define _MPI_COMMONS_H


#include "containers/array_list.h"

#include "bioformats/fastq/fastq_read.h"

//--------------------------------------------------------------------

#define END_OF_FILE_MSG     0
#define END_OF_WORKER_MSG   1
#define GET_FASTQ_BATCH_MSG 3
#define FASTQ_BATCH_MSG     4
#define SAM_BATCH_MSG       5


//--------------------------------------------------------------------

//void *mpi_receiver(void *input);
//int mpi_sender(void *data);


//--------------------------------------------------------------------

char *pack_fastq(array_list_t *reads, int batch_size, int *msg_size);
array_list_t *unpack_fastq(int msg_size, char *msg_data);

//--------------------------------------------------------------------

char *pack_sam(array_list_t *reads, int batch_size, int *msg_size);
array_list_t *unpack_sam(int msg_size, char *msg_data);

//--------------------------------------------------------------------
//--------------------------------------------------------------------
#endif // _MPI_COMMONS_H
