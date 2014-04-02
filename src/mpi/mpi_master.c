#include "mpi_master.h"

//--------------------------------------------------------------------

void master_fastq_reader(options_t *options, list_t *in_list);
void master_send_recv(options_t *options, list_t *in_list, list_t *out_list);
void master_sam_writer(options_t *options, list_t *out_list);

//--------------------------------------------------------------------
// MPI master
//--------------------------------------------------------------------

void mpi_master(options_t *options) {

  int bam_format = 0;
  int num_threads = options->num_cpu_threads;

  list_t in_list;
  list_init("in", 1, 10, &in_list);

  list_t out_list;
  list_init("out", 1, 10, &out_list);

  int tid;
  #pragma omp parallel num_threads(3) private(tid)
  {
    tid = omp_get_thread_num();
    if (tid == 0) {
      // fastq reader
      printf("\tmaster: fastq reader (%i), begin...\n", tid);      
      master_fastq_reader(options, &in_list);
      printf("\tmaster:...end of fastq reader (%i). Done !!\n", tid);      
    } else if (tid == 1) {
      // sam writer
      printf("\tmaster: sam writer (%i)\n", tid);      
      master_sam_writer(options, &out_list);
      printf("\tmaster:...end of sam writer (%i). Done !!\n", tid);      
    } else if (tid == 2) {
      // master receiver
      printf("\tmaster: send, recv (%i)\n", tid);
      master_send_recv(options, &in_list, &out_list);
      printf("\tmaster:...end of send, recv (%i). Done !!\n", tid);      
    }
  }
}

//--------------------------------------------------------------------

void master_fastq_reader(options_t *options, list_t *in_list) {

  // preparing input FastQ file
  fastq_batch_reader_input_t reader_input;

  char *fq_list1 = options->in_filename, *fq_list2 = options->in_filename2;
  char token[2] = ",";
  char *ptr;
  array_list_t *files_fq1 = array_list_new(50,
					   1.25f,
					   COLLECTION_MODE_ASYNCHRONIZED);

  array_list_t *files_fq2 = array_list_new(50,
					   1.25f,
					   COLLECTION_MODE_ASYNCHRONIZED);
  int num_files1 = 0, num_files2 = 0;

  if (fq_list1) {
    ptr = strtok(fq_list1, token);    // Primera llamada => Primer token
    array_list_insert(strdup(ptr), files_fq1);
    while( (ptr = strtok( NULL, token)) != NULL ) {    // Posteriores llamadas
      array_list_insert(strdup(ptr), files_fq1);
    }
    num_files1 = array_list_size(files_fq1);
  }
  
  if (fq_list2) {
   ptr = strtok(fq_list2, token);    // Primera llamada => Primer token
   array_list_insert(strdup(ptr), files_fq2);
   while( (ptr = strtok( NULL, token)) != NULL ) {    // Posteriores llamadas
     array_list_insert(strdup(ptr), files_fq2);
   }    
   num_files2 = array_list_size(files_fq2);
  }
 
  
  if (fq_list2 && (num_files1 != num_files2)) {
    LOG_FATAL("Diferent number of files in paired-end/mate-pair mode");
  }
  
  list_item_t *item;
  size_t num_reads;
  array_list_t *reads;

  char *file1, *file2;
  for (int f = 0; f < num_files1; f++) {
    // open FastQ files
    file1 = array_list_get(f, files_fq1);
    if (num_files2) {
      file2 = array_list_get(f, files_fq2);
    } else {
      file2 = NULL;
    }
    
    fastq_batch_reader_input_init(file1, file2, 
				  options->pair_mode, 
				  options->batch_size, 
				  NULL, options->gzip, 
				  &reader_input);
    
    if (options->pair_mode == SINGLE_END_MODE) {
      if (options->gzip) {
	reader_input.fq_gzip_file1 = fastq_gzopen(file1);
      } else {
	reader_input.fq_file1 = fastq_fopen(file1);
      }
    } else {
      if (options->gzip) {
	reader_input.fq_gzip_file1 = fastq_gzopen(file1);
	reader_input.fq_gzip_file2 = fastq_gzopen(file2);	
      } else {
	reader_input.fq_file1 = fastq_fopen(file1);
	reader_input.fq_file2 = fastq_fopen(file2);
      }
    }

    // main loop
    while (1) {
      reads = array_list_new(reader_input.batch_size, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

      if (reader_input.gzip) {
	// Gzip FastQ file
	if (reader_input.flags == SINGLE_END_MODE) {
	  fastq_gzread_bytes_se(reads, reader_input.batch_size, reader_input.fq_gzip_file1);
	} else {
	  fastq_gzread_bytes_pe(reads, reader_input.batch_size, reader_input.fq_gzip_file1, reader_input.fq_gzip_file2);
	}
      } else {
	// FastQ file
	if (reader_input.flags == SINGLE_END_MODE) {
	  fastq_fread_bytes_se(reads, reader_input.batch_size, reader_input.fq_file1);
	} else {
	  fastq_fread_bytes_aligner_pe(reads, reader_input.batch_size, 
				       reader_input.fq_file1, reader_input.fq_file2);
	}
      }

      num_reads = array_list_size(reads);
      printf("\t\t\t\tmaster_fastq_reader: %i reads from FastQ file\n", num_reads);
      if (num_reads == 0) {
	array_list_free(reads, (void *) fastq_read_free);
	break;
      } else {
	// insert this batch to the corresponding list
	item = list_item_new(0, 0, (void *)reads);
	list_insert_item(item, in_list);
	printf("\t\t\t\tmaster_fastq_reader: inserting item, num. items = %i\n", array_list_size(in_list));
      }
    }
    //printf("--------> reader finished\n");

    // close FastQ files
    if (options->gzip) {
      if (options->pair_mode == SINGLE_END_MODE) {
	fastq_gzclose(reader_input.fq_gzip_file1);
      } else {
	fastq_gzclose(reader_input.fq_gzip_file1);
	fastq_gzclose(reader_input.fq_gzip_file2);
      }
    } else {
      if (options->pair_mode == SINGLE_END_MODE) {
	fastq_fclose(reader_input.fq_file1);
      } else {
	fastq_fclose(reader_input.fq_file1);
	fastq_fclose(reader_input.fq_file2);
      }
    }
  }

  // free memory
  array_list_free(files_fq1, (void *) free);
  array_list_free(files_fq2, (void *) free);

  // decrease the number of writer for this list
  list_decr_writers(in_list);
}

//--------------------------------------------------------------------

void master_send_recv(options_t *options, list_t *in_list, list_t *out_list) {
 
  int skip;
  list_item_t *item;

  int msg_size;
  char *msg_data;
  array_list_t *reads;

  MPI_Status status;

  int np;
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  while (1) {
    // probe for an incoming message from any process
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    printf("\t\t\t\tmaster_send_recv: probe, source = %i, tag = %i\n", status.MPI_SOURCE, status.MPI_TAG);
    if (status.MPI_TAG == GET_FASTQ_BATCH_MSG) {
      // you must consumn this message
      MPI_Recv(&skip, 1, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);

      item = list_remove_item(in_list);
      if (item == NULL) {
	assert(list_get_writers(in_list) <= 0);
	printf("\t\t\t\tmaster_send_recv: send END_OF_FILE_MSG, to process %i\n", status.MPI_SOURCE);
	MPI_Send(&skip, 1, MPI_INT, status.MPI_SOURCE, END_OF_FILE_MSG, MPI_COMM_WORLD);
      } else {
	reads = (array_list_t *) item->data_p;
	printf("\t\t\t\tmaster_send_recv: pack FastQ and send %i reads\n", array_list_size(reads));
	msg_data = pack_fastq(reads, options->batch_size, &msg_size);
	MPI_Send(msg_data, msg_size, MPI_CHAR, status.MPI_SOURCE, FASTQ_BATCH_MSG, MPI_COMM_WORLD);
	printf("\t\t\t\tmaster_send_recv: -------> sent %i bytes (%s)\n", msg_size, msg_data);

	// free memory
	// we must free here the msg_data !!!
	array_list_free(reads, (void *) fastq_read_free);
	list_item_free(item);
      }
    } else if (status.MPI_TAG == SAM_BATCH_MSG) {
      // when probe returns, the status object has the size and other
      // attributes of the incoming message. Get the size of the message
      MPI_Get_count(&status, MPI_CHAR, &msg_size);
      printf("\t\t\tmaster_send_recv: received %i bytes\n", msg_size);

      // allocate a buffer just big enough to hold the incoming numbers
      msg_data = (char *) calloc(msg_size + 1, sizeof(char));
      
      // Now receive the message with the allocated buffer
      MPI_Recv(msg_data, msg_size, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      printf("\t\t\tmaster_send_recv: dynamically received %d numbers from %i ---> %s\n", msg_size, status.MPI_SOURCE, msg_data);
	     

      //      reads = unpack_sam(msg_size, msg_data);
      
      item = list_item_new(0, 0, (void *) NULL);
      list_insert_item(item, out_list);

    } else if (status.MPI_TAG == END_OF_WORKER_MSG) {
      np--;
      printf("\t\t\tmaster_send_recv: received END_OF_WORKER_MSG (num. left workers = %i)\n", (np - 1));
      exit(-1);

      if (np == 1) break;
      break;
      //new_wf_batch = NULL;
    } else {
      printf("Error: master received an unknown message (%i)\n", status.MPI_TAG);
      break;
    }
  }

  list_decr_writers(out_list);
}

//--------------------------------------------------------------------

void master_sam_writer(options_t *options, list_t *out_list) {

  void *batch;
  list_item_t *item;
  while ((item = list_remove_item(out_list)) != NULL) {
    
    // batch = item->data_p;
    //sa_sam_writer(batch);
    
    list_item_free(item);
  }
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
