#include "mpi_worker.h"

//--------------------------------------------------------------------

void worker_receiver(options_t *options, sa_index3_t *sa_index, list_t *in_list);
void worker_mapper(options_t *options, sa_index3_t *sa_index, list_t *in_list, list_t *out_list);
void worker_sender(options_t *options, sa_index3_t *sa_index, list_t *out_list);

//--------------------------------------------------------------------
// MPI worker
//--------------------------------------------------------------------

void mpi_worker(options_t *options) {
  sa_index3_t *sa_index = sa_index3_new(options->bwt_dirname);
  printf("sa_index loaded!!\n");

  int bam_format = 0;
  int num_threads = options->num_cpu_threads;

  list_t in_list;
  list_init("in", 1, 10, &in_list);

  list_t out_list;
  list_init("out", num_threads - 2, 10, &out_list);

  assert(num_threads > 2);

  int tid;
  #pragma omp parallel num_threads(num_threads) private(tid)
  {
    tid = omp_get_thread_num();
    if (tid == 0) {
      // worker receiver
      printf("\tworker: receiver (%i), begin...\n", tid);      
      worker_receiver(options, sa_index, &in_list);
      printf("\tworker:...end of receiver (%i). Done !!\n", tid);      
    } else if (tid == 1) {
      // worker sender
      printf("\tworker: sender (%i), begin...\n", tid);      
      worker_sender(options, sa_index, &out_list);
      printf("\tworker:...end of sender (%i). Done !!\n", tid);      
    } else {
      // worker mappers
      printf("\tworker: mapper (%i), begin...\n", tid);
      worker_mapper(options, sa_index, &in_list, &out_list);
      printf("\tworker:...end of mapper (%i). Done !!\n", tid);      
    }
  }

  if (sa_index) sa_index3_free(sa_index);
}

//--------------------------------------------------------------------

void worker_receiver(options_t *options, sa_index3_t *sa_index, list_t *in_list) {

  int skip, msg_size, eof = 0;
  list_item_t *item;
  char *msg_data;

  array_list_t *reads;
  sa_mapping_batch_t *sa_mapping_batch;
  sa_wf_batch_t *wf_batch;

  MPI_Status status;

  printf("\t\t\t\tworker_receiver, send GET_FASTQ_BATCH_MSG\n");
  MPI_Send(&skip, 1, MPI_INT, 0, GET_FASTQ_BATCH_MSG, MPI_COMM_WORLD);

  while (1) {
    // probe for an incoming message from master process
    MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    printf("\t\t\tworker_receiver, probe, source = %i, tag = %i\n", status.MPI_SOURCE, status.MPI_TAG);
    if (status.MPI_TAG == FASTQ_BATCH_MSG) {
      printf("\t\t\tworker_receiver: received FASTQ_BATCH_MSG\n");
      
      // when probe returns, the status object has the size and other
      // attributes of the incoming message. Get the size of the message
      MPI_Get_count(&status, MPI_CHAR, &msg_size);
      printf("\t\t\tworker_receiver: received %i bytes\n", msg_size);

      // allocate a buffer just big enough to hold the incoming numbers
      msg_data = (char *) calloc(msg_size + 1, sizeof(char));
      
      // Now receive the message with the allocated buffer
      MPI_Recv(msg_data, msg_size, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      printf("1 dynamically received %d numbers from master\n", msg_size);

      reads = unpack_fastq(msg_size, msg_data);
      free(msg_data);
      
      sa_mapping_batch = sa_mapping_batch_new(reads);
      sa_mapping_batch->bam_format = 0;
      wf_batch = sa_wf_batch_new(options, sa_index, NULL, sa_mapping_batch, NULL);
      
      item = list_item_new(0, 0, (void *) wf_batch);
      list_insert_item(item, in_list);

      printf("************************************** sending GET_FASTQ_BATCH_MSG\n");
      MPI_Send(&skip, 1, MPI_INT, 0, GET_FASTQ_BATCH_MSG, MPI_COMM_WORLD);

      // can we avoid this free ?
      // free(msg_data);
    } else if (status.MPI_TAG == END_OF_FILE_MSG) {
      printf("\t\t\tworker_receiver: END_OF_FILE_MSG\n");
      break;
    } else {
      printf("Error: master received an unknown message (%i)\n", status.MPI_TAG);
      abort();
    }
  }

  // decrease the number of writer for this list
  list_decr_writers(in_list);
}

//--------------------------------------------------------------------

void worker_mapper(options_t *options, sa_index3_t *sa_index, 
		   list_t *in_list, list_t *out_list) {

  void *batch;
  list_item_t *item, *new_item;
  while ((item = list_remove_item(in_list)) != NULL) {
    printf("\t\t\tworker_mapper: processing batch\n");
        
    batch = item->data_p;

    if (options->pair_mode == SINGLE_END_MODE) {
      sa_single_mapper(batch);
    } else {
      sa_pair_mapper(batch);
    }
    
    // insert this batch to the corresponding list
    new_item = list_item_new(0, 0, batch);
    list_insert_item(new_item, out_list);

    // free memory
    list_item_free(item);
  }
  printf("\t\t\tworker_mapper: Done!!\n");
  //  printf("----------> mapper finished\n");
  list_decr_writers(out_list);
}


//--------------------------------------------------------------------

void worker_sender(options_t *options, sa_index3_t *sa_index, list_t *out_list) {
  void *batch;

  int skip;
  list_item_t *item;

  int msg_size;
  char *msg_data;

  while ((item = list_remove_item(out_list)) != NULL) {

    printf("\t\t\tworker_sender: sending SAM...\n");

    batch = item->data_p;

    //reads = (array_list_t *) item->data_p;
    printf("\t\t\t\tworker_sender: pack SAM and send\n");
    msg_data = pack_sam(((sa_wf_batch_t *)batch)->mapping_batch, sa_index->genome, options->batch_size, &msg_size);
    MPI_Send(msg_data, msg_size, MPI_CHAR, 0, SAM_BATCH_MSG, MPI_COMM_WORLD);
    printf("\t\t\t\tworker_sender: -------> sent %i bytes\n", msg_size);

    // free memory
    list_item_free(item);
    free(msg_data);
  }

  printf("\t\t\tworker_sender: sending END_OF_WORKER...\n");
  MPI_Send(&skip, 1, MPI_INT, 0, END_OF_WORKER_MSG, MPI_COMM_WORLD);
  printf("\t\t\tworker_sender: Done!!\n");
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------





  /*
  //--------------------------------------------------------------------------------------
  // workflow management
  //
  sa_wf_batch_t *wf_batch = sa_wf_batch_new(options, (void *)sa_index, NULL, NULL);
  sa_wf_input_t *wf_input = sa_wf_input_new(bam_format, NULL, wf_batch);
  
  // create and initialize workflow
  workflow_t *wf = workflow_new();
  
  workflow_stage_function_t stage_functions[1];
  char *stage_labels[1] = {"SA mapper"};
  if (options->pair_mode == SINGLE_END_MODE) {
    stage_functions[0] = sa_single_mapper;
  } else {
    stage_functions[0] = sa_pair_mapper;
  }
  workflow_set_stages(1, &stage_functions, stage_labels, wf);
  
  // optional producer and consumer functions
  workflow_set_producer(mpi_receiver, "MPI receiver", wf);
  workflow_set_consumer(mpi_sender, "MPI sender", wf);
  
  printf("----------------------------------------------\n");
  printf("Starting mapping...\n");
  gettimeofday(&start, NULL);
  workflow_run_with(num_threads, wf_input, wf);
  gettimeofday(&stop, NULL);
  printf("End of mapping in %0.2f min. Done!!\n", 
	 ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f)/60.0f);  
  
  printf("----------------------------------------------\n");
  workflow_display_timing(wf);
  printf("----------------------------------------------\n");
  
  
  printf("----------------------------------------------\n");
  printf("Num. reads         : %lu\nNum. mapped reads  : %lu (%0.2f %%)\nNum. unmapped reads: %lu (%0.2f %%)\n",
	 num_mapped_reads + num_unmapped_reads,
	 num_mapped_reads, 100.0f * num_mapped_reads / (num_mapped_reads + num_unmapped_reads),
	 num_unmapped_reads, 100.0f * num_unmapped_reads / (num_mapped_reads + num_unmapped_reads));
  printf("----------------------------------------------\n");
  
  
  #ifdef _TIMING
  char func_name[1024];
  double total_func_times = 0;
  for (int i = 0; i < NUM_TIMING; i++) {
    if (i != FUNC_SEARCH_PREFIX && i != FUNC_SEARCH_SA 
	&& i < FUNC_INIT_CALS_FROM_SUFFIXES || i > FUNC_CAL_MNG_INSERT) {
      total_func_times += func_times[i];
    }
  }
  printf("Timing in seconds:\n");
  for (int i = 0; i < NUM_TIMING; i++) {
    if (i == FUNC_SEARCH_PREFIX || i == FUNC_SEARCH_SA ||
	(i >= FUNC_INIT_CALS_FROM_SUFFIXES && i <= FUNC_CAL_MNG_INSERT)) {
      printf("\t");
    }
    printf("\t%0.2f %%\t%0.4f\tof %0.4f\t%s\n", 
	   100.0 * func_times[i] / total_func_times, func_times[i], total_func_times, func_names[i]);
  }
  #endif
  
  // free memory
  sa_wf_input_free(wf_input);
  sa_wf_batch_free(wf_batch);
  workflow_free(wf);
}
  */
