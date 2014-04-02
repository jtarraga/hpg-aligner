#include "mpi_commons.h"


//--------------------------------------------------------------------

char *pack_fastq(array_list_t *reads, int batch_size, int *msg_size) {
  printf("pack_fastq: num. reads = %i, batch size = %i\n", array_list_size(reads), batch_size);

  int data_size = 0;

  fastq_read_t *read;
  int num_reads = array_list_size(reads);

  char *data = (char *) malloc((batch_size + (num_reads * 2 * sizeof(int))) * sizeof(char));

  int id_len, char_count = 0;
  char *char_p = &data[(num_reads * 2 * sizeof(int)) + 1];
  
  int *int_p = (int *) data;

  *int_p = num_reads;
  int_p++;
  for (int i = 0; i < num_reads; i++) {
    read = array_list_get(i, reads);
    id_len = strlen(read->id); 

    *int_p = id_len;
    int_p++;
    *int_p = read->length; 
    int_p++;
    
    printf("\tpack: char_p = %i\n", char_p);
    strcpy(char_p, read->id);
    printf("\t\tpack: char_p = %s\n", char_p);
    char_p += id_len;
    strcpy(char_p, read->sequence);
    char_p += read->length;
    strcpy(char_p, read->quality);
    char_p += read->length;

    printf("(id, seq, qual) = (%s, %s, %s)\n", read->id, read->sequence, read->quality);    
  }
  data_size = abs(char_p - data);
  printf("data_size = %i\n", data_size);

  //  unpack_fastq(data_size, data);
  {
    printf("unpacking fastq:\n");
    int len;
    int *int_p = (int *) data;

    num_reads = *int_p;
    int_p++;
    printf("num_reads = %i\n", num_reads);

    char *char_p = &data[(num_reads * 2 * sizeof(int)) + 1];
    printf("\tunpack: char_p = %i\n", char_p);
    printf("char_p = %s\n", &char_p[1]);

    for (int i = 0; i < num_reads; i++) {
      printf("\t\tread %i\n", i);
      // id length
      len = *int_p;
      printf("\t\t\tid length = %i\n", len);
      int_p++;
      
      // sequence length
      len = *int_p;
      printf("\t\t\tseq length = %i\n", len);
      int_p++;           
    }

    array_list_t *reads = array_list_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    for (int i = 0; i < num_reads; i++) {
      //read = array_list_get(i, reads);
      //printf("(id, seq, qual) = (%s, %s, %s)\n", read->id, read->sequence, read->quality);    
    }
  }



  exit(-1);

  sprintf(data, "Hello world - FastQ");
  data_size = strlen(data);

  *msg_size = data_size;
  return data;
}

//--------------------------------------------------------------------

array_list_t *unpack_fastq(int msg_size, char *msg_data) {
  array_list_t *reads;
  printf("unpack_fastq: msg size = %i, msg data = %x\n", msg_size, msg_data);
  printf("\tunpack_fastq_reads: msg data = %s\n", msg_data);

  return reads;
}

//--------------------------------------------------------------------

char *pack_sam(array_list_t *reads, int batch_size, int *msg_size) {
  printf("pack_sam: batch size = %i\n", batch_size);
  int data_size;
  char *data = (char *) calloc(1000, sizeof(char));

  sprintf(data, "Hello world - SAM");
  data_size = strlen(data);

  *msg_size = data_size;
  return data;
}

//--------------------------------------------------------------------

array_list_t *unpack_sam(int msg_size, char *msg_data) {
  array_list_t *reads;
  printf("unpack_fastq_reads: msg size = %i, msg data = %x\n", msg_size, msg_data);
  printf("\tunpack_fastq_reads: msg data = %s\n", msg_data);

  return reads;
}

/*
//====================================================================
// WORKER
//====================================================================

//--------------------------------------------------------------------
// mpi_worker_receiver
//--------------------------------------------------------------------

void *mpi_worker_receiver(void *input) {
  sa_wf_input_t *wf_input = (sa_wf_input_t *) input;
  
  sa_wf_batch_t *new_wf_batch;
  sa_wf_batch_t *curr_wf_batch = wf_input->wf_batch;
  
  fastq_batch_reader_input_t *fq_reader_input = wf_input->fq_reader_input;
  array_list_t *reads = array_list_new(fq_reader_input->batch_size, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  // probe for an incoming message from any process
  MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

  //    printf("master: probe, source = %i, tag = %i\n", status.MPI_SOURCE, status.MPI_TAG);
  if (status.MPI_TAG == EOF) {
    new_wf_batch = NULL;
  } else if (status.MPI_TAG == FASTQ_BATCH) {
    int msg_size;
    char *mgs_data;

    // when probe returns, the status object has the size and other
    // attributes of the incoming message. Get the size of the message
    MPI_Get_count(&status, MPI_INT, &msg_size);
 
    // allocate a buffer just big enough to hold the incoming numbers
    msg_data = (char *) malloc(sizeof(char) * msg_size);
 
    // Now receive the message with the allocated buffer
    MPI_Recv(msg_data, msg_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("1 dynamically received %d numbers from 0.\n", msg_size);

    sa_mapping_batch_t *sa_mapping_batch = sa_mapping_batch_new(reads);
    sa_mapping_batch->bam_format = wf_input->bam_format;

    new_wf_batch = sa_wf_batch_new(curr_wf_batch->options,
				   curr_wf_batch->sa_index,
				   curr_wf_batch->writer_input, 
				   sa_mapping_batch);

    free(msg_data);
  } else {
    printf("Unknow MPI_TAG (%i)\n", status.MPI_TAG);
    abort(-1);
  }

  return new_wf_batch;
}

//--------------------------------------------------------------------
// mpi_worker_sender
//--------------------------------------------------------------------

int mpi_worker_sender(void *data) {
  sa_wf_batch_t *wf_batch = (sa_wf_batch_t *) data;
  
  sa_mapping_batch_t *mapping_batch = (sa_mapping_batch_t *) wf_batch->mapping_batch;
  if (mapping_batch == NULL) {
    printf("mpi_sender: error, NULL mapping batch\n");
    return 0;
  }

  //  for (int i = 0; i < NUM_COUNTERS; i++) {
  //    counters[i] += mapping_batch->counters[i];
  //  }

  #ifdef _TIMING
  for (int i = 0; i < NUM_TIMING; i++) {
    func_times[i] += mapping_batch->func_times[i];
  }
  #endif

  int num_mismatches, num_cigar_ops;
  size_t flag, pnext = 0, tlen = 0;
  char *cigar_string, *cigar_M_string, *rnext = "*";

  fastq_read_t *read;
  array_list_t *read_list = mapping_batch->fq_reads;

  array_list_t *mapping_list;
  FILE *out_file = (FILE *) wf_batch->writer_input->bam_file;

  sa_genome3_t *genome = wf_batch->sa_index->genome;

  size_t num_reads, num_mappings;
  num_reads = mapping_batch->num_reads;

  seed_cal_t *cal;
  for (size_t i = 0; i < num_reads; i++) {
    read = (fastq_read_t *) array_list_get(i, read_list);
    mapping_list = mapping_batch->mapping_lists[i];
    num_mappings = array_list_size(mapping_list);
#ifdef _VERBOSE
    if (num_mappings > 1) {
      num_dup_reads++;
      num_total_dup_reads += num_mappings;
    }
#endif
    
    if (num_mappings > 0) {
      num_mapped_reads++;
      for (size_t j = 0; j < num_mappings; j++) {
	cal = (seed_cal_t *) array_list_get(j, mapping_list);
	
	flag = (cal->strand ? 16 : 0);
	cigar_string = cigar_to_string(&cal->cigar);
	cigar_M_string = cigar_to_M_string(&num_mismatches, &num_cigar_ops, &cal->cigar);
	fprintf(out_file, "%s\t%i\t%s\t%lu\t%i\t%s\t%s\t%lu\t%lu\t%s\t%s\tNH:i:%i\tNM:i:%i\tXC:Z:%s\n", 
		read->id,
		flag,
		genome->chrom_names[cal->chromosome_id],
		cal->start + 1,
		cal->AS,
		cigar_M_string,
		rnext,
		pnext,
		tlen,
		read->sequence,
		read->quality,
		num_mappings,
		num_mismatches,
		cigar_string
		);
	free(cigar_M_string);
	free(cigar_string);
	seed_cal_free(cal);	 
      }
    } else {
      num_unmapped_reads++;
      fprintf(out_file, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", 
	      read->id,
	      read->sequence,
	      read->quality
	      );
    }
    
      array_list_free(mapping_list, (void *) NULL);
  }

  // free memory
  sa_mapping_batch_free(mapping_batch);

  if (wf_batch) sa_wf_batch_free(wf_batch);

  return 0;
}
*/


