#include "mp_mapper.h"

//--------------------------------------------------------------------
// main functions
//--------------------------------------------------------------------

void mp_fastq_reader(int id, int np, options_t *options, sa_index3_t *sa_index, 
		      list_t *in_list);
void mp_read_mapper(int id, options_t *options, sa_index3_t *sa_index, 
		     list_t *in_list, list_t *out_list);
void mp_sam_writer(int id, options_t *options, sa_index3_t *sa_index, 
		    list_t *out_list);
void mp_sam_concat(int id, int np, options_t *options, sa_index3_t *sa_index);

//--------------------------------------------------------------------
// Read functions
//--------------------------------------------------------------------

size_t mp_fastq_fread_bytes_se(int id, int np, array_list_t *reads, size_t bytes, 
				fastq_file_t *fq_file);
size_t mp_fastq_fread_bytes_aligner_pe(int id, int np, array_list_t *reads, size_t bytes, 
					fastq_file_t *fq_file1, fastq_file_t *fq_file2);

//--------------------------------------------------------------------
// wrtite functions
//--------------------------------------------------------------------

void write_mapping_batch(FILE *out_file, sa_mapping_batch_t *mapping_batch,
			 sa_genome3_t *genome);

//--------------------------------------------------------------------
// MP mapper
//--------------------------------------------------------------------

void mp_mapper(int id, int np, options_t *options) {

  sa_index3_t *sa_index = sa_index3_new(options->bwt_dirname);

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
      // fastq reader
      mp_fastq_reader(id, np, options, sa_index, &in_list);
    } else if (tid == 1) {
      // sam writer
      mp_sam_writer(id, options, sa_index, &out_list);
      /*
    } else if (tid == 2 && id == 0) {
      list_decr_writers(&out_list);
      mp_sam_concat(id, np, options, sa_index);
      */
    } else {
      // master receiver
      mp_read_mapper(id, options, sa_index, &in_list, &out_list);
    }
  }

  if (sa_index) sa_index3_free(sa_index);
}

//--------------------------------------------------------------------

void mp_fastq_reader(int id, int np, options_t *options, sa_index3_t *sa_index, 
		      list_t *in_list) {
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

  sa_mapping_batch_t *sa_mapping_batch;
  sa_wf_batch_t *wf_batch;

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
	  //	  mp_fastq_gzread_bytes_se(id, reads, reader_input.batch_size, reader_input.fq_gzip_file1);
	} else {
	  //	  mp_fastq_gzread_bytes_pe(id, reads, reader_input.batch_size, reader_input.fq_gzip_file1, reader_input.fq_gzip_file2);
	}
      } else {
	// FastQ file
	if (reader_input.flags == SINGLE_END_MODE) {
	  mp_fastq_fread_bytes_se(id, np, reads, reader_input.batch_size, reader_input.fq_file1);
	} else {
	  mp_fastq_fread_bytes_aligner_pe(id, np, reads, reader_input.batch_size, 
					   reader_input.fq_file1, reader_input.fq_file2);
	}
      }

      num_reads = array_list_size(reads);
      //printf("\t\t\t\tmaster_fastq_reader: %i reads from FastQ file\n", num_reads);
      if (num_reads == 0) {
	array_list_free(reads, (void *) fastq_read_free);
	break;
      } else {
	// insert this batch to the corresponding list
	sa_mapping_batch = sa_mapping_batch_new(reads);
	sa_mapping_batch->bam_format = 0;
	wf_batch = sa_wf_batch_new(options, sa_index, NULL, sa_mapping_batch, NULL);
      
	item = list_item_new(0, 0, (void *) wf_batch);
	list_insert_item(item, in_list);
	//printf("\t\t\t\tmaster_fastq_reader: inserting item, num. items = %i\n", array_list_size(in_list));
      }
    }
    ////printf("--------> reader finished\n");

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

void mp_read_mapper(int id, options_t *options, sa_index3_t *sa_index, 
		     list_t *in_list, list_t *out_list) {
  void *batch;
  list_item_t *item, *new_item;
  while ((item = list_remove_item(in_list)) != NULL) {
    //printf("\t\t\tworker_mapper: processing batch\n");

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
  //printf("\t\t\tworker_mapper: Done!!\n");
  //printf("----------> mapper finished\n");
  list_decr_writers(out_list);
}

//--------------------------------------------------------------------

void mp_sam_writer(int id, options_t *options, sa_index3_t *sa_index, list_t *out_list) {

  int max_batch_counter = 100;
  int file_counter = 0, batch_counter = 0;

  // setting output name
  int bam_format = 0;
  int len = 200;
  if (options->prefix_name) {
    len += strlen(options->prefix_name);
  }
  if (options->output_name) {
    len += strlen(options->output_name);
  }
  char filename[len], out_filename[len], out_filename_done[len];
  filename[0] = 0;
  strcat(filename, (options->output_name ? options->output_name : "."));
  strcat(filename, "/");
  if (options->prefix_name) {
    strcat(filename, options->prefix_name);
    strcat(filename, "_");
  }

  strcat(filename, "out.sam");

  sprintf(out_filename, "%s.%i.%i", filename, id, file_counter);

  // open file
  FILE *out_file = fopen(out_filename, "w");    
  //write_sam_header(genome, out_file);

  // main loop
  sa_wf_batch_t *wf_batch;
  sa_mapping_batch_t *mapping_batch;

  list_item_t *item;
  while ((item = list_remove_item(out_list)) != NULL) {
    
    wf_batch = (sa_wf_batch_t *) item->data_p;
    mapping_batch = wf_batch->mapping_batch;
    //printf("\t\t\master_sam_writer: writing SAM...\n");

    write_mapping_batch(out_file, mapping_batch, sa_index->genome);

    batch_counter++;
    if (batch_counter == max_batch_counter) {
      batch_counter = 0;
      fclose(out_file);

      sprintf(out_filename_done, "%s__subfile_done__", out_filename);
      rename(out_filename, out_filename_done);
      //      if (!rename(out_filename, out_filename_done)) {
      //	printf("Error (%i) rename from %s to %s: %s\n", errno, out_filename, out_filename_done, strerror(errno));
      //      }

      file_counter++;
      sprintf(out_filename, "%s.%i.%i", filename, id, file_counter);

      // open file
      out_file = fopen(out_filename, "w");    
    }
    
    // free memory
    sa_mapping_batch_free(mapping_batch);
    sa_wf_batch_free(wf_batch);
    list_item_free(item);
  }

  // free memory and close file
  fclose(out_file);
  
  sprintf(out_filename_done, "%s__subfile_done__", out_filename);
  rename(out_filename, out_filename_done);  

  sprintf(out_filename_done, "%s__worker_done__", out_filename);
  out_file = fopen(out_filename_done, "w");
  fclose(out_file);
  
  //printf("\t\t\tmaster_sam_writer: Done!!\n");
}

//--------------------------------------------------------------------

void mp_sam_concat(int id, int np, options_t *options, sa_index3_t *sa_index) {

  // setting output name
  int bam_format = 0;
  int len = 200;
  if (options->prefix_name) {
    len += strlen(options->prefix_name);
  }
  if (options->output_name) {
    len += strlen(options->output_name);
  }
  char dirname[len], out_filename[len], filename_aux[len];
  dirname[0] = 0;
  strcat(dirname, (options->output_name ? options->output_name : "."));
  strcat(dirname, "/");

  strcpy(out_filename, dirname);
  if (options->prefix_name) {
    strcat(out_filename, options->prefix_name);
    strcat(out_filename, "_");
  }

  strcat(out_filename, "out.sam");

  int max_bytes = 2048, read_bytes;
  char data[max_bytes];

  // open file
  FILE *f, *out_file = fopen(out_filename, "w");    
  write_sam_header(sa_index->genome, out_file);

  DIR *dir;
  struct dirent *ent;

  int num_workers = np;

  while (num_workers > 0) {
    if ((dir = opendir(dirname)) != NULL) {
      /* print all the files and directories within directory */
      while ((ent = readdir (dir)) != NULL) {
	if (strstr(ent->d_name, "__worker_done__")) {
	  sprintf(filename_aux, "%s/%s", dirname, ent->d_name);
	  // decrease num workers and delete file
	  --num_workers;
	  remove(filename_aux);
	} else if (strstr(ent->d_name, "__subfile_done__")) {
	  sprintf(filename_aux, "%s/%s", dirname, ent->d_name);
	  // concat file and delete it
	  f = fopen(filename_aux, "r");
	  while ((read_bytes = fread(data, 1, max_bytes, f))) {
	    fwrite(data, 1, read_bytes, out_file);
	  }
	  fclose(f);
	  remove(filename_aux);	  
	}
      }
      closedir (dir);
    } else {
      /* could not open directory */
      perror ("could not open directory");
    }
    sleep(2);
  }

  // close file
  fclose(out_file);
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

void write_mapping_batch(FILE *out_file, sa_mapping_batch_t *mapping_batch,
			 sa_genome3_t *genome) {
  if (mapping_batch == NULL) {
    //printf("bam_writer1: error, NULL mapping batch\n");
    return;
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

  size_t num_reads, num_mappings;
  num_reads = mapping_batch->num_reads;


  if (mapping_batch->options->pair_mode != SINGLE_END_MODE) {

    // PAIR MODE
    alignment_t *alig;
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
	//	num_mapped_reads++;
	for (size_t j = 0; j < num_mappings; j++) {
	  alig = (alignment_t *) array_list_get(j, mapping_list);

	  flag = 0;
	  if (alig->is_paired_end)                              flag += BAM_FPAIRED;
	  if (alig->is_paired_end_mapped)                       flag += BAM_FPROPER_PAIR;
	  if (!alig->is_seq_mapped)                             flag += BAM_FUNMAP;   
	  if ((!alig->is_mate_mapped) && (alig->is_paired_end)) flag += BAM_FMUNMAP;
	  if (alig->seq_strand)                                 flag += BAM_FREVERSE;
	  if (alig->mate_strand)                                flag += BAM_FMREVERSE;
	  if (alig->pair_num == 1)	                        flag += BAM_FREAD1;
	  if (alig->pair_num == 2)                              flag += BAM_FREAD2;
	  if (alig->secondary_alignment)                        flag += BAM_FSECONDARY;
	  if (alig->fails_quality_check)                        flag += BAM_FQCFAIL;
	  if (alig->pc_optical_duplicate)                       flag += BAM_FDUP;

	  fprintf(out_file, "%s\t%i\t%s\t%lu\t%i\t%s\t%s\t%lu\t%i\t%s\t%s\t%s\n", 
		  read->id,
		  flag,
		  genome->chrom_names[alig->chromosome],
		  alig->position + 1,
		  alig->map_quality,
		  alig->cigar,
		  genome->chrom_names[alig->mate_chromosome],
		  alig->mate_position,
		  alig->template_length,
		  read->sequence,
		  read->quality,
		  (alig->optional_fields ? alig->optional_fields : "")
		  );
	}
	alignment_free(alig);	 
      } else {
	//	num_unmapped_reads++;
	fprintf(out_file, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", 
		read->id,
		read->sequence,
		read->quality
		);
      }
      
      array_list_free(mapping_list, (void *) NULL);
    }
  } else {

    // SINGLE MODE
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
	//	num_mapped_reads++;
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
	//	num_unmapped_reads++;
	fprintf(out_file, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", 
		read->id,
		read->sequence,
		read->quality
		);
      }
      
      array_list_free(mapping_list, (void *) NULL);
    }
  }
}

//--------------------------------------------------------------------

size_t mp_fastq_fread_bytes_se(int id, int np, array_list_t *reads, size_t bytes, fastq_file_t *fq_file) {
  size_t accumulated_size = 0;
  char header1[MAX_READ_ID_LENGTH];
  char sequence[MAX_READ_SEQUENCE_LENGTH];
  char header2[MAX_READ_ID_LENGTH];
  char qualities[MAX_READ_SEQUENCE_LENGTH];
  int header_length, sequence_length, quality_length;
  fastq_read_t *read;

  size_t counter = 0;

  while (accumulated_size < bytes && fgets(header1, MAX_READ_ID_LENGTH, fq_file->fd) != NULL) {
    fgets(sequence, MAX_READ_SEQUENCE_LENGTH, fq_file->fd);
    fgets(header2, MAX_READ_ID_LENGTH, fq_file->fd);
    fgets(qualities, MAX_READ_SEQUENCE_LENGTH, fq_file->fd);

    if (counter % np  == id) {

      header_length = strlen(header1);
      sequence_length = strlen(sequence);
      quality_length = strlen(qualities);
    
      // '\n' char is removed, but '\0' is left
      chomp_at(header1, header_length - 1);
      chomp_at(sequence, sequence_length - 1);
      chomp_at(qualities, quality_length - 1);
      
      read = fastq_read_new(header1, sequence, qualities);
      array_list_insert(read, reads);

      accumulated_size += header_length + sequence_length + quality_length;
    }

    counter++;
  }
  
  return accumulated_size;
}

//--------------------------------------------------------------------

size_t mp_fastq_fread_bytes_aligner_pe(int id, int np, array_list_t *reads, size_t bytes, 
					fastq_file_t *fq_file1, fastq_file_t *fq_file2) {
  size_t accumulated_size = 0;
  char header1[MAX_READ_ID_LENGTH];
  char header2[MAX_READ_ID_LENGTH];
  char read_separator[MAX_READ_ID_LENGTH];
  char sequence1[MAX_READ_SEQUENCE_LENGTH];
  char sequence2[MAX_READ_SEQUENCE_LENGTH];
  char qualities1[MAX_READ_SEQUENCE_LENGTH];
  char qualities2[MAX_READ_SEQUENCE_LENGTH];
  int header_length1, sequence_length1, quality_length1;
  int header_length2, sequence_length2, quality_length2;
  fastq_read_t *read1, *read2;
  
  while (accumulated_size < bytes && fgets(header1, MAX_READ_ID_LENGTH, fq_file1->fd) != NULL) {
    fgets(sequence1, MAX_READ_SEQUENCE_LENGTH, fq_file1->fd);
    fgets(read_separator, MAX_READ_ID_LENGTH, fq_file1->fd);
    fgets(qualities1, MAX_READ_SEQUENCE_LENGTH, fq_file1->fd);
    
    header_length1 = strlen(header1);
    sequence_length1 = strlen(sequence1);
    quality_length1 = strlen(qualities1);
    
    // '\n' char is removed, but '\0' is left
    chomp_at(header1, header_length1 - 1);
    chomp_at(sequence1, sequence_length1 - 1);
    chomp_at(qualities1, quality_length1 - 1);
    
    // second file
    fgets(header2, MAX_READ_ID_LENGTH, fq_file2->fd);
    fgets(sequence2, MAX_READ_SEQUENCE_LENGTH, fq_file2->fd);
    fgets(read_separator, MAX_READ_ID_LENGTH, fq_file2->fd);
    fgets(qualities2, MAX_READ_SEQUENCE_LENGTH, fq_file2->fd);
    
    header_length2 = strlen(header2);
    sequence_length2 = strlen(sequence2);
    quality_length2 = strlen(qualities2);
    
    // '\n' char is removed, but '\0' is left
    chomp_at(header2, header_length2 - 1);
    chomp_at(sequence2, sequence_length2 - 1);
    chomp_at(qualities2, quality_length2 - 1);
    
    read1 = fastq_read_new(header1, sequence1, qualities1);
    read2 = fastq_read_new(header2, sequence2, qualities2);
    
    array_list_insert(read1, reads);
    array_list_insert(read2, reads);
    
    accumulated_size += header_length1 + sequence_length1 + quality_length1 + 
      header_length2 + sequence_length2 + quality_length2;
  }
  
  return accumulated_size;
}



//--------------------------------------------------------------------
//--------------------------------------------------------------------
