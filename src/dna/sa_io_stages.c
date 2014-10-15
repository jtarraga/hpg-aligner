#include "sa_io_stages.h"

//====================================================================
// PRODUCER
//====================================================================

//--------------------------------------------------------------------
// sa fq reader
//--------------------------------------------------------------------

void *sa_fq_reader(void *input) {
  sa_wf_input_t *wf_input = (sa_wf_input_t *) input;
  
  sa_wf_batch_t *new_wf_batch = NULL;
  sa_wf_batch_t *curr_wf_batch = wf_input->wf_batch;
  
  fastq_batch_reader_input_t *fq_reader_input = wf_input->fq_reader_input;
  array_list_t *reads = array_list_new(fq_reader_input->batch_size, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  if (fq_reader_input->gzip) {
    // Gzip fastq file
    if (fq_reader_input->flags == SINGLE_END_MODE) {
      fastq_gzread_bytes_se(reads, fq_reader_input->batch_size, fq_reader_input->fq_gzip_file1);
    } else {
      fastq_gzread_bytes_pe(reads, fq_reader_input->batch_size, fq_reader_input->fq_gzip_file1, fq_reader_input->fq_gzip_file2);
    }
  } else {
    // Fastq file
    if (fq_reader_input->flags == SINGLE_END_MODE) {
      fastq_fread_bytes_se(reads, fq_reader_input->batch_size, fq_reader_input->fq_file1);
    } else {
      fastq_fread_bytes_aligner_pe(reads, fq_reader_input->batch_size, 
				   fq_reader_input->fq_file1, fq_reader_input->fq_file2);
    }
  }
  
  size_t num_reads = array_list_size(reads);
  
  if (num_reads == 0) {
    array_list_free(reads, (void *)fastq_read_free);
  } else {
    sa_mapping_batch_t *sa_mapping_batch = sa_mapping_batch_new(reads);
    sa_mapping_batch->bam_format = wf_input->bam_format;

    new_wf_batch = sa_wf_batch_new(curr_wf_batch->options,
				   curr_wf_batch->sa_index,
				   curr_wf_batch->writer_input, 
				   sa_mapping_batch,
				   NULL);
  }

  return new_wf_batch;

}

//====================================================================
// CONSUMER
//====================================================================

#ifdef _VERBOSE
int num_dup_reads = 0;
int num_total_dup_reads = 0;
#endif

size_t num_mapped_reads = 0;
size_t num_unmapped_reads = 0;
size_t num_total_mappings = 0;
size_t num_multihit_reads = 0;
size_t num_unmapped_reads_by_invalid_cal = 0;
size_t num_unmapped_reads_by_cigar_length = 0;

//--------------------------------------------------------------------
// SAM writer
//--------------------------------------------------------------------

void *write_sam_header(sa_genome3_t *genome, FILE *f) {
  fprintf(f, "@PG\tID:HPG-Aligner\tVN:%s\n", HPG_ALIGNER_VERSION);
  for (int i = 0; i < genome->num_seqs; i++) {
    fprintf(f, "@SQ\tSN:%s\tLN:%lu\n", genome->seq_names[i], genome->seq_lengths[i]);
  }
}

//--------------------------------------------------------------------

int sa_sam_writer(void *data) {
  sa_wf_batch_t *wf_batch = (sa_wf_batch_t *) data;
  
  sa_mapping_batch_t *mapping_batch = (sa_mapping_batch_t *) wf_batch->mapping_batch;
  if (mapping_batch == NULL) {
    printf("bam_writer1: error, NULL mapping batch\n");
    return 0;
  }

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

  array_list_t *mapping_list, *mate_list;
  FILE *out_file = (FILE *) wf_batch->writer_input->bam_file;

  sa_genome3_t *genome = wf_batch->sa_index->genome;

  size_t num_reads, num_mappings, num_mate_mappings;
  num_reads = mapping_batch->num_reads;

  if (mapping_batch->options->pair_mode != SINGLE_END_MODE) {
    // PAIR MODE
    int len;
    char *sequence, *quality;

    char *seq, *opt_fields;
    alignment_t *alig;
  

    for (size_t i = 0; i < num_reads; i++) {
      read = (fastq_read_t *) array_list_get(i, read_list);

      mapping_list = mapping_batch->mapping_lists[i];
      num_mappings = array_list_size(mapping_list);
      num_total_mappings += num_mappings;

      #ifdef _VERBOSE
      if (num_mappings > 1) {
	num_dup_reads++;
	num_total_dup_reads += num_mappings;
      }
      #endif
      
      if (num_mappings > 0) {
	num_mapped_reads++;
	if (num_mappings > 1) {
	  num_multihit_reads++;
	}
	for (size_t j = 0; j < num_mappings; j++) {
	  alig = (alignment_t *) array_list_get(j, mapping_list);

	  if (alig->optional_fields) {
	    opt_fields = alig->optional_fields;
	  } else {
	    opt_fields = NULL;
	  }

	  flag = 0;
	  if (alig->is_paired_end)                              flag += BAM_FPAIRED;
	  if (alig->is_paired_end_mapped)                       flag += BAM_FPROPER_PAIR;
	  if (!alig->is_seq_mapped)                             flag += BAM_FUNMAP;   
	  if ((!alig->is_mate_mapped) && (alig->is_paired_end)) flag += BAM_FMUNMAP;
	  if (alig->mate_strand)                                flag += BAM_FMREVERSE;
	  if (alig->pair_num == 1)	                        flag += BAM_FREAD1;
	  if (alig->pair_num == 2)                              flag += BAM_FREAD2;
	  if (alig->secondary_alignment)                        flag += BAM_FSECONDARY;
	  if (alig->fails_quality_check)                        flag += BAM_FQCFAIL;
	  if (alig->pc_optical_duplicate)                       flag += BAM_FDUP;
	  if (alig->seq_strand)                                 flag += BAM_FREVERSE;

	  fprintf(out_file, "%s\t%lu\t%s\t%i\t%i\t%s\t%s\t%i\t%i\t%s\t%s\t%s\n", 
		  read->id,
		  flag,
		  genome->seq_names[alig->chromosome],
		  alig->position + 1,
		  (num_mappings > 1 ? 0 : alig->mapq),
		  alig->cigar,
		  (alig->chromosome == alig->mate_chromosome ? "=" : genome->seq_names[alig->mate_chromosome]),
		  alig->mate_position + 1,
		  alig->template_length,
		  alig->sequence,
		  alig->quality,
		  (opt_fields == NULL ? "" : opt_fields)
		  );

	  // free memory
	  alignment_free(alig);	 
	} // end for num_mappings
      } else {
	num_unmapped_reads++;

	if (read->adapter) {
	  len = read->length + abs(read->adapter_length);
	  sequence = (char *) malloc(len + 1);
	  quality = (char *) malloc(len + 1);

	  if (read->adapter_length < 0) {
	    strcpy(quality, read->adapter_quality);
	    strcat(quality, read->quality);
	  } else {
	    strcpy(quality, read->quality);
	    strcat(quality, read->adapter_quality);
	  }
	  
	  if ((read->adapter_strand == 0 && read->adapter_length < 0) || 
	      (read->adapter_strand == 1 && read->adapter_length > 0)) {
	    strcpy(sequence, read->adapter);
	    strcat(sequence, read->sequence);
	  } else {
	    strcpy(sequence, read->sequence);
	    strcat(sequence, read->adapter);
	  }

	  sequence[len] = 0; 
	  quality[len] = 0; 
	} else {
	  sequence = read->sequence;
	  quality = read->quality;
	}

	fprintf(out_file, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", 
		read->id,
		sequence,
		quality
		);

	if (read->adapter) {
	  free(sequence);
	  free(quality);
	}
      }
      array_list_free(mapping_list, (void *) NULL);
    }
  } else {
    // SINGLE MODE
    int len, mapq;
    char *seq;
    seed_cal_t *cal;

    cigar_t *cigar;
    char *sequence, *revcomp, *quality;

    for (size_t i = 0; i < num_reads; i++) {
      read = (fastq_read_t *) array_list_get(i, read_list);
      mapping_list = mapping_batch->mapping_lists[i];
      num_mappings = array_list_size(mapping_list);
      num_total_mappings += num_mappings;

      #ifdef _VERBOSE
      if (num_mappings > 1) {
	num_dup_reads++;
	num_total_dup_reads += num_mappings;
      }
      #endif
      
      if (num_mappings > 0) {
	num_mapped_reads++;
	if (num_mappings > 1) {
	  num_multihit_reads++;
	}

	for (size_t j = 0; j < num_mappings; j++) {
	  cal = (seed_cal_t *) array_list_get(j, mapping_list);
	  
	  if (read->adapter) {
	    // sequences and cigar
	    len = read->length + abs(read->adapter_length);
	    sequence = (char *) malloc(len + 1);
	    revcomp = (char *) malloc(len + 1);
	    quality = (char *) malloc(len + 1);
	    cigar = cigar_new_empty();

	    if (read->adapter_length < 0) {
	      strcpy(quality, read->adapter_quality);
	      strcat(quality, read->quality);
	    } else {
	      strcpy(quality, read->quality);
	      strcat(quality, read->adapter_quality);
	    }
	    
	    if ( (cal->strand == 1 && 
		  ((read->adapter_strand == 0 && read->adapter_length > 0) || 
		   (read->adapter_strand == 1 && read->adapter_length < 0)))
		 ||
		 (cal->strand == 0 && 
		  ((read->adapter_strand == 0 && read->adapter_length < 0) ||
		   (read->adapter_strand == 1 && read->adapter_length > 0))) ) {
	      strcpy(sequence, read->adapter);
	      strcat(sequence, read->sequence);
	      strcpy(revcomp, read->adapter_revcomp);
	      strcat(revcomp, read->revcomp);
	      
	      cigar_append_op(abs(read->adapter_length), 'S', cigar);
	      cigar_concat(&cal->cigar, cigar);
	    } else {
	      strcpy(sequence, read->sequence);
	      strcat(sequence, read->adapter);
	      strcpy(revcomp, read->revcomp);
	      strcat(revcomp, read->adapter_revcomp);
	      
	      cigar_concat(&cal->cigar, cigar);
	      cigar_append_op(read->adapter_length, 'S', cigar);
	    }
	    sequence[len] = 0; 
	    revcomp[len] = 0; 
	    quality[len] = 0; 
	  } else {
	    // sequences and cigar
	    sequence = read->sequence;
	    revcomp = read->revcomp;
	    quality = read->quality;
	    cigar = &cal->cigar;
	  }

	  if (cal->strand) {
	    flag = 16;
	    seq = revcomp;
	  } else {
	    flag = 0;
	    seq = sequence;
	  }

	  cigar_string = cigar_to_string(cigar);
	  cigar_M_string = cigar_to_M_string(&num_mismatches, &num_cigar_ops, cigar);
	  if (num_mappings > 1) {
	    cal->mapq = 0;
	  }
	  fprintf(out_file, "%s\t%lu\t%s\t%lu\t%i\t%s\t%s\t%lu\t%lu\t%s\t%s\tAS:i:%i\tNM:i:%i\n", 
		  read->id,
		  flag,
		  genome->seq_names[cal->chromosome_id],
		  cal->start + 1,
		  (num_mappings == 1 ? cal->mapq : 0),
		  cigar_M_string,
		  rnext,
		  pnext,
		  tlen,
		  seq,
		  quality,
		  (int) cal->score,
		  num_mismatches
		  );

	  // free memory
	  free(cigar_M_string);
	  free(cigar_string);
	  seed_cal_free(cal);	 
	  if (read->adapter) {
	    free(sequence);
	    free(revcomp);
	    free(quality);
	    cigar_free(cigar);
	  }
	}
      } else {
	num_unmapped_reads++;

	if (read->adapter) {
	  // sequences and cigar
	  len = read->length + abs(read->adapter_length);
	  sequence = (char *) malloc(len + 1);
	  quality = (char *) malloc(len + 1);

	  if (read->adapter_length < 0) {
	    strcpy(quality, read->adapter_quality);
	    strcat(quality, read->quality);
	  } else {
	    strcpy(quality, read->quality);
	    strcat(quality, read->adapter_quality);
	  }
	  
	  if ((read->adapter_strand == 0 && read->adapter_length < 0) || 
	      (read->adapter_strand == 1 && read->adapter_length > 0)) {
	    strcpy(sequence, read->adapter);
	    strcat(sequence, read->sequence);
	  } else {
	    strcpy(sequence, read->sequence);
	    strcat(sequence, read->adapter);
	  }

	  sequence[len] = 0; 
	  quality[len] = 0; 
	} else {
	  // sequences
	  sequence = read->sequence;
	  quality = read->quality;
	}
	
	fprintf(out_file, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", 
		read->id,
		sequence,
		quality
		);

	if (read->adapter) {
	  free(sequence);
	  free(quality);
	}
      }
      
      array_list_free(mapping_list, (void *) NULL);
    } // end for num_reads
  }

  // free memory
  sa_mapping_batch_free(mapping_batch);

  if (wf_batch) sa_wf_batch_free(wf_batch);

  return 0;
}


//--------------------------------------------------------------------
// BAM writer
//--------------------------------------------------------------------

bam_header_t *create_bam_header(sa_genome3_t *genome) {

  bam_header_t *bam_header = (bam_header_t *) calloc(1, sizeof(bam_header_t));

  int num_seqs = genome->num_seqs;
  int target, num_targets = genome->num_refs;

  bam_header->n_targets = num_targets;
  bam_header->target_name = (char **) calloc(num_targets, sizeof(char *));
  bam_header->target_len = (uint32_t*) calloc(num_targets, sizeof(uint32_t));
  target = 0;
  for (int i = 0; i < num_seqs; i++) {
    if (IS_REF_SEQ(genome->seq_types, i)) {
      bam_header->target_name[target] = strdup(genome->seq_names[i]);
      bam_header->target_len[target] = genome->seq_lengths[i];
      target++;
    }
  }

  char pg[1024];
  sprintf(pg, "@PG\tID:HPG-Aligner\tVN:%s\n", HPG_ALIGNER_VERSION);
  bam_header->text = strdup(pg);
  bam_header->l_text = strlen(bam_header->text);

  return bam_header;
}

//--------------------------------------------------------------------

int sa_bam_writer(void *data) {
  sa_wf_batch_t *wf_batch = (sa_wf_batch_t *) data;
  
  sa_mapping_batch_t *mapping_batch = (sa_mapping_batch_t *) wf_batch->mapping_batch;
  if (mapping_batch == NULL) {
    printf("bam_writer1: error, NULL mapping batch\n");
    return 0;
  }

  #ifdef _TIMING
  for (int i = 0; i < NUM_TIMING; i++) {
    func_times[i] += mapping_batch->func_times[i];
  }
  #endif

  int flag, len;
  char *sequence, *quality;

  fastq_read_t *read;
  array_list_t *read_list = mapping_batch->fq_reads;

  bam1_t *bam1;
  alignment_t *alig;
  array_list_t *mapping_list;
  bam_file_t *out_file = wf_batch->writer_input->bam_file;

  sa_genome3_t *genome = wf_batch->sa_index->genome;

  size_t num_reads, num_mappings, num_mate_mappings;
  num_reads = mapping_batch->num_reads;
  for (size_t i = 0; i < num_reads; i++) {
    read = (fastq_read_t *) array_list_get(i, read_list);
    mapping_list = mapping_batch->mapping_lists[i];
    num_mappings = array_list_size(mapping_list);
    num_total_mappings += num_mappings;

    #ifdef _VERBOSE
    if (num_mappings > 1) {
      num_dup_reads++;
      num_total_dup_reads += num_mappings;
    }
    #endif

    if (num_mappings > 0) {
      num_mapped_reads++;
      if (num_mappings > 1) {
	num_multihit_reads++;
      }
      for (size_t j = 0; j < num_mappings; j++) {
	alig = (alignment_t *) array_list_get(j, mapping_list);

	// update alignment
	if (num_mappings > 1) {
	  alig->map_quality = 0;
	} else {
	  alig->map_quality = alig->mapq;
	}

	bam1 = convert_to_bam(alig, 33);
	bam_fwrite(bam1, out_file);
	bam_destroy1(bam1);
	alignment_free(alig);
      }
    } else {
      num_unmapped_reads++;

      if (read->adapter) {
	// sequences and cigar
	len = read->length + abs(read->adapter_length);
	sequence = (char *) malloc(len + 1);
	quality = (char *) malloc(len + 1);

	if (read->adapter_length < 0) {
	  strcpy(quality, read->adapter_quality);
	  strcat(quality, read->quality);
	} else {
	  strcpy(quality, read->quality);
	  strcat(quality, read->adapter_quality);
	}
	
	if ((read->adapter_strand == 0 && read->adapter_length < 0) || 
	    (read->adapter_strand == 1 && read->adapter_length > 0)) {
	  strcpy(sequence, read->adapter);
	  strcat(sequence, read->sequence);
	} else {
	  strcpy(sequence, read->sequence);
	  strcat(sequence, read->adapter);
	}
	sequence[len] = 0; 
	quality[len] = 0; 
      } else {
	// sequences
	sequence = read->sequence;
	quality = read->quality;
      }
      
      alig = alignment_new();       
      alignment_init_single_end(strdup(read->id), sequence, quality,
				0, -1, -1, "", 0, 0, 0, 0, 0, NULL, alig);
      
      bam1 = convert_to_bam(alig, 33);
      bam_fwrite(bam1, out_file);
        
      // free memory
      bam_destroy1(bam1);
      alig->sequence = NULL;
      alig->quality = NULL;
      alig->cigar = NULL;
      alignment_free(alig);
      if (read->adapter) {
	free(sequence);
	free(quality);
      }
    }
    array_list_free(mapping_list, (void *) NULL);
  }

  // free memory
  sa_mapping_batch_free(mapping_batch);

  if (wf_batch) sa_wf_batch_free(wf_batch);

  return 0;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
