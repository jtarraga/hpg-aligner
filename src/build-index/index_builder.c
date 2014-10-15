#include "index_builder.h"


//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

index_options_t *index_options_new() {
  index_options_t *options = (index_options_t *)malloc(sizeof(index_options_t));
  
  options->version = 0;
  options->help = 0;
  options->bs_index = 0;
  options->index_ratio = 0;

  options->ref_genome = NULL;
  options->alt_filename = NULL;
  options->scaffold_filename = NULL;
  options->index_filename = NULL;

  return options;
}


void index_options_free(index_options_t *options) {
  if (options) {
    if (options->ref_genome) { free(options->ref_genome); }
    if (options->alt_filename) { free(options->alt_filename); }
    if (options->scaffold_filename) { free(options->scaffold_filename); }
    if (options->index_filename) { free(options->index_filename); }
    free(options);
  }
}


void** argtable_index_options_new(int mode) {
  int num_options = NUM_INDEX_OPTIONS;
  if (mode == BWT_INDEX) { 
    num_options += NUM_INDEX_BWT_OPTIONS; 
  }
    
  // NUM_OPTIONS +1 to allocate end structure
  void **argtable = (void**)malloc((num_options + 1) * sizeof(void*));	

  int count = 0;
  argtable[count++] = arg_lit0("h", "help", "Help option");
  argtable[count++] = arg_file0("i", "index", NULL, "Index directory name");
  argtable[count++] = arg_file0("g", "ref-genome", NULL, "Reference genome filename");
  argtable[count++] = arg_file0("a", "alternative-map", NULL, "Alternative mapping filename. This two-columns file contains the alternative sequence names with their corresponding chromosome names (only for SA index)");
  argtable[count++] = arg_file0("s", "scaffold-names", NULL, "Scaffold names filename. This one-column file contains the names of the scaffold sequences (only for SA index)");
  argtable[count++] = arg_lit0("v", "version", "Display HPG Aligner version");

  if (mode == BWT_INDEX) {
    argtable[count++] = arg_int0("r", "index-ratio", NULL, "BWT index compression ratio");
  }

  argtable[num_options] = arg_end(count);
  
  return argtable;
}

void argtable_index_options_free(void **argtable, int num_options) {
  if(argtable != NULL) {
    arg_freetable(argtable, num_options);        // struct end must also be freed
    free(argtable);
  }
}

index_options_t *read_CLI_index_options(void **argtable, index_options_t *options, int mode) {        

  int count = -1;
  if (((struct arg_int*)argtable[++count])->count) { options->help = ((struct arg_int*)argtable[count])->count; }
  if (((struct arg_file*)argtable[++count])->count) { options->index_filename = strdup(*(((struct arg_file*)argtable[count])->filename)); }
  if (((struct arg_file*)argtable[++count])->count) { options->ref_genome = strdup(*(((struct arg_file*)argtable[count])->filename)); }
  if (((struct arg_file*)argtable[++count])->count) { options->alt_filename = strdup(*(((struct arg_file*)argtable[count])->filename)); }
  if (((struct arg_file*)argtable[++count])->count) { options->scaffold_filename = strdup(*(((struct arg_file*)argtable[count])->filename)); }
  if (((struct arg_int*)argtable[++count])->count) { options->version = ((struct arg_int*)argtable[count])->count; }
  if (mode == BWT_INDEX) {
    if (((struct arg_int*)argtable[++count])->count) { options->index_ratio = *(((struct arg_int*)argtable[count])->ival); }
  }

  return options;
}


void usage_index(void **argtable, int mode) {
  printf("Usage:\nhpg-aligner {build-sa-index | build-bwt-index}");

  arg_print_syntaxv(stdout, argtable, "\n");
  arg_print_glossary(stdout, argtable, "%-50s\t%s\n");

  exit(0);
}


index_options_t *parse_index_options(int argc, char **argv) {
  int mode, num_options = NUM_INDEX_OPTIONS;

  if (strcmp(argv[0], "build-bwt-index") == 0) {
    mode = BWT_INDEX;
    num_options += NUM_INDEX_BWT_OPTIONS;
  } else if (strcmp(argv[0], "build-sa-index") == 0) {
    mode = SA_INDEX;
  } 

  void **argtable = argtable_index_options_new(mode);
  index_options_t *options = index_options_new();
  if (argc < 2) {
    usage_index(argtable, mode);
    exit(-1);
  } else {
    int num_errors = arg_parse(argc, argv, argtable);

    if (num_errors > 0) {
      fprintf(stdout, "Errors:\n");
      // struct end is always allocated in the last position               
      arg_print_errors(stdout, argtable[num_options], "hpg-aligner");
      usage_index(argtable, mode);
      exit(-1);
    }else {
      options = read_CLI_index_options(argtable, options, mode);
      if(options->help) {
        usage_index(argtable, mode);
        argtable_index_options_free(argtable, num_options);
        index_options_free(options);
        exit(0);
      }
      if (options->version) {
	display_version();
        argtable_index_options_free(argtable, num_options);
        index_options_free(options);
        exit(0);
      }
    }
  }
  
  argtable_index_options_free(argtable, num_options);
  
  return options;
}

void validate_index_options(index_options_t *options, int mode) {
  if (!exists(options->ref_genome)) {
    LOG_FATAL("Reference genome does not exist.\n");
  }
  
  if (!exists(options->index_filename)) {
    LOG_FATAL("Index directory does not exist.\n");
  }

  if (options->alt_filename && !exists(options->alt_filename)) {
    LOG_FATAL("Alternative mapping file does not exist.\n");
  }

  if (options->scaffold_filename && !exists(options->scaffold_filename)) {
    LOG_FATAL("Scaffold names file does not exist.\n");
  }
  
  if (mode == BWT_INDEX && options->index_ratio <= 0) {
    LOG_FATAL("Invalid BWT index ratio. It must be greater than 0.\n");
  }

  if (mode == BWT_INDEX && options->alt_filename) {
    LOG_FATAL("BWT index does not support alternative mapping.\n");
  }

  if (mode == BWT_INDEX && options->scaffold_filename) {
    LOG_FATAL("BWT index does not support scaffolds.\n");
  }
}


//------------------------------------------------------------------------------------

void run_index_builder(int argc, char **argv, char *mode_str) {
  int mode;
  if (!strcmp(mode_str, "build-bwt-index")) {
    mode = BWT_INDEX;
  } else if (!strcmp(mode_str, "build-sa-index")) {
    mode = SA_INDEX;
  }

  index_options_t *options = parse_index_options(argc, argv);

  argtable_index_options_new(mode);
  validate_index_options(options, mode);
  
  if (mode == SA_INDEX) {
    const uint prefix_value = 18;
    char binary_filename[strlen(options->index_filename) + 128];
    sprintf(binary_filename, "%s/dna_compression.bin", options->index_filename);
    printf("Generating SA Index...\n");
    sa_index3_build_k18_ex(options->ref_genome, options->alt_filename, 
			   options->scaffold_filename, prefix_value, 
			   options->index_filename);
    generate_codes(binary_filename, options->ref_genome);
    printf("SA Index generated!\n");

    LOG_DEBUG("Compressing reference genome...\n");
    generate_codes(binary_filename, options->ref_genome);
    LOG_DEBUG("...done !\n");

  } else {
    char binary_filename[strlen(options->index_filename) + 128];
    sprintf(binary_filename, "%s/dna_compression.bin", options->index_filename);
    
    LOG_DEBUG("Compressing reference genome...\n");
    generate_codes(binary_filename, options->ref_genome);
    LOG_DEBUG("...done !\n");
    
    LOG_DEBUG("Building BWT index...\n");
    bwt_generate_index_files(options->ref_genome, options->index_filename, options->index_ratio, false, "ACGT");
    LOG_DEBUG("...done !\n");    
  }

  exit(0);
}


//const uint prefix_value = 18;
//run_index_builder_sa(options->genome_filename, prefix_value, options->bwt_dirname);

//------------------------------------------------------------------------------------

void run_index_builder_bs(char *genome_filename, char *bwt_dirname, 
			  int bwt_ratio, bool duplicate_strand, char *bases) {
      /** **************************************************************************	*
       * 										*
       * Generates the genome transform from the input and builds the index		*
       * 										*
       * The genome transformed are stored in the directory give by the user,		*
       * and the index are stored in subfolders				       		*
       * 										*
       * ***************************************************************************	*/
  /*
      run_index_builder(options->genome_filename, options->bwt_dirname, options->index_ratio, false, "ACGT");

      // generate binary code for original genome
      char binary_filename[strlen(options->bwt_dirname) + 128];
      sprintf(binary_filename, "%s/dna_compression.bin", options->bwt_dirname);
      generate_codes(binary_filename, options->genome_filename);

      char bs_dir1[256];
      sprintf(bs_dir1, "%s/AGT_index", options->bwt_dirname);
      //if (is_directory(bs_dir1) == 0) {
      create_directory(bs_dir1);
      //}

      LOG_DEBUG("Generation of AGT index\n");
      char genome_1[256];
      sprintf(genome_1, "%s/AGT_genome.fa", options->bwt_dirname);
      char gen1[256];
      sprintf(gen1, "sed 's/C/T/g' %s > %s",options->genome_filename, genome_1);
      system(gen1);

      run_index_builder(genome_1, bs_dir1, options->index_ratio, false, "AGT");
      LOG_DEBUG("AGT index Done !!\n");

      LOG_DEBUG("Generation of ACT index\n");
      char bs_dir2[256];
      sprintf(bs_dir2, "%s/ACT_index", options->bwt_dirname);
      //if (is_directory(bs_dir2) == 0) {
      create_directory(bs_dir2);
      //}

      char genome_2[256];
      sprintf(genome_2, "%s/ACT_genome.fa", options->bwt_dirname);
      char gen2[256];
      sprintf(gen2, "sed 's/G/A/g' %s > %s",options->genome_filename, genome_2);
      system(gen2);

      run_index_builder(genome_2, bs_dir2, options->index_ratio, false, "ACT");
      LOG_DEBUG("ACT index Done !!\n");
  */
}
