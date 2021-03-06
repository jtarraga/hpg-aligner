#include "dna/dna_aligner.h"
#include "rna/rna_aligner.h"
#include "bs/bs_aligner.h"
#include "build-index/index_builder.h"


double emboss_matrix_t = 0.0f, emboss_tracking_t = 0.0f;
double sse_matrix_t = 0.0f, sse_tracking_t = 0.0f;
double sse1_matrix_t = 0.0f, sse1_tracking_t = 0.0f;
double avx_matrix_t = 0.0f, avx_tracking_t = 0.0f;
double avx1_matrix_t = 0.0f, avx1_tracking_t = 0.0f;

//--------------------------------------------------------------------
// constants
//--------------------------------------------------------------------

#define OPTIONS 		       29
#define MIN_ARGC  			5
#define NUM_SECTIONS_STATISTICS 	5
#define NUM_SECTIONS_STATISTICS_SB     21

//--------------------------------------------------------------------
// global variables for timing and capacity meausures
//--------------------------------------------------------------------

basic_statistics_t *basic_st;

pthread_mutex_t mutex_sp;

FILE *fd_log;
size_t junction_id;

size_t total_reads = 0;
size_t reads_no_map = 0;

size_t total_sw = 0;

unsigned char mute;

double time_write = 0;
double time_free = 0;
double time_free_batch = 0;

double time_timer0 = 0;
double time_timer1 = 0;
double time_timer2 = 0;
double time_timer3 = 0;

double time_read_fq   = 0;
double time_read_fq_process   = 0;
double time_read_alig = 0;
double time_read_proc = 0;


char convert_ASCII[128];


size_t search_calls = 0;
size_t insert_calls = 0;
double time_search = 0.0;
double time_insert = 0.0;
pthread_mutex_t mutex_calls;

size_t fd_read_bytes = 0;
size_t fd_total_bytes = 0;

size_t total_reads_ph2 = 0;
size_t reads_ph2 = 0;

int redirect_stdout = 0;
int gziped_fileds = 0;

st_bwt_t st_bwt;
int w1_end;
int w2_end;
int w3_end;

size_t total_reads_w2, total_reads_w3;
size_t reads_w2, reads_w3;

//--------------------------------------------------------------------
// main parameters support
//--------------------------------------------------------------------
int main(int argc, char* argv[]) {
  redirect_stdout = 0;
  gziped_fileds = 0;

  if (!isatty(fileno(stdout))) {
    redirect_stdout = 1;
  }

  st_bwt.multi_alig = 0;
  st_bwt.single_alig = 0;
  st_bwt.total_reads = 0;
  st_bwt.map_bwt = 0;

  st_bwt.map_w1 = 0;
  st_bwt.map_w2 = 0;
  st_bwt.map_w3 = 0;

  st_bwt.tot_sj = 0;
  st_bwt.dif_sj = 0;
  st_bwt.cannonical_sj = 0;
  st_bwt.semi_cannonical_sj = 0;

  pthread_mutex_init(&mutex_sp, NULL);

  //memset(tot_cals, 0, sizeof(int)*50);
  //const char HEADER_FILE[1024] = "Human_NCBI37.hbam\0";

  basic_st = basic_statistics_new();

  // init logs, after parsing the command-line
  // logs will be re-set according to the command-line
  log_level = LOG_FATAL_LEVEL;
  log_verbose = 1;
  log_file = NULL;

  if (argc <= 1) {
    LOG_FATAL("Missing command.\nValid commands are:\n\tdna: to map DNA sequences\n\trna: to map RNA sequences\n\tbuild-sa-index: to create the genome SA index (suffix array).\n\tbuild-bwt-index: to create the genome BWT index.\nUse -h or --help to display hpg-aligner options.\n");
  }

  if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
    usage_cli(DNA_MODE);
  }

  char *command = argv[1];  

  // We need to consume command: {dna | rna | bs | build-index}
  argc -= 1;
  argv += 1;
  
  if(strcmp(command, "dna") != 0 && 
     strcmp(command, "rna") != 0 &&
     strcmp(command, "bs" ) != 0 && 
     strcmp(command, "build-sa-index") != 0 &&
     strcmp(command, "build-bwt-index") != 0) {
    LOG_FATAL("Command unknown.\nValid commands are:\n\tdna: to map DNA sequences\n\trna: to map RNA sequences\n\tbs: to map BS sequences\n\tbuild-sa-index: to create the genome sa index.\n\tbuild-bwt-index: to create the genome bwt index.\nUse -h or --help to display hpg-aligner options.\n");

  }

  if (!strcmp(command, "build-bwt-index") || 
      !strcmp(command, "build-sa-index")) {
      run_index_builder(argc, argv, command);
  }

  // parsing options
  options_t *options = parse_options(argc, argv);

  if (options->adapter) {
    options->adapter_revcomp = strdup(options->adapter);
    seq_reverse_complementary(options->adapter_revcomp, strlen(options->adapter_revcomp));
  }

  // now, we can set logs according to the command-line
  init_log_custom(options->log_level, 1, "hpg-aligner.log", "w");
  LOG_DEBUG_F("Command Mode: %s\n", command);

  //convert ASCII fill 
  
  convert_ASCII['a'] = 'T';
  convert_ASCII['A'] = 'T';

  convert_ASCII['c'] = 'G';
  convert_ASCII['C'] = 'G';

  convert_ASCII['g'] = 'C';
  convert_ASCII['G'] = 'C';

  convert_ASCII['t'] = 'a';
  convert_ASCII['T'] = 'A';

  convert_ASCII['n'] = 'N';
  convert_ASCII['N'] = 'N';

  
  if(strcmp(command, "dna") == 0) { 
    // DNA command
    validate_options(options);
    dna_aligner(options);
  } else if (strcmp(command, "rna") == 0)  { 
    // RNA command
    validate_options(options);
    rna_aligner(options);
  }

  //} else if (strcmp(command, "bs") == 0) {
  // BS commnad
  //printf("Run BS mode...");
  //run_bs_aligner(genome1, genome2, genome, bwt_index1, bwt_index2,
  //                   bwt_optarg, cal_optarg, pair_mng, report_optarg, options);
  //}   

  /*
  if (options->fast_mode) {  
    size_t mapped_reads = total_reads - unmapped_reads;
    printf("TOTAL READS    : %lu\n", total_reads);
    printf("    UNMAPPED READS : %lu (%f%%)\n", unmapped_reads, ((float)unmapped_reads * 100)/((float)total_reads));
    printf("    MAPPED READS   : %lu (%f%%)\n", mapped_reads, ((float)mapped_reads * 100)/((float)total_reads));
    printf("-------------------------------\n");
    printf("    CORRECT READS    : %lu (%f%%)\n", correct_reads,  ((float)correct_reads * 100)/((float)total_reads));
    printf("    INCORRECT READS  : %lu (%f%%)\n", mapped_reads - correct_reads,  ((float)(mapped_reads - correct_reads) * 100)/((float)total_reads));
    

    //printf("Seeds with one Error %i (%0.2f)\n", seeds_1err, seeds_1err*100/total_reads);
  }
  */

  options_free(options);

  return 0;

}


//--------------------------------------------------------------------
