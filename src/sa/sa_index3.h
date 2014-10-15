#ifndef SA_INDEX3_H
#define SA_INDEX3_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <assert.h>

#include "sa_tools.h"

//--------------------------------------------------------------------------------------

#define max_uint 4294967295

#define MAX_GENOME_LENGTH 4294967295
#define MAX_NUM_SEQUENCES      16384

#define UNKNOWN_CHROM  15555

//--------------------------------------------------------------------------------------

#define CHROM_SEQ_TYPE    0
#define ALT_SEQ_TYPE      1
#define DECOY_SEQ_TYPE    2
#define SCAFFOLD_SEQ_TYPE 3

#define GET_SEQ_TYPE_NAME(type) ((type) == CHROM_SEQ_TYPE ? "CHROM" : ((type) == ALT_SEQ_TYPE ? "ALT" : ((type) == DECOY_SEQ_TYPE ? "DECOY" : ((type) == SCAFFOLD_SEQ_TYPE ? "SCAFFOLD" : "UNKNOWN"))))

#define IS_REF_SEQ(types, i) ((types)[i] == CHROM_SEQ_TYPE || (types)[i] == SCAFFOLD_SEQ_TYPE)
#define IS_CHROM_SEQ(types, i) ((types)[i] == CHROM_SEQ_TYPE)
#define IS_ALT_SEQ(types, i) ((types)[i] == ALT_SEQ_TYPE)
#define IS_SCAFFOLD_SEQ(types, i) ((types)[i] == SCAFFOLD_SEQ_TYPE)

//--------------------------------------------------------------------------------------

typedef struct alt_names {
  size_t size;
  char **alt_names;
  char **chrom_names;
} alt_names_t;

alt_names_t *alt_names_new(char *alt_filename);
void alt_names_free(alt_names_t *p);

int alt_names_exists(char *alt_name, alt_names_t *p);
char *alt_names_get_chrom_name(char *alt_name, alt_names_t *p);
char *alt_names_display(alt_names_t *p);

//--------------------------------------------------------------------------------------

typedef struct sa_genome3 {
  size_t length;
  size_t num_refs;
  size_t num_seqs;
  size_t num_A;
  size_t num_C;
  size_t num_G;
  size_t num_N;
  size_t num_T;
  int *seq_types;       // CHROM, SCAFFOLD or HAP
  size_t *seq_chroms;   // used by HAP sequences (chromosome)
  size_t *seq_starts;   // used by HAP sequences (start point)
  size_t *seq_ends;    // used by HAP sequences (end point)
  size_t *left_flanks;   // used by HAP sequences (left flank)
  size_t *right_flanks;    // used by HAP sequences (right flank)
  size_t *seq_lengths;
  size_t *seq_offsets;
  char **seq_names;
  char *S;
} sa_genome3_t;

//--------------------------------------------------------------------------------------

sa_genome3_t *read_genome3(char *genome_filename);
sa_genome3_t *read_genome3_ex(char *genome_filename, char *alt_filename, 
			      char *scaffold_filename);

//--------------------------------------------------------------------------------------

static inline sa_genome3_t *sa_genome3_new(size_t length, size_t num_seqs,
					   size_t *seq_lengths,
					   int *seq_types,
					   size_t *seq_chroms,
					   size_t *seq_starts,
					   size_t *seq_ends,
					   char **seq_names, char *S) {
  sa_genome3_t *p = (sa_genome3_t *) calloc(1, sizeof(sa_genome3_t));
  p->length = length;
  p->num_seqs = num_seqs;
  p->num_refs = 0;
  p->seq_lengths = seq_lengths;
  if (num_seqs && seq_lengths) {
    p->seq_offsets = (size_t *) calloc(num_seqs, sizeof(size_t));
    size_t offset = 0;
    for (size_t i = 0; i < num_seqs; i++) {
      if (IS_REF_SEQ(seq_types, i)) {
	p->num_refs++;
      }
      p->seq_offsets[i] = offset;
      offset += seq_lengths[i];
    }
  } else {
    p->seq_offsets = NULL;
  }
  p->seq_names = seq_names;
  p->S = S;

  p->seq_types = seq_types;
  p->seq_chroms = seq_chroms;
  p->seq_starts = seq_starts;
  p->seq_ends = seq_ends;

  // calculate flanks
  size_t *left_flanks = (size_t *) calloc(num_seqs, sizeof(size_t));
  size_t *right_flanks = (size_t *) calloc(num_seqs, sizeof(size_t));
  size_t flank_size;
  char *alt_seq, *chrom_seq;
  for (size_t i = 0; i < num_seqs; i++) {
    if (IS_ALT_SEQ(seq_types, i)) {
      // calculate left flank
      flank_size = 0;
      alt_seq = &S[p->seq_offsets[i]];
      chrom_seq = &S[p->seq_offsets[seq_chroms[i]] + seq_starts[i]];
      while (*alt_seq == *chrom_seq) {
	alt_seq++;
	chrom_seq++;
	flank_size++;
      }
      left_flanks[i] = flank_size;

      // calculate right flank
      flank_size = 0;
      alt_seq = &S[p->seq_offsets[i] + seq_lengths[i] - 1];
      chrom_seq = &S[p->seq_offsets[seq_chroms[i]] + seq_lengths[seq_chroms[i]] - seq_ends[i] - 1];
      while (*alt_seq == *chrom_seq) {
	alt_seq--;
	chrom_seq--;
	flank_size++;
      }
      right_flanks[i] = flank_size;
    }
  }
  p->left_flanks = left_flanks;
  p->right_flanks = right_flanks;

  return p;
}

//--------------------------------------------------------------------------------------

static inline void sa_genome3_free(sa_genome3_t *p) {
  if (p) {
    if (p->seq_lengths) free(p->seq_lengths);
    if (p->seq_offsets) free(p->seq_offsets);
    if (p->seq_names) {
      for (unsigned short int i = 0; i < p->num_seqs; i++) {
	free(p->seq_names[i]);
      }
      free(p->seq_names);
    }
    if (p->S) free(p->S);
    if (p->seq_types) free(p->seq_types);
    if (p->seq_chroms) free(p->seq_chroms);
    if (p->seq_starts) free(p->seq_starts);
    if (p->seq_ends) free(p->seq_ends);
    if (p->left_flanks) free(p->left_flanks);
    if (p->right_flanks) free(p->right_flanks);
    free(p);
  }
}

//--------------------------------------------------------------------------------------

static inline char *sa_genome_get_sequence(unsigned short int chrom, size_t start, size_t end, sa_genome3_t *p) {
  size_t pos, len = end - start + 1;
  char *seq = (char *) malloc((len + 1) * sizeof(char));
  for (size_t i = 0, pos = start + p->seq_offsets[chrom]; i < len; i++, pos++) {
    seq[i] = p->S[pos];
  }
  seq[len] = 0;
  return seq;
}

//--------------------------------------------------------------------------------------

static inline void sa_genome3_set_nt_counters(size_t num_A, size_t num_C, size_t num_G,
					      size_t num_N, size_t num_T, sa_genome3_t *p) {
  if (p) {
    p->num_A = num_A;
    p->num_C = num_C;
    p->num_G = num_G;
    p->num_N = num_N;
    p->num_T = num_T;
  }
}

//--------------------------------------------------------------------------------------

static inline void sa_genome3_display(sa_genome3_t *p) {
  if (!p) return;

  printf("Genome length: %lu\n", p->length);
  printf("Number of refs.: %lu\n", p->num_refs);
  printf("Number of sequences: %lu\n", p->num_seqs);
  for (unsigned short int i = 0; i < p->num_seqs; i++) {
    printf("%u\t%s\t%s\t%lu\t%lu\t%u\t%s\t%lu\t%lu\t%lu\t%lu\n", 
	   i, GET_SEQ_TYPE_NAME(p->seq_types[i]), p->seq_names[i],
	   p->seq_lengths[i], p->seq_offsets[i],
	   p->seq_chroms[i], (IS_ALT_SEQ(p->seq_types, i) ? p->seq_names[p->seq_chroms[i]] : ""), 
	   p->seq_starts[i], p->seq_ends[i],
	   p->left_flanks[i], p->right_flanks[i]);
  }
  printf("Nucleotide counters:\n");
  printf("\tNumber of A: %lu\n", p->num_A);
  printf("\tNumber of C: %lu\n", p->num_C);
  printf("\tNumber of G: %lu\n", p->num_G);
  printf("\tNumber of T: %lu\n", p->num_T);
  printf("\tNumber of N: %lu\n", p->num_N);
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

typedef struct sa_index3 {
  uint k_value;
  uint num_suffixes;
  uint prefix_length;
  uint A_items; // JA_items = A_items
  uint IA_items;
  unsigned short int *CHROM;
  uint *PRE;
  uint *SA;
  uint *A;
  uint *IA;
  unsigned char *JA;
  sa_genome3_t *genome;
} sa_index3_t;

//--------------------------------------------------------------------------------------

void sa_index3_build(char *genome_filename, uint k_value, char *sa_index_dirname);
void sa_index3_build_k18(char *genome_filename, uint k_value, char *sa_index_dirname);
void sa_index3_build_k18_ex(char *genome_filename, char *alt_filename, char *scaffold_filename,
			    uint k_value, char *sa_index_dirname);

//--------------------------------------------------------------------------------------

sa_index3_t *sa_index3_parallel_new(char *sa_index_dirname, int num_threads);
sa_index3_t *sa_index3_new(char *sa_index_dirname);
void sa_index3_free(sa_index3_t *sa_index);

//--------------------------------------------------------------------------------------

static inline void sa_index3_display(sa_index3_t *p) {
  if (!p) return;

  printf("Num. suffixes          : %i\n", p->num_suffixes);
  printf("Prefix table length    : %i\n", p->prefix_length);
  printf("Prefix length (k-value): %i\n", p->k_value);
  printf("A length (JA length)   : %i\n", p->A_items);
  printf("AI length              : %i\n", p->IA_items);

  if (p->genome) sa_genome3_display(p->genome);

}

//--------------------------------------------------------------------------------------

static inline void display_prefix(char *S, int len) {
  for (int i = 0; i < len; i++) {
    printf("%c", S[i]);
  }
}



//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

#endif // SA_INDEX3_H
