#ifndef _MPI_WORKER_H
#define _MPI_WORKER_H

#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

#include "commons/workflow_scheduler.h"
#include "bioformats/bam/bam_file.h"

#include "options.h"
#include "batch_writer.h"

#include "tools/bam/aligner/alig.h"
#include "tools/bam/recalibrate/bam_recal_library.h"

#include "sa/sa_index3.h"
#include "sa/sa_search.h"

#include "dna/sa_dna_commons.h"
#include "dna/doscadfun.h"
#include "dna/sa_io_stages.h"
#include "dna/sa_mapper_stage.h"

#include "mpi/mpi_commons.h"

//--------------------------------------------------------------------

void mpi_worker(options_t *options);

//--------------------------------------------------------------------
//--------------------------------------------------------------------
#endif // _MPI_WORKER_H
