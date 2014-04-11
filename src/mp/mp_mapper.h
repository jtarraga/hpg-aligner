#ifndef _MP_MAPPER_H
#define _MP_MAPPER_H

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

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

//--------------------------------------------------------------------

void mp_mapper(int id, int np, options_t *options);

//--------------------------------------------------------------------
//--------------------------------------------------------------------
#endif // _MP_MAPPER_H
