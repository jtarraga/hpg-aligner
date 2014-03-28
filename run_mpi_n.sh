mpirun --np 12 -machinefile mm_host \
       ./bin/hpg-aligner-mpi dna -f ~/datasets/chr21_100nt_r0.01.bwa.read1.fastq -i ~/datasets/sa/index/chrom_20_21_22_k18 --cpu-threads 12 && samtools view -bS hpg_aligner_output/out.sam > hpg_aligner_output/out.bam
