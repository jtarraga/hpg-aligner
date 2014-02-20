#!/bin/bash -l 
#$ -N hpg-aligner-dna
#$ -l h_vmem=32G
##$ -pe ompi 1 
#$ -pe ompi-rr 1 
#$ -wd /home/jtarraga/appl/bioinfo-c/hpg-aligner
#$ -o $HOME/queue/stdout
#$ -e $HOME/queue/stderr
#$ -l clinicq=true
##$ -q clinic.q@clinic10,clinic.q@clinic12
#$ -q clinic.q@clinic12

. /etc/profile.d/modules.sh
module load openmpi/1.6.3

echo "Got $NSLOTS processors."
export OMP_NUM_THREADS=12
#mpirun -np $NSLOTS /bin/hostname
#mpirun -np $NSLOTS `pwd`/mpi_program/mpi_hello

# Bind to socket and assign each process to one of them, when launching 6 threads (MPI * OMP) in more than 1 MPI process
#time mpirun --np $NSLOTS --bysocket --bind-to-socket -x OMP_NUM_THREADS --mca orte_base_help_aggregate 0 --mca btl_sm_use_knem 0 /home/jtarraga/appl/bioinfo-c/hpg-aligner/bin/hpg-aligner dna -f /home/jtarraga/datasets/chr21_100nt_r0.01.bwa.read1.fastq -i /home/jtarraga/datasets/sa/index/chrom_20_21_22_k18 --cpu-threads $OMP_NUM_THREADS

# Bind to socket, when launching (less than) 6 threads (MPI * OMP)
#time mpirun --np $NSLOTS --bind-to-socket -x OMP_NUM_THREADS --mca orte_base_help_aggregate 0 --mca btl_sm_use_knem 0 /home/jtarraga/appl/bioinfo-c/hpg-aligner/bin/hpg-aligner dna -f /home/jtarraga/datasets/chr21_100nt_r0.01.bwa.read1.fastq -i /home/jtarraga/datasets/sa/index/chrom_20_21_22_k18 --cpu-threads $OMP_NUM_THREADS

# When launching more than 6 threads in total (MPI * OMP)
#time mpirun --np $NSLOTS --num-cores 12 -x $OMP_NUM_THREADS --mca orte_base_help_aggregate 0 --mca btl_sm_use_knem 0 /home/jtarraga/appl/bioinfo-c/hpg-aligner/bin/hpg-aligner dna -f /home/jtarraga/datasets/chr21_100nt_r0.01.bwa.read1.fastq -i /tmp/sa-index-chr20-21-22/ --cpu-threads $OMP_NUM_THREADS

time mpirun --np $NSLOTS --bysocket --num-cores 12 -x $OMP_NUM_THREADS --mca orte_base_help_aggregate 0 --mca btl_sm_use_knem 0 /home/jtarraga/appl/bioinfo-c/hpg-aligner/bin/hpg-aligner dna -f /home/jtarraga/datasets/chr21_100nt_r0.01.bwa.read1.fastq -i /tmp/sa-index-chr20-21-22/ --cpu-threads $OMP_NUM_THREADS


