#!/usr/bin/sh

#PBS -N bwa
#PBS -P as_lkc_sunny.wong
#PBS -l select=1:ncpus=30:mem=100g
#PBS -l walltime=240:00:00
#PBS -q qintel_wfly

source /home/han_xiao/xiaohan/Software/miniconda3/bin/activate
conda activate spatial_gus

ncore=30
OUTDIR="/home/han_xiao/xiaohan/Others/Sunny/03.result/04-Spatial_GUS/03.DNA_bwaAdj"

mkdir -p $OUTDIR && cd $OUTDIR

# Define input file paths
FQ_IN="/home/han_xiao/xiaohan/Others/Sunny/01.data/04-Spatial_GUS/D04707D6.classified.fq"
FA_REF="/home/han_xiao/xiaohan/Others/Sunny/01.data/04-Spatial_GUS/GUS_dna.fna"

# 1. Replace spaces with underscores in FASTQ headers (prevents information truncation in the BAM file)
sed 's/ /_/g' ${FQ_IN} > D04707D6.fixed.fq

# 2. Build the index for the GUS gene reference sequence
bwa index ${FA_REF}

# 3. Align reads with BWA-MEM (assuming single-end reads)
# -t 8: use 30 threads
bwa mem -t 8 -k 13 -A 1 -B 2 -O 5 -T 15 -t $ncore ${FA_REF} D04707D6.fixed.fq > mapped_GUS.sam

# 4. Convert SAM to BAM, sort, and build the index
samtools view -@ $ncore -Sb mapped_GUS.sam | samtools sort -@ $ncore -o mapped_GUS.sorted.bam
samtools index mapped_GUS.sorted.bam

# (Optional) Remove intermediate files
# rm D04707D6.fixed.fq mapped_GUS.sam