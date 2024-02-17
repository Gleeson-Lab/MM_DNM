#!/bin/bash
#$ -cwd
#PBS -l nodes=1:ppn=6
#PBS -l walltime=400:00:00
#PBS -q home

module load GenomeAnalysisTK
echo $CHR


mkdir -p $OutPath/"tmpGT_"$INTERVAL
gatk --java-options "-Xms32G -Xmx32G -XX:ParallelGCThreads=4" GenotypeGVCFs \
--variant $OutPath/$Output".interval"$INTERVAL".vcf.gz" \
--reference $REF \
--intervals /projects/ps-gleesonlab8/User/hiyoothere/NTD/WorkFlow/2.HC_Combine/intervals/interval_$INTERVAL'.bed' \
--tmp-dir=$OutPath/"tmpGT_"$INTERVAL \
--use-jdk-inflater \
--use-jdk-deflater \
--output /oasis/tscc/scratch/sak015/NTD/$Output".interval"$INTERVAL".genotyped.vcf.gz" 

#--variant $OutPath/$Output".interval"$INTERVAL".vcf.gz" \
