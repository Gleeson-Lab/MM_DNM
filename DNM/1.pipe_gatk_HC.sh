#!/bin/bash
#$ -cwd
#PBS -l nodes=1:ppn=1
#PBS -l walltime=46:00:00

module load GenomeAnalysisTK

gatk --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=4" HaplotypeCaller \
-R $REF \
-I $DataPath/$Data \
-O $OutPath/$IDX'.HC.gvcf.gz' \
-ERC GVCF


