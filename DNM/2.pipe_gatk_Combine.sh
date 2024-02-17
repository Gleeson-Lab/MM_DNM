#!/bin/bash
#$ -cwd
#PBS -l nodes=1:ppn=6
#PBS -l walltime=128:00:00
#PBS -q home
#PBS -l mem=120gb

module load GenomeAnalysisTK
echo $CHR

mkdir -p $OutPath"tmp_SSC1_"$CHR

gatk --java-options "-Xmx90G -Xms90G -XX:ParallelGCThreads=4" CombineGVCFs \
--variant $VARIANT \
--reference $REF \
--intervals chr$CHR \
--tmp-dir $OutPath"tmp_SSC1_"$CHR \
--output $OutPath/$Output".chr"$CHR".vcf.gz"


