#!/bin/bash
#$ -cwd
#PBS -l nodes=1:ppn=4
#PBS -l walltime=128:00:00

java -Xmx4g -jar /home/dantakli/bin/picard-2.20.7/picard.jar MergeVcfs \
I=/projects/ps-gleesonlab8/User/hiyoothere/NTD/WorkFlow/3.With_Combined/forMerge.list \
O=/oasis/tscc/scratch/sak015/NTD_All_re/NTD_All_re.genotyped.vcf.gz
