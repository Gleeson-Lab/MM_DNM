#!/bin/bash
#$ -cwd
#PBS -l nodes=1:ppn=2
#PBS -l walltime=128:00:00

module load bcftools
bcftools norm -m -any $OutPath/$Output".interval"$INTERVAL".genotyped.vcf.gz" -O z -o $OutPath/$Output".interval"$INTERVAL".genotyped.split.vcf.gz"

module load samtools
tabix -p vcf $OutPath/$Output".interval"$INTERVAL".genotyped.split.vcf.gz"
