#!/bin/bash
#$ -cwd
#PBS -l nodes=1:ppn=16
#PBS -l walltime=128:00:00

module load GenomeAnalysisTK
DIR=/projects/ps-gleesonlab8/User/hiyoothere/NTD/WorkFlow/3.With_Combined
REF=$DIR/resources/Homo_sapiens_assembly38.fasta
dbsnp=$DIR/resources/dbsnp_146.hg38.vcf.gz
mills=$DIR/resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
axiom=$DIR/resources/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz

gatk VariantRecalibrator \
   -R $REF \
   -V /oasis/tscc/scratch/sak015/NTD_All_re/NTD_All_re.genotyped.vcf.gz \
   --resource mills,known=false,training=true,truth=true,prior=12.0:$mills \
   --resource axiom,known=false,training=true,truth=false,prior=10.0:$axiom \
   --resource dbsnp,known=true,training=false,truth=false,prior=2.0:$dbsnp \
   -an QD  -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
   -mode INDEL \
   -O /oasis/tscc/scratch/hiyoothere/NTD/4.Analysis/1.HC/0.gvcf/5.VQSR/indel_recal/joint.indels.recal \
   --tranches-file /oasis/tscc/scratch/hiyoothere/NTD/4.Analysis/1.HC/0.gvcf/5.VQSR/indel_recal/joint.indels.tranches \
   --rscript-file /oasis/tscc/scratch/hiyoothere/NTD/4.Analysis/1.HC/0.gvcf/5.VQSR/indel_recal/joint.indels.output.plots.R


gatk ApplyVQSR \
   -R $REF \
   -V /oasis/tscc/scratch/hiyoothere/NTD/4.Analysis/1.HC/0.gvcf/5.VQSR/NTD_All_re.genotyped.recal_snps.vcf.gz \
   -O /oasis/tscc/scratch/hiyoothere/NTD/4.Analysis/1.HC/0.gvcf/5.VQSR/NTD_All_re.genotyped.VQSR.vcf.gz \
   --truth-sensitivity-filter-level 99.0 \
   --tranches-file /oasis/tscc/scratch/hiyoothere/NTD/4.Analysis/1.HC/0.gvcf/5.VQSR/indel_recal/joint.indels.tranches \
   --recal-file /oasis/tscc/scratch/hiyoothere/NTD/4.Analysis/1.HC/0.gvcf/5.VQSR/indel_recal/joint.indels.recal \
   -mode INDEL

