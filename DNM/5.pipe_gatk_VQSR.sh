#!/bin/bash
#$ -cwd
#PBS -l nodes=1:ppn=16
#PBS -l walltime=128:00:00

module load GenomeAnalysisTK
DIR=/projects/ps-gleesonlab8/User/hiyoothere/NTD/WorkFlow/3.With_Combined
REF=$DIR/resources/Homo_sapiens_assembly38.fasta
hapmap=$DIR/resources/hapmap_3.3.hg38.vcf.gz
omni=$DIR/resources/1000G_omni2.5.hg38.vcf.gz
t000G=$DIR/resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz
dbsnp=$DIR/resources/dbsnp_146.hg38.vcf.gz


gatk VariantRecalibrator \
   -R $REF \
   -V /oasis/tscc/scratch/sak015/NTD_All_re/NTD_All_re.genotyped.vcf.gz \
   --resource hapmap,known=false,training=true,truth=true,prior=15.0:$hapmap \
   --resource omni,known=false,training=true,truth=false,prior=12.0:$omni \
   --resource t000G,known=false,training=true,truth=false,prior=10.0:$t000G \
   --resource dbsnp,known=true,training=false,truth=false,prior=2.0:$dbsnp \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
   -mode SNP \
   -O /oasis/tscc/scratch/hiyoothere/NTD/4.Analysis/1.HC/0.gvcf/5.VQSR/snp_recal/joint.snps.recal \
   --tranches-file /oasis/tscc/scratch/hiyoothere/NTD/4.Analysis/1.HC/0.gvcf/5.VQSR/snp_recal/joint.snps.tranches \
   --rscript-file /oasis/tscc/scratch/hiyoothere/NTD/4.Analysis/1.HC/0.gvcf/5.VQSR/snp_recal/joint.snps.output.plots.R


 gatk ApplyVQSR \
   -R $REF \
   -V /oasis/tscc/scratch/sak015/NTD_All_re/NTD_All_re.genotyped.vcf.gz \
   -O /oasis/tscc/scratch/hiyoothere/NTD/4.Analysis/1.HC/0.gvcf/5.VQSR/NTD_All_re.genotyped.recal_snps.vcf.gz \
   --truth-sensitivity-filter-level 99.7 \
   --tranches-file /oasis/tscc/scratch/hiyoothere/NTD/4.Analysis/1.HC/0.gvcf/5.VQSR/snp_recal/joint.snps.tranches \
   --recal-file /oasis/tscc/scratch/hiyoothere/NTD/4.Analysis/1.HC/0.gvcf/5.VQSR/snp_recal/joint.snps.recal \
   -mode SNP
