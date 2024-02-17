#!/bin/bash
#$ -cwd
#PBS -l nodes=1:ppn=16
#PBS -l walltime=128:00:00

module load GenomeAnalysisTK

module load samtools


DIR=/projects/ps-gleesonlab8/User/hiyoothere/NTD/WorkFlow/3.With_Combined
hapmap=$DIR/resources/hapmap_3.3.hg38.vcf.gz

gatk --java-options "-Xmx4g" CalculateGenotypePosteriors \
   -V  /oasis/tscc/scratch/hiyoothere/NTD/4.Analysis/1.HC/0.gvcf/5.VQSR/NTD_All_re.genotyped.VQSR.afheader.vcf.gz \
   -O /oasis/tscc/scratch/hiyoothere/NTD/4.Analysis/1.HC/0.gvcf/6.CGP/NTD_All_re.genotyped.VQSR.CGP.chr$CHR'.vcf.gz' \
   -supporting $hapmap \
   -ped /projects/ps-gleesonlab8/User/hiyoothere/NTD/0.Data/0.SampleInfo/20220919_NTD_WES_Trios.ped \
   -L chr$CHR
   
#Collect De novo
/home/hiyoothere/miniconda3/envs/snakemake/bin/gatk --java-options "-Xmx4g" VariantAnnotator \
-L chr$CHR \
--annotation PossibleDeNovo \
--annotation ClippingRankSumTest \
--pedigree /projects/ps-gleesonlab8/User/hiyoothere/NTD/0.Data/0.SampleInfo/20220919_NTD_WES_Trios.ped \
-V /oasis/tscc/scratch/hiyoothere/NTD/4.Analysis/1.HC/0.gvcf/6.CGP/NTD_All_re.genotyped.VQSR.CGP.chr$CHR'.vcf.gz' \
-O /oasis/tscc/scratch/hiyoothere/NTD/4.Analysis/1.HC/0.gvcf/7.Ann/NTD_All_re.genotyped.VQSR.CGP.ann.chr$CHR'.vcf.gz' 

##VEP
DIR=/projects/ps-gleesonlab8/User/hiyoothere/NTD/WorkFlow/3.With_Combined
VEPDIR=/projects/ps-gleesonlab7/gleeson3/resources/vep_106
FASTA=/home/hiyoothere/resource/Homo_sapiens_assembly38.fasta


OUTPUT1=/oasis/tscc/scratch/hiyoothere/NTD/4.Analysis/1.HC/0.gvcf/7.Ann/NTD_All_re.genotyped.VQSR.CGP.ann.chr$CHR'.vcf.gz'
OUTPUT2=/oasis/tscc/scratch/hiyoothere/NTD/4.Analysis/1.HC/0.gvcf/7.Ann/NTD_All_re.genotyped.VQSR.CGP.ann.chr$CHR'.vep.vcf'

# --compress_output bgzip

perl $VEPDIR/ensembl-vep/vep \
    -v -assembly "GRCh38" --everything --terms "SO" \
    --cache \
    --canonical \
    --tsl \
    --dir_cache /projects/ps-gleesonlab7/gleeson3/resources/vep_106 \
    --offline \
    --dir $VEPDIR \
     -s homo_sapiens \
    --format "vcf" \
    --dir_plugins /projects/ps-gleesonlab7/gleeson3/resources/vep_110/plugins \
    --plugin dbNSFP,/projects/ps-gleesonlab7/gleeson3/resources/dbnsfp/dbNSFP4.1a.gz,CADD_phred,MetaSVM_pred,MetaSVM_rankscore,Polyphen2_HDIV_score,SIFT_pred,MPC,MTR,pLI,pLI_values,gnomAD_exomes_AC,gnomAD_exomes_AF,gnomAD_exomes_AN,SplicAI,SpliceRegion,DisGeNET,Phenotypes,Mastermind \
    --fasta $FASTA \
    --vcf \
    -i $OUTPUT1 \
    --force_overwrite \
    -o $OUTPUT2

cat /oasis/tscc/scratch/hiyoothere/NTD/4.Analysis/1.HC/0.gvcf/7.Ann/NTD_All_re.genotyped.VQSR.CGP.ann.chr$CHR'.vep.vcf' | egrep "#|DeNovo" > /oasis/tscc/scratch/hiyoothere/NTD/4.Analysis/1.HC/0.gvcf/8.DNM/NTD_All_re.genotyped.VQSR.CGP.ann.chr$CHR'.vep.DNM.vcf'






