#!/bin/bash

# run truffle
truffle --vcf ../../0.data/4.joint_vcf/NTD_WES_HC_all.VQSR.CGP.DNV.vcf.gz --cpu 12 --maf 0.01
# extract columns ID1, ID2, IBD1 from truffle output to tsv format
sed -e "s/[ \t]\+/\t/g" truffle.ibd | sed 's/^\t//g' | cut -f 1,2,8 > truffle_maf0.01_IBD1.tsv
# run kinship analysis script to remove kinship failures and contaminations
python truffle_kinship_filter.py -ped ../..//0.data/2.ped/1.Before_Truffle/All_sampleqc.ped  -truffle truffle_maf0.01_IBD1.tsv -IBD1 0.75 -contam 5 -out ../../0.data/2.ped/2.After_Truffle/All_sampleqc_truffle_corrected.ped
