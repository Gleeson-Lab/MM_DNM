#!/usr/bin/python

import gzip


ped = open("/projects/ps-gleesonlab8/User/arzoo/GMKF_full_ped/complete_peds/20230216_GMKF_NTD_WGS_trios_fixed.ped", 'r')
f_samples = open("/projects/ps-gleesonlab8/User/hiyoothere/NTD/30.WGS/SNV/sample.samplelist", 'r')
l_samples = []
for line in f_samples:
    sample = line.strip()
    l_samples.append(sample)
print (len(l_samples))

                     
#Get trio information from pedigree file
def get_ped():
    ped = open("/projects/ps-gleesonlab8/User/arzoo/GMKF_full_ped/complete_peds/20230216_GMKF_NTD_WGS_trios_fixed.ped", 'r')
    #dic_ped = {'pro': {'mom':'extid_mom', 'dad':'extid_dad', 'family':'id_fam'}}
    h_ped = ["Family", "Sample", "Father", "Mother", "affected", "Sex"]
    for line in ped:
        if line.split('\t')[2] != '0': #if the sample is proband
            d_ped = dict(zip(h_ped, line.strip().split('\t')))
            dic_ped[d_ped['Sample']] = {"dad":d_ped['Father'], "mom":d_ped['Mother'], "family":d_ped["Family"], "sex":d_ped['Sex']}
    return dic_ped

dic_ped = {}
get_ped()
    
def parse_info(line, pro_5, l_header):
    s = line.strip().split("\t")
    info = s[7]
    variant = info.split("CSQ=")[0] #variant information
    variant_feature = variant.split(';') #variant feature
    info_latter = info.split(';')
    QD = 'NA'
    MQ = 'NA'
    for each in variant_feature:
        if 'MQ' in each and "Rank" not in each and "RAW" not in each:
            MQ = each.lstrip("MQ=")
        elif 'QD' in each:
            QD = each.lstrip("QD=")
    #print (info)
    vep = info.split("CSQ=")[1] #Annotaion
    #print ("vep", vep)
    vep = vep.split(';')[0]
    l_ann = vep.split(',')
    vep = l_ann[0]
    l_gene = []
    l_conseq = []
    l_green_header = []
    l_green_value = []
    l_splice_value = []
    gene_state = 0
    sp_state = 0
    gr_state = 0
    for each_info in info_latter:
        if 'green' in each_info:
            l_green_value.append(each_info.split('=')[1])
            gr_state +=1
        if 'SpliceAI' in each_info:
            spliceai = each_info.split('=')[1]
            l_splice_value.append('\t'.join(spliceai.split('|')))
            sp_state +=1
    if len(l_green_value) == 6 and l_green_value != ['.\t.\t.\t.\t.\t.\t.']:
            l_green_value.insert(4, "")
    #        print (l_green_value, len(l_green_value))
    #print (len(l_green_value))
    if sp_state == 0:
        l_splice_value.append(".\t.\t.\t.\t.\t.\t.\t.\t.\t.") #10
    if gr_state == 0:
        l_green_value.append(".\t.\t.\t.\t.\t.\t.") #7
    #print (l_green_value)
    
    for each_ann in l_ann:
        tmp_gene = each_ann.split('|')[3]
        if tmp_gene != "":
            gene_state = 1
               # print (tmp_gene)
            tmp_cons = each_ann.split('|')[1]
            if tmp_gene not in l_gene:
                l_gene.append(tmp_gene)
                l_conseq.append(tmp_cons)
            if "YES" in each_ann:
                vep = each_ann
                break
    if len(l_splice_value) != 7:
        print (l_splice_value, l_green_value)
    if gene_state == 0:
        vep = l_ann[0]
            
            
    #print (l_gene)
    l_vep = vep.split("|")
    d_vep = dict(zip(l_header, l_vep))
    l_format = s[9:]
    d_format = dict(zip(l_h_format, l_format))
    format_seq = s[8].split(':')
    
    locus_1 = s[0] + ':' + s[1]
    ref_2 = s[3]
    alt_3 = s[4]
    deno_info = variant.split(';')
    #deno_info =  variant2.split(';')
    
    for tmp_ann in deno_info[-5:]:
        if "DeNovo" in tmp_ann:
            deno_ann = tmp_ann.split('=')[0]
    #print (pro_5)
    family_4 = get_ped()[pro_5]['family']
    dad_6 = get_ped()[pro_5]['dad']
    mom_7 = get_ped()[pro_5]['mom']
    pro_sex_8 = get_ped()[pro_5]['sex']
    filters = s[6]
    format_pro = d_format[pro_5].split(':')
    format_dad = d_format[dad_6].split(':')
    format_mom = d_format[mom_7].split(':')
    d_format_pro = dict(zip(format_seq, format_pro))
    d_format_dad = dict(zip(format_seq, format_dad))
    d_format_mom = dict(zip(format_seq, format_mom))
    #
    pro_DP_9 = d_format_pro['DP']
    pro_ref_10 = d_format_pro['AD'].split(',')[0]
    pro_alt_11 = d_format_pro['AD'].split(',')[1]
    pro_GT_12 = d_format_pro['GT']
    pro_GQ_13 = d_format_pro['GQ']
    dad_DP_14 = d_format_dad['DP']
    dad_ref_15 = d_format_dad['AD'].split(',')[0]
    dad_alt_16 = d_format_dad['AD'].split(',')[1]
    dad_GT_17 = d_format_dad['GT']
    dad_GQ_18 = d_format_dad['GQ']
    mom_DP_19 = d_format_mom['DP']
    mom_ref_20 = d_format_mom['AD'].split(',')[0]
    mom_alt_21 = d_format_mom['AD'].split(',')[1]
    mom_GT_22 = d_format_mom['GT']
    mom_GQ_23 = d_format_mom['GQ']
    p_denovo_24 = '.'
    confidence_25 = deno_ann
    #From annotation
    l_annotation = []
    for e_ann in l_header:
        ann_value = d_vep[e_ann]
        l_annotation.append(ann_value)
   
    #polyphen_score_38 = highest_pp
    l_new_header = [locus_1, ref_2, alt_3, family_4, pro_5, dad_6, mom_7,pro_sex_8, pro_DP_9,pro_ref_10,
                   pro_alt_11, pro_GT_12,pro_GQ_13,dad_DP_14,dad_ref_15,
                           dad_alt_16,dad_GT_17, dad_GQ_18, mom_DP_19, mom_ref_20,
                          mom_alt_21, mom_GT_22, mom_GQ_23, p_denovo_24, confidence_25]
    l_new_header.extend(l_annotation)
    l_new_header.extend(l_green_value)
    l_new_header.extend(l_splice_value)
    outline = '\t'.join(l_new_header) + '\n'
    
    return outline

cnt = 0
l_header = []
f_vcf = gzip.open('/projects/ps-gleesonlab8/User/hiyoothere/NTD/30.WGS/SNV/joint.multisplit.VQSR.CGP.ann.vep.DNM.spliceAI.annotated.vcf.gz', 'rt')
out = open("/projects/ps-gleesonlab8/User/hiyoothere/NTD/30.WGS/SNV/joint.multisplit.VQSR.CGP.ann.vep.DNM.spliceAI.annotated.den_input", "w")

h_input = ['locus', 'ref', 'alt', 'family', 'pro', 'dad', 'mom','pro_sex', 'pro_DP','pro_ref', 
                           'pro_alt', 'pro_GT','pro_GQ','dad_DP','dad_ref',
                           'dad_alt','dad_GT', 'dad_GQ', 'mom_DP', 'mom_ref',
                           'mom_alt', 'mom_GT', 'mom_GQ', 'p_denovo', 'confidence']


for line in f_vcf:
    if line[0] == '#':
        if "Consequence annotations from Ensembl VEP." in line:
            header = line.strip().split()[-1][:-2]
            l_header = header.split('|')
            h_input.extend(l_header)
            h_input.extend(["greendb_id", "greendb_stdtype", "greendb_dbsource", "greendb_genes", "greendb_constraint", "greendb_level", "greendb_more_support"])
            h_input.extend(["ALLELE", "SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL"])
            h_input.extend(['occurence', 'probands'])
            out.write('\t'.join(h_input) + '\n')
           # print (len(l_header))
            print (l_header)
        elif "#CHROM" in line:
            l_h_format = line.strip().split('\t')[9:]
            
            
    else:
        #print (line)
        s = line.strip().split("\t")
        info = s[7]
        if "DeNovo" in info and "CSQ" in info:
            #print (info)
            variant = info.split(";CSQ=")[0]
            pro_5 = variant.split("DeNovo=")[1]
            #print (pro_5)
            if ',' in pro_5:
                l_pro = pro_5.split(',')
                for each_pro in l_pro:
                    if each_pro in l_samples:
                        out.write(parse_info(line, each_pro, l_header).rstrip('\n') + '\t' + str(len(l_pro)) + '\t' + '/'.join(l_pro) + '\n' ) 
                        cnt +=1
                    
            else:
                ######## why only sample detected among non-final sample list
                if pro_5 in l_samples:
                #print (pro_5, info)
                    out.write(parse_info(line, pro_5, l_header).rstrip('\n') + '\t' + str(1) + '\t' + pro_5 + '\n' ) 
                    cnt += 1
            
        else:
            pass


f_vcf.close()
print (cnt)

