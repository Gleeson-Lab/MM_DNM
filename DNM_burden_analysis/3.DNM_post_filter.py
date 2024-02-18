#!/usr/bin/python

## <1> PF1: Remove if DNM are not 0/0 in any other parents
#1. For FLNC like variants- to remove variants in bad alignment by filtering mother/father observed variants

#1-a extract locus for extraction from gvcf

f_final = open("/projects/ps-gleesonlab8/User/hiyoothere/NTD/5.Analysis/DNM/Filteration/NTD_WES_DNM.sampFil.HeuFil.den_input", 'r')
out_locus = open("/projects/ps-gleesonlab8/User/hiyoothere/NTD/5.Analysis/DNM/Filteration/NTD_WES_DNM.sampFil.HeuFil.locus", 'w')


for line in f_final:
    if 'pro_DP' not in line:
        locus = line.split('\t')[0]
        CHROM = locus.split(':')[0]
        start = locus.split(':')[1]
        out_locus.write(CHROM + '\t' + start + '\t' + start + '\n'

### bcftools view -R NTD_WES_DNM.sampFil.HeuFil.locus NTD_All_re.genotyped.VQSR.CGP.ann.vep.DNM.vcf.gz > NTD_WES_DNM.sampFil.HeuFil.locus.extracted.vcf

#2-c list up if a dnm candidate was shown from any mother or father 
f = open("/projects/ps-gleesonlab8/User/hiyoothere/NTD/5.Analysis/DNM/Filteration/NTD_WES_DNM.sampFil.HeuFil.locus.extracted.vcf", "r")

l_out_align = []


for line in f:
    if line[0] == '#':
        header = f.readline().strip().split('\t')
    else:
        d = dict(zip(header, line.strip().split('\t')))
        l = []
        for each in d:
            if '0/1' in d[each] and each != 'INFO':
                if int(d[each].split(':')[2]) >= 12: #If dp >= 12 
                    #print (each, d[each])
                    l.append(each)
        dic_pro_member = {}
        #print (l)
        if len(l) > 1:
            locus = line.split('\t')[0] + ':' + line.split('\t')[1]
            ped = open("/projects/ps-gleesonlab8/User/hiyoothere/NTD/0.Data/0.SampleInfo/20220919_NTD_WES_Trios.truffle.wH.ped", 'r' )
            h_ped = ped.readline().strip().split('\t')
            for line in ped:
                if line.split('\t')[2] != '0': #if the sample is proband
                    d_ped = dict(zip(h_ped, line.strip().split('\t')))
                Family = d_ped['Family']
                dic_pro_member[Family] = {'locus':locus, 'member':[]}
                if d_ped['Sample'] in l:
                    dic_pro_member[Family]['member'].append('proband')
                if d_ped['Father'] in l:
                    dic_pro_member[Family]['member'].append('father')
                if d_ped['Mother'] in l:
                    dic_pro_member[Family]['member'].append('mother')
        
        for family in dic_pro_member:
            members = dic_pro_member[family]['member']
            locus = dic_pro_member[family]['locus']
            if 'father' in members and 'proband' not in members:
                l_out_align.append(locus)
            elif 'mother' in members and 'proband' not in members:
                l_out_align.append(locus)
  

l_out_align = list(set(l_out_align))
print (len(l_out_align))
#print (l_out_align)

f_final = open("/projects/ps-gleesonlab8/User/hiyoothere/NTD/5.Analysis/DNM/Filteration/NTD_WES_DNM.sampFil.HeuFil.den_input", 'r')


out_filter1 = open("/projects/ps-gleesonlab8/User/hiyoothere/NTD/5.Analysis/DNM/Filteration/NTD_WES_DNM.sampFil.HeuFil.PF1.den_input", 'w')
cnt = 0
for line in f_final:
    if 'pro_DP' in line:
        out_filter1.write(line)
    else:
        locus = line.split('\t')[0]
        #print (locus)
        
        if locus not in l_out_align:
            cnt +=1
            out_filter1.write(line)
out_filter1.close()            
print (cnt)

## <2> PF2: Clustered event of one sample

#2. remove homopolymer

f_pf1 = open("/projects/ps-gleesonlab8/User/hiyoothere/NTD/5.Analysis/DNM/Filteration/NTD_WES_DNM.sampFil.HeuFil.PF1.den_input", 'r')
out_pf2 = open("/projects/ps-gleesonlab8/User/hiyoothere/NTD/5.Analysis/DNM/Filteration/NTD_WES_DNM.sampFil.HeuFil.PF2.den_input", 'w')
#out_pf2 = open("/projects/ps-gleesonlab8/User/hiyoothere/NTD/5.Analysis/DNM/Filteration/NTD_WES_DNM.sampFil.HeuFil.PF2.den_input", 'a')
dic_pro_locus = {}

l_homo = []
for line in f_pf1:
    if 'pro_DP' in line:
        pass
    else:
        s = line.split('\t')
        locus = s[0]
        pro = s[4]
        if pro not in dic_pro_locus:
            dic_pro_locus[pro] = [locus]
        else:
            dic_pro_locus[pro].append(locus)

for pro in dic_pro_locus:
    if len(dic_pro_locus[pro]) > 1:
        l = dic_pro_locus[pro]
        #print (l)
        base_chr = l[0].split(':')[0]
        base_pos = int(l[0].split(':')[1])
        for i in range(1,len(l)):
            if l[i].split(':')[0] != base_chr:
                base_chr = l[i].split(':')[0]
                base_pos = int(l[i].split(':')[1])
            else:
                if int(l[i].split(':')[1]) - int(base_pos) < 10 :
                    #print ("HOMOPOLYMER:", l)
                    l_homo.append(base_chr + ':' + str(base_pos))
                    base_chr = l[i].split(':')[0]
                    base_pos = l[i].split(':')[1]
                    l_homo.append(base_chr + ':' + str(base_pos))
                else:
                    base_chr = l[i].split(':')[0]
                    base_pos = l[i].split(':')[1]
f_pf1.close()
l_homo = list(set(l_homo))
print (len(l_homo))
if 'chr2:226796947' in l_homo:
    print ()
cnt = 0
#print (l_homo)
f_pf1 = open("/projects/ps-gleesonlab8/User/hiyoothere/NTD/5.Analysis/DNM/Filteration/NTD_WES_DNM.sampFil.HeuFil.PF1.den_input", 'r')
for line in f_pf1:
    if 'pro_DP' in line:
        out_pf2.write(line)
    else:
        if line.split('\t')[0] not in l_homo:
            out_pf2.write(line)
            #if line.split('\t')[4] == '6672029182':
             #   print (line)

        else:
            cnt +=1
print ("Removed PF2:", cnt)



