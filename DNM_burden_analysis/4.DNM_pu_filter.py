#!/usr/bin/python


## <1> PF1: Remove if DNM are not 0/0 in any other parents

dir_pu = "/projects/ps-gleesonlab8/User/hiyoothere/NTD/6.manual_check/PU_check/output/"

def get_AD(fname, alt):
    f = open(dir_pu + fname, 'r')
    state = 'invalid'
    for line in f:
        if line[0] != '#':
            s = line.strip().split('\t')
            l_alt = s[4].split(',')
            AD = s[-1].split(':')[-1]
            for i in range(len(l_alt)):
                if l_alt[i] == alt:
                    if AD.split(',')[i+1] != '0':
                        state = 'valid'
                    else:
                        state = 'invalid'
    return (state)


f = open("/projects/ps-gleesonlab8/User/hiyoothere/NTD/5.Analysis/DNM/Filteration/NTD_WES_DNM.sampFil.PF2.Insp.den_input", 'r')
out = open("/projects/ps-gleesonlab8/User/hiyoothere/NTD/5.Analysis/DNM/Filteration/NTD_WES_DNM.sampFil.PF2.Insp.pu.den_input", 'w')

cnt_valid = 0
cnt_invalid = 0
cnt_indel = 0

for line in f:
    s = line.strip().split('\t')
    if s[0] == 'locus':
        out.write(line)
    else:
        locus = s[0]
        sample = s[4]
        fname = sample + '_' + locus + '.pu.vcf'
        ref = s[1]
        alt = s[2]
        if len(ref) ==1 and len(alt) == 1: #only for snp
           # print (fname, alt)
            v_state = get_AD(fname, alt)

            if v_state == 'valid':
                out.write(line)
                cnt_valid += 1
            elif v_state == 'invalid':
                cnt_invalid += 1
            else:
                print (sample, locus)
        else:
            out.write(line)
            cnt_indel +=1
print (cnt_valid)
print (cnt_invalid)
print (cnt_indel)


## Check if alt exist in parents of the DNM proband & Binom test with pu allele counts
from scipy import stats
def get_alt_cnt(line):
    s = line.strip().split('\t')
    alt = s[4].split(',')[0]
    AD = s[-1].split(':')[-1]
    alt_cnt = AD.split(',')[1]
    return ([alt, alt_cnt])

def get_ref_cnt(line):
    s = line.strip().split('\t')
    AD = s[-1].split(':')[-1]
    alt_cnt = AD.split(',')[0]
    return ([alt_cnt])

f = open("NTD_WES_DNM.sampFil.PF2.Insp.pu.den_input", 'r')
out = open("NTD_WES_DNM.sampFil.PF2.Insp.pu.Binom0.05.v2.den_input", 'w')
pu_dir="/projects/ps-gleesonlab8/User/hiyoothere/NTD/6.manual_check/PU_check/output/"
pu_dir_parents="/projects/ps-gleesonlab8/User/hiyoothere/NTD/6.manual_check/PU_check/output_parents/"
dic_parents = {}
f_ped = open("/projects/ps-gleesonlab8/User/hiyoothere/NTD/5.Analysis/DNM/ped/NTD_WES_DNM.FixedSample.ped", 'r')
for line in f_ped:
        if 'family' not in line:
                s = line.split('\t')
                pro  = s[1]
                fa = s[2]
                mo = s[3]
                if fa != '0':
                        dic_parents[pro] = [fa, mo]


cnt_valid = 0
cnt_invalid = 0
cnt_indel = 0

def get_alt_cnt(line):
    s = line.strip().split('\t')
    alt = s[4].split(',')[0]
    AD = s[-1].split(':')[-1]
    alt_cnt = AD.split(',')[1]
    return ([alt, alt_cnt])

for line in f:
    s = line.strip().split('\t')
    if s[0] == 'locus':
        out.write(line)
        pass
    else:
        locus = s[0]
        sample = s[4]
        fname = sample + '_' + locus + '.pu.vcf'
        state = 'until_valid'
        #print (locus, sample)
        f_pro = open(pu_dir + fname, 'r')
        rl = f_pro.readlines()[-1]
        pro_alt = get_alt_cnt(rl)
        pro_alt_base = s[2]
        pro_ref_base = s[1]
        if len(pro_alt_base) == 1 and len(pro_ref_base) == 1: #only for SNVs
            pu_ref_cnt = int(get_ref_cnt(rl)[0])
            pu_alt_cnt = int(pro_alt[1])
            if p_binom >= 0.05:
                cnt_valid += 1
                out.write(line)
        else: #If indel
            out.write(line)
            f_fa.close()
            f_mo.close()
            f_pro.close()



