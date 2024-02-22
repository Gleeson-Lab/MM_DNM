import argparse

parser = argparse.ArgumentParser(description='extract de novo SVs from a survivor trio vcf (9 genotype columns, 3 tools for each individual, individuals should be in the order: proband, father, mother)')
parser.add_argument('-i', required=True, help='trio vcf, uncompressed')
parser.add_argument('-p', required=True, help='ped file')
parser.add_argument('-o', required=True, help='output de novo vcf, uncompressed')
args=parser.parse_args()

ref_genotypes = ['0/0','0','.','./.']
alt_genotypes = ['0/1','1/1','1','./1']

with open(args.i,'r') as input_file:
    with open(args.o,'w') as output_file:
        for line in input_file:
            if line[0:2] == '##':
                output_file.write(line)
            elif line[0] == '#':
                output_file.write(line)
                order = line[:-1].split('\t')[-3:]
                with open(args.p,'r') as ped:
                    for line in ped:
                        info = line[:-1].split('\t')
                        if info[1] in order and info[5] == '2':
                            proband = order.index(info[1])
                            break
            else:
                if 'DENOVO_FILTER=PASS' not in line[:-1].split('\t')[7]:
                    continue
                elif line[:-1].split('\t')[6] not in ['PASS','.']:
                    continue
                genotypes = [i.split(':')[0] for i in line[:-1].split('\t')[-3:]]
                is_denovo = True
                for i in range(len(genotypes)):
                    if i == proband and genotypes[i] in ref_genotypes:
                        is_denovo = False
                        break
                    elif i != proband and genotypes[i] in alt_genotypes:
                        is_denovo = False
                        break
                if is_denovo:
                    output_file.write(line)
                    print(line)
