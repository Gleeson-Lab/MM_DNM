import argparse

parser = argparse.ArgumentParser(description='filter records that have 0/0 accross all samples')
parser.add_argument('-i', required=True, help='input vcf, uncompressed')
parser.add_argument('-o', required=True, help='output vcf, uncompressed')
args=parser.parse_args()

alt_genotypes = ['0/1','1/1']

with open(args.i,'r') as input_file:
    with open(args.o,'w') as output_file:
        for line in input_file:
            if line[0] == '#':
                output_file.write(line)
            else:
                GTs = [i.split(':')[0] for i in line[:-1].split('\t')[9:]]
                for i in GTs:
                    if '1' in i:
                        output_file.write(line)
                        break
