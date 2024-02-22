import argparse

parser = argparse.ArgumentParser(description='extract SVs shorter than 50,000 bp from survivor trio vcf')
parser.add_argument('-i', required=True, help='trio vcf, uncompressed')
parser.add_argument('-o', required=True, help='output vcf, uncompressed')
args=parser.parse_args()

with open(args.i,'r') as input_file:
    with open(args.o,'w') as output_file:
        for line in input_file:
            if line[0] == '#':
                output_file.write(line)
            else:
                info = line[:-1].split('\t')[7]
                for i in info.split(';'):
                    if 'SVLEN' in i:
                        svlen = int(i.split('=')[1])
                        if svlen < 50000 and svlen > -50000:
                            output_file.write(line)

