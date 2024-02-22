import argparse

parser = argparse.ArgumentParser(description='fix sample names in header after survivor merge')
parser.add_argument('-i', required=True, help='survivor merged samples list')
parser.add_argument('-o', required=True, help='output fixed samples list')
args=parser.parse_args()

with open(args.i,'r') as input_file:
    with open(args.o,'w') as output_file:
        for line in input_file:
            if '/' in line:
                line = line[:-1].split('/')[-1].split('.')[0] + '\n'
            output_file.write(line)
