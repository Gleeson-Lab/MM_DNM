import argparse
import numpy as np

parser = argparse.ArgumentParser(description='converts intersected bedg file into bedg file with minimum coverage per position')
parser.add_argument('-i', required=True, help='path to input bedg')
parser.add_argument('-o', required=True, help='path to output bedg')
args=parser.parse_args()

with open(args.i) as bedg:
    with open(args.o, 'w') as output:
        """
        for each row in the input bedg find the minimum coverage value 
        and output a new bedg with the min coverage as the 4th column
        """
        for line in bedg:
            info = line[:-1].split('\t')
            if info[1] == 'start':
                continue
            nums = [int(i) for i in info[3:]]
            m = min(nums)
            new_line = info[:3]
            new_line.append(str(m))
            output.write('\t'.join(new_line) + '\n')
