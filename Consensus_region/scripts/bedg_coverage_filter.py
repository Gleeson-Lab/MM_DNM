import argparse
import numpy as np
import gzip

parser = argparse.ArgumentParser(description='outputs lines where <p> percent of feilds are greater than or equal to <f> coverage')
parser.add_argument('-i', required=True, help='path to input gzipped bedg')
parser.add_argument('-f', required=True, help='coverage value')
parser.add_argument('-p', required=True, help='percentage of values that must be above <f>')
parser.add_argument('-o', required=True, help='path to output bedg')
args=parser.parse_args()

with gzip.open(args.i, 'rt') as bedg:
    with open(args.o, 'w') as output:
        """
        for each row in bedg file calculates the percentage of coverages
        that are >= <f> and keeps the row if the percentage >= <p>
        """
        for line in bedg:
            if line[:-1].split('\t')[1] == 'start':
                output.write(line)
                continue
            info = line[:-1].split('\t')[3:]
            feilds = len(info)
            passed = 0
            for i in info:
                if int(i) >= int(args.f):
                    passed += 1
            if passed/feilds >= float(args.p):
                output.write(line)
