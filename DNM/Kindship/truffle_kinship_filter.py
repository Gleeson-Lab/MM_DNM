#!/usr/bin/env python
# coding: utf-8

# import packages
import pandas as pd
import argparse

# add agruements
parser = argparse.ArgumentParser(description='create a new ped file removing families with kinship failure and contamination using IBD1 values calulated by truffle')
parser.add_argument('-truffle', required=True, help='path to truffle results in tsv format with columns ID1, ID2, and IBD1')
parser.add_argument('-ped', required=True, help='path to ped file (no header)')
parser.add_argument('-IBD1', type=float, required=True, help='IBD1 threshold for relatedness, remove families where parent-child relationships are below <IBD1> threshold')
parser.add_argument('-contam', type=int, required=True, help='contamination cutoff, remove families where samples have at least <contam> relationships outside of the family')
parser.add_argument('-out', required=True, help='output path for corrected ped file')
args=parser.parse_args()

ped_families = {}
truffle_IBD1_related = []
truffle_IBD1 = {}
ped_pairs = []
ped_samples = []
families = set()
children = []
fathers = []
mothers = []
problematic_families = []

# read in ped and gather lists for mothers, fathers, children, and parent-child pairs
with open(args.ped) as ped:
    for line in ped:
        info = line[:-1].split('\t')
        ped_samples.append(info[1])
        ped_families[info[1]] = info[0]
        families.add(info[0])
        if info[5] == '2' and info[2] != '0' and info[3] != '0':
            children.append(info[1])
            fathers.append(info[2])
            mothers.append(info[3])
            ped_pair1 = sorted([info[1],info[2]])
            ped_pair2 = sorted([info[1],info[3]])
            ped_pairs.append(ped_pair1)
            ped_pairs.append(ped_pair2)

# read tsv with truffle IBD1 output and gather all pairs where IBD1 is greater than 0.75 into truffle_pairs 
# gather all pairs of IBD1 values in truffle_IBD1
with open(args.truffle) as truffle:
    for line in truffle:
        truffle_pair = sorted(line[:-1].split('\t')[:2])
        if truffle_pair[0] not in ped_samples or truffle_pair[1] not in ped_samples:
            continue
        IBD1 = float(line[:-1].split('\t')[2])
        if IBD1 >= args.IBD1:
            truffle_IBD1_related.append(truffle_pair)
        truffle_IBD1[str(truffle_pair)] = IBD1

        
        
# 1. filter out families with kinship failures (parent-child pairs below the IBD1 threshold)
# create a data frame with all the information gathered for kinship failures
data = {
    "sample1": [],
    "sample2": [],
    "family1": [],
    "family2": [],
    "relation1": [],
    "relation2": [],
    'IBD1': []
}
for i in ped_pairs:
    if i not in truffle_IBD1_related:
        if ped_families[i[0]] > ped_families[i[1]]:
            i.reverse()
        sample1 = i[0]
        sample2 = i[1]
        family1 = ped_families[sample1]
        family2 = ped_families[sample2]
        if sample1 in children:
            relation1 = "proband"
        elif sample1 in fathers:
            relation1 = "father"
        elif sample1 in mothers:
            relation1 = "mother"
        if sample2 in children:
            relation2 = "proband"
        elif sample2 in fathers:
            relation2 = "father"
        elif sample2 in mothers:
            relation2 = "mother"
        data["sample1"].append(sample1)
        data["sample2"].append(sample2)
        data["family1"].append(family1)
        data["family2"].append(family2)
        data["relation1"].append(relation1)
        data["relation2"].append(relation2)
        i = sorted(i)
        if str(i) in truffle_IBD1.keys():
            data['IBD1'].append(truffle_IBD1[str(i)])
        else:
            data['IBD1'].append(-1.0)

df = pd.DataFrame(data=data)
df = df.sort_values(by='family1')

# add all kinship failures to the problematic families list
problematic_families = sorted(list(set(df['family1'])))

# print families with kinship failures
print("families with kinship failures")
print(sorted(set(df['family1'])))



# 2. filter out contaminated families (where samples have at least <contam> relationships outside of the family)
# create a data frame with all the information gathered for potentially contaminated samples
data2 = {
    "sample1": [],
    "sample2": [],
    "family1": [],
    "family2": [],
    "relation1": [],
    "relation2": [],
    'IBD1': []
}

for i in truffle_IBD1_related:
    if i not in ped_pairs:
        i = sorted(i)
        sample1 = i[0]
        sample2 = i[1]
        family1 = ped_families[sample1]
        family2 = ped_families[sample2]
        if family1 in problematic_families or family2 in problematic_families:
            continue
        if sample1 in children:
            relation1 = "proband"
        elif sample1 in fathers:
            relation1 = "father"
        elif sample1 in mothers:
            relation1 = "mother"
        if sample2 in children:
            relation2 = "proband"
        elif sample2 in fathers:
            relation2 = "father"
        elif sample2 in mothers:
            relation2 = "mother"
        data2["sample1"].append(sample1)
        data2["sample2"].append(sample2)
        data2["family1"].append(family1)
        data2["family2"].append(family2)
        data2["relation1"].append(relation1)
        data2["relation2"].append(relation2)
        if str(i) in truffle_IBD1.keys():
            data2['IBD1'].append(truffle_IBD1[str(i)])
        else:
            data2['IBD1'].append(-1.0)

df2 = pd.DataFrame(data=data2)
df2 = df2.sort_values(by='family1')

# count number of relationships outside of family per sample and create a sorted dataframe 
problems_outside_dict = {
    'samples': [],
    'values': [],
    'family': []
}
problems_outside = {}
problems_outside_samples = list(df2['sample1']) + list(df2['sample2'])
for i in problems_outside_samples:
    if i in problems_outside:
        problems_outside[i] += 1
    else:
        problems_outside[i] = 1
for i in problems_outside.keys():
    problems_outside_dict['samples'].append(i)
    problems_outside_dict['family'].append(ped_families[i])
    problems_outside_dict['values'].append(problems_outside[i])
problems_outside_df = pd.DataFrame(data=problems_outside_dict).sort_values(by='values',ascending=False)

# print contaminated samples
print("families with " + str(args.contam) + " or more relationships to other families")
print(set(problems_outside_df[problems_outside_df['values'] >= args.contam]['family']))

# update problematic_families list with contaminated samples
problematic_families += list(problems_outside_df[problems_outside_df['values'] >= args.contam]['family'])
problematic_families = sorted(set(problematic_families))



# write a new ped with kinship failures and contaminated samples removed
with open(args.ped) as ped:
    with open(args.out,'w') as ped_out:
        for line in ped:
            info = line[:-1].split('\t')
            if info[0] in problematic_families:
                continue
            else:
                ped_out.write(line)

# write problematic families to a file in the same directory as output
truffle_error_path = "/".join(args.out.split('/')[:-1]) + "/truffle_problem_families.txt"
with open(truffle_error_path,'w') as truffle_problems:
    for i in problematic_families:
        truffle_problems.write(i + '\n')
