"""
Check the relatedness of the set of samples against the specified PED file.
Specifically we check for cryptic relatedness, relatedness in the expected range
for child/parent pairs, and sufficient relatedness between siblings/replicates.
"""

import os
from common import email_user_with_error
from ped import Pedigree

RELATEDNESS_CHECK_FAILURE_SUBJECT = "Failed relatedness check."


def format_sample_list(list_of_sample_pairs):
    return f"({'), ('.join([', '.join(s) for s in list_of_sample_pairs])})\n"


# checked so as to not email more than once
already_failed = os.path.isfile(snakemake.params.failure_file)
if snakemake.log:
    log_fh = open(str(snakemake.log), "w")
else:
    import sys

    log_fh = sys.stdout
try:
    relatedness_values = {}
    with open(snakemake.input["pairs"]) as relatedness_fh:
        header = relatedness_fh.readline().strip().split("\t")
        for line in relatedness_fh:
            d = dict(zip(header, line.strip().split("\t")))
            sample_a, sample_b = d["#sample_a"], d["sample_b"]
            if sample_a > sample_b:
                sample_a, sample_b = sample_b, sample_a
            relatedness_values[(sample_a, sample_b)] = float(d["relatedness"])
    with open(snakemake.input["ped"]) as ped_fh:
        ped = Pedigree(ped_fh)
    valid_count = 0
    child_parent_errors = []
    # in case of siblings we warn the user if > max first degree relatedness
    # rather than raise an exception because the "siblings" can be twins or
    # replicates
    sibling_warnings = []
    sibling_errors = []
    cryptic_relatedness_errors = []
    # in case of unexpected relatedness within a family we warn the user instead
    # of raising an exception
    intra_family_cryptic_relatedness_warnings = []
    (
        child_parent_pairs,
        sibling_pairs,
        other_intra_family_pairs,
    ) = ped.get_intra_family_relationships()
    for sample_a, sample_b in child_parent_pairs:
        if sample_a > sample_b:
            sample_a, sample_b = sample_b, sample_a
        if (
            snakemake.params.first_degree_relatedness_min
            <= relatedness_values[(sample_a, sample_b)]
            <= snakemake.params.first_degree_relatedness_max
        ):
            valid_count += 1
        else:
            child_parent_errors.append((sample_a, sample_b))
    for sample_a, sample_b in sibling_pairs:
        if sample_a > sample_b:
            sample_a, sample_b = sample_b, sample_a
        r = relatedness_values[(sample_a, sample_b)]
        if r < snakemake.params.first_degree_relatedness_min:
            sibling_errors.append((sample_a, sample_b))
        elif r > snakemake.params.first_degree_relatedness_max:
            sibling_warnings.append((sample_a, sample_b))
        else:
            valid_count += 1
    for sample_a, sample_b in other_intra_family_pairs:
        if sample_a > sample_b:
            sample_a, sample_b = sample_b, sample_a
        if (
            relatedness_values[(sample_a, sample_b)]
            > snakemake.params.intra_family_cryptic_relatedness_threshold
        ):
            intra_family_cryptic_relatedness_warnings.append((sample_a, sample_b))
        else:
            valid_count += 1
    for sample_a, sample_b in ped.get_inter_family_relationships():
        if sample_a > sample_b:
            sample_a, sample_b = sample_b, sample_a
        if (
            relatedness_values[(sample_a, sample_b)]
            > snakemake.params.cryptic_relatedness_threshold
        ):
            cryptic_relatedness_errors.append((sample_a, sample_b))
        else:
            valid_count += 1
    error = False
    if (
        child_parent_errors
        or sibling_errors
        or cryptic_relatedness_errors
        or sibling_warnings
        or intra_family_cryptic_relatedness_warnings
    ):
        messages = []
        if child_parent_errors:
            #error = True
            messages.append(
                "The following pairs of parent-child had outside "
                "the tolerable range of relatedness "
                #f"({snakemake.params.first_degree_relatedness_min}-"
                f"{snakemake.params.first_degree_relatedness_max}): "
                f"{format_sample_list(child_parent_errors)}"
            )
        if sibling_errors:
            #error = True
            messages.append(
                "The following pairs of siblings had less than the minimum "
                "tolerable relatedness "
                f"({snakemake.params.first_degree_relatedness_min}): "
                f"{format_sample_list(sibling_errors)}"
            )
        if cryptic_relatedness_errors:
            #error = True
            messages.append(
                "The following pairs of unrelated individuals had greater than "
                "the maximum tolerable relatedness "
                f"({snakemake.params.cryptic_relatedness_threshold}): "
                f"{format_sample_list(cryptic_relatedness_errors)}"
            )
        if sibling_warnings:
            messages.append(
                "The following pairs of siblings had greater than the maximum "
                "tolerable relatedness "
                f"({snakemake.params.first_degree_relatedness_max}) "
                "(perhaps they're replicates/twins?): "
                f"{format_sample_list(sibling_warnings)}"
            )
    else:
        messages = ["No errors detected."]
    num_samples = ped.get_num_samples()
    num_pairs = int(num_samples * (num_samples - 1) / 2)
    messages.append(f"{valid_count}/{num_pairs} pairs of samples were within tolerable range.")
    message = "\n".join(messages) + "\n"
    log_fh.write(message)
finally:
    if snakemake.log:
        log_fh.close()
if error:
    with open(snakemake.params.failure_file, "w") as failure_fh:
        failure_fh.write(message)
    if not already_failed:
        email_user_with_error(
            snakemake.params.email,
            subject=RELATEDNESS_CHECK_FAILURE_SUBJECT,
            message=message,
        )
    raise ValueError(RELATEDNESS_CHECK_FAILURE_SUBJECT)
# success
if already_failed:
    os.remove(snakemake.params.failure_file)

samples = []
def format_tsv(pair,error_type):
    samples.append(pair[0])
    samples.append(pair[1])
    return '\t'.join([error_type,pair[0],pair[1]]) + '\n'

open(snakemake.output.tsv,'w') as out_file:
    for i in child_parent_errors:
        out_file.write(format_tsv(i,'child-parent')
    for i in sibling_errors:
        out_file.write(format_tsv(i,'sibling_error')
    for i in sibling_warnings:
        out_file.write(format_tsv(i,'sibling_waring')
    for i in cryptic_relatedness_errors:
        out_file.write(format_tsv(i,'interfamily')

open(snakemake.output.list,'w') as out_list:
    samples = list(set(samples))
    for i in samples:
        out_list.write(i + '\n')
