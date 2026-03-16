#!/usr/bin/env python3

import sys
import argparse
import pandas as pd
from csv import DictReader

def parse_som_output(som_output, output_counts):
    """
    Parses the som.py output file to extract TP, FP, FN counts.

    Parameters:
    - som_output: Path to the som.py output file.
    - output_counts: Path to the output counts file.
    """
    # Initialize counts
    counts_dict = {}

    # Read the som.py output
    with open(som_output, 'r') as f:
        reader = DictReader(f)
        for row in reader:
            variant_type = row['type'].upper()
            counts_dict[variant_type] = {"TP": row['tp'], "FP": row['fp'], "FN": row['fn']}
    # Write the counts to the output_counts file
    with open(output_counts, 'w') as out_f:
        out_f.write(f"Variant\tTP\tFP\tFN\n")
        for variant_type, counts in counts_dict.items():
            out_f.write(f"{variant_type}\t{counts['TP']}\t{counts['FP']}\t{counts['FN']}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse som.py output to extract TP, FP, FN counts.")
    parser.add_argument("som_output", help="Path to som.py output file.")
    parser.add_argument("output_counts", help="Path to output counts file.")
    args = parser.parse_args()

    parse_som_output(args.som_output, args.output_counts)
