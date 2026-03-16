#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import argparse
import sys

def plot_results(counts_file, output_file, variant, title):
    """
    Generates a stacked horizontal bar chart from aggregated counts,
    including a summary median bar at the bottom.

    Parameters:
    - counts_file: Path to the aggregated_counts.csv file.
    - output_file: Path to save the generated plot (e.g., results_plot.png).
    - variant: Variant type to filter on (e.g., SNVS, INDELS).
    - title: Title to display on the plot.
    """
    # Read the aggregated counts CSV
    try:
        df = pd.read_csv(counts_file)
    except Exception as e:
        print(f"Error reading {counts_file}: {e}")
        sys.exit(1)

    # Check if necessary columns exist
    required_columns = {'Sample ID', 'Variant', 'TP', 'FP', 'FN'}
    if not required_columns.issubset(df.columns):
        print(f"Input CSV must contain the following columns: {required_columns}")
        sys.exit(1)

    # Filter by variant and set Sample ID as index
    df = df[df['Variant'] == variant].set_index('Sample ID')

    # Rename columns: TP -> 'Both', FP -> 'Nanopore Only', FN -> 'Illumina Only'
    df = df.rename(columns={'TP': 'Both', 'FP': 'Nanopore Only', 'FN': 'Illumina Only'})

    # Keep only the three relevant columns
    df = df[['Both', 'Illumina Only', 'Nanopore Only']]

    # Drop rows where all values are zero
    df = df[(df != 0).any(axis=1)]

    # Sort Sample IDs alphabetically
    df = df.sort_index()

    # Store original counts for labeling
    original_counts = df.copy()

    # Compute median across samples
    median_vals = original_counts.median().to_frame().T
    median_vals.index = ['Median']

    # Append median counts to original_counts
    original_counts = pd.concat([original_counts, median_vals])

    # Normalize to proportions (each row sums to 1)
    prop_df = original_counts.div(original_counts.sum(axis=1), axis=0)

    # Extract proportions for plotting
    plot_df = prop_df.loc[original_counts.index]

    # Plotting
    ax = plot_df.plot(
        kind='barh',
        stacked=True,
        figsize=(18, 8),
        color=['green', 'orange', 'blue']
    )

    plt.title(f"Variant Comparison Results: {title}")
    plt.xlabel("Proportion of Variants")
    plt.ylabel("Sample ID")
    ax.invert_yaxis()

    plt.legend(title="Sequencing Platform", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    # Annotate with original counts, but skip zeroes
    for i, container in enumerate(ax.containers):
        counts = original_counts.iloc[:, i]
        labels = [f'{int(val)}' if val > 0 else '' for val in counts]
        ax.bar_label(container, labels=labels, label_type='center', fontsize=8, color='white')

    # Save the plot
    try:
        plt.savefig(output_file, dpi=600, bbox_inches='tight')
        print(f"Plot successfully saved to {output_file}")
    except Exception as e:
        print(f"Error saving plot to {output_file}: {e}")
        sys.exit(1)

    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Plot ONT, Both, NGS counts from aggregated_counts.csv with summary median bar"
    )
    parser.add_argument("--counts", type=str, required=True, help="Path to aggregated_counts.csv")
    parser.add_argument("--output", type=str, required=True, help="Path to output plot image (e.g., results_plot.png)")
    parser.add_argument("--variant", type=str, default="SNVS", help="Variant: SNVS or INDELS or RECORDS (i.e. ALL)")
    parser.add_argument("--title", type=str, default="TITLE", help="Enter title text")
    args = parser.parse_args()

    plot_results(args.counts, args.output, args.variant, args.title)

if __name__ == "__main__":
    main()
