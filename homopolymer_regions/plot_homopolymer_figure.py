#!/usr/bin/env python3
"""
Standalone script to generate median recall plots from TSV input data.

Input TSV format should have columns:
- region_type: e.g., 'all_regions', 'non_homopolymer', 'homopolymer_any', etc.
- display_label: Human-readable label for the region
- variant_type: 'SNV' or 'INDEL'
- recall: Median recall value (0-1)

Usage:
    python plot_median_recall.py input_data.tsv
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse


def plot_grouped_bars(combined_data, ax, title, region_order):
    """
    Create grouped bar chart comparing SNV vs INDEL recall values.

    Args:
        combined_data: DataFrame with columns [region_type, display_label, variant_type, recall]
        ax: Matplotlib axis object
        title: Plot title
        region_order: List defining the order of regions on x-axis
    """
    if combined_data.empty:
        ax.set_title(title, fontsize=14, fontweight='bold', loc='left')
        ax.text(0.5, 0.5, 'No data available', ha='center', va='center', transform=ax.transAxes)
        return

    regions = [r for r in region_order if r in combined_data['region_type'].values]
    n_regions = len(regions)

    if n_regions == 0:
        ax.set_title(title, fontsize=18, fontweight='bold', loc='left')
        ax.text(0.5, 0.5, 'No data available', ha='center', va='center', transform=ax.transAxes)
        return

    x = np.arange(n_regions)
    width = 0.35

    snv_data = combined_data[combined_data['variant_type'] == 'SNV'].copy()
    indel_data = combined_data[combined_data['variant_type'] == 'INDEL'].copy()

    snv_values = []
    indel_values = []
    labels = []

    for region in regions:
        snv_val = snv_data[snv_data['region_type'] == region]['recall'].values
        indel_val = indel_data[indel_data['region_type'] == region]['recall'].values
        snv_values.append(snv_val[0] if len(snv_val) > 0 else 0)
        indel_values.append(indel_val[0] if len(indel_val) > 0 else 0)

        label = snv_data[snv_data['region_type'] == region]['display_label'].values
        if len(label) == 0:
            label = indel_data[indel_data['region_type'] == region]['display_label'].values
        labels.append(label[0] if len(label) > 0 else region)

    bars1 = ax.bar(x - width/2, snv_values, width, label='SNV', color='steelblue', alpha=0.8)
    bars2 = ax.bar(x + width/2, indel_values, width, label='INDEL', color='darkorange', alpha=0.8)

    ax.set_ylabel('Median Recall', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold', loc='left')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=14)
    ax.set_ylim(0, 1)
    ax.legend(fontsize=12)

    # Add value labels on bars
    for i, (snv_val, indel_val) in enumerate(zip(snv_values, indel_values)):
        if snv_val > 0:
            ax.text(i - width/2, snv_val + 0.01, f'{snv_val:.3f}', 
                    ha='center', va='bottom', fontsize=12, fontweight='bold')
        if indel_val > 0:
            ax.text(i + width/2, indel_val + 0.01, f'{indel_val:.3f}', 
                    ha='center', va='bottom', fontsize=12, fontweight='bold')


def main():
    parser = argparse.ArgumentParser(
        description='Generate median recall plots from TSV input',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('input_tsv', help='Input TSV file with columns: region_type, display_label, variant_type, recall')
    parser.add_argument('-o', '--output', default='cohort_median_recall_snv_indel_plot.pdf', 
                        help='Output PDF filename (default: cohort_median_recall_snv_indel_plot.pdf)')
    parser.add_argument('--dpi', type=int, default=600, help='Output DPI (default: 300)')
    parser.add_argument('--figsize', nargs=2, type=float, default=[16, 9], 
                        help='Figure size in inches (width height) (default: 16 6)')

    args = parser.parse_args()

    # Read input data
    print(f"Reading data from {args.input_tsv}")
    df = pd.read_csv(args.input_tsv, sep='\t')

    # Validate required columns
    required_cols = ['region_type', 'display_label', 'variant_type', 'recall']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"Error: Missing required columns: {missing_cols}")
        print(f"Found columns: {list(df.columns)}")
        sys.exit(1)

    print(f"Loaded {len(df)} rows")
    print(f"Variant types: {df['variant_type'].unique()}")
    print(f"Region types: {df['region_type'].unique()}")

    # Define region orders for two panels
    region_order_a = ['all_regions', 'non_homopolymer', 'homopolymer_any']
    region_order_b = ['homopolymer_4to6', 'homopolymer_7to11', 'homopolymer_gt11']

    # Split data for two panels
    panel_a_data = df[df['region_type'].isin(region_order_a)].copy()
    panel_b_data = df[df['region_type'].isin(region_order_b)].copy()

    # Create two-panel plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=tuple(args.figsize))

    # Panel A: General regions
    plot_grouped_bars(panel_a_data, ax1, 'A) Inherited Variants by Region', region_order_a)

    # Panel B: Homopolymer stratification
    plot_grouped_bars(panel_b_data, ax2, 'B) Homopolymer Region Length', region_order_b)

    plt.tight_layout()
    plt.subplots_adjust(top=0.9)
    plt.savefig(args.output, dpi=args.dpi, bbox_inches='tight')
    png_output = args.output.replace('.pdf', '.png')
    plt.savefig(png_output, dpi=args.dpi, bbox_inches='tight')
    plt.close()

    print(f"Plot saved as {args.output}")
    print(f"Plot saved as {png_output}")

    # Save panel data as CSV
    panel_a_csv = args.output.replace('.pdf', '_panel_a_data.csv')
    panel_b_csv = args.output.replace('.pdf', '_panel_b_data.csv')
    panel_a_data.to_csv(panel_a_csv, index=False)
    panel_b_data.to_csv(panel_b_csv, index=False)
    print(f"Panel A data saved as {panel_a_csv}")
    print(f"Panel B data saved as {panel_b_csv}")


if __name__ == '__main__':
    main()
