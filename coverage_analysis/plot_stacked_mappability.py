import argparse
import pandas as pd
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(description="Stacked bar chart of mappable/unmappable gene regions")
    parser.add_argument("--input", required=True, help="TSV from compare_summaries_all_columns.py")
    parser.add_argument("--top", type=int, default=10, help="Top N genes to plot")
    parser.add_argument("--direction", choices=["illumina", "nanopore"], required=True, help="Highlight where other platform has higher unmappable")
    parser.add_argument("--output", default="stacked_discrepant_plot.png", help="Output image file")
    parser.add_argument("--gene-list", help="TSV file with one gene name per row (optional)")
    args = parser.parse_args()

    # Load data
    df = pd.read_csv(args.input, sep="\t")

    # Filter by gene list if provided
    if args.gene_list:
        gene_list_df = pd.read_csv(args.gene_list, sep="\t")
        genes_to_keep = gene_list_df["Gene"].str.strip().tolist()
        df = df[df["GENE"].isin(genes_to_keep)]
        print(f"Filtered to {len(df)} genes from gene list")


    # Calculate TOTAL_DARK for Illumina and Nanopore
    df["TOTAL_DARK_Illumina"] = df["%MAPQ0_ONLY_MEDIAN_Illumina"] + df["%ZERO_MEDIAN_Illumina"]
    df["TOTAL_DARK_Nanopore"] = df["%MAPQ0_ONLY_MEDIAN_Nanopore"] + df["%ZERO_MEDIAN_Nanopore"]
    df["delta_TOTAL_DARK"] = df["TOTAL_DARK_Illumina"] - df["TOTAL_DARK_Nanopore"]

    # Filter top hits based on TOTAL_DARK
    if args.direction == "illumina":
        top_df = df[df["delta_TOTAL_DARK"] > 0].nlargest(args.top, "delta_TOTAL_DARK")
    else:
        top_df = df[df["delta_TOTAL_DARK"] < 0].nsmallest(args.top, "delta_TOTAL_DARK")

    # Sort in reverse so worst-performing gene is on top
    top_df = top_df.sort_values(by="delta_TOTAL_DARK", ascending=False)

    # Calculate mappable/unmappable for stacked bars
    genes = top_df["GENE"]
    illum_mapq0 = top_df["%MAPQ0_ONLY_MEDIAN_Illumina"]
    illum_zero = top_df["%ZERO_MEDIAN_Illumina"]
    illum_mappable = 100 - (illum_mapq0 + illum_zero)
    nano_mapq0 = top_df["%MAPQ0_ONLY_MEDIAN_Nanopore"]
    nano_zero = top_df["%ZERO_MEDIAN_Nanopore"]
    nano_mappable = 100 - (nano_mapq0 + nano_zero)

    y_pos = range(len(genes))
    bar_height = 0.35

    fig, ax = plt.subplots(figsize=(12, max(6, len(genes) * 0.6)))

    # Illumina bars
    illum_bars_mappable = ax.barh(y_pos, illum_mappable, height=bar_height, color="dodgerblue", edgecolor='black', label="Illumina: Mappable")
    illum_bars_mapq0 = ax.barh(y_pos, illum_mapq0, left=illum_mappable, height=bar_height, 
                               color="dodgerblue", alpha=0.4, edgecolor='black', label="Illumina: Dark by Mapping Quality")
    illum_bars_zero = ax.barh(y_pos, illum_zero, left=illum_mappable + illum_mapq0, height=bar_height, 
                              color="dodgerblue", hatch="xx", alpha=0.2, edgecolor='black', label="Illumina: Dark by Depth")

    # Nanopore bars
    nano_y_pos = [y + bar_height for y in y_pos]
    nano_bars_mappable = ax.barh(nano_y_pos, nano_mappable, height=bar_height, color="orange", edgecolor='black', label="Nanopore: Mappable")
    nano_bars_mapq0 = ax.barh(nano_y_pos, nano_mapq0, left=nano_mappable, height=bar_height, 
                              color="orange", alpha=0.4, edgecolor='black', label="Nanopore: Dark by Mapping Qaulity")
    nano_bars_zero = ax.barh(nano_y_pos, nano_zero, left=nano_mappable + nano_mapq0, height=bar_height, 
                             color="orange", hatch="xx", alpha=0.2, edgecolor='black', label="Nanopore: Dark by Depth")

    # Aesthetics
    ax.set_yticks([y + bar_height / 2 for y in y_pos])
    ax.set_yticklabels(genes, style='italic')
    ax.set_xlim(0, 100)
    ax.set_xlabel("Percent of Gene Bases")
    ax.set_title(f"Top {args.top} Genes with Greater % Unmappable Bases ({args.direction} underperforms)")
    ax.invert_yaxis()

    # Move legend outside right
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=True)

    plt.tight_layout()
    plt.savefig(args.output, dpi=600, bbox_inches="tight", 
            facecolor='white', edgecolor='none')

    print(f"Saved plot to {args.output}")

if __name__ == "__main__":
    main()