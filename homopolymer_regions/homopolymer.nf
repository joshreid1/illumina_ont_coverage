#!/usr/bin/env nextflow

/*
 * Consolidated Multi-Family Inherited Variant Analysis and Som.py Comparison Pipeline
 * Nextflow DSL2 Implementation - SNV and INDEL Processing
 */

nextflow.enable.dsl = 2

/*
 * Define pipeline parameters
 */
// Pipeline 1 parameters
params.input_vcf = "joint_germline_pass.vcf.gz"
params.sample_csv = "sample_ids.csv"
params.stratifications_dir = "giab/release/genome-stratifications/v3.1/GRCh38/LowComplexity/"
params.ref_fasta = "Homo_sapiens_assembly38.fasta"

// Pipeline 2 parameters  
params.lrs_csv = 'affected_rows_snv_high_coverage.csv'
params.ref_ont = 'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'
params.container = 'hap.py.sif'

// Common parameters
params.outdir = "results"
params.force = false
params.run_sompy = false  // Set to true to run som.py analysis

/*
 * Homopolymer stratification files
 */
def homopolymer_bed = [
    "GRCh38_SimpleRepeat_homopolymer_4to6_slop5.bed.gz",
    "GRCh38_SimpleRepeat_homopolymer_7to11_slop5.bed.gz", 
    "GRCh38_SimpleRepeat_homopolymer_gt11_slop5.bed.gz"
]

// Help message
def helpMessage() {
    log.info"""
    Consolidated Multi-Family Inherited Variant Analysis and Som.py Comparison Pipeline
    ==================================================================================

    Usage:
        nextflow run consolidated_pipeline.nf [options]

    Required Parameters (Pipeline 1):
        --input_vcf         Input VCF file with joint-called family data
        --sample_csv        CSV file with family information (PROBAND_ID,SIBLING_ID,MATERNAL_ID,PATERNAL_ID)
        --stratifications_dir   GIAB stratification directory

    Optional Parameters (Pipeline 2):
        --lrs_csv           CSV file with NGS samples and LRS VCF paths (NGS_ID,LRS_SNV_VCF)
        --run_sompy         Enable som.py comparison analysis [default: false]
        --container         Som.py container path [default: ${params.container}]
        --ref_ont           ONT reference FASTA [default: ${params.ref_ont}]

    Common Parameters:
        --ref_fasta         Reference FASTA file [default: ${params.ref_fasta}]
        --outdir            Output directory [default: ${params.outdir}]
        --force             Force restart, overwrite existing files [default: false]

    Examples:
        # Run only variant analysis
        nextflow run consolidated_pipeline.nf \\
            --input_vcf variants.vcf.gz \\
            --sample_csv families.csv \\
            --stratifications_dir /path/to/GRCh38/

        # Run both variant analysis and som.py comparison
        nextflow run consolidated_pipeline.nf \\
            --input_vcf variants.vcf.gz \\
            --sample_csv families.csv \\
            --stratifications_dir /path/to/GRCh38/ \\
            --lrs_csv lrs_samples.csv \\
            --run_sompy true
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

'''
// Validate required parameters
if (!params.input_vcf || !params.sample_csv || !params.stratifications_dir) {
    error "Missing required parameters for Pipeline 1. Use --help for usage information."
}
'''

if (params.run_sompy && !params.lrs_csv) {
    error "Missing --lrs_csv parameter required for som.py analysis. Use --help for usage information."
}

/*
 * PIPELINE 1 PROCESSES: Variant Analysis
 */

process PARSE_FAMILIES {
    tag "Parsing families"

    input:
    path sample_csv

    output:
    path "families.txt"

    script:
    """
    #!/bin/bash
    tail -n +2 ${sample_csv} | while IFS=',' read -r proband_id sibling_id maternal_id paternal_id; do
        if [[ -z "\$proband_id" ]]; then
            continue
        fi
        proband_id=\$(echo "\$proband_id" | xargs)
        maternal_id=\$(echo "\$maternal_id" | xargs)
        paternal_id=\$(echo "\$paternal_id" | xargs)
        proband_vcf="\${proband_id}_\${proband_id}"
        maternal_vcf="\${maternal_id}_\${maternal_id}"
        paternal_vcf="\${paternal_id}_\${paternal_id}"
        echo "\$proband_id,\$maternal_id,\$paternal_id,\$proband_vcf,\$maternal_vcf,\$paternal_vcf"
    done > families.txt
    """
}

process SUBSET_TRIO {
    tag "${family_meta.proband_id}"

    input:
    path input_vcf
    path input_vcf_index
    val family_meta

    output:
    tuple val(family_meta), path("${family_meta.proband_id}_trio.vcf.gz"), path("${family_meta.proband_id}_trio.vcf.gz.tbi")

    script:
    """
    trio_samples="${family_meta.proband_vcf},${family_meta.maternal_vcf},${family_meta.paternal_vcf}"
    bcftools view -s "\$trio_samples" -O z -o "${family_meta.proband_id}_trio.vcf.gz" ${input_vcf}
    tabix "${family_meta.proband_id}_trio.vcf.gz"
    """
}

process NORMALIZE_VARIANTS {
    tag "${family_meta.proband_id}"

    input:
    tuple val(family_meta), path(vcf), path(vcf_index)
    path ref_fasta

    output:
    tuple val(family_meta), path("${family_meta.proband_id}_trio_AC_norm.vcf.gz"), path("${family_meta.proband_id}_trio_AC_norm.vcf.gz.tbi")

    script:
    """
    bcftools norm -f ${ref_fasta} -m -any -O z -o "${family_meta.proband_id}_trio_AC_norm.vcf.gz" ${vcf}
    tabix "${family_meta.proband_id}_trio_AC_norm.vcf.gz"
    """
}

process FILTER_AC {
    tag "${family_meta.proband_id}"

    input:
    tuple val(family_meta), path(norm_vcf), path(norm_vcf_index)

    output:
    tuple val(family_meta), path("${family_meta.proband_id}_trio_AC_final.vcf.gz"), path("${family_meta.proband_id}_trio_AC_final.vcf.gz.tbi")

    script:
    """
    bcftools view -i 'AC>=0' -O z -o "${family_meta.proband_id}_trio_AC_final.vcf.gz" ${norm_vcf}
    tabix "${family_meta.proband_id}_trio_AC_final.vcf.gz"
    """
}

process MISSING2REF {
    tag "${family_meta.proband_id}"

    input:
    tuple val(family_meta), path(final_AC_vcf), path(final_AC_vcf_index)

    output:
    tuple val(family_meta), path("${family_meta.proband_id}_trio_AC_missing2ref.vcf.gz"), path("${family_meta.proband_id}_trio_AC_missing2ref.vcf.gz.tbi")

    script:
    """
    bcftools +missing2ref ${final_AC_vcf} -O z -o "${family_meta.proband_id}_trio_AC_missing2ref.vcf.gz"
    tabix "${family_meta.proband_id}_trio_AC_missing2ref.vcf.gz"
    """
}

process INHERITANCE_FILTER {
    tag "${family_meta.proband_id}"
    publishDir "${params.outdir}/${family_meta.proband_id}", mode: 'copy'

    module 'gatk'

    input:
    tuple val(family_meta), path(missing2ref_vcf), path(missing2ref_vcf_index)

    output:
    tuple val(family_meta), path("${family_meta.proband_id}_inherited_variants.vcf.gz"), path("${family_meta.proband_id}_inherited_variants.vcf.gz.tbi")

    script:
    """
    proband_idx=\$(bcftools query -l ${missing2ref_vcf} | grep -n "^${family_meta.proband_vcf}\$" | cut -d: -f1)
    maternal_idx=\$(bcftools query -l ${missing2ref_vcf} | grep -n "^${family_meta.maternal_vcf}\$" | cut -d: -f1)
    paternal_idx=\$(bcftools query -l ${missing2ref_vcf} | grep -n "^${family_meta.paternal_vcf}\$" | cut -d: -f1)

    proband_idx=\$((proband_idx - 1))
    maternal_idx=\$((maternal_idx - 1))
    paternal_idx=\$((paternal_idx - 1))

    inheritance_expr="GT[\${proband_idx}]='alt' && (GT[\${maternal_idx}]='alt' || GT[\${paternal_idx}]='alt')"

    cat > "${family_meta.proband_id}_inheritance_filter_info.txt" << EOF
Family: ${family_meta.proband_id}
Proband: ${family_meta.proband_id} -> ${family_meta.proband_vcf} (index: \${proband_idx})
Mother:  ${family_meta.maternal_id} -> ${family_meta.maternal_vcf} (index: \${maternal_idx})
Father:  ${family_meta.paternal_id} -> ${family_meta.paternal_vcf} (index: \${paternal_idx})
Filter: \${inheritance_expr}
EOF

    bcftools view -i "\${inheritance_expr} && FORMAT/AD[\${proband_idx}:1]/FORMAT/DP[\${proband_idx}] >= 0.20 && FORMAT/GQ[\${proband_idx}] >= 20 && FORMAT/DP[\${proband_idx}] >= 5" \
        -O z -o "${family_meta.proband_id}_inherited_variants.vcf.gz" \
        ${missing2ref_vcf}

    tabix "${family_meta.proband_id}_inherited_variants.vcf.gz"
    """
}

process FILTER_VARIANTS {
    tag "${family_meta.proband_id}"
    publishDir "${params.outdir}/${family_meta.proband_id}", mode: 'copy'

	module 'gatk'

    input:
    tuple val(family_meta), path(inherited_vcf), path(inherited_vcf_index)

    output:
    tuple val(family_meta), path("${family_meta.proband_id}_inherited_snvs.vcf.gz"), path("${family_meta.proband_id}_inherited_snvs.vcf.gz.tbi"), path("${family_meta.proband_id}_inherited_indels_filtered.vcf.gz"), path("${family_meta.proband_id}_inherited_indels_filtered.vcf.gz.tbi")

    script:
    """
    # Filter SNVs
    bcftools view -v snps -O z -o "${family_meta.proband_id}_inherited_snvs.vcf.gz" ${inherited_vcf}
    tabix "${family_meta.proband_id}_inherited_snvs.vcf.gz"

    # Filter INDELs
    bcftools view -v indels -O z -o "${family_meta.proband_id}_inherited_indels.vcf.gz" ${inherited_vcf}
    tabix "${family_meta.proband_id}_inherited_indels.vcf.gz"

	gatk SelectVariants \
	-V "${family_meta.proband_id}_inherited_indels.vcf.gz" \
	--max-indel-size 50 \
	-O "${family_meta.proband_id}_inherited_indels_filtered.vcf"

	bgzip "${family_meta.proband_id}_inherited_indels_filtered.vcf"
	tabix -f "${family_meta.proband_id}_inherited_indels_filtered.vcf.gz"
    """
}

process CREATE_HP_UNION {
    tag "Creating HP union"

    module 'bedtools'

    input:
    path stratifications_dir

    output:
    path "all_homopolymers_merged.bed"

    script:
    def hp_files = homopolymer_bed.join(' ')
    """
    for hp_file in ${hp_files}; do
        hp_path="${stratifications_dir}/\${hp_file}"
        if [ -f "\$hp_path" ]; then
            if [[ "\$hp_path" == *.gz ]]; then
                zcat "\$hp_path" >> all_homopolymers_union.bed
            else
                cat "\$hp_path" >> all_homopolymers_union.bed
            fi
        fi
    done

    sort -k1,1 -k2,2n all_homopolymers_union.bed > all_homopolymers_sorted.bed
    bedtools merge -i all_homopolymers_sorted.bed > all_homopolymers_merged.bed
    """
}

process HP_INTERSECT_SNV {
    tag "${family_meta.proband_id} SNV"
    publishDir "${params.outdir}/${family_meta.proband_id}", mode: 'copy'

    module 'bedtools'
    cpus 2
    memory '32 GB'
    time '2 h'

    input:
    tuple val(family_meta), path(snv_vcf), path(snv_vcf_index)
    path stratifications_dir
    path hp_union_bed

    output:
    tuple val(family_meta), path("*_snv*.vcf.gz"), path("*_snv*.vcf.gz.tbi")

    script:
    // Reverse order for prioritizing largest homopolymers first
    def hp_order = homopolymer_bed.reverse().join(' ')
    """
    mkdir -p homopolymers_snv

    # Create union (any homopolymer) for snvs
    bedtools intersect -a ${snv_vcf} -b ${hp_union_bed} -header -wa | \
    bcftools view -O z -o "${family_meta.proband_id}_inherited_snv_homopolymer_any.vcf.gz"

    if [ -s "${family_meta.proband_id}_inherited_snv_homopolymer_any.vcf.gz" ]; then
        tabix "${family_meta.proband_id}_inherited_snv_homopolymer_any.vcf.gz"
    fi

    # Create non-homopolymer snv VCF as variants NOT in homopolymer_any
    bedtools intersect -a  ${snv_vcf} -b "${family_meta.proband_id}_inherited_snv_homopolymer_any.vcf.gz" -header -v | 
    bcftools view -O z -o "${family_meta.proband_id}_inherited_snv_non_homopolymer.vcf.gz"
   
    if [ -s "${family_meta.proband_id}_inherited_snv_non_homopolymer.vcf.gz" ]; then
        tabix "${family_meta.proband_id}_inherited_snv_non_homopolymer.vcf.gz"
    else
        rm -f "${family_meta.proband_id}_inherited_snv_non_homopolymer.vcf.gz"
    fi

    prev_vcfs=""

    for hp_file in ${hp_order}; do
        hp_path="${stratifications_dir}/\${hp_file}"
        if [ ! -f "\$hp_path" ]; then
            continue
        fi

        category=\$(echo "\$hp_file" | sed 's/GRCh38_SimpleRepeat_homopolymer_//' | sed 's/_slop5\\.bed\\.gz\$//')
        output_vcf="${family_meta.proband_id}_inherited_snv_homopolymer_\${category}.vcf.gz"

        bedtools intersect -a ${snv_vcf} -b "\$hp_path" -header -wa | \
          bcftools view -O z -o tmp_\${category}.vcf.gz

        tabix tmp_\${category}.vcf.gz

        if [ -n "\$prev_vcfs" ]; then
            bcftools isec -C -w 1 -O z -o "\$output_vcf" tmp_\${category}.vcf.gz \$prev_vcfs
            rm -f tmp_\${category}.vcf.gz
        else
            mv tmp_\${category}.vcf.gz "\$output_vcf"
        fi

        if [ -s "\$output_vcf" ]; then
            tabix "\$output_vcf"
            prev_vcfs="\$prev_vcfs \$output_vcf"
        else
            rm -f "\$output_vcf"
        fi
    done

    # Copy full snvs VCF for all regions
    cp ${snv_vcf} "${family_meta.proband_id}_inherited_snv_all_regions.vcf.gz"
    cp ${snv_vcf_index} "${family_meta.proband_id}_inherited_snv_all_regions.vcf.gz.tbi"
    """
}


process HP_INTERSECT_INDEL {
    tag "${family_meta.proband_id} INDEL"
    publishDir "${params.outdir}/${family_meta.proband_id}", mode: 'copy'

    module 'bedtools'
    cpus 2
    memory '32 GB'
    time '2 h'

    input:
    tuple val(family_meta), path(indel_vcf), path(indel_vcf_index)
    path stratifications_dir
    path hp_union_bed

    output:
    tuple val(family_meta), path("*_indel*.vcf.gz"), path("*_indel*.vcf.gz.tbi")

    script:
    // Reverse order for prioritizing largest homopolymers first
    def hp_order = homopolymer_bed.reverse().join(' ')
    """
    mkdir -p homopolymers_indel

    # Create union (any homopolymer) for INDELs
    bedtools intersect -a ${indel_vcf} -b ${hp_union_bed} -header -wa | \
    bcftools view -O z -o "${family_meta.proband_id}_inherited_indel_homopolymer_any.vcf.gz"

    if [ -s "${family_meta.proband_id}_inherited_indel_homopolymer_any.vcf.gz" ]; then
        tabix "${family_meta.proband_id}_inherited_indel_homopolymer_any.vcf.gz"
    fi

    # Create non-homopolymer INDEL VCF as variants NOT in homopolymer_any
    bedtools intersect -a  ${indel_vcf} -b "${family_meta.proband_id}_inherited_indel_homopolymer_any.vcf.gz" -header -v | 
    bcftools view -O z -o "${family_meta.proband_id}_inherited_indel_non_homopolymer.vcf.gz"
   
    if [ -s "${family_meta.proband_id}_inherited_indel_non_homopolymer.vcf.gz" ]; then
        tabix "${family_meta.proband_id}_inherited_indel_non_homopolymer.vcf.gz"
    else
        rm -f "${family_meta.proband_id}_inherited_indel_non_homopolymer.vcf.gz"
    fi

    prev_vcfs=""

    for hp_file in ${hp_order}; do
        hp_path="${stratifications_dir}/\${hp_file}"
        if [ ! -f "\$hp_path" ]; then
            continue
        fi

        category=\$(echo "\$hp_file" | sed 's/GRCh38_SimpleRepeat_homopolymer_//' | sed 's/_slop5\\.bed\\.gz\$//')
        output_vcf="${family_meta.proband_id}_inherited_indel_homopolymer_\${category}.vcf.gz"

        bedtools intersect -a ${indel_vcf} -b "\$hp_path" -header -wa | \
          bcftools view -O z -o tmp_\${category}.vcf.gz

        tabix tmp_\${category}.vcf.gz

        if [ -n "\$prev_vcfs" ]; then
            bcftools isec -C -w 1 -O z -o "\$output_vcf" tmp_\${category}.vcf.gz \$prev_vcfs
            rm -f tmp_\${category}.vcf.gz
        else
            mv tmp_\${category}.vcf.gz "\$output_vcf"
        fi

        if [ -s "\$output_vcf" ]; then
            tabix "\$output_vcf"
            prev_vcfs="\$prev_vcfs \$output_vcf"
        else
            rm -f "\$output_vcf"
        fi
    done

    # Copy full INDELs VCF for all regions
    cp ${indel_vcf} "${family_meta.proband_id}_inherited_indel_all_regions.vcf.gz"
    cp ${indel_vcf_index} "${family_meta.proband_id}_inherited_indel_all_regions.vcf.gz.tbi"
    """
}

/*
 * PIPELINE 2 PROCESSES: Som.py Comparison Analysis
 */

process SPLIT_MULTIALLELIC {
    tag "${ngs_id}"

    module 'gatk'

    input:
    tuple val(ngs_id), path(lrs_vcf)

    output:
    tuple val(ngs_id), path("${ngs_id}_normalized.vcf.gz")

    script:
    """
    bcftools norm -m -any -f ${params.ref_ont} ${lrs_vcf} \
    | bcftools view -i "FORMAT/AF >= 0.20 && FORMAT/GQ >= 20 && FORMAT/DP >= 5" -Oz -o ${ngs_id}_LRS_filtered.vcf.gz
    tabix ${ngs_id}_LRS_filtered.vcf.gz
    
    gatk SelectVariants -V ${ngs_id}_LRS_filtered.vcf.gz --max-indel-size 50 --output ${ngs_id}_normalized.vcf
    bgzip ${ngs_id}_normalized.vcf
     
    tabix ${ngs_id}_normalized.vcf.gz
    """
}

process RUN_SOMPY {
    tag "${ngs_id}-${variant_file.simpleName}"
    publishDir "${params.outdir}/som_py_outputs", mode: 'copy'

    memory '32 GB'

    input:
    tuple val(ngs_id), path(normalized_vcf), path(variant_file)

    output:
    path("${ngs_id}_${variant_basename}.*csv")

    container "${params.container}"

    script:
    variant_basename = variant_file.simpleName
    """
    /opt/hap.py/bin/som.py \
      -o ${ngs_id}_${variant_basename} \
      -r ${params.ref_fasta} \
      --normalize-truth --normalize-query \
      ${variant_file} ${normalized_vcf}
    """
}

process PLOT_MEDIAN_RECALL {
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path(stats_files)
    
    output:
    path("cohort_median_recall_snv_indel_plot.pdf")
    path("panel_a_data.csv")
    path("panel_b_data.csv")
    
    script:
    """
#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import re
import numpy as np

# Gather all stats files (exclude any gt20 category)
files = glob.glob('*.stats.csv')
print(f"Found {len(files)} stats files")

all_dfs = []
for file in files:
    basename = file.replace('.stats.csv', '')
    parts = basename.split('_')
    if len(parts) >= 3:
        comparison_type = '_'.join(parts[2:])
    else:
        comparison_type = basename

    df = pd.read_csv(file, skiprows=[1])
    df['comparison_type'] = comparison_type
    df['sample_id'] = parts[0] if parts else 'unknown'
    print(comparison_type)
    all_dfs.append(df)

if not all_dfs:
    print("No stats files found!")
    exit(1)

# Combine all results
df = pd.concat(all_dfs, ignore_index=True)

print(f"Total records: {len(df)}")
print("Types found:", df['type'].unique())
print("Comparison types found:", df['comparison_type'].unique())

# Adjust labels: relabel gt11 to >=12bp, remove gt20 category
label_mapping = {
    'inherited_snv_all_regions': 'All regions',
    'inherited_snv_non_homopolymer': 'Non-homopolymer', 
    'inherited_snv_homopolymer_any': 'Homopolymer',
    'inherited_snv_homopolymer_4to6': '4-6bp',
    'inherited_snv_homopolymer_7to11': '7-11bp',
    'inherited_snv_homopolymer_gt11': '≥12bp',
    'inherited_indel_all_regions': 'All regions',
    'inherited_indel_non_homopolymer': 'Non-homopolymer',
    'inherited_indel_homopolymer_any': 'Homopolymer',
    'inherited_indel_homopolymer_4to6': '4-6bp',
    'inherited_indel_homopolymer_7to11': '7-11bp',
    'inherited_indel_homopolymer_gt11': '≥12bp'
}

# Extract variant type and region type
def extract_region_type(x):
    if 'all_regions' in x:
        return 'all_regions'
    elif 'non_homopolymer' in x:
        return 'non_homopolymer'
    elif 'homopolymer_any' in x:
        return 'homopolymer_any'
    elif 'homopolymer_4to6' in x:
        return 'homopolymer_4to6'
    elif 'homopolymer_7to11' in x:
        return 'homopolymer_7to11'
    elif 'homopolymer_gt11' in x:
        return 'homopolymer_gt11'
    else:
        return 'unknown'

df['variant_type'] = df['comparison_type'].apply(lambda x: 'SNV' if '_snv_' in x else 'INDEL' if '_indel_' in x else 'Unknown')
df['region_type'] = df['comparison_type'].apply(extract_region_type)

# Apply display label mapping
df['display_label'] = df['comparison_type'].map(label_mapping)
df = df[df['display_label'].notna()].copy()
df = df[df['region_type'] != 'unknown'].copy()

# Separate SNVs and INDELs
df_snv = df[(df['variant_type'] == 'SNV')].copy()
df_indel = df[(df['variant_type'] == 'INDEL')].copy()

print(f"SNV records: {len(df_snv)}")
print(f"INDEL records: {len(df_indel)}")

# Define region orders (exclude gt20)
region_order_a = ['all_regions', 'non_homopolymer', 'homopolymer_any']
region_order_b = ['homopolymer_4to6', 'homopolymer_7to11', 'homopolymer_gt11']

def prepare_grouped_data(df_snv, df_indel, region_order):
    # Compute median recall and median total.truth per group for SNV and INDEL
    snv_medians = (
        df_snv
        .groupby(['region_type', 'display_label'])
        .agg({'recall': 'median', 'total.truth': 'median'})
        .reset_index()
    )
    indel_medians = (
        df_indel
        .groupby(['region_type', 'display_label'])
        .agg({'recall': 'median', 'total.truth': 'median'})
        .reset_index()
    )

    # Filter by region order
    snv_medians = snv_medians[snv_medians['region_type'].isin(region_order)].copy()
    indel_medians = indel_medians[indel_medians['region_type'].isin(region_order)].copy()
    
    # Add variant type
    snv_medians['variant_type'] = 'SNV'
    indel_medians['variant_type'] = 'INDEL'
    
    # Combine
    combined = pd.concat([snv_medians, indel_medians], ignore_index=True)
    
    # Sort by region order
    combined['order'] = combined['region_type'].map({v: i for i, v in enumerate(region_order)})
    combined = combined.sort_values(['order', 'variant_type'])

    print(combined)
    return combined

# Create two-panel plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

def plot_grouped_bars(combined_data, ax, title, region_order):
    if combined_data.empty:
        ax.set_title(title, fontsize=14, fontweight='bold', loc='left')
        ax.text(0.5, 0.5, 'No data available', ha='center', va='center', transform=ax.transAxes)
        return
    
    regions = [r for r in region_order if r in combined_data['region_type'].values]
    n_regions = len(regions)
    if n_regions == 0:
        ax.set_title(title, fontsize=14, fontweight='bold', loc='left')
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
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=11)
    ax.set_ylim(0, 1)
    ax.legend()
    for i, (snv_val, indel_val) in enumerate(zip(snv_values, indel_values)):
        if snv_val > 0:
            ax.text(i - width/2, snv_val + 0.01, f'{snv_val:.3f}', 
                    ha='center', va='bottom', fontsize=9, fontweight='bold')
        if indel_val > 0:
            ax.text(i + width/2, indel_val + 0.01, f'{indel_val:.3f}', 
                    ha='center', va='bottom', fontsize=9, fontweight='bold')

# Panel A: General
panel_a_data = prepare_grouped_data(
    df_snv[df_snv['region_type'].isin(region_order_a)],
    df_indel[df_indel['region_type'].isin(region_order_a)],
    region_order_a
)
plot_grouped_bars(panel_a_data, ax1, 'A) Inherited Variants by Region', region_order_a)

# Panel B: Homopolymer, max bin now >=12bp
panel_b_data = prepare_grouped_data(
    df_snv[df_snv['region_type'].isin(region_order_b)],
    df_indel[df_indel['region_type'].isin(region_order_b)],
    region_order_b
)
plot_grouped_bars(panel_b_data, ax2, 'B) Homopolymer Region Length', region_order_b)

plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig("cohort_median_recall_snv_indel_plot.pdf", dpi=300, bbox_inches='tight')
plt.close()

panel_a_data.to_csv('panel_a_data.csv', index=False)
panel_b_data.to_csv('panel_b_data.csv', index=False)

print("Two-panel grouped SNV vs INDEL plot saved as cohort_median_recall_snv_indel_plot.pdf")
    """
}

/*
 * MAIN WORKFLOW
 */
workflow {
    // Input validation
    if (!file(params.input_vcf).exists()) {
        error "Input VCF file not found: ${params.input_vcf}"
    }
    if (!file(params.sample_csv).exists()) {
        error "Sample CSV file not found: ${params.sample_csv}"
    }
    if (!file(params.stratifications_dir).exists()) {
        error "Stratifications directory not found: ${params.stratifications_dir}"
    }
    if (!file(params.ref_fasta).exists()) {
        error "Reference FASTA file not found: ${params.ref_fasta}"
    }

    // Setup input paths
    input_vcf = params.input_vcf
    input_vcf_index = params.input_vcf + ".tbi"
    sample_csv = params.sample_csv
    stratifications_dir = params.stratifications_dir
    ref_fasta = params.ref_fasta

    // PIPELINE 1: Variant Analysis
    PARSE_FAMILIES(sample_csv)

    families_ch = PARSE_FAMILIES.out
        .splitText()
        .map { line -> 
            def fields = line.trim().split(',')
            return [
                proband_id: fields[0],
                maternal_id: fields[1], 
                paternal_id: fields[2],
                proband_vcf: fields[3],
                maternal_vcf: fields[4],
                paternal_vcf: fields[5],
            ]
        }

    // Process families through variant analysis pipeline
    SUBSET_TRIO(input_vcf, input_vcf_index, families_ch)
    NORMALIZE_VARIANTS(SUBSET_TRIO.out, ref_fasta)
    FILTER_AC(NORMALIZE_VARIANTS.out)
    MISSING2REF(FILTER_AC.out)
    INHERITANCE_FILTER(MISSING2REF.out)
    
    // Filter variants into SNV and INDEL channels after inheritance filtering
    FILTER_VARIANTS(INHERITANCE_FILTER.out)
    
    // Extract SNV and INDEL channels
    snv_channel = FILTER_VARIANTS.out.map { family_meta, snv_vcf, snv_index, indel_vcf, indel_index ->
        [family_meta, snv_vcf, snv_index]
    }
    
    indel_channel = FILTER_VARIANTS.out.map { family_meta, snv_vcf, snv_index, indel_vcf, indel_index ->
        [family_meta, indel_vcf, indel_index]
    }

    CREATE_HP_UNION(stratifications_dir)
    
    // Process SNVs and INDELs through homopolymer intersection
    snv_intersect_results = HP_INTERSECT_SNV(snv_channel, stratifications_dir, CREATE_HP_UNION.out)
    indel_intersect_results = HP_INTERSECT_INDEL(indel_channel, stratifications_dir, CREATE_HP_UNION.out)

    // PIPELINE 2: Som.py Comparison Analysis - Simplified approach
    if (params.run_sompy) {
        if (!file(params.lrs_csv).exists()) {
            error "LRS CSV file not found: ${params.lrs_csv}"
        }

        // Simple completion check - just use one of the completion channels
        pipeline1_done = indel_intersect_results
            .collect()
            .map { results -> 
                log.info "Pipeline 1 completed. Starting Pipeline 2..."
                return "DONE"
            }

        // Create LRS samples channel
        ngs_samples = Channel
            .fromPath(params.lrs_csv)
            .splitCsv(header:true)
            .filter { row -> row.NGS_ID && row.LRS_SNV_VCF }
            .map { row -> [ row.NGS_ID, file(row.LRS_SNV_VCF) ] }

        // Split multiallelic variants
        normalized_vcfs_raw = SPLIT_MULTIALLELIC(ngs_samples)
        
        // Wait for Pipeline 1 completion
        normalized_vcfs = normalized_vcfs_raw
            .combine(pipeline1_done)
            .map { ngs_id, normalized_vcf, done_signal -> 
                [ngs_id, normalized_vcf]
            }
        
        // Create comparison pairs for both SNVs and INDELs
        comparison_pairs = normalized_vcfs
            .flatMap { ngs_id, normalized_vcf ->
                def variant_dir = file("${params.outdir}/${ngs_id}")
                if (!variant_dir.exists()) {
                    log.warn "Results directory not found for ${ngs_id}: ${variant_dir}"
                    return []
                }
                
                def all_pairs = []
                
                // Add SNV comparisons
                def snv_pattern = "${params.outdir}/${ngs_id}/*inherited*snv*.vcf.gz"
                def snv_files = file(snv_pattern, type: 'file')
                if (snv_files instanceof List) {
                    snv_files.each { snv_file ->
                        all_pairs.add([ngs_id, normalized_vcf, snv_file])
                    }
                } else if (snv_files.exists()) {
                    all_pairs.add([ngs_id, normalized_vcf, snv_files])
                }
                
                // Add INDEL comparisons
                def indel_pattern = "${params.outdir}/${ngs_id}/*inherited*indel*.vcf.gz"
                def indel_files = file(indel_pattern, type: 'file')
                if (indel_files instanceof List) {
                    indel_files.each { indel_file ->
                        all_pairs.add([ngs_id, normalized_vcf, indel_file])
                    }
                } else if (indel_files.exists()) {
                    all_pairs.add([ngs_id, normalized_vcf, indel_files])
                }
                
                if (all_pairs.isEmpty()) {
                    log.warn "No inherited variant files found for ${ngs_id}"
                }
                
                return all_pairs
            }
        
        // Run som.py comparisons
        RUN_SOMPY(comparison_pairs)

        // Generate plots
        stats_files = RUN_SOMPY.out
            .flatten()
            .filter { it.name.endsWith('.stats.csv') }
            .collect()

        PLOT_MEDIAN_RECALL(stats_files)
    }
}


/*
 * Workflow completion handler
 */
workflow.onComplete {
    log.info """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    Exit status : ${workflow.exitStatus}
    Error report: ${(workflow.errorReport ?: 'None')}
    """.stripIndent()
}

