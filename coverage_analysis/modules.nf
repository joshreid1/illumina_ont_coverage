// ============================================================================
// MODULES.NF - Process Definitions for Coverage Analysis Workflow  
// ============================================================================

process GencodePreprocessing {
	
	input:
	path gff_path
	path db_path
	val gene_list

	output:
	path "qualifying_genes.txt"

	script:
	"""
	#!/usr/bin/env python3
	import gffutils

	# Load GTF DB
	db = gffutils.FeatureDB('${db_path}', keep_order=True)

	# Build gene list from GTF, filtering for protein_coding genes, excluding chrM and chrY, and ensuring canonical transcript
	all_gene_dict = {}
	for g in db.features_of_type("gene"):
		if (g.attributes.get("gene_type", [None])[0] == "protein_coding" and 
			g.chrom not in ["chrM", "chrY"]):
			# Check if the gene has a canonical transcript
			transcripts = list(db.children(g, featuretype='transcript'))
			canonical = [t for t in transcripts if 'Ensembl_canonical' in t.attributes.get('tag', [])]
			if canonical:
				all_gene_dict[g.attributes.get("gene_name", [g.id])[0]] = g

	# Filter if gene list provided
	if '${gene_list}' != 'null' and '${gene_list}' != '':
		with open('${gene_list}') as f:
			input_genes = set(line.strip() for line in f if line.strip())
		gene_dict = {g: all_gene_dict[g] for g in input_genes if g in all_gene_dict}
	else:
		gene_dict = all_gene_dict

	# Write qualifying genes to a file
	with open('qualifying_genes.txt', 'w') as f:
		for gene in gene_dict:
			f.write(f"{gene}\\n")
	
	print(f"Found {len(gene_dict)} qualifying protein-coding genes with canonical transcripts")
	"""
}

process SplitGeneList {

    publishDir "${params.output_dir}/gene_list", mode: 'copy'
	
	input:
	path gene_file

	output:
	path "genes_chunk_*"

	script:
	"""
	#!/usr/bin/env bash
	split -l ${params.genes_per_job} ${gene_file} genes_chunk_
	echo "Split gene list into chunks of ${params.genes_per_job} genes each"
	"""
}

process CoverageAnalysis_v3 {
    tag "${platform}_${cram_file.baseName}_${gene_chunk.baseName}"

    //clusterOptions '--qos=bonus'

    cpus = 2
    memory = { 64 * task.attempt + ' GB' }
    time = { 12 * task.attempt + ' h'}

    publishDir "${params.output_dir}/sample_coverage", mode: 'copy'
    
    input:
    tuple val(platform), path(cram_file), path(cram_index), path(reference_fasta), path(gene_chunk)
    path db_path
    
    output:
    tuple val(platform), val(cram_file.baseName), path("*_exon_coverage.tsv"), emit: exon_coverage
    tuple val(platform), val(cram_file.baseName), path("*_intron_coverage.tsv"), emit: intron_coverage
    tuple val(platform), val(cram_file.baseName), path("*_gene_coverage.tsv"), emit: gene_coverage
    
    script:
    def sample_id = "${platform}_${cram_file.baseName.replaceAll(/\.(cram|bam)$/, '')}"
    def chunk_id = gene_chunk.baseName
    """
#!/usr/bin/env python3
import gffutils
import pysam
import time
import sys
import os
import numpy as np


def timestamp(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def find_gene_by_name(db, gene_symbol):
    for gene in db.features_of_type('gene'):
        if 'gene_name' in gene.attributes:
            if (gene_symbol in gene.attributes['gene_name'] and
                gene.chrom not in ["chrM", "chrY"]):

                gene_type = gene.attributes.get("gene_type",
                                                gene.attributes.get("gene_biotype", [None]))[0]
                if gene_type != "protein_coding":
                    continue

                transcripts = list(db.children(gene, featuretype='transcript'))
                canonical = [t for t in transcripts
                             if 'Ensembl_canonical' in t.attributes.get('tag', [])]
                if canonical:
                    return gene, canonical[0]
    return None, None


def analyze_feature_per_base(feature, bam):
    chrom, start, end = feature.chrom, feature.start, feature.end
    length = end - start
    depths = np.zeros(length, dtype=np.int32)
    any_depths = np.zeros(length, dtype=np.int32)
    mapq0_only = 0

    for pcol in bam.pileup(chrom, start, end,
                           truncate=True, stepper="samtools", min_base_quality=0):
        idx = pcol.reference_pos - start
        total_pr = 0
        non0_count = 0
        for pr in pcol.pileups:
            if not pr.alignment.is_unmapped:
                total_pr += 1
                if pr.alignment.mapping_quality > 0:
                    non0_count += 1
        any_depths[idx] = total_pr
        depths[idx] = non0_count
        if total_pr > 0 and non0_count == 0:
            mapq0_only += 1

    zero = np.count_nonzero(any_depths == 0)
    nonzero_vals = depths[depths > 0]

    if nonzero_vals.size > 0:
        mean_cov = np.mean(nonzero_vals)
        median_cov = np.median(nonzero_vals)
        std_cov = np.std(nonzero_vals, ddof=1)
    else:
        mean_cov = median_cov = std_cov = 0.0

    return length, zero, mapq0_only, mean_cov, median_cov, std_cov


def analyze_intron_region(chrom, start, end, bam):
    # ---Analyze coverage for an intron region defined by coordinates.---
    class IntronFeature:
        def __init__(self, chrom, start, end):
            self.chrom = chrom
            self.start = start
            self.end = end
    
    intron = IntronFeature(chrom, start, end)
    return analyze_feature_per_base(intron, bam)


# -----------------------
# Main script parameters
# -----------------------
platform = "${platform}"
reference_fasta = "${reference_fasta}"
sample_id = "${sample_id}"
chunk_id = "${chunk_id}"
num_threads = ${task.cpus}

timestamp(f"Processing {platform} sample: {sample_id}, chunk: {chunk_id}")
timestamp(f"Using reference: {reference_fasta}")

timestamp("Loading gene list...")
with open('${gene_chunk}') as f:
    genes = [line.strip() for line in f if line.strip()]

timestamp(f"Processing {len(genes)} genes from chunk {chunk_id}...")

timestamp("Loading database...")
db = gffutils.FeatureDB('${db_path}', keep_order=True)

timestamp(f"Opening {platform} CRAM/BAM file with reference...")
bam = pysam.AlignmentFile('${cram_file}', "rc",
                          threads=num_threads, reference_filename=reference_fasta)

# Output files
exon_output_file = f"{sample_id}_{chunk_id}_exon_coverage.tsv"
intron_output_file = f"{sample_id}_{chunk_id}_intron_coverage.tsv"
gene_output_file = f"{sample_id}_{chunk_id}_gene_coverage.tsv"

timestamp("Starting analysis...")
genes_processed = 0
genes_skipped = 0

with open(exon_output_file, 'w') as exon_out, \
     open(intron_output_file, 'w') as intron_out, \
     open(gene_output_file, 'w') as gene_out:
    
    # Headers
    exon_header = ["GENE", "TRANSCRIPT", "EXON_NUMBER", "CHROM", "START", "END",
                   "LENGTH", "ZERO_BASES", "MAPQ0_ONLY_BASES", "%ZERO", "%MAPQ0_ONLY",
                   "MEAN_COV", "MEDIAN_COV", "STD_COV", "CODING_STATUS"]
    exon_out.write("\\t".join(exon_header) + "\\n")

    intron_header = ["GENE", "TRANSCRIPT", "INTRON_NUMBER", "CHROM", "START", "END",
                     "LENGTH", "ZERO_BASES", "MAPQ0_ONLY_BASES", "%ZERO", "%MAPQ0_ONLY",
                     "MEAN_COV", "MEDIAN_COV", "STD_COV"]
    intron_out.write("\\t".join(intron_header) + "\\n")

    gene_header = ["GENE", "TRANSCRIPT", "CHROM", "START", "END",
                   "LENGTH", "ZERO_BASES", "MAPQ0_ONLY_BASES", "%ZERO", "%MAPQ0_ONLY",
                   "MEAN_COV", "MEDIAN_COV", "STD_COV"]
    gene_out.write("\\t".join(gene_header) + "\\n")

    progress_threshold = len(genes) * 0.1 if len(genes) > 10 else 1
    next_prog = progress_threshold

    for i, gene_name in enumerate(genes, start=1):
        if i >= next_prog:
            pct = i / len(genes) * 100
            timestamp(f"Progress: {pct:.0f}% ({i}/{len(genes)})")
            next_prog += progress_threshold

        timestamp(f"Processing gene: {gene_name}")

        gene_feat, canonical_tx = find_gene_by_name(db, gene_name)
        if gene_feat is None or canonical_tx is None:
            timestamp(f"WARNING: Gene {gene_name} not found or missing canonical transcript, skipping...")
            genes_skipped += 1
            continue

        try:
            # --- Get exons ---
            exon_features = list(db.children(canonical_tx, featuretype='exon', order_by='start'))
            if not exon_features:
                timestamp(f"WARNING: No exon features found for {gene_name}, skipping...")
                genes_skipped += 1
                continue

            # Determine transcript strand
            transcript_strand = canonical_tx.strand
            timestamp(f"  Found {len(exon_features)} exon features on strand {transcript_strand}")

            # Collect CDS and UTR features for annotation
            cds_features = list(db.children(canonical_tx, featuretype='CDS', order_by='start'))
            utr_features = list(db.children(canonical_tx, featuretype='UTR', order_by='start'))
            
            # Gene-level analysis (whole transcript span)
            gene_L, gene_zero, gene_mapq0, gene_mean, gene_median, gene_std = analyze_feature_per_base(canonical_tx, bam)
            gene_pct_zero = (gene_zero / gene_L * 100) if gene_L else 0
            gene_pct_mapq0 = (gene_mapq0 / gene_L * 100) if gene_L else 0

            gene_out.write(
                f"{gene_name}\\t{canonical_tx.id}\\t{canonical_tx.chrom}\\t{canonical_tx.start}\\t{canonical_tx.end}\\t"
                f"{gene_L}\\t{gene_zero}\\t{gene_mapq0}\\t{gene_pct_zero:.2f}\\t{gene_pct_mapq0:.2f}\\t"
                f"{gene_mean:.2f}\\t{gene_median:.2f}\\t{gene_std:.2f}\\n"
            )

            # Exon-level analysis
            for exon in exon_features:
                exon_num = exon.attributes.get('exon_number', ['Unknown'])[0]

                # Determine coding/UTR status
                coding_status = "noncoding"
                for cds in cds_features:
                    if exon.start <= cds.end and exon.end >= cds.start:
                        coding_status = "coding"
                        break
                if coding_status == "noncoding":
                    for utr in utr_features:
                        if exon.start <= utr.end and exon.end >= utr.start:
                            coding_status = "UTR"
                            break
                
                exon_L, exon_zero, exon_mapq0, exon_mean, exon_median, exon_std = analyze_feature_per_base(exon, bam)
                exon_pct_zero = (exon_zero / exon_L * 100) if exon_L else 0
                exon_pct_mapq0 = (exon_mapq0 / exon_L * 100) if exon_L else 0

                exon_out.write(
                    f"{gene_name}\\t{canonical_tx.id}\\t{exon_num}\\t{exon.chrom}\\t{exon.start}\\t{exon.end}\\t"
                    f"{exon_L}\\t{exon_zero}\\t{exon_mapq0}\\t{exon_pct_zero:.2f}\\t{exon_pct_mapq0:.2f}\\t"
                    f"{exon_mean:.2f}\\t{exon_median:.2f}\\t{exon_std:.2f}\\t{coding_status}\\n"
                )

            # Intron-level analysis
            # Calculate introns as intervals between consecutive exons
            # Account for strand: on negative strand, exons are numbered in reverse transcriptional order
            if len(exon_features) > 1:
                # Create list of exon pairs with proper intron numbering
                if transcript_strand == '-':
                    # For negative strand: intron 1 is between highest numbered exons
                    # Need to reverse the genomic order for biological numbering
                    exon_pairs = []
                    for idx in range(len(exon_features) - 1):
                        # Intron number counts backward from the biological perspective
                        intron_num = len(exon_features) - idx - 1
                        exon_pairs.append((intron_num, exon_features[idx], exon_features[idx + 1]))
                else:
                    # For positive strand: intron 1 is between exon 1 and exon 2
                    exon_pairs = []
                    for idx in range(len(exon_features) - 1):
                        intron_num = idx + 1
                        exon_pairs.append((intron_num, exon_features[idx], exon_features[idx + 1]))
                
                for intron_num, exon_current, exon_next in exon_pairs:
                    # Intron spans from end of current exon to start of next exon
                    intron_start = exon_current.end + 1
                    intron_end = exon_next.start - 1
                    
                    # Skip if invalid intron (shouldn't happen with proper GTF)
                    if intron_end < intron_start:
                        timestamp(f"WARNING: Invalid intron coordinates for {gene_name} intron {intron_num}, skipping...")
                        continue
                    
                    intron_L, intron_zero, intron_mapq0, intron_mean, intron_median, intron_std = \
                        analyze_intron_region(canonical_tx.chrom, intron_start, intron_end, bam)
                    intron_pct_zero = (intron_zero / intron_L * 100) if intron_L else 0
                    intron_pct_mapq0 = (intron_mapq0 / intron_L * 100) if intron_L else 0

                    intron_out.write(
                        f"{gene_name}\\t{canonical_tx.id}\\t{intron_num}\\t{canonical_tx.chrom}\\t{intron_start}\\t{intron_end}\\t"
                        f"{intron_L}\\t{intron_zero}\\t{intron_mapq0}\\t{intron_pct_zero:.2f}\\t{intron_pct_mapq0:.2f}\\t"
                        f"{intron_mean:.2f}\\t{intron_median:.2f}\\t{intron_std:.2f}\\n"
                    )
                
                timestamp(f"  Analyzed {len(exon_features) - 1} introns")

            genes_processed += 1
            timestamp(f"Completed gene: {gene_name} ({len(exon_features)} exons, {max(0, len(exon_features) - 1)} introns)")

        except Exception as e:
            timestamp(f"ERROR processing gene {gene_name}: {str(e)}")
            genes_skipped += 1
            continue

timestamp("Analysis complete!")
timestamp(f"Summary: Processed {genes_processed} genes, skipped {genes_skipped} genes")
timestamp(f"Outputs: {exon_output_file} (exon-level), {intron_output_file} (intron-level), and {gene_output_file} (gene-level)")
bam.close()
    """
}

process ConcatenateResults {
    tag "${platform}"

    //clusterOptions '--qos=bonus'
    
    input:
    tuple val(platform), val(sample_id), path(exon_files), path(intron_files), path(gene_files)
    
    output:
    tuple val(platform), path("${platform}_${sample_id}_exon_coverage_combined.tsv"), emit: exon_combined
    tuple val(platform), path("${platform}_${sample_id}_intron_coverage_combined.tsv"), emit: intron_combined
    tuple val(platform), path("${platform}_${sample_id}_gene_coverage_combined.tsv"), emit: gene_combined
    
    script:
    """
#!/usr/bin/env python3
import os
import glob


def concatenate_files(file_list, output_file, file_type):
    print(f"Concatenating {len(file_list)} {file_type} files into {output_file}")
    
    if not file_list:
        print(f"Warning: No {file_type} files found for sample ${sample_id}")
        # Create empty file with header
        if file_type == "Exon":
            header = ["GENE","TRANSCRIPT","EXON_NUMBER","CHROM","START","END",
                        "LENGTH","ZERO_BASES","MAPQ0_ONLY_BASES","%ZERO","%MAPQ0_ONLY",
                        "MEAN_COV","MEDIAN_COV","STD_COV","CODING_STATUS"]
        elif file_type == "Intron":
            header = ["GENE","TRANSCRIPT","INTRON_NUMBER","CHROM","START","END",
                        "LENGTH","ZERO_BASES","MAPQ0_ONLY_BASES","%ZERO","%MAPQ0_ONLY",
                        "MEAN_COV","MEDIAN_COV","STD_COV"]
        else:  # Gene
            header = ["GENE","TRANSCRIPT","CHROM","START","END",
                        "LENGTH","ZERO_BASES","MAPQ0_ONLY_BASES","%ZERO","%MAPQ0_ONLY",
                        "MEAN_COV","MEDIAN_COV","STD_COV"]
        
        with open(output_file, 'w') as out:
            out.write("\\t".join(header) + "\\n")
        return
    
    # Sort files to ensure consistent ordering
    file_list = sorted(file_list)
    
    with open(output_file, 'w') as out:
        header_written = False
        total_lines = 0
        
        for file_path in file_list:
            print(f"  Processing: {file_path}")
            try:
                with open(file_path, 'r') as f:
                    lines = f.readlines()
                    
                if not lines:
                    continue
                
                # Write header from first file only
                if not header_written:
                    out.write(lines[0])  # Header line
                    header_written = True
                
                # Write data lines (skip header)
                data_lines = lines[1:] if len(lines) > 1 else []
                for line in data_lines:
                    out.write(line)
                    total_lines += 1
                    
            except Exception as e:
                print(f"Error processing {file_path}: {e}")
                continue
        
        print(f"  Total data lines written: {total_lines}")


# Get file lists
exon_files = "${exon_files}".split()
intron_files = "${intron_files}".split()
gene_files = "${gene_files}".split()


print(f"Platform: ${platform}, Sample: ${sample_id}")
print(f"Exon files to concatenate: {len(exon_files)}")
print(f"Intron files to concatenate: {len(intron_files)}")
print(f"Gene files to concatenate: {len(gene_files)}")


# Concatenate Exon files
concatenate_files(exon_files, "${platform}_${sample_id}_exon_coverage_combined.tsv", "Exon")


# Concatenate Intron files
concatenate_files(intron_files, "${platform}_${sample_id}_intron_coverage_combined.tsv", "Intron")


# Concatenate Gene files  
concatenate_files(gene_files, "${platform}_${sample_id}_gene_coverage_combined.tsv", "Gene")


print("Concatenation complete!")
    """
}

process SummarizeCoverage {
    tag "${platform}"
    
    cpus = 4  // Increased for parallel processing
    memory = { 128 * task.attempt + ' GB' } 
    time = { 8 * task.attempt + ' h'}  
    
    //publishDir "${params.output_dir}/${platform}/cohort_summary", mode: 'copy'
    
    input:
    tuple val(platform), path(exon_files)
    tuple val(platform), path(intron_files)
    tuple val(platform), path(gene_files)
    
    output:
    tuple val(platform), path("${platform.toLowerCase()}_cohort_exon_coverage_summary_v2.tsv"), emit: exon_summary
    tuple val(platform), path("${platform.toLowerCase()}_cohort_intron_coverage_summary_v2.tsv"), emit: intron_summary
    tuple val(platform), path("${platform.toLowerCase()}_cohort_gene_coverage_summary_v2.tsv"), emit: gene_summary
    
    script:
    """
#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os
from typing import List, Dict, Optional



def optimize_dtypes(df: pd.DataFrame) -> pd.DataFrame:
    #Optimize DataFrame data types for memory efficiency.
    # Convert string columns to categorical if they have repeated values
    for col in df.select_dtypes(include=['object']).columns:
        if df[col].nunique() < len(df) * 0.5:  # If less than 50% unique values
            df[col] = df[col].astype('category')
    
    # Downcast numeric columns
    for col in df.select_dtypes(include=['int64']).columns:
        df[col] = pd.to_numeric(df[col], downcast='integer')
    
    for col in df.select_dtypes(include=['float64']).columns:
        df[col] = pd.to_numeric(df[col], downcast='float')
    
    return df



def read_and_process_file(filepath: str, feature_type: str) -> Optional[pd.DataFrame]:
    #Read and preprocess a single coverage file.
    try:
        # Read with optimized settings
        dtype_spec = None
        if feature_type == "Exon":
            dtype_spec = {'EXON_NUMBER': 'Int32'}
        elif feature_type == "Intron":
            dtype_spec = {'INTRON_NUMBER': 'Int32'}
        
        df = pd.read_csv(filepath, sep="\\t", dtype=dtype_spec)
        
        # Extract sample name more efficiently
        sample_name = os.path.basename(filepath)
        sample_name = (sample_name.replace('_exon_coverage_combined.tsv', '')
                                  .replace('_intron_coverage_combined.tsv', '')
                                  .replace('_gene_coverage_combined.tsv', ''))
        df['SAMPLE'] = sample_name
        
        # Optimize data types
        df = optimize_dtypes(df)
        
        # Create feature key
        if feature_type == "Exon":
            df['FEATURE_KEY'] = df['GENE'].astype(str) + '_' + df['TRANSCRIPT'].astype(str) + '_EXON' + df['EXON_NUMBER'].astype(str)
        elif feature_type == "Intron":
            df['FEATURE_KEY'] = df['GENE'].astype(str) + '_' + df['TRANSCRIPT'].astype(str) + '_INTRON' + df['INTRON_NUMBER'].astype(str)
        else:  # Gene
            df['FEATURE_KEY'] = df['GENE'].astype(str) + '_' + df['TRANSCRIPT'].astype(str)
        
        return df
        
    except Exception as e:
        print(f"    Error reading {filepath}: {e}", flush=True)
        return None



def calculate_summary_stats(df: pd.DataFrame, feature_type: str) -> pd.DataFrame:
    #Calculate summary statistics efficiently using vectorized operations.
    
    # Define columns based on feature type
    if feature_type == "Exon":
        group_cols = ['GENE', 'TRANSCRIPT', 'EXON_NUMBER', 'CHROM', 'START', 'END', 'LENGTH']
    elif feature_type == "Intron":
        group_cols = ['GENE', 'TRANSCRIPT', 'INTRON_NUMBER', 'CHROM', 'START', 'END', 'LENGTH']
    else:  # Gene
        group_cols = ['GENE', 'TRANSCRIPT', 'CHROM', 'START', 'END', 'LENGTH']
    
    # Numeric columns for statistics
    numeric_cols = ['ZERO_BASES', 'MAPQ0_ONLY_BASES', '%ZERO', '%MAPQ0_ONLY', 'MEAN_COV', 'MEDIAN_COV']
    available_numeric_cols = [col for col in numeric_cols if col in df.columns]
    
    # Get positional information (first occurrence of each feature)
    feature_info = df.groupby('FEATURE_KEY')[group_cols].first().reset_index()
    
    # Calculate all statistics at once using agg
    agg_funcs = {
        col: ['median', 'mean', 'min', 'max', 'std'] for col in available_numeric_cols
    }
    
    # Perform aggregation
    stats_df = df.groupby('FEATURE_KEY')[available_numeric_cols].agg(agg_funcs)
    
    # Flatten column names
    stats_df.columns = [f"{col[0].upper()}_{col[1].upper()}" for col in stats_df.columns]
    stats_df = stats_df.reset_index()
    
    # Round values appropriately
    for col in stats_df.columns:
        if col != 'FEATURE_KEY':
            if 'MEDIAN' in col and ('ZERO_BASES' in col or 'MAPQ0_ONLY_BASES' in col or 'COV' in col):
                stats_df[col] = stats_df[col].round().astype('Int64')
            else:
                stats_df[col] = stats_df[col].round(2)
    
    # Merge with positional information
    result_df = pd.merge(feature_info, stats_df, on='FEATURE_KEY', how='inner')
    
    # Remove the FEATURE_KEY column as it's not needed in output
    result_df = result_df.drop('FEATURE_KEY', axis=1)
    
    # Sort the results
    if feature_type == "Exon":
        result_df = result_df.sort_values(['GENE', 'TRANSCRIPT', 'EXON_NUMBER'])
    elif feature_type == "Intron":
        result_df = result_df.sort_values(['GENE', 'TRANSCRIPT', 'INTRON_NUMBER'])
    else:  # Gene
        result_df = result_df.sort_values(['GENE', 'TRANSCRIPT'])
    
    return result_df


def create_empty_output(feature_type: str, filename: str):
    #Create empty output file with appropriate headers.
    base_headers = ['GENE', 'TRANSCRIPT']
    if feature_type == "Exon":
        base_headers.append('EXON_NUMBER')
    elif feature_type == "Intron":
        base_headers.append('INTRON_NUMBER')
    base_headers.extend(['CHROM', 'START', 'END', 'LENGTH'])
    
    stat_headers = []
    for metric in ['ZERO_BASES', 'MAPQ0_ONLY_BASES', '%ZERO', '%MAPQ0_ONLY', 'MEAN_COV', 'MEDIAN_COV']:
        for stat in ['MEDIAN', 'MEAN', 'MIN', 'MAX', 'STD']:
            stat_headers.append(f"{metric}_{stat}")
    
    all_headers = base_headers + stat_headers
    
    with open(filename, 'w') as f:
        f.write("\\t".join(all_headers) + "\\n")


def summarize_coverage_files(files: List[str], output_name: str, feature_type: str):
    #Main function to summarize coverage files with optimizations.
    print(f"→ Processing {len(files)} {feature_type} files for ${platform} cohort summary", flush=True)
    
    if not files:
        print(f"No {feature_type} files found for ${platform}", flush=True)
        create_empty_output(feature_type, output_name)
        return
    
    # Sort files for consistent processing
    files = sorted(files)
    print(f"Files to process: {[os.path.basename(f) for f in files]}", flush=True)
    
    # Read and process all files
    dataframes = []
    for filepath in files:
        print(f"  Reading: {os.path.basename(filepath)}", flush=True)
        df = read_and_process_file(filepath, feature_type)
        if df is not None:
            dataframes.append(df)
    
    if not dataframes:
        print(f"No valid {feature_type} files could be processed for ${platform}", flush=True)
        create_empty_output(feature_type, output_name)
        return
    
    # Combine all dataframes efficiently
    print("Combining all sample data...", flush=True)
    combined_df = pd.concat(dataframes, ignore_index=True, copy=False)
    
    # Clear individual dataframes to free memory
    del dataframes
    
    print(f"Combined data shape: {combined_df.shape}", flush=True)
    print(f"Samples: {combined_df['SAMPLE'].nunique()}", flush=True)
    print(f"Unique features: {combined_df['FEATURE_KEY'].nunique()}", flush=True)
    
    # Calculate summary statistics
    print("Calculating summary statistics...", flush=True)
    summary_df = calculate_summary_stats(combined_df, feature_type)
    
    # Save with optimized settings
    summary_df.to_csv(output_name, sep="\\t", index=False, float_format='%.2f')
    
    print(f"  ✔ Summary written: {output_name}", flush=True)
    print(f"    Features summarized: {len(summary_df)}", flush=True)
    print(f"    Total columns: {len(summary_df.columns)}", flush=True)
    
    # Print examples for verification
    print("\\nFirst few features:", flush=True)
    for i, row in summary_df.head().iterrows():
        if feature_type == "Exon":
            print(f"  {row['GENE']} {row['TRANSCRIPT']} exon{row['EXON_NUMBER']}: MEAN_COV_MEDIAN={row.get('MEAN_COV_MEDIAN', 'N/A')}", flush=True)
        elif feature_type == "Intron":
            print(f"  {row['GENE']} {row['TRANSCRIPT']} intron{row['INTRON_NUMBER']}: MEAN_COV_MEDIAN={row.get('MEAN_COV_MEDIAN', 'N/A')}", flush=True)
        else:
            print(f"  {row['GENE']} {row['TRANSCRIPT']}: MEAN_COV_MEDIAN={row.get('MEAN_COV_MEDIAN', 'N/A')}", flush=True)


# Process Exon files
exon_files = "${exon_files}".split()
exon_files = [f for f in exon_files if f and f != "null"]

if exon_files:
    print(f"Found {len(exon_files)} Exon files for ${platform}", flush=True)
    summarize_coverage_files(exon_files, "${platform.toLowerCase()}_cohort_exon_coverage_summary_v2.tsv", "Exon")
else:
    print(f"No Exon files found for ${platform}")
    create_empty_output("Exon", "${platform.toLowerCase()}_cohort_exon_coverage_summary_v2.tsv")


# Process Intron files
intron_files = "${intron_files}".split()
intron_files = [f for f in intron_files if f and f != "null"]

if intron_files:
    print(f"Found {len(intron_files)} Intron files for ${platform}", flush=True)
    summarize_coverage_files(intron_files, "${platform.toLowerCase()}_cohort_intron_coverage_summary_v2.tsv", "Intron")
else:
    print(f"No Intron files found for ${platform}", flush=True)
    create_empty_output("Intron", "${platform.toLowerCase()}_cohort_intron_coverage_summary_v2.tsv")


# Process Gene files
gene_files = "${gene_files}".split()
gene_files = [f for f in gene_files if f and f != "null"]

if gene_files:
    print(f"Found {len(gene_files)} Gene files for ${platform}", flush=True)
    summarize_coverage_files(gene_files, "${platform.toLowerCase()}_cohort_gene_coverage_summary_v2.tsv", "Gene")
else:
    print(f"No Gene files found for ${platform}", flush=True)
    create_empty_output("Gene", "${platform.toLowerCase()}_cohort_gene_coverage_summary_v2.tsv")


print(f"${platform} cohort summary analysis complete!", flush=True)
    """
}

process GnomadCoverageAnalysis {
    tag "${gene_chunk.baseName}"


    cpus = 2
    memory = { 128 * task.attempt + ' GB' }
    time = { 12 * task.attempt + ' h'}

    publishDir "${params.output_dir}/gnomad", mode: 'copy'

    //clusterOptions '--qos=bonus'
    
    input:
    tuple path(gene_chunk), path(exome_coverage), path(exome_index), path(genome_coverage), path(genome_index)
    path db_path
    
    output:
    path("*_exon_coverage.tsv"), emit: exon_coverage
    path("*_intron_coverage.tsv"), emit: intron_coverage

    script:
    def chunk_id = gene_chunk.baseName
    """
#!/usr/bin/env python3
import gffutils
import pysam
import time
import sys
import os
import numpy as np

def timestamp(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

def normalize_chromosome(chrom):
    # Convert chromosome format to match gnomAD files (chr1, chr2, etc.).
    chrom_str = str(chrom).strip()
    
    # If it already starts with 'chr', return as is
    if chrom_str.startswith('chr'):
        return chrom_str
    
    # Add 'chr' prefix for standard chromosomes
    if chrom_str in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 
                     '11', '12', '13', '14', '15', '16', '17', '18', '19', 
                     '20', '21', '22', 'X', 'Y', 'MT', 'M']:
        return f'chr{chrom_str}'
    
    # Handle mitochondrial chromosome
    if chrom_str in ['MT', 'M']:
        return 'chrM'
    
    # For any other format, add chr prefix
    return f'chr{chrom_str}'

def find_gene_by_name(db, gene_symbol):
    for gene in db.features_of_type('gene'):
        if 'gene_name' in gene.attributes:
            if (gene_symbol in gene.attributes['gene_name'] and
                gene.chrom not in ["chrM", "chrY"]):

                gene_type = gene.attributes.get("gene_type",
                                                gene.attributes.get("gene_biotype", [None]))[0]
                if gene_type != "protein_coding":
                    continue


                transcripts = list(db.children(gene, featuretype='transcript'))
                canonical = [t for t in transcripts
                             if 'Ensembl_canonical' in t.attributes.get('tag', [])]
                if canonical:
                    return gene, canonical[0]
    return None, None

def extract_coverage_from_gnomad(chrom, start, end, tabix_file):
    # Normalize chromosome format for gnomAD compatibility
    query_chrom = normalize_chromosome(chrom)

    # Adjust coordinates (1-based to 0-based for tabix)
    start_0based = start - 1

    timestamp(f"    Querying region: {chrom}:{start}-{end} (normalized: {query_chrom}:{start}-{end})")

    means = []
    medians = []
    zero_positions = 0

    try:
        # Query the tabix file for overlapping regions
        for record in tabix_file.fetch(query_chrom, start_0based, end):
            fields = record.strip().split('\t')
            if len(fields) >= 5:
                # gnomAD BED format: chrom, start, end, mean, median_approx, ...
                rec_start = int(fields[1])
                rec_end = int(fields[2])
                mean_cov = float(fields[3])
                median_cov = float(fields[4])

                # Check if this position overlaps with our feature
                if rec_start < end and rec_end > start_0based:
                    means.append(mean_cov)
                    medians.append(median_cov)

                    # Count zero positions (where median coverage is 0)
                    if median_cov == 0.0:
                        zero_positions += 1

    except Exception as e:
        timestamp(f"    Warning: Could not query {query_chrom}:{start}-{end}: {e}")
        return 0, 0.0, 0.0, 0

    if means and medians:
        length = len(means)
        mean_of_means = np.mean(means)
        median_of_medians = np.median(medians)
        return length, mean_of_means, median_of_medians, zero_positions
    else:
        return 0, 0.0, 0.0, 0

# -----------------------
# Main script parameters
# -----------------------
chunk_id = "${chunk_id}"

timestamp(f"Processing chunk: {chunk_id}")

timestamp("Loading gene list...")
with open('${gene_chunk}') as f:
    genes = [line.strip() for line in f if line.strip()]

timestamp(f"Processing {len(genes)} genes from chunk {chunk_id}...")

timestamp("Loading database...")
db = gffutils.FeatureDB('${db_path}', keep_order=True)

timestamp("Opening gnomAD coverage files...")
exome_tabix = pysam.TabixFile('${exome_coverage}')
genome_tabix = pysam.TabixFile('${genome_coverage}')

# Output files
exon_output_file = f"{chunk_id}_exon_coverage.tsv"
intron_output_file = f"{chunk_id}_intron_coverage.tsv"

timestamp("Starting analysis...")
genes_processed = 0
genes_skipped = 0

with open(exon_output_file, 'w') as exon_out, open(intron_output_file, 'w') as intron_out:
    # Headers
    exon_header = ["GENE", "TRANSCRIPT", "EXON_NUMBER", "CHROM", "START", "END",
                   "LENGTH", "EXOME_MEAN_COV", "EXOME_MEDIAN_COV", "GENOME_MEAN_COV", "GENOME_MEDIAN_COV", 
                   "EXOME_ZERO_POS", "GENOME_ZERO_POS", "EXOME_ZERO_PERC", "GENOME_ZERO_PERC", "CODING_STATUS"]
    exon_out.write("\\t".join(exon_header) + "\\n")

    intron_header = ["GENE", "TRANSCRIPT", "INTRON_NUMBER", "CHROM", "START", "END",
                     "LENGTH", "EXOME_MEAN_COV", "EXOME_MEDIAN_COV", "GENOME_MEAN_COV", "GENOME_MEDIAN_COV",
                     "EXOME_ZERO_POS", "GENOME_ZERO_POS", "EXOME_ZERO_PERC", "GENOME_ZERO_PERC"]
    intron_out.write("\\t".join(intron_header) + "\\n")

    progress_threshold = len(genes) * 0.1 if len(genes) > 10 else 1
    next_prog = progress_threshold

    for i, gene_name in enumerate(genes, start=1):
        if i >= next_prog:
            pct = i / len(genes) * 100
            timestamp(f"Progress: {pct:.0f}% ({i}/{len(genes)})")
            next_prog += progress_threshold

        timestamp(f"Processing gene: {gene_name}")

        gene_feat, canonical_tx = find_gene_by_name(db, gene_name)
        if gene_feat is None or canonical_tx is None:
            timestamp(f"WARNING: Gene {gene_name} not found or missing canonical transcript, skipping...")
            genes_skipped += 1
            continue

        try:
            # --- Get exons ---
            exon_features = list(db.children(canonical_tx, featuretype='exon', order_by='start'))
            if not exon_features:
                timestamp(f"WARNING: No exon features found for {gene_name}, skipping...")
                genes_skipped += 1
                continue

            # Determine transcript strand
            transcript_strand = canonical_tx.strand
            timestamp(f"  Found {len(exon_features)} exon features on strand {transcript_strand}")

            # Collect CDS and UTR features for annotation
            cds_features = list(db.children(canonical_tx, featuretype='CDS', order_by='start'))
            utr_features = list(db.children(canonical_tx, featuretype='UTR', order_by='start'))
            
            # Use normalized chromosome for output consistency
            output_chrom = normalize_chromosome(canonical_tx.chrom)
            
            # Exon-level analysis for both exome and genome
            for exon in exon_features:
                exon_num = exon.attributes.get('exon_number', ['Unknown'])[0]

                # Determine coding/UTR status
                coding_status = "noncoding"
                for cds in cds_features:
                    if exon.start <= cds.end and exon.end >= cds.start:
                        coding_status = "coding"
                        break
                if coding_status == "noncoding":
                    for utr in utr_features:
                        if exon.start <= utr.end and exon.end >= utr.start:
                            coding_status = "UTR"
                            break
                
                # Extract coverage from both exome and genome data
                timestamp(f"  Processing exon {exon_num}: {exon.chrom}:{exon.start}-{exon.end}")
                timestamp(f"    Extracting exome coverage...")
                exome_length, exome_mean, exome_median, exome_zero_pos = extract_coverage_from_gnomad(
                    exon.chrom, exon.start, exon.end, exome_tabix)
                timestamp(f"    Extracting genome coverage...")
                genome_length, genome_mean, genome_median, genome_zero_pos = extract_coverage_from_gnomad(
                    exon.chrom, exon.start, exon.end, genome_tabix)
                
                # Use the maximum length found
                max_length = max(exome_length, genome_length)
                
                # Calculate percentages
                if max_length > 0:
                    exome_zero_perc = (exome_zero_pos / max_length) * 100
                    genome_zero_perc = (genome_zero_pos / max_length) * 100
                else:
                    exome_zero_perc = 0.0
                    genome_zero_perc = 0.0
                
                exon_out.write(
                    f"{gene_name}\\t{canonical_tx.id}\\t{exon_num}\\t{output_chrom}\\t{exon.start}\\t{exon.end}\\t"
                    f"{max_length}\\t{exome_mean:.2f}\\t{exome_median:.2f}\\t{genome_mean:.2f}\\t{genome_median:.2f}\\t"
                    f"{exome_zero_pos}\\t{genome_zero_pos}\\t{exome_zero_perc:.2f}\\t{genome_zero_perc:.2f}\\t{coding_status}\\n"
                )
                
                timestamp(f"    Exome - Mean: {exome_mean:.2f}, Median: {exome_median:.2f}, Zero: {exome_zero_pos}")
                timestamp(f"    Genome - Mean: {genome_mean:.2f}, Median: {genome_median:.2f}, Zero: {genome_zero_pos}")

            # Intron-level analysis
            if len(exon_features) > 1:
                # Create list of exon pairs with proper intron numbering based on strand
                if transcript_strand == '-':
                    exon_pairs = []
                    for idx in range(len(exon_features) - 1):
                        intron_num = len(exon_features) - idx - 1
                        exon_pairs.append((intron_num, exon_features[idx], exon_features[idx + 1]))
                else:
                    exon_pairs = []
                    for idx in range(len(exon_features) - 1):
                        intron_num = idx + 1
                        exon_pairs.append((intron_num, exon_features[idx], exon_features[idx + 1]))
                
                for intron_num, exon_current, exon_next in exon_pairs:
                    # Intron spans from end of current exon to start of next exon
                    intron_start = exon_current.end + 1
                    intron_end = exon_next.start - 1
                    
                    # Skip if invalid intron
                    if intron_end < intron_start:
                        timestamp(f"WARNING: Invalid intron coordinates for {gene_name} intron {intron_num}, skipping...")
                        continue
                    
                    timestamp(f"  Processing intron {intron_num}: {canonical_tx.chrom}:{intron_start}-{intron_end}")
                    timestamp(f"    Extracting exome coverage...")
                    exome_length, exome_mean, exome_median, exome_zero_pos = extract_coverage_from_gnomad(
                        canonical_tx.chrom, intron_start, intron_end, exome_tabix)
                    timestamp(f"    Extracting genome coverage...")
                    genome_length, genome_mean, genome_median, genome_zero_pos = extract_coverage_from_gnomad(
                        canonical_tx.chrom, intron_start, intron_end, genome_tabix)
                    
                    max_length = max(exome_length, genome_length)
                    
                    if max_length > 0:
                        exome_zero_perc = (exome_zero_pos / max_length) * 100
                        genome_zero_perc = (genome_zero_pos / max_length) * 100
                    else:
                        exome_zero_perc = 0.0
                        genome_zero_perc = 0.0
                    
                    intron_out.write(
                        f"{gene_name}\\t{canonical_tx.id}\\t{intron_num}\\t{output_chrom}\\t{intron_start}\\t{intron_end}\\t"
                        f"{max_length}\\t{exome_mean:.2f}\\t{exome_median:.2f}\\t{genome_mean:.2f}\\t{genome_median:.2f}\\t"
                        f"{exome_zero_pos}\\t{genome_zero_pos}\\t{exome_zero_perc:.2f}\\t{genome_zero_perc:.2f}\\n"
                    )
                    
                    timestamp(f"    Exome - Mean: {exome_mean:.2f}, Median: {exome_median:.2f}, Zero: {exome_zero_pos}")
                    timestamp(f"    Genome - Mean: {genome_mean:.2f}, Median: {genome_median:.2f}, Zero: {genome_zero_pos}")
                
                timestamp(f"  Analyzed {len(exon_features) - 1} introns")

            genes_processed += 1
            timestamp(f"Completed gene: {gene_name} ({len(exon_features)} exons, {max(0, len(exon_features) - 1)} introns)")

        except Exception as e:
            timestamp(f"ERROR processing gene {gene_name}: {str(e)}")
            genes_skipped += 1
            continue

timestamp("Analysis complete!")
timestamp(f"Summary: Processed {genes_processed} genes, skipped {genes_skipped} genes")
timestamp(f"Outputs: {exon_output_file} (exon-level), {intron_output_file} (intron-level)")

# Close tabix files
exome_tabix.close()
genome_tabix.close()
    """
}

process MergeGnomadCoverage {
    tag "${feature_type}"

    //publishDir "${params.output_dir}/gnomad_coverage", mode: 'copy'

    input:
    path(tsv_files)
    val(feature_type)

    output:
    path "gnomad_${feature_type}_coverage.tsv", emit: merged_coverage

    script:
    """
    # Get header from first file
    head -n 1 ${tsv_files[0]} > gnomad_${feature_type}_coverage.tsv

    # Append data from all files (skip headers)
    for file in ${tsv_files}; do
        tail -n +2 "\$file" >> gnomad_${feature_type}_coverage.tsv
    done

    echo "Merged ${tsv_files.size()} ${feature_type} coverage files"
    """
}

process MergeCohortData {
    tag "merge_${feature_type}_data"

    cpus = 2
    memory = '16 GB'
    time = '2 h'

    publishDir "${params.output_dir}/merged_coverage", mode: 'copy'

    input:
    path gnomad_output
    tuple val(platform1), path(cohort_output1)
    tuple val(platform2), path(cohort_output2)
    val feature_type

    output:
    path("gnomad_cohort_${feature_type}_coverage_merged.tsv"), emit: merged_coverage

    script:
    def feature_col = feature_type == 'exon' ? 'EXON_NUMBER' : 'INTRON_NUMBER'
    """
#!/usr/bin/env python3
import pandas as pd
import numpy as np

def merge_coverage_data():
    print("="*80, flush=True)
    print(f"Merging ${feature_type} coverage data", flush=True)
    print("="*80, flush=True)

    # Load gnomAD data
    print("\\nLoading gnomAD data: ${gnomad_output}", flush=True)
    gnomad_df = pd.read_csv("${gnomad_output}", sep="\\t")
    print(f"  Shape: {gnomad_df.shape}", flush=True)
    print(f"  Columns: {list(gnomad_df.columns)}", flush=True)

    # Load platform data
    print(f"\\nLoading ${platform1} data: ${cohort_output1}", flush=True)
    platform1_df = pd.read_csv("${cohort_output1}", sep="\\t")
    print(f"  Shape: {platform1_df.shape}", flush=True)
    print(f"  Columns: {list(platform1_df.columns)}", flush=True)

    print(f"\\nLoading ${platform2} data: ${cohort_output2}", flush=True)
    platform2_df = pd.read_csv("${cohort_output2}", sep="\\t")
    print(f"  Shape: {platform2_df.shape}", flush=True)
    print(f"  Columns: {list(platform2_df.columns)}", flush=True)

    # Define merge columns
    merge_cols = ['GENE', 'TRANSCRIPT', '${feature_col}', 'CHROM', 'START', 'END']

    # Columns to keep from cohort data (these will get platform suffixes)
    cohort_cols_to_keep = [
        '%ZERO_MEDIAN', '%ZERO_MIN', '%ZERO_MAX',
        '%MAPQ0_ONLY_MEDIAN', '%MAPQ0_ONLY_MIN', '%MAPQ0_ONLY_MAX',
        'MEDIAN_COV_MEDIAN', 'MEDIAN_COV_MIN', 'MEDIAN_COV_MAX'
    ]

    # Rename gnomAD columns
    print("\\nRenaming gnomAD columns...", flush=True)
    gnomad_df = gnomad_df.rename(columns={
        'EXOME_MEAN_COV': 'GNOMAD_EXOME_MEAN_COV',
        'EXOME_MEDIAN_COV': 'GNOMAD_EXOME_MEDIAN_COV',
        'GENOME_MEAN_COV': 'GNOMAD_GENOME_MEAN_COV',
        'GENOME_MEDIAN_COV': 'GNOMAD_GENOME_MEDIAN_COV',
        'EXOME_ZERO_POS': 'GNOMAD_EXOME_ZERO_POS',
        'GENOME_ZERO_POS': 'GNOMAD_GENOME_ZERO_POS',
        'EXOME_ZERO_PERC': 'GNOMAD_EXOME_ZERO_PERC',
        'GENOME_ZERO_PERC': 'GNOMAD_GENOME_ZERO_PERC'
    })

    # Prepare platform 1 data
    print("\\nPreparing platform-specific columns...", flush=True)

    # Columns to select from platform1: merge cols + LENGTH + CODING_STATUS (if exon) + cohort stats
    platform1_cols = merge_cols.copy()
    if 'LENGTH' in platform1_df.columns:
        platform1_cols.append('LENGTH')
    if "${feature_type}" == "exon" and 'CODING_STATUS' in platform1_df.columns:
        platform1_cols.append('CODING_STATUS')
    platform1_cols.extend(cohort_cols_to_keep)

    # Select columns
    platform1_subset = platform1_df[platform1_cols].copy()

    # Rename only the cohort stats columns (not LENGTH or CODING_STATUS)
    rename_dict1 = {col: f"{col}_${platform1}" for col in cohort_cols_to_keep}
    platform1_subset = platform1_subset.rename(columns=rename_dict1)

    # Prepare platform 2 data
    platform2_cols = merge_cols + cohort_cols_to_keep
    platform2_subset = platform2_df[platform2_cols].copy()

    # Rename cohort stats columns
    rename_dict2 = {col: f"{col}_${platform2}" for col in cohort_cols_to_keep}
    platform2_subset = platform2_subset.rename(columns=rename_dict2)

    # Perform merges
    print("\\nMerging datasets...", flush=True)

    # First merge: gnomAD + platform1
    merged_df = pd.merge(gnomad_df, platform1_subset, on=merge_cols, how='outer')
    print(f"  After gnomAD + ${platform1}: {merged_df.shape}", flush=True)

    # Second merge: add platform2
    merged_df = pd.merge(merged_df, platform2_subset, on=merge_cols, how='outer')
    print(f"  After adding ${platform2}: {merged_df.shape}", flush=True)

    # Calculate total dark regions
    print("\\nCalculating derived metrics...", flush=True)
    for platform in ['${platform1}', '${platform2}']:
        zero_col = f"%ZERO_MEDIAN_{platform}"
        mapq0_col = f"%MAPQ0_ONLY_MEDIAN_{platform}"
        total_dark_col = f"TOTAL_DARK_{platform}"

        if zero_col in merged_df.columns and mapq0_col in merged_df.columns:
            merged_df[total_dark_col] = (
                merged_df[zero_col].fillna(0) + merged_df[mapq0_col].fillna(0)
            )
            print(f"  Created {total_dark_col}", flush=True)

    # Calculate total dark difference
    merged_df['TOTAL_DARK_DIFF'] = (
        merged_df['TOTAL_DARK_${platform1}'].fillna(0) - 
        merged_df['TOTAL_DARK_${platform2}'].fillna(0)
    )
    print("  Created TOTAL_DARK_DIFF", flush=True)

    # Define column order
    base_cols = ['GENE', 'TRANSCRIPT', '${feature_col}']
    if "${feature_type}" == "exon":
        base_cols.append('CODING_STATUS')
    base_cols.extend(['CHROM', 'START', 'END', 'LENGTH'])

    final_column_order = base_cols + [
        'TOTAL_DARK_${platform1}', 'TOTAL_DARK_${platform2}', 'TOTAL_DARK_DIFF',
        'GNOMAD_EXOME_MEAN_COV', 'GNOMAD_EXOME_MEDIAN_COV',
        'GNOMAD_GENOME_MEAN_COV', 'GNOMAD_GENOME_MEDIAN_COV',
        'GNOMAD_EXOME_ZERO_POS', 'GNOMAD_GENOME_ZERO_POS',
        'GNOMAD_EXOME_ZERO_PERC', 'GNOMAD_GENOME_ZERO_PERC',
        '%ZERO_MEDIAN_${platform1}', '%ZERO_MIN_${platform1}', '%ZERO_MAX_${platform1}',
        '%MAPQ0_ONLY_MEDIAN_${platform1}', '%MAPQ0_ONLY_MIN_${platform1}', 
        '%MAPQ0_ONLY_MAX_${platform1}',
        'MEDIAN_COV_MEDIAN_${platform1}', 'MEDIAN_COV_MIN_${platform1}', 
        'MEDIAN_COV_MAX_${platform1}',
        '%ZERO_MEDIAN_${platform2}', '%ZERO_MIN_${platform2}', '%ZERO_MAX_${platform2}',
        '%MAPQ0_ONLY_MEDIAN_${platform2}', '%MAPQ0_ONLY_MIN_${platform2}', 
        '%MAPQ0_ONLY_MAX_${platform2}',
        'MEDIAN_COV_MEDIAN_${platform2}', 'MEDIAN_COV_MIN_${platform2}', 
        'MEDIAN_COV_MAX_${platform2}'
    ]

    # Reorder columns (keep only those that exist)
    available_columns = [col for col in final_column_order if col in merged_df.columns]
    merged_df = merged_df[available_columns]

    # Sort by genomic coordinates
    print("\\nSorting by genomic coordinates...", flush=True)
    merged_df = merged_df.sort_values([
        'CHROM', 'START', 'END', 'GENE', 'TRANSCRIPT', '${feature_col}'
    ])

    # Save results
    output_file = "gnomad_cohort_${feature_type}_coverage_merged.tsv"
    merged_df.to_csv(output_file, sep="\\t", index=False, float_format='%.2f')

    print("\\n" + "="*80, flush=True)
    print(f"✓ Merged ${feature_type} coverage written: {output_file}", flush=True)
    print(f"  Final shape: {merged_df.shape}", flush=True)
    print(f"  Columns: {list(merged_df.columns[:10])}...", flush=True)
    print("="*80, flush=True)

# Execute merge
try:
    merge_coverage_data()
except Exception as e:
    print(f"ERROR: {e}", flush=True)
    import traceback
    traceback.print_exc()
    raise
    """
}


process MergeGeneLevelData {
    tag "gene_level_merge"
    
    cpus = 2
    memory = '8 GB'
    time = '1 h'
    
    publishDir "${params.output_dir}/merged_coverage", mode: 'copy'
    
    input:
    tuple val(platform1), path(gene_output1)
    tuple val(platform2), path(gene_output2)
    
    output:
    path("gene_level_coverage_merged.tsv"), emit: merged_coverage
    
    script:
    """
#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os


def merge_gene_coverage():
    print("Loading input files...", flush=True)
    
    # Load platform 1 data
    print("  Reading ${platform1} gene file: ${gene_output1}", flush=True)
    platform1_df = pd.read_csv("${gene_output1}", sep="\\t")
    print(f"    ${platform1} shape: {platform1_df.shape}", flush=True)
    
    # Load platform 2 data
    print("  Reading ${platform2} gene file: ${gene_output2}", flush=True)
    platform2_df = pd.read_csv("${gene_output2}", sep="\\t")
    print(f"    ${platform2} shape: {platform2_df.shape}", flush=True)
    
    # Define merge columns
    merge_cols = ['GENE', 'TRANSCRIPT', 'CHROM', 'START', 'END']
    
    # Define columns to retain from each platform
    cohort_cols_to_keep = [
        '%ZERO_MEDIAN', '%ZERO_MIN', '%ZERO_MAX',
        '%MAPQ0_ONLY_MEDIAN', '%MAPQ0_ONLY_MIN', '%MAPQ0_ONLY_MAX',
        'MEDIAN_COV_MEDIAN', 'MEDIAN_COV_MIN', 'MEDIAN_COV_MAX'
    ]
    
    # Select and rename platform 1 columns
    print(f"Processing ${platform1} columns...", flush=True)
    platform1_subset = platform1_df[merge_cols + ['LENGTH'] + cohort_cols_to_keep].copy()
    platform1_rename_dict = {col: f"{col}_${platform1}" for col in cohort_cols_to_keep}
    platform1_subset = platform1_subset.rename(columns=platform1_rename_dict)
    
    # Select and rename platform 2 columns
    print(f"Processing ${platform2} columns...", flush=True)
    platform2_subset = platform2_df[merge_cols + cohort_cols_to_keep].copy()
    platform2_rename_dict = {col: f"{col}_${platform2}" for col in cohort_cols_to_keep}
    platform2_subset = platform2_subset.rename(columns=platform2_rename_dict)
    
    # Perform merge
    print("Merging datasets...", flush=True)
    merged_df = pd.merge(
        platform1_subset,
        platform2_subset,
        on=merge_cols,
        how='outer',
        suffixes=('', '_dup')
    )
    print(f"  After merge: {merged_df.shape}", flush=True)
    
    # Calculate total dark regions for each platform
    print("Calculating derived columns...", flush=True)
    
    # Platform 1 total dark
    platform1_zero_col = f"%ZERO_MEDIAN_${platform1}"
    platform1_mapq0_col = f"%MAPQ0_ONLY_MEDIAN_${platform1}"
    platform1_total_dark_col = f"TOTAL_DARK_${platform1}"
    
    if platform1_zero_col in merged_df.columns and platform1_mapq0_col in merged_df.columns:
        merged_df[platform1_total_dark_col] = (
            merged_df[platform1_zero_col].fillna(0) + merged_df[platform1_mapq0_col].fillna(0)
        )
        print(f"  Created {platform1_total_dark_col}", flush=True)
    
    # Platform 2 total dark
    platform2_zero_col = f"%ZERO_MEDIAN_${platform2}"
    platform2_mapq0_col = f"%MAPQ0_ONLY_MEDIAN_${platform2}"
    platform2_total_dark_col = f"TOTAL_DARK_${platform2}"
    
    if platform2_zero_col in merged_df.columns and platform2_mapq0_col in merged_df.columns:
        merged_df[platform2_total_dark_col] = (
            merged_df[platform2_zero_col].fillna(0) + merged_df[platform2_mapq0_col].fillna(0)
        )
        print(f"  Created {platform2_total_dark_col}", flush=True)
    
    # Calculate total dark difference
    if platform1_total_dark_col in merged_df.columns and platform2_total_dark_col in merged_df.columns:
        merged_df['TOTAL_DARK_DIFF'] = (
            merged_df[platform1_total_dark_col].fillna(0) - merged_df[platform2_total_dark_col].fillna(0)
        )
        print("  Created TOTAL_DARK_DIFF", flush=True)
    
    # Define the exact column order
    print("Reordering columns...", flush=True)
    final_column_order = [
        'GENE', 'TRANSCRIPT', 'CHROM', 'START', 'END', 'LENGTH',
        'TOTAL_DARK_${platform1}', 'TOTAL_DARK_${platform2}', 'TOTAL_DARK_DIFF',
        '%ZERO_MEDIAN_${platform1}', '%ZERO_MIN_${platform1}', '%ZERO_MAX_${platform1}',
        '%MAPQ0_ONLY_MEDIAN_${platform1}', '%MAPQ0_ONLY_MIN_${platform1}', '%MAPQ0_ONLY_MAX_${platform1}',
        'MEDIAN_COV_MEDIAN_${platform1}', 'MEDIAN_COV_MIN_${platform1}', 'MEDIAN_COV_MAX_${platform1}',
        '%ZERO_MEDIAN_${platform2}', '%ZERO_MIN_${platform2}', '%ZERO_MAX_${platform2}',
        '%MAPQ0_ONLY_MEDIAN_${platform2}', '%MAPQ0_ONLY_MIN_${platform2}', '%MAPQ0_ONLY_MAX_${platform2}',
        'MEDIAN_COV_MEDIAN_${platform2}', 'MEDIAN_COV_MIN_${platform2}', 'MEDIAN_COV_MAX_${platform2}'
    ]
    
    # Keep only columns that exist
    available_columns = [col for col in final_column_order if col in merged_df.columns]
    merged_df = merged_df[available_columns]
    
    print(f"  Kept {len(available_columns)} columns", flush=True)
    
    # Sort by genomic coordinates
    print("Sorting results by genomic coordinates...", flush=True)
    merged_df = merged_df.sort_values(['CHROM', 'START', 'END', 'GENE', 'TRANSCRIPT'])
    
    # Save merged results
    output_file = "gene_level_coverage_merged.tsv"
    merged_df.to_csv(output_file, sep="\\t", index=False, float_format='%.2f')
    
    print(f"✔ Merged coverage data written to: {output_file}", flush=True)
    print(f"  Final shape: {merged_df.shape}", flush=True)
    
    # Show sample of results
    print("\\nFirst few rows preview:", flush=True)
    preview_cols = ['GENE', 'TRANSCRIPT', 'TOTAL_DARK_${platform1}', 'TOTAL_DARK_${platform2}', 'TOTAL_DARK_DIFF']
    preview_available = [col for col in preview_cols if col in merged_df.columns]
    if preview_available:
        print(merged_df[preview_available].head().to_string(index=False))
    
    return merged_df


# Execute the merge
try:
    result_df = merge_gene_coverage()
    print("\\nMerge completed successfully!", flush=True)
except Exception as e:
    print(f"Error during merge: {e}", flush=True)
    raise
    """
}

process PlotStackedMappability {
    tag "mappability_plot"
    
    cpus = 1
    memory = '4 GB'
    time = '30 min'
    
    publishDir "${params.output_dir}/plots", mode: 'copy'
    
    input:
    path merged_gene_data
    
    output:
    path("*.png"), emit: plot
    
    script:
    """
    
    python ${project_dir}/coverage_analysis/plot_stacked_mappability.py \
        --input ${merged_gene_data} \
        --direction illumina \
        --top 20 \
        --gene-list EpilepsyGenes_NamesOnly_v2025-03.tsv
    """
}
