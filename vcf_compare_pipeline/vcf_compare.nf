#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define input parameters with defaults if not provided
params.output_dir    = params.output_dir ?: 'results'
params.samples_csv   = params.samples_csv   ?: 'samples.csv'
params.plot_script   = 'plot_results.py'
params.region_bed_file = params.region_bed_file ?: 'regions.bed'
params.ngs_fasta     = params.ngs_fasta ?: 'ngs_reference.fasta'
params.lrs_fasta     = params.lrs_fasta ?: 'lrs_reference.fasta'
params.vep_cache     = params.vep_cache ?: '/path/to/vep_cache'

// ------------------
// Workflow Definition
// ------------------

//SNV
include { filter_variants } from './modules.nf'

include { gnomad as gnomad_ngs } from './modules.nf'
include { gnomad as gnomad_lrs } from './modules.nf'

include { vep_ngs } from './modules.nf'
include { vep_lrs } from './modules.nf'

include { CADD_Run_Container as CADD_Run_Container_ngs } from './modules.nf'
include { CADD_Run_Container as CADD_Run_Container_lrs } from './modules.nf'

include { Process_CADD as Process_CADD_ngs } from './modules.nf'
include { Process_CADD as Process_CADD_lrs } from './modules.nf'

include { strict_filter_variants } from './modules.nf'
include { g4e_filter_variants } from './modules.nf'

include { run_som } from './modules.nf'
include { run_hap } from './modules.nf'
include { parse_som_output } from './modules.nf'
include { collect_counts } from './modules.nf'
include { plot_results } from './modules.nf'

include { vcf_plot_all } from './modules.nf'
include { vcf_plot_filter } from './modules.nf'
include { vcf_plot_g4e } from './modules.nf'

//SV
include { filter_sv_variants } from './modules.nf'
include { vep_ngs_sv } from './modules.nf'
include { vep_lrs_sv } from './modules.nf'
include { strict_filter_sv_variants } from './modules.nf'
include { g4e_filter_sv_variants } from './modules.nf'
include { run_truvari } from './modules.nf'
include { parse_truvari_output } from './modules.nf'
include { collect_truvari_counts } from './modules.nf'
include { plot_sv_results } from './modules.nf'
include { vcf_sv_plot_all } from './modules.nf'
include { vcf_sv_plot_filter } from './modules.nf'
include { vcf_sv_plot_g4e } from './modules.nf'

//Summary
include { plot_summary_results } from './modules.nf'

workflow {

	//Single nucleotide variants/indels
	// Step 1: Create output directory
	file(params.output_dir).mkdirs()

	// Step 2: Read the samples CSV and create a channel of tuples: (sample_id, ngs_vcf, lrs_vcf)
	samples_ch = Channel.fromPath(params.snv_samples_csv)
				  .splitCsv(header: true)
				  .map { row -> tuple(row['Sample ID'], file(row['NGS_SNV_VCF']), file(row['LRS_SNV_VCF'])) }

	// Step 3: Filter variants and generate intermediate files with data type in filename
	filtered_ch = filter_variants(samples_ch)

	// Step 4: Split filtered files into separate channels for NGS and LRS with the data type included in the tuple
	filtered_ngs_ch = filtered_ch.map { sample_id, ngs_vcf, lrs_vcf -> tuple(sample_id, "NGS", ngs_vcf) }
	filtered_lrs_ch = filtered_ch.map { sample_id, ngs_vcf, lrs_vcf -> tuple(sample_id, "LRS", lrs_vcf) }

	// Then run them on each channel
	processed_ngs_ch = gnomad_ngs(filtered_ngs_ch) | vep_ngs
	processed_lrs_ch = gnomad_lrs(filtered_lrs_ch) | vep_lrs


	// Pair processed NGS and LRS files based on sample_id
	paired_vcfs_ch = processed_ngs_ch.join(processed_lrs_ch, by: 0)
		.map { sample_id, ngs_type, ngs_vcf, lrs_type, lrs_vcf -> 
			tuple(sample_id, ngs_vcf, lrs_vcf)
		}


	strict_ch = strict_filter_variants(paired_vcfs_ch) 
	g4e_ch = g4e_filter_variants(strict_ch)

	vcf_plot_all(paired_vcfs_ch)
	vcf_plot_filter(strict_ch)
	vcf_plot_g4e(g4e_ch)



	//Structural variants
	sv_samples_ch = Channel.fromPath(params.sv_samples_csv)
				  .splitCsv(header: true)
				  .map { row -> tuple(row['Sample ID'], file(row['NGS_SV_VCF']), file(row['LRS_SV_VCF'])) }

	// Step 3: Filter variants and generate intermediate files with data type in filename
	filtered_sv_ch = filter_sv_variants(sv_samples_ch)

	// Step 4: Split filtered files into separate channels for NGS and LRS with the data type included in the tuple
	filtered_ngs_sv_ch = filtered_sv_ch.map { sample_id, ngs_vcf, lrs_vcf -> tuple(sample_id, "NGS", ngs_vcf, "${ngs_vcf}.tbi") }
	filtered_lrs_sv_ch = filtered_sv_ch.map { sample_id, ngs_vcf, lrs_vcf -> tuple(sample_id, "LRS", lrs_vcf, "${lrs_vcf}.tbi") }

	// Then run them on each channel
	processed_ngs_sv_ch = vep_ngs_sv(filtered_ngs_sv_ch) 
	processed_lrs_sv_ch = vep_lrs_sv(filtered_lrs_sv_ch) 


	// Pair processed NGS and LRS files based on sample_id
	paired_sv_vcfs_ch = processed_ngs_sv_ch.join(processed_lrs_sv_ch, by: 0)
		.map { sample_id, ngs_type, ngs_vcf, lrs_type, lrs_vcf -> 
			tuple(sample_id, ngs_vcf, lrs_vcf)
		}


	strict_sv_ch = strict_filter_sv_variants(paired_sv_vcfs_ch) 
	g4e_sv_ch = g4e_filter_sv_variants(strict_sv_ch)

	vcf_sv_plot_all(paired_sv_vcfs_ch)
	vcf_sv_plot_filter(strict_sv_ch)
	vcf_sv_plot_g4e(g4e_sv_ch)

	plot_summary_results(vcf_plot_all.out.collected_counts,
							vcf_plot_filter.out.collected_counts,
							vcf_plot_g4e.out.collected_counts,
							vcf_sv_plot_all.out.collected_counts,
							vcf_sv_plot_filter.out.collected_counts,
							vcf_sv_plot_g4e.out.collected_counts)
}


