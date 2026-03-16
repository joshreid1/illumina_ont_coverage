// ============================================================================
// Coverage Analysis Workflow
// ============================================================================

nextflow.enable.dsl = 2

// ============================================================================
// PARAMETERS
// ============================================================================

params.gff_path = "gencode.v43.annotation.gtf.gz"
params.db_path = "gencode.v43.annotation.db"
params.gene_list = ""
params.output_dir = "results"
params.genes_per_job = 1000 // Number of genes to process per job

// Input parameters for both platforms with separate reference files
params.illumina_cram_list = "ngs_cram_list.txt"
params.nanopore_cram_list = "lrs_cram_list.txt"

params.illumina_ref_fasta = "Homo_sapiens_assembly38.fasta"
params.nanopore_ref_fasta = "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

// gnomAD coverage files
params.gnomad_exome_coverage = "gnomad.exomes.r2.1.1.coverage.summary.tsv.bgz"
params.gnomad_genome_coverage = "gnomad.genomes.r2.1.1.coverage.summary.tsv.bgz"


// ============================================================================
// IMPORT PROCESSES FROM MODULES
// ============================================================================

include { GencodePreprocessing } from './modules.nf'
include { SplitGeneList } from './modules.nf'
include { CoverageAnalysis_v3 } from './modules.nf'
include { ConcatenateResults } from './modules.nf'
include { SummarizeCoverage } from './modules.nf'
include { GnomadCoverageAnalysis } from './modules.nf'
include { MergeGeneLevelData } from './modules.nf'
include { PlotStackedMappability } from './modules.nf'

// Import MergeGnomadCoverage with aliases for exon and intron
include { MergeGnomadCoverage as MergeGnomadCoverage_Exon } from './modules.nf'
include { MergeGnomadCoverage as MergeGnomadCoverage_Intron } from './modules.nf'

// Import MergeCohortData with aliases for exon and intron
include { MergeCohortData as MergeExonData } from './modules.nf'
include { MergeCohortData as MergeIntronData } from './modules.nf'


// ============================================================================
// WORKFLOW
// ============================================================================

// Add these parameters at the top of your nextflow.config or script
params.use_existing_chunks = true
params.existing_chunks_dir = ""
params.use_existing_gnomad_chunks = true
params.existing_gnomad_chunks_dir = ""
params.skip_coverage_analysis = true
params.skip_gnomad_analysis = true

workflow {
    // ========================================
    // Input Validation
    // ========================================
    if (!params.skip_coverage_analysis && !params.use_existing_chunks) {
        if (!params.illumina_cram_list && !params.nanopore_cram_list) {
            error "At least one platform CRAM list must be provided (--illumina_cram_list or --nanopore_cram_list)"
        }

        if (params.illumina_cram_list && !params.illumina_ref_fasta) {
            error "Illumina reference FASTA is required when processing Illumina samples (--illumina_ref_fasta)"
        }

        if (params.nanopore_cram_list && !params.nanopore_ref_fasta) {
            error "Nanopore reference FASTA is required when processing Nanopore samples (--nanopore_ref_fasta)"
        }
    }

    // ========================================
    // Gene List Preparation (needed for fresh runs only)
    // ========================================
    if (!params.skip_coverage_analysis && !params.use_existing_chunks) {
        qualifying_genes = GencodePreprocessing(
            params.gff_path,
            params.db_path,
            params.gene_list
        )
        gene_chunks = SplitGeneList(qualifying_genes).flatten()
    }

    // ========================================
    // Coverage Analysis - Run Fresh or Load Existing Chunks
    // ========================================
    if (params.use_existing_chunks && params.existing_chunks_dir) {
        log.info "Loading existing chunk-level coverage files from: ${params.existing_chunks_dir}"

        // Load chunk-level files (platform, sample, chunk_file)
        // Expected format: {Platform}_{Sample}_genes_chunk_{XX}_{type}_coverage.tsv

        exon_chunks = Channel
            .fromPath("${params.existing_chunks_dir}/*_genes_chunk_*_exon_coverage.tsv")
            .map { file ->
                def basename = file.baseName
                // Parse: Platform_Sample_genes_chunk_XX_exon_coverage
                def parts = basename.split('_')
                def platform = parts[0]
                // Find where genes_chunk starts
                def genes_chunk_idx = parts.findIndexOf { it == 'genes' }
                // Everything between platform and genes_chunk is the sample
                def sample = parts[1..(genes_chunk_idx-1)].join('_')

                tuple(platform, sample, file)
            }

        intron_chunks = Channel
            .fromPath("${params.existing_chunks_dir}/*_genes_chunk_*_intron_coverage.tsv")
            .map { file ->
                def basename = file.baseName
                def parts = basename.split('_')
                def platform = parts[0]
                def genes_chunk_idx = parts.findIndexOf { it == 'genes' }
                def sample = parts[1..(genes_chunk_idx-1)].join('_')

                tuple(platform, sample, file)
            }

        gene_chunks = Channel
            .fromPath("${params.existing_chunks_dir}/*_genes_chunk_*_gene_coverage.tsv")
            .map { file ->
                def basename = file.baseName
                def parts = basename.split('_')
                def platform = parts[0]
                def genes_chunk_idx = parts.findIndexOf { it == 'genes' }
                def sample = parts[1..(genes_chunk_idx-1)].join('_')

                tuple(platform, sample, file)
            }

        // Group by platform and sample
        exon_by_sample = exon_chunks.groupTuple(by: [0, 1])
        intron_by_sample = intron_chunks.groupTuple(by: [0, 1])
        gene_by_sample = gene_chunks.groupTuple(by: [0, 1])

        // Join all three types for each sample
        sample_results = exon_by_sample
            .join(intron_by_sample, by: [0, 1])
            .join(gene_by_sample, by: [0, 1])
            .map { platform, sample, exon_files, intron_files, gene_files ->
                tuple(platform, sample, exon_files, intron_files, gene_files)
            }

        log.info "Loaded chunk-level files, proceeding to concatenation..."

    } else if (!params.skip_coverage_analysis) {
        log.info "Running full CRAM coverage analysis pipeline"

        // ========================================
        // Platform-Specific CRAM Channel Setup
        // ========================================
        illumina_crams = params.illumina_cram_list ? 
            Channel.fromPath(params.illumina_cram_list)
                .splitText()
                .map { it.trim() }
                .filter { it }
                .map { cram_path ->
                    def cram_file = file(cram_path)
                    def index_path = cram_path.endsWith('.cram') ? "${cram_path}.crai" : 
                                     cram_path.endsWith('.bam') ? "${cram_path}.bai" : 
                                     error("Unsupported file format: ${cram_path}")
                    def index_file = file(index_path)

                    if (!index_file.exists()) {
                        error "Index file not found: ${index_path}"
                    }

                    return tuple('Illumina', cram_file, index_file, file(params.illumina_ref_fasta))
                } : Channel.empty()

        nanopore_crams = params.nanopore_cram_list ? 
            Channel.fromPath(params.nanopore_cram_list)
                .splitText()
                .map { it.trim() }
                .filter { it }
                .map { cram_path ->
                    def cram_file = file(cram_path)
                    def index_path = cram_path.endsWith('.cram') ? "${cram_path}.crai" : 
                                     cram_path.endsWith('.bam') ? "${cram_path}.bai" : 
                                     error("Unsupported file format: ${cram_path}")
                    def index_file = file(index_path)

                    if (!index_file.exists()) {
                        error "Index file not found: ${index_path}"
                    }

                    return tuple('Nanopore', cram_file, index_file, file(params.nanopore_ref_fasta))
                } : Channel.empty()

        // ========================================
        // Coverage Analysis (Sample Level)
        // ========================================
        all_crams = illumina_crams.mix(nanopore_crams)

        cram_gene_combinations = all_crams
            .combine(gene_chunks)
            .map { platform, cram, index, ref, gene_chunk -> 
                tuple(platform, cram, index, ref, gene_chunk)
            }

        sample_coverage_results = CoverageAnalysis_v3(
            cram_gene_combinations,
            params.db_path
        )

        // ========================================
        // Group Results by Sample
        // ========================================
        exon_by_sample = sample_coverage_results.exon_coverage.groupTuple(by: [0, 1])
        intron_by_sample = sample_coverage_results.intron_coverage.groupTuple(by: [0, 1])
        gene_by_sample = sample_coverage_results.gene_coverage.groupTuple(by: [0, 1])

        sample_results = exon_by_sample
            .join(intron_by_sample, by: [0, 1])
            .join(gene_by_sample, by: [0, 1])
            .map { platform, sample, exon_files, intron_files, gene_files ->
                tuple(platform, sample, exon_files, intron_files, gene_files)
            }
    } else {
        error "Either run coverage analysis or provide existing chunk files with --use_existing_chunks"
    }

    // ========================================
    // Concatenate Chunks by Sample
    // ========================================
    concatenated_results = ConcatenateResults(sample_results)

    // ========================================
    // Cohort-Level Summary Statistics
    // ========================================
    exon_by_platform = concatenated_results.exon_combined.groupTuple(by: 0)
    intron_by_platform = concatenated_results.intron_combined.groupTuple(by: 0)
    gene_by_platform = concatenated_results.gene_combined.groupTuple(by: 0)

    cohort_summaries = SummarizeCoverage(
        exon_by_platform,
        intron_by_platform,
        gene_by_platform
    )

    // ========================================
    // gnomAD Coverage Analysis - Run or Load Existing Chunks
    // ========================================
    if (params.use_existing_gnomad_chunks && params.existing_gnomad_chunks_dir) {
        log.info "Loading existing gnomAD chunk files from: ${params.existing_gnomad_chunks_dir}"

        // Load chunk-level gnomAD files
        // Expected format: genes_chunk_{XX}_{type}_coverage.tsv
        gnomad_exon_chunks = Channel
            .fromPath("${params.existing_gnomad_chunks_dir}/genes_chunk_*_exon_coverage.tsv")
            .collect()

        gnomad_intron_chunks = Channel
            .fromPath("${params.existing_gnomad_chunks_dir}/genes_chunk_*_intron_coverage.tsv")
            .collect()

        log.info "Loaded gnomAD chunk files, proceeding to merge..."

    } else if (!params.skip_gnomad_analysis) {
        log.info "Running gnomAD coverage analysis"

        exome_files = Channel.of([
            file('gnomad.exomes.v4.0.coverage.summary.sorted.bed.gz'), 
            file('gnomad.exomes.v4.0.coverage.summary.sorted.bed.gz.tbi')
        ])

        genome_files = Channel.of([
            file('gnomad.genomes.r3.0.1.coverage.summary.sorted.bed.gz'), 
            file('gnomad.genomes.r3.0.1.coverage.summary.sorted.bed.gz.tbi')
        ])

        gnomad_gene_chunks = gene_chunks
            .combine(exome_files)
            .combine(genome_files)

        gnomad_coverage_results = GnomadCoverageAnalysis(
            gnomad_gene_chunks,
            params.db_path
        )

        gnomad_exon_chunks = gnomad_coverage_results.exon_coverage.collect()
        gnomad_intron_chunks = gnomad_coverage_results.intron_coverage.collect()

    } else {
        error "Either run gnomAD analysis or provide existing chunk files with --use_existing_gnomad_chunks"
    }

    // ========================================
    // Merge gnomAD Chunks
    // ========================================
    gnomad_exon_merged = MergeGnomadCoverage_Exon(
        gnomad_exon_chunks,
        'exon'
    )

    gnomad_intron_merged = MergeGnomadCoverage_Intron(
        gnomad_intron_chunks,
        'intron'
    )

    // ========================================
    // Merge gnomAD with Cohort Data - Exons
    // ========================================
    cohort_summaries.exon_summary
        .branch { platform, file ->
            nanopore: platform == 'Nanopore'
                return [platform, file]
            illumina: platform == 'Illumina'
                return [platform, file]
        }
        .set { exon_summaries_branched }

    exon_summaries_branched.nanopore.view { "Nanopore exon summary: $it" }
    exon_summaries_branched.illumina.view { "Illumina exon summary: $it" }

    exon_merged_output = MergeExonData(
        gnomad_exon_merged.merged_coverage,
        exon_summaries_branched.nanopore,
        exon_summaries_branched.illumina,
        'exon'
    )

    // ========================================
    // Merge gnomAD with Cohort Data - Introns
    // ========================================
    cohort_summaries.intron_summary
        .branch { platform, file ->
            nanopore: platform == 'Nanopore'
                return [platform, file]
            illumina: platform == 'Illumina'
                return [platform, file]
        }
        .set { intron_summaries_branched }

    intron_summaries_branched.nanopore.view { "Nanopore intron summary: $it" }
    intron_summaries_branched.illumina.view { "Illumina intron summary: $it" }

    intron_merged_output = MergeIntronData(
        gnomad_intron_merged.merged_coverage,
        intron_summaries_branched.nanopore,
        intron_summaries_branched.illumina,
        'intron'
    )

    // ========================================
    // Merge Gene-Level Data (No gnomAD)
    // ========================================
    cohort_summaries.gene_summary
        .branch { platform, file ->
            nanopore: platform == 'Nanopore'
                return [platform, file]
            illumina: platform == 'Illumina'
                return [platform, file]
        }
        .set { gene_summaries_branched }

    gene_merged_output = MergeGeneLevelData(
        gene_summaries_branched.nanopore,
        gene_summaries_branched.illumina
    )

    // ========================================
    // Visualization
    // ========================================
    PlotStackedMappability(gene_merged_output.merged_coverage)
}