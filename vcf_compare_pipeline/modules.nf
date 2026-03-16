/////////////////////
//SINGLE NUCLEOTIDE VARIANTS AND INDELS
/////////////////////

//Process 1: Filter Variants
process filter_variants {
	module 'bcftools'
	module 'gatk'

	memory = '8 GB'
	 
	input:
		tuple val(sample_id), path(ngs_vcf), path(lrs_vcf)

	output:
		tuple val(sample_id), path("${sample_id}_NGS_size_filtered.vcf.gz"), path("${sample_id}_LRS_size_filtered.vcf.gz")

	script:
	"""
	# Index the input VCFs
	tabix ${ngs_vcf}
	tabix ${lrs_vcf}

	##Filters: Allele frequency >= 0.20, Genotype quality >= 20, Depth >= 5

	### Process NGS VCF ###
	bcftools norm -c w --fasta-ref ${params.ngs_fasta} --multiallelics '-' ${ngs_vcf} \
	| bcftools view -i "FORMAT/AD[0:1]/FORMAT/DP >= 0.20 && FORMAT/GQ >= 20 && FORMAT/DP >= 5" -Oz -o ${sample_id}_NGS_filtered.vcf.gz
	tabix ${sample_id}_NGS_filtered.vcf.gz
	gatk SelectVariants -V ${sample_id}_NGS_filtered.vcf.gz --max-indel-size 50 --output ${sample_id}_NGS_size_filtered.vcf
	bgzip ${sample_id}_NGS_size_filtered.vcf


	### Process LRS VCF ###
	bcftools norm -c w --fasta-ref ${params.lrs_fasta} --multiallelics '-' ${lrs_vcf} \
	| bcftools view -i "FORMAT/AF >= 0.20 && FORMAT/GQ >= 20 && FORMAT/DP >= 5" -Oz -o ${sample_id}_LRS_filtered.vcf.gz
	tabix ${sample_id}_LRS_filtered.vcf.gz
	gatk SelectVariants -V ${sample_id}_LRS_filtered.vcf.gz --max-indel-size 50 --output ${sample_id}_LRS_size_filtered.vcf
	bgzip ${sample_id}_LRS_size_filtered.vcf
	"""
}

// Process 2: Annotate with gnomAD
process gnomad {
	module 'bcftools/1.20'

	cpus = 4
	memory = '16 GB'

	input:
		tuple val(sample_id), val(data_type), path(vcf)
	
	output:
		tuple val(sample_id), val(data_type), path("${sample_id}_${data_type}_gnomad.vcf.gz"), path("${sample_id}_${data_type}_gnomad.vcf.gz.tbi")

	script:
	def vcfanno = "/path/to/vcfanno"
	"""
	# Set vcfanno path as a variable
	VCFANNO="${vcfanno}"

	# First, sort the VCF and write an index
	#bcftools sort --temp-dir tmp ${vcf} -O z -o ${sample_id}_${data_type}_filtered.vcf.gz --write-index=tbi 

	# Annotate using vcfanno with two steps
	#\$VCFANNO -p ${task.cpus} ${projectDir}/vcf_compare_pipeline/gnomad_v4.toml ${sample_id}_${data_type}_filtered.vcf.gz > ${sample_id}_${data_type}_gnomad_v1.vcf
	#\$VCFANNO -p ${task.cpus} ${projectDir}/vcf_compare_pipeline/gnomad_v4_postprocess.toml ${sample_id}_${data_type}_gnomad_v1.vcf > ${sample_id}_${data_type}_gnomad.vcf

	\$VCFANNO -p ${task.cpus} ${projectDir}/vcf_compare_pipeline/gnomad_v4.toml ${vcf} > ${sample_id}_${data_type}_gnomad_v1.vcf
	\$VCFANNO -p ${task.cpus} ${projectDir}/vcf_compare_pipeline/gnomad_v4_postprocess.toml ${vcf} > ${sample_id}_${data_type}_gnomad.vcf

	bgzip ${sample_id}_${data_type}_gnomad.vcf
	tabix ${sample_id}_${data_type}_gnomad.vcf.gz 
	"""
}

// Process 3: Variant Effect Predictor (VEP)
process vep_ngs {
	module 'ensembl-vep/112'

	cpus = 4
	memory = '16 GB'
	
	input:
		tuple val(sample_id), val(data_type), path(vcf), path(vcf_index)

	output:
		tuple val(sample_id), val(data_type), path("${sample_id}_${data_type}_vep.vcf")

	shell:
	'''
	output_file=!{sample_id}_!{data_type}_vep.vcf
	vep --cache --dir !{params.vep_cache} --cache_version 104 --assembly GRCh38 \
	-i !{vcf} -o $output_file --format vcf --vcf --symbol --terms SO --tsl --hgvs \
	--fasta !{params.ngs_fasta} --offline --sift b --polyphen b --ccds --hgvs --hgvsg --symbol \
	--numbers --protein --variant_class --pick_allele_gene --no_stats --fork 4
	'''
}

process vep_lrs {
	module 'ensembl-vep/112'

	cpus = 4
	memory = '16 GB'
	
	input:
		tuple val(sample_id), val(data_type), path(vcf), path(vcf_index)

	output:
		tuple val(sample_id), val(data_type), path("${sample_id}_${data_type}_vep.vcf")

	shell:
	'''
	output_file=!{sample_id}_!{data_type}_vep.vcf
	vep --cache --dir !{params.vep_cache} --cache_version 104 --assembly GRCh38 \
	-i !{vcf} -o $output_file --format vcf --vcf --symbol --terms SO --tsl --hgvs \
	--fasta !{params.lrs_fasta} --offline --sift b --polyphen b --ccds --hgvs --hgvsg --symbol \
	--numbers --protein --variant_class --pick_allele_gene --no_stats --fork 4
	'''
}

//Process 4: Strict filtering of variants based on VEP annotations
process strict_filter_variants {

	input:
		tuple val(sample_id), path(ngs_vcf), path(lrs_vcf)
   

	output:
		tuple val(sample_id), path("${sample_id}_ngs_filtered.vcf"), path("${sample_id}_lrs_filtered.vcf")

	shell:
	'''
	grep -E "^#|MODERATE|HIGH" !{ngs_vcf} > !{sample_id}_ngs_filtered.vcf

	grep -E "^#|MODERATE|HIGH" !{lrs_vcf} > !{sample_id}_lrs_filtered.vcf
	'''
}

//Process 5: Filter variants to those in Genes4Epilepsy list
process g4e_filter_variants {

	input:
		tuple val(sample_id), path(ngs_vcf), path(lrs_vcf)
   

	output:
		tuple val(sample_id), path("${sample_id}_ngs_g4e.vcf"), path("${sample_id}_lrs_g4e.vcf")

	shell:
	'''
	#NGS
	bgzip -c !{ngs_vcf} > ngs.vcf.gz
	tabix ngs.vcf.gz

	bcftools view --regions-file ${projectDir}/vcf_compare_pipeline/EpilepsyGenes_v2024-09.bed \
	 ngs.vcf.gz -O v -o !{sample_id}_ngs_g4e.vcf

	#LRS
	bgzip -c !{lrs_vcf} > lrs.vcf.gz
	tabix lrs.vcf.gz

	bcftools view --regions-file ${projectDir}/vcf_compare_pipeline/EpilepsyGenes_v2024-09.bed \
	 lrs.vcf.gz -O v -o !{sample_id}_lrs_g4e.vcf
	'''
}


process CADD_Run_Container {
	cpus = 1
	memory = { 10 * task.attempt + ' GB' }
	time = { 8 * task.attempt + ' h'}

	container 'cadd-scoring_latest.sif'
	containerOptions '-B /path/to/CADD-scripts-master/data/annotations:/CADD-scripts/data/annotations --writable-tmpfs'
	// Requires download of CADD-scripts and annotation data to a local directory, and bind mounting that directory into the container as shown above. See CADD documentation for details.


	input:
		tuple val(sample_id), val(data_type), path(vcf)
		
	output:
		tuple val(sample_id), val(data_type), path(vcf), path("${sample_id}_${data_type}_cadd_run.tsv.gz")
		
	shell:
	'''
	#Link local software
	ln -s /path/to/snakemake/8.11.3/bin/snakemake /usr/local/bin/snakemake
	ln -s /path/to/bcftools/1.20/bin/bcftools /usr/local/bin/bcftools

	#Current version requires symlink snakemake and bcftools into the container

	# Run CADD
	if [[ $(bcftools query -f '%ALT\n' "!{vcf}" | uniq) == "*" ]]; then
		cp "!{vcf}" "!{sample_id}_!{data_type}.cadd_run.vcf"
	else
		zgrep -v "^#" "!{vcf}" | cut -f 1-5 | sed 's/^chr//' > "!{sample_id}_!{data_type}_cadd_run.vcf"
		/CADD-scripts/CADD.sh "!{sample_id}_!{data_type}_cadd_run.vcf"
	fi
	'''
}

process Process_CADD {
	cpus = 4
	memory = { 10 * task.attempt + ' GB' }
	time = { 1 * task.attempt + ' h'}
		
	input:
		tuple val(sample_id), val(data_type), path(vcf), path(cadd_tsv) 
		
	output:
		tuple val(sample_id), val(data_type), path("${sample_id}_${data_type}.cadd.vcf")
	
	shell:
	'''
	cp "$(readlink -f !{cadd_tsv})" ./cadd_tmp

	tabix -f -b 2 -e 2 -s 1 ./cadd_tmp

	cadd_path=$(realpath ./cadd_tmp)
	echo '[[annotation]]' > cadd.toml
	echo "file= '${cadd_path}'" >> cadd.toml
	echo 'names=["CADD_Score"]' >> cadd.toml
	echo 'ops=["mean"]' >> cadd.toml
	echo 'columns=[6]' >> cadd.toml

	vcfanno -p ${task.cpus} cadd.toml !{vcf} > !{sample_id}_!{data_type}_cadd.vcf 
	'''
}

process run_som {
	tag "Running som.py for sample ${sample_id}"

	container 'hap.py.sif'
	
	input:
		tuple val(sample_id), path(ngs_vcf), path(lrs_vcf)

	output:
		tuple val(sample_id), path("som_output_${sample_id}_NGS_vs_LRS.txt.stats.csv")

	script:
	"""
	/opt/hap.py/bin/som.py -o som_output_${sample_id}_NGS_vs_LRS.txt \
	  -r ${params.ngs_fasta} \
	  --normalize-truth --normalize-query ${ngs_vcf} ${lrs_vcf}
	"""
}

process run_hap {
	tag "Running hap.py for sample ${sample_id}"

	cpus = 4
	memory = { 10 * task.attempt + ' GB' }
	time = { 2 * task.attempt + ' h'}

	container 'hap.py.sif'
	
	input:
		tuple val(sample_id), path(ngs_vcf), path(lrs_vcf)

	output:
		tuple val(sample_id), path("hap_output_${sample_id}_NGS_vs_LRS.txt.summary.csv")

	script:
	"""
	/opt/hap.py/bin/hap.py -o hap_output_${sample_id}_NGS_vs_LRS.txt \
	  -r ${params.ngs_fasta}  \
	  ${ngs_vcf} ${lrs_vcf}
	"""
}

process parse_som_output {
	tag "Parsing som.py output for sample ${sample_id}"

	input:
		tuple val(sample_id), path(som_output)

	output:
		tuple val(sample_id), path("counts_${sample_id}_NGS_vs_LRS.txt")

	script:
	"""
	$NXF_PYTHON ${params.projectDir}/vcf_compare_pipeline/scripts/parse_som_output.py ${som_output} counts_${sample_id}_NGS_vs_LRS.txt
	"""
}

process collect_counts {
	tag "Collecting counts for plotting"

	input:
		path counts_files

	output:
		path "aggregated_counts.csv"

	script:
	"""
	echo "Sample ID,Variant,TP,FP,FN" > aggregated_counts.csv

	for file in ${counts_files}; do
		if [[ -r "\$file" ]]; then
			sample_id=\$(basename "\$file" | sed 's/^counts_\\(.*\\)_NGS_vs_LRS\\.txt/\\1/')
			tail -n +2 "\$file" | while read -r variant tp fp fn; do
				if [[ -n "\$variant" && -n "\$tp" && -n "\$fp" && -n "\$fn" ]]; then
					echo "\${sample_id},\${variant},\${tp},\${fp},\${fn}" >> aggregated_counts.csv
				else
					echo "Warning: Malformed line in \$file: \$variant \$tp \$fp \$fn"
				fi
			done
		else
			echo "Error: Cannot read file \$file"
		fi
	done
	"""
}

process plot_results {
	tag "Plotting results"

	input:
		path aggregated_counts
		val type

	script:
	"""
	which python
	##SNV
	$NXF_PYTHON ${params.plot_script} --counts ${aggregated_counts} --variant SNVS --output ${params.output_dir}/${type}_snvs_results_plot.png --title "${type} SNVs (HaplotypeCaller vs. Clair3)"

	##INDELS
	$NXF_PYTHON ${params.plot_script} --counts ${aggregated_counts} --variant INDELS --output ${params.output_dir}/${type}_indels_results_plot.png --title "${type} Indels (HaplotypeCaller vs. Clair3)"
	"""
}

process plot_sv_results {
	tag "Plotting results"

	input:
		path aggregated_counts
		val type

	script:
	"""
	##Structural Variants
	$NXF_PYTHON ${params.plot_script} --counts ${aggregated_counts} --variant RECORDS --output ${params.output_dir}/${type}_sv_results_plot.png --title "${type} Structural Variants (Manta vs. Sniffles2)"
	"""
}


workflow vcf_plot_all {
	take:
	input_tuple

	main:
	som_results_ch = run_som(input_tuple) 
	parsed_som_ch = parse_som_output(som_results_ch)
	counts_files = parsed_som_ch[0].map { it[1] }.collect()
	collected_counts = collect_counts(counts_files)
	plot_results(collected_counts, 'All')

	emit: collected_counts
}

workflow vcf_plot_filter {
	take:
	input_tuple

	main:
	som_results_ch = run_som(input_tuple) 
	parsed_som_ch = parse_som_output(som_results_ch)
	counts_files = parsed_som_ch[0].map { it[1] }.collect()
	collected_counts = collect_counts(counts_files)
	plot_results(collected_counts, 'VEP_HighModerate')

	emit: collected_counts
}

workflow vcf_plot_g4e {
	take:
	input_tuple

	main:
	som_results_ch = run_som(input_tuple) 
	parsed_som_ch = parse_som_output(som_results_ch)
	counts_files = parsed_som_ch[0].map { it[1] }.collect()
	collected_counts = collect_counts(counts_files)
	plot_results(collected_counts, 'Genes4Epilepsy')

	emit: collected_counts
}


/////////////////////
//STRUCTURAL VARIANTS
/////////////////////

// Process 1: Filter Variants 
process filter_sv_variants {
	module 'bcftools'
	module 'gatk'

	memory = '8 GB' 

	input:
		tuple val(sample_id), path(ngs_vcf), path(lrs_vcf)

	output:
		tuple val(sample_id), path("${sample_id}_NGS_size_filtered.vcf.gz"), path("${sample_id}_LRS_size_filtered.vcf.gz")

	script:
	"""
	# Index the input VCFs
	tabix ${ngs_vcf}
	tabix ${lrs_vcf}

	##Filters: Genotype quality >= 20, Depth >= 5

	### Process NGS VCF ###
	bcftools view ${ngs_vcf} -i 'FILTER="PASS" && FORMAT/GQ >= 20' | bcftools view -e 'INFO/SVTYPE=="BND"' -Oz -o ${sample_id}_NGS_filtered.vcf.gz
	tabix ${sample_id}_NGS_filtered.vcf.gz
	gatk SelectVariants -V ${sample_id}_NGS_filtered.vcf.gz --min-indel-size 51 --output ${sample_id}_NGS_size_filtered.vcf
	bgzip ${sample_id}_NGS_size_filtered.vcf

	### Process LRS VCF ###
	bcftools view ${lrs_vcf} -i 'FILTER="PASS" && FORMAT/GQ >= 20' -Oz -o ${sample_id}_LRS_filtered.vcf.gz
	tabix ${sample_id}_LRS_filtered.vcf.gz
	gatk SelectVariants -V ${sample_id}_LRS_filtered.vcf.gz --min-indel-size 51 --output ${sample_id}_LRS_size_filtered.vcf
	bgzip ${sample_id}_LRS_size_filtered.vcf

	"""
}

process vep_ngs_sv {
	module 'ensembl-vep/112'

	cpus = 4
	memory = '16 GB'
	
	input:
		tuple val(sample_id), val(data_type), path(vcf), path(vcf_index)

	output:
		tuple val(sample_id), val(data_type), path("${sample_id}_${data_type}_vep.vcf")

	shell:
	'''
	output_file=!{sample_id}_!{data_type}_vep.vcf
	vep --cache --dir !{params.vep_cache} --cache_version 104 --assembly GRCh38 \
	-i !{vcf} -o $output_file --format vcf --vcf --symbol --terms SO --tsl --hgvs \
	--fasta !{params.ngs_fasta} --offline --sift b --polyphen b --ccds --hgvs --hgvsg --symbol \
	--numbers --protein --variant_class --pick_allele_gene --no_stats --fork 4
	'''
}

process vep_lrs_sv {
	module 'ensembl-vep/112'

	cpus = 4
	memory = '16 GB'
	
	input:
		tuple val(sample_id), val(data_type), path(vcf), path(vcf_index)

	output:
		tuple val(sample_id), val(data_type), path("${sample_id}_${data_type}_vep.vcf")

	shell:
	'''
	output_file=!{sample_id}_!{data_type}_vep.vcf
	vep --cache --dir !{params.vep_cache} --cache_version 104 --assembly GRCh38 \
	-i !{vcf} -o $output_file --format vcf --vcf --symbol --terms SO --tsl --hgvs \
	--fasta !{params.lrs_fasta} --offline --sift b --polyphen b --ccds --hgvs --hgvsg --symbol \
	--numbers --protein --variant_class --pick_allele_gene --no_stats --fork 4
	'''
}

process strict_filter_sv_variants {

	input:
		tuple val(sample_id), path(ngs_vcf), path(lrs_vcf)
   

	output:
		tuple val(sample_id), path("${sample_id}_ngs_filtered.vcf"), path("${sample_id}_lrs_filtered.vcf")

	shell:
	'''
	grep -E "^#|MODERATE|HIGH" !{ngs_vcf} > !{sample_id}_ngs_filtered.vcf

	grep -E "^#|MODERATE|HIGH" !{lrs_vcf} > !{sample_id}_lrs_filtered.vcf
	'''
}


process g4e_filter_sv_variants {

	input:
		tuple val(sample_id), path(ngs_vcf), path(lrs_vcf)
   

	output:
		tuple val(sample_id), path("${sample_id}_ngs_g4e.vcf"), path("${sample_id}_lrs_g4e.vcf")

	shell:
	'''
	#NGS#
	bgzip -c !{ngs_vcf} > ngs.vcf.gz
	tabix ngs.vcf.gz

	bcftools view --regions-file ${params.projectDir}/EpilepsyGenes_v2024-09.bed \
	 ngs.vcf.gz -O v -o !{sample_id}_ngs_g4e.vcf

	#LRS#
	bgzip -c !{lrs_vcf} > lrs.vcf.gz
	tabix lrs.vcf.gz

	bcftools view --regions-file ${params.projectDir}/EpilepsyGenes_v2024-09.bed \
	 lrs.vcf.gz -O v -o !{sample_id}_lrs_g4e.vcf
	'''
}



// Process 6: Run truvari to compare NGS and LRS variants
process run_truvari {
	tag "Running truvari for sample ${sample_id}"
   
	input:
		tuple val(sample_id), path(ngs_vcf), path(lrs_vcf)

	output:
		tuple val(sample_id), path("output/summary.json")

	script:
	"""
	# Remove problematic config header lines
	cat ${ngs_vcf} | \
	# 2. Filter out unwanted ##contig lines
	awk 'BEGIN{IGNORECASE=1}
		/^##contig/ && /_decoy|_alt|_random|chrEBV|chrUn|HLA/ { next }
		{ print }' > "${ngs_vcf}_edit.vcf"

	bgzip  "${ngs_vcf}_edit.vcf"
	bgzip ${lrs_vcf}

	tabix  "${ngs_vcf}_edit.vcf.gz"
	tabix "${lrs_vcf}.gz"

	truvari bench -b "${ngs_vcf}_edit.vcf.gz" -c "${lrs_vcf}.gz"  -o output --sizemin 50 --sizefilt 30 --sizemax 5000000
	"""
}

// Process 7: Parse som.py output
process parse_truvari_output {
	tag "Parsing truari output for sample ${sample_id}"

	input:
		tuple val(sample_id), path(truvari_output)

	output:
		tuple val(sample_id), path("counts_${sample_id}_NGS_vs_LRS.txt")

	script:
	"""
	$NXF_PYTHON ${projectDir}/vcf_compare_pipeline/scripts/parse_truvari_output.py ${truvari_output} counts_${sample_id}_NGS_vs_LRS.txt
	"""
}

// Process 8: Collect counts from all parsed outputs
process collect_truvari_counts {
	tag "Collecting counts for plotting"

	input:
		path counts_files

	output:
		path "aggregated_counts.csv"

	script:
	"""
	echo "Sample ID,Variant,TP,FP,FN" > aggregated_counts.csv

	for file in ${counts_files}; do
		if [[ -r "\$file" ]]; then
			sample_id=\$(basename "\$file" | sed 's/^counts_\\(.*\\)_NGS_vs_LRS\\.txt/\\1/')
			tail -n +2 "\$file" | while read -r variant tp fp fn; do
				if [[ -n "\$variant" && -n "\$tp" && -n "\$fp" && -n "\$fn" ]]; then
					echo "\${sample_id},\${variant},\${tp},\${fp},\${fn}" >> aggregated_counts.csv
				else
					echo "Warning: Malformed line in \$file: \$variant \$tp \$fp \$fn"
				fi
			done
		else
			echo "Error: Cannot read file \$file"
		fi
	done
	"""
}

workflow vcf_sv_plot_all {
	take:
	input_tuple

	main:
	truvari_results_ch = run_truvari(input_tuple)

	parsed_truvari_ch = parse_truvari_output(truvari_results_ch)
	counts_files = parsed_truvari_ch[0].map { it[1] }.collect()
	collected_counts = collect_counts(counts_files)
	plot_sv_results(collected_counts, 'All')

	emit: collected_counts
}


workflow vcf_sv_plot_filter {
	take:
	input_tuple

	main:
	truvari_results_ch = run_truvari(input_tuple)

	parsed_truvari_ch = parse_truvari_output(truvari_results_ch)
	counts_files = parsed_truvari_ch[0].map { it[1] }.collect()
	collected_counts = collect_counts(counts_files)
	plot_sv_results(collected_counts, 'VEP_HighModerate')

	emit: collected_counts
}



workflow vcf_sv_plot_g4e {
	take:
	input_tuple

	main:
	truvari_results_ch = run_truvari(input_tuple)

	parsed_truvari_ch = parse_truvari_output(truvari_results_ch)
	counts_files = parsed_truvari_ch[0].map { it[1] }.collect()
	collected_counts = collect_counts(counts_files)
	plot_sv_results(collected_counts, 'Genes4Epilepsy')

	emit: collected_counts
}


// Process 10: Plot summary results
process plot_summary_results {
	tag "Plotting summary results"

	publishDir "${params.output_dir}", mode: 'copy'

	input:
		path snv_indel_all_aggregated_counts, stageAs: 'snv_indel_all_aggregated_counts.csv'
		path snv_indel_filter_aggregated_counts, stageAs: 'snv_indel_filter_aggregated_counts.csv'
		path snv_indel_g4e_aggregated_counts, stageAs: 'snv_indel_g4e_aggregated_counts.csv'
		path sv_all_aggregated_counts, stageAs: 'sv_all_aggregated_counts.csv'
		path sv_filter_aggregated_counts, stageAs: 'sv_filter_aggregated_counts.csv'
		path sv_g4e_aggregated_counts, stageAs: 'sv_g4e_aggregated_counts.csv'

	output:
		path "combined_results_plot.png"

	script:
	"""
	plot_script="${projectDir}/vcf_compare_pipeline/scripts/plot_summary_results.py"
	$NXF_PYTHON \$plot_script \
	--snv-indel-all ${snv_indel_all_aggregated_counts} \
	--snv-indel-filter ${snv_indel_filter_aggregated_counts} \
	--snv-indel-g4e ${snv_indel_g4e_aggregated_counts} \
	--sv-all ${sv_all_aggregated_counts} \
	--sv-filter ${sv_filter_aggregated_counts} \
	--sv-g4e ${sv_g4e_aggregated_counts} \
	--output combined_results_plot.png \
	--title "Comprehensive Variant Analysis"
	"""
}