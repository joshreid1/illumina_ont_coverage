#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import argparse
import sys

def calculate_median_for_file(counts_file, variant):
	"""
	Calculate median values for a specific variant type from a counts file.


	Parameters:
	- counts_file: Path to the aggregated_counts.csv file.
	- variant: Variant type to filter on (e.g., SNVS, INDELS).


	Returns:
	- Series with median values for ['Both', 'Illumina Only', 'Nanopore Only']
	"""
	try:
		df = pd.read_csv(counts_file)
	except Exception as e:
		print(f"Error reading {counts_file}: {e}")
		return None



	# Check if necessary columns exist
	required_columns = {'Sample ID', 'Variant', 'TP', 'FP', 'FN'}
	if not required_columns.issubset(df.columns):
		print(f"File {counts_file} must contain columns: {required_columns}")
		return None



	# Filter by variant
	df = df[df['Variant'] == variant]


	if df.empty:
		print(f"No data found for variant {variant} in {counts_file}")
		return None



	# Rename columns
	df = df.rename(columns={'TP': 'Both', 'FP': 'Nanopore Only', 'FN': 'Illumina Only'})


	# Keep only the three relevant columns
	df = df[['Both', 'Illumina Only', 'Nanopore Only']]


	# Drop rows where all values are zero
	df = df[(df != 0).any(axis=1)]


	# Calculate median
	median_vals = df.median()


	return median_vals



def print_summary_statistics(original_counts):
	"""
	Print summary statistics to standard output.


	Parameters:
	- original_counts: DataFrame with raw median counts
	"""
	print("\n" + "="*80)
	print("SUMMARY STATISTICS")
	print("="*80)


	print("\n RAW MEDIAN COUNTS:")
	print("-" * 50)
	# Format the counts table
	counts_display = original_counts.round(1)
	print(counts_display.to_string())


	# Calculate proportions correctly - row-wise division
	row_sums = original_counts.sum(axis=1)
	prop_df = original_counts.div(row_sums, axis=0)


	print("\n📈 PROPORTIONS (as percentages):")
	print("-" * 50)
	# Convert to percentages and format
	prop_percentage = (prop_df * 100).round(1)
	print(prop_percentage.to_string())


	print("\n📋 SUMMARY BY VARIANT TYPE:")
	print("-" * 50)


	# Group by variant type for summary
	for variant_type in ['SNV', 'INDEL', 'SV']:
		variant_rows = [idx for idx in original_counts.index if idx.startswith(variant_type)]
		if variant_rows:
			print(f"\n{variant_type}:")
			variant_data = original_counts.loc[variant_rows]
			variant_props = prop_df.loc[variant_rows]


			for idx in variant_rows:
				condition = idx.split(' - ')[1]  # Extract condition (All, Filter, G4E)
				counts = variant_data.loc[idx]
				props = variant_props.loc[idx]
				total = counts.sum()


				print(f"  {condition}:")
				print(f"    Total variants: {total:.1f}")
				print(f"    Both: {counts['Both']:.1f} ({props['Both']*100:.1f}%)")
				print(f"    Illumina Only: {counts['Illumina Only']:.1f} ({props['Illumina Only']*100:.1f}%)")
				print(f"    Nanopore Only: {counts['Nanopore Only']:.1f} ({props['Nanopore Only']*100:.1f}%)")


	print("\n CROSS-CONDITION COMPARISON:")
	print("-" * 50)


	# Compare conditions across variant types
	conditions = ['All', 'Filter', 'G4E']
	for condition in conditions:
		condition_rows = [idx for idx in original_counts.index if idx.endswith(f' - {condition}')]
		if condition_rows:
			print(f"\n{condition} Condition:")
			condition_data = original_counts.loc[condition_rows]
			condition_props = prop_df.loc[condition_rows]


			for idx in condition_rows:
				variant = idx.split(' - ')[0]  # Extract variant type
				counts = condition_data.loc[idx]
				props = condition_props.loc[idx]
				total = counts.sum()


				print(f"  {variant}: Total={total:.1f}, Both={props['Both']*100:.1f}%, Illumina Only={props['Illumina Only']*100:.1f}%, Nanopore Only={props['Nanopore Only']*100:.1f}%")



def format_number_with_commas(num):
	"""Format numbers with commas for values >= 1000"""
	if num >= 1000:
		return f"{int(num):,}"
	else:
		return f"{int(num)}"



def plot_multiple_results(file_configs, output_file, title):
	"""
	Generates a stacked horizontal bar chart from multiple aggregated counts files,
	showing only median values grouped by variant type.
	Parameters:
	- file_configs: List of tuples (file_path, file_label, variant_types)
	- output_file: Path to save the generated plot
	- title: Title to display on the plot
	"""


	# Group by variant types
	variant_groups = {'SNV': [], 'INDEL': [], 'SV': []}


	for file_path, file_label, is_sv in file_configs:
		if is_sv:
			# For SV files, use 'RECORDS' as variant type
			median_vals = calculate_median_for_file(file_path, 'RECORDS')
			if median_vals is not None:
				variant_groups['SV'].append((file_label, median_vals))
		else:
			# For SNV/INDEL files, process both variant types
			for variant_type in ['SNVS', 'INDELS']:
				median_vals = calculate_median_for_file(file_path, variant_type)
				if median_vals is not None:
					group_key = 'SNV' if variant_type == 'SNVS' else 'INDEL'
					variant_groups[group_key].append((file_label, median_vals))


	# Build the final dataframe for plotting WITHOUT initial spacer row
	all_data = []
	all_labels = []
	group_label_positions = {}  # Track positions for group labels
	current_pos = 0


	variant_names = {
		'SNV': 'Single Nucleotide Variants',
		'INDEL': 'Insertions/Deletions', 
		'SV': 'Structural Variants'
	}


	for variant_type in ['SNV', 'INDEL', 'SV']:
		if variant_groups[variant_type]:  # Only process if group has data
			# Add group title in spacer row and track its position
			all_data.append([0, 0, 0])
			all_labels.append('')  # Empty label for spacer row
			group_label_positions[variant_names[variant_type]] = current_pos
			current_pos += 1


			# Add data rows for this group
			for file_label, median_vals in variant_groups[variant_type]:
				# Map label names
				display_label = file_label
				if file_label == 'Filter':
					display_label = 'Filtered'
				elif file_label == 'G4E':
					display_label = 'G4E'


				all_data.append(median_vals)
				all_labels.append(display_label)
				current_pos += 1


	if not all_data:
		print("No valid data found in any files")
		sys.exit(1)


	# Create DataFrame with explicit column ordering
	column_order = ['Both', 'Illumina Only', 'Nanopore Only']


	# Convert all_data to proper format for DataFrame creation
	data_rows = []
	for row in all_data:
		if isinstance(row, list):
			data_rows.append(row)
		else:
			# Convert pandas Series to list in correct order
			data_rows.append([row[col] for col in column_order])


	df = pd.DataFrame(data_rows, index=all_labels, columns=column_order)


	# Store original counts for labeling and summary (excluding spacer rows)
	original_counts = df[df.sum(axis=1) > 0].copy()  # Remove spacer rows for summary


	# Calculate proportions correctly - each row should sum to 1
	prop_df = df.div(df.sum(axis=1), axis=0).fillna(0)  # Fill NaN with 0 for spacer rows


	# Print summary statistics to stdout
	print_summary_statistics(original_counts)


	# Plotting - increased figure size for larger fonts and better spacing
	fig, ax = plt.subplots(figsize=(22, 18))  # Increased from (19.8, 14)
	

	# Create the plot with reduced bar width for more vertical spacing
	bars = prop_df.plot(
		kind='barh',
		stacked=True,
		ax=ax,
		color=['green', 'orange', 'blue'],
		width=0.4  # Reduced from 0.6 to create more vertical spacing
	)

	# Remove main plot title
	plt.xlabel("Proportion of Variants", fontsize=24, fontweight='bold')  # Increased from 18
	plt.ylabel("")  # Remove Y axis main label


	# Increase X and Y axis tick font sizes
	ax.tick_params(axis='x', labelsize=20)  # Increased from 14
	ax.tick_params(axis='y', labelsize=20)  # Increased from 14


	# Determine which rows have data vs spacer rows
	rows_with_data = df.sum(axis=1) > 0


	# Prepare yticklabels with proper renaming for sub-labels
	yticklabels = []
	for idx, has_data in zip(df.index, rows_with_data):
		if has_data:
			if idx == 'All':
				yticklabels.append('A')
			elif idx == 'Filtered':
				yticklabels.append('B')
			elif idx == 'G4E':
				yticklabels.append('C')
			else:
				yticklabels.append(idx)
		else:
			yticklabels.append('')  # Empty for spacer rows


	ax.set_yticklabels(yticklabels)


	# Set font weight bold for all visible labels
	for label in ax.get_yticklabels():
		if label.get_text() != '':
			label.set_fontweight('bold')

	# Add group labels as text annotations with larger font
	for group_name, position in group_label_positions.items():
		# Add less offset for positioning to move labels up closer to bars
		vertical_offset = 0.25 if position == 0 else 0.25  # Less offset for first group
		ax.text(
			0.5, position + vertical_offset, group_name,
			ha='center', va='center', 
			fontsize=26, fontweight='bold',
			color='black'
		)

	# Hide tick marks only for spacer rows (empty labels)
	for i, label in enumerate(ax.get_yticklabels()):
		if label.get_text() == '':
			ax.tick_params(axis='y', which='major', length=0)


	# Invert y-axis so SNV appears at top
	ax.invert_yaxis()


	# Maximize bar width by setting xlim to exactly 0-1 with no extra margins
	ax.set_xlim(0, 1)


	# Remove extra padding around the plot
	ax.margins(x=0)


	# Move legend to bottom with larger font
	plt.legend(title="Sequencing Platform", 
			loc='upper center', 
			bbox_to_anchor=(0.5, -0.06),  # Adjusted for larger figure
			ncol=3, 
			fontsize=24,  # Increased from 16
			title_fontsize=24)  # Increased from 16


	# Add "n = XXXXXX" text below each bar with larger font
	for i, (idx, row) in enumerate(df.iterrows()):
		if rows_with_data.iloc[i]:  # Only for rows with actual data
			total_count = row.sum()
			if total_count > 0:  # Only show text for non-zero totals
				# Position text below the bar
				ax.text(0.5, i + 0.35, f"n = {format_number_with_commas(total_count)}", 
					   ha='center', va='center', 
					   fontsize=22, color='black')  # Increased from 13


	plt.tight_layout(rect=[0, 0.04, 1, 0.98])  # [left, bottom, right, top] in figure coordinates
	
	# Get current ylim and adjust to reduce top spacing
	current_ylim = ax.get_ylim()
	ax.set_ylim(current_ylim[0], current_ylim[1] * 0.95)  # Reduce top by 5%


	# Save the plot
	try:
		plt.savefig(output_file, dpi=600, bbox_inches='tight')
		print(f"\n Plot successfully saved to {output_file}")
	except Exception as e:
		print(f" Error saving plot to {output_file}: {e}")
		sys.exit(1)


	plt.close()



def main():
	parser = argparse.ArgumentParser(
		description="Plot median counts from multiple aggregated_counts.csv files"
	)
	parser.add_argument("--snv-indel-all", type=str, required=True, 
					   help="Path to SNV/INDEL all aggregated_counts.csv")
	parser.add_argument("--snv-indel-filter", type=str, required=True, 
					   help="Path to SNV/INDEL filter aggregated_counts.csv")
	parser.add_argument("--snv-indel-g4e", type=str, required=True, 
					   help="Path to SNV/INDEL G4E aggregated_counts.csv")
	parser.add_argument("--sv-all", type=str, required=True, 
					   help="Path to SV all aggregated_counts.csv")
	parser.add_argument("--sv-filter", type=str, required=True, 
					   help="Path to SV filter aggregated_counts.csv")
	parser.add_argument("--sv-g4e", type=str, required=True, 
					   help="Path to SV G4E aggregated_counts.csv")
	parser.add_argument("--output", type=str, required=True, 
					   help="Path to output plot image (e.g., results_plot.png)")
	parser.add_argument("--title", type=str, default="Multi-file Comparison", 
					   help="Title for the plot")


	args = parser.parse_args()


	# Configure files: (file_path, label, is_sv)
	file_configs = [
		(args.snv_indel_all, "All", False),
		(args.snv_indel_filter, "Filter", False), 
		(args.snv_indel_g4e, "G4E", False),
		(args.sv_all, "All", True),
		(args.sv_filter, "Filter", True),
		(args.sv_g4e, "G4E", True)
	]


	plot_multiple_results(file_configs, args.output, args.title)



if __name__ == "__main__":
	main()
