#!/bin/bash

# Define input files
vcf_file="input_vcf"
metadata_file="input_metadata"
output_dir="output_path"
output_file="$output_dir/$(basename "$vcf_file" .vcf.gz)_summary.txt"
mkdir -p "$output_dir"

# to index the vcf and make the following commands faster
#bcftools index "$vcf_file"

#echo "vcf indexed"

# Generate stats file
#bcftools stats "$vcf_file" > "$output_dir/$(basename "$vcf_file" .vcf.gz)stats.txt"

#echo "bcftools stats completed"

# to generate a summarized version of the vcf file, the format of the output can be changed according to the need.
#output_summarized_vcf="$output_dir/$(basename "$vcf_file" .vcf.gz)_summarized_vcf.txt"
#echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER" > "$output_summarized_vcf"
#bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\n' "$vcf_file" >> "$output_summarized_vcf"

#echo "vcf summarized"

# to analyze and report how many samples have genotype information for each position of the genome in the vcf file
#output_genotyped_samples="$output_dir/$(basename "$vcf_file" .vcf.gz)_genotyped_samples_for_position.txt"
#echo "Genotyped_Samples_for_Positions_Count"
#echo -e "CHROM\tPOS\tGenotyped_Samples_Count" > "$output_genotyped_samples"
#bcftools query -f '%CHROM\t%POS[\t%GT]\n' "$vcf_file" | \
#awk 'BEGIN {FS="\t"; OFS="\t"} {count = 0; for (i = 3; i <= NF; i++) if ($i != ".") count++; print $1, $2, count}' >> #"$output_genotyped_samples"

#echo "count of genotyped samples for position complete"

# Define stats file variable
stats_file="$output_dir/$(basename "$vcf_file" .vcf.gz)stats.txt"

# Initialize output file and write the header
echo -e "Metric\tCount" > "$output_file"

# Extract sample names from VCF header
samples=($(bcftools query -l "$vcf_file"))

echo "names extracted"

# Count number of samples
num_samples=${#samples[@]}
echo -e "Number of Samples\t$num_samples" >> "$output_file"

# Parse metadata file
declare -A populations super_pops genders

while IFS=$'\t' read -r sample pop super_pop gender; do
    	# Encode the population name
    	encoded_pop=$(printf "%s" "$pop" | base64)

    	if [[ " ${samples[@]} " =~ " $sample " ]]; then
        	# Decode the population name when accessing the associative array
        	populations["$encoded_pop"]=$((populations["$encoded_pop"] + 1))
        	super_pops["$super_pop"]=$((super_pops["$super_pop"] + 1))
		genders["$gender"]=$((genders["$gender"] + 1))
	fi

done < "$metadata_file"

# Count number of samples from each population
echo -e "\nNumber of Samples from Each Population" >> "$output_file"
for encoded_pop in "${!populations[@]}"; do
    	# Decode the population name
    	pop=$(printf "%s" "$encoded_pop" | base64 --decode)
    	echo -e "$pop\t${populations[$encoded_pop]}" >> "$output_file"
done

# Count number of samples from each super population
echo -e "\nNumber of Samples from Each Super Population" >> "$output_file"
for super_pop in "${!super_pops[@]}"; do
    	echo -e "$super_pop\t${super_pops[$super_pop]}" >> "$output_file"
done

# Print the counts for each gender category
echo -e "\nNumber of Samples by Gender" >> "$output_file"
for gender in "${!genders[@]}"; do
    	echo -e "$gender: ${genders[$gender]}" >> "$output_file"
done

# Function to extract and debug the specified statistic
extract_stat() {
    	local stat_name=$1
    	local pattern=$2
    	local extracted_value=$(awk -v pattern="$pattern" '$1 == "SN" && $0 ~ pattern {print $NF}' "$stats_file")
    	echo -e "Debug: Extracting $stat_name"
    	echo -e "Debug: Pattern: $pattern"
    	echo -e "Debug: $stat_name extracted:\n$extracted_value"
    	if [ -z "$extracted_value" ]; then
        	echo -e "Debug: No $stat_name extracted!"
    	else
        	echo -e "$stat_name:\t$extracted_value" >> "$output_file"
    	fi
}

# Extract and debug the number of no-ALTs
extract_stat "Number of no-ALTs" "number of no-ALTs:"

# Extract and debug the number of SNPs
extract_stat "Number of SNPs" "number of SNPs:"

# Extract and debug the number of MNPs
extract_stat "Number of MNPs" "number of MNPs:"

# Extract and debug the number of indels
extract_stat "Number of indels" "number of indels:"

# Extract and debug the number of others
extract_stat "Number of others" "number of others:"

# Extract and debug the number of multiallelic sites
extract_stat "Number of multiallelic sites" "number of multiallelic sites:"

# Extract and debug the number of multiallelic SNP sites
extract_stat "Number of multiallelic SNP sites" "number of multiallelic SNP sites:"

# Function to calculate average read depth
calculate_average_read_depth() {
	# Initialize variables
    	sum_read_depths=0
    	total_sites=0

    	# Iterate through each line in the DP section
    	while IFS=$'\t' read -r header id bin genotypes fraction_genotypes sites fraction_sites; do
        	# Skip the header line
        	if [[ "$header" != "DP" ]]; then
            		continue
        	fi
        
        	# Extract the lower bound of the bin
        	lower_bound=$(echo "$bin" | awk -F'-' '{print $1}')

        	# Check if lower_bound is a valid number
        	if [[ ! "$lower_bound" =~ ^[0-9]+$ ]]; then
            		continue
        	fi
        
        	# Accumulate the read depths and total sites
        	sum_read_depths=$((sum_read_depths + (lower_bound * sites)))
        	total_sites=$((total_sites + sites))
    	done < <(grep -E "^DP" "$stats_file")

    	# Calculate the average read depth
    	if (( total_sites != 0 )); then
        	average_read_depth=$(awk "BEGIN {printf \"%.2f\", $sum_read_depths / $total_sites}")
        	echo "Average Read Depth: $average_read_depth" >> "$output_file"
    	else
        	echo "Total number of sites is zero. Unable to calculate average read depth." >> "$output_file"
    	fi
}

# Calculate the average read depth
calculate_average_read_depth

echo "Summary complete. Output written to $output_file."