#!/bin/bash

# Define input and output directories
input_dir="input_path"
output_dir="output_path"
genotyped_samples_dir="output:path_for_final_output"

# Create output directories if they don't exist
mkdir -p "$output_dir"
mkdir -p "$genotyped_samples_dir"

# Loop through each VCF file in the input directory
for input_vcf in "$input_dir"/*.vcf.gz; do
	# Extract the filename without extension
    	filename=$(basename "$input_vcf" .vcf.gz)

    	# Define the output VCF file path
    	output_vcf="$output_dir/${filename}_newbiallelic_snps2.vcf.gz"
    	genotyped_output_vcf="${output_vcf%.vcf.gz}_genotyped.vcf.gz"

    	# Step 1: Filter for biallelic SNPs
    	bcftools view -m2 -M2 -v snps "$input_vcf" -Oz -o "$output_vcf"
    	
    	# Index the compressed VCF file
    	tabix -p vcf "$output_vcf"
    
    	# Step 2: Identify samples with only missing genotypes in a single pass
    	non_missing_samples_file="$genotyped_samples_dir/${filename}_genotyped_samples.txt"
    	bcftools query -f '[%GT\n]' "$output_vcf" | awk '
    	BEGIN {
        	OFS="\t";
    	}
    	{
        	for (i = 1; i <= NF; i++) {
            		if ($i != "./.") {
                		non_missing[i] = 1;
            		}
        	}
    	}
    	END {
        	for (i in non_missing) {
            		print i;
        	}
    	}
    	' | bcftools query -l "$output_vcf" | awk '{print $1}' > "$non_missing_samples_file"

    	# Step 3: Remove samples with only missing genotypes
    	if [ -s "$non_missing_samples_file" ]; then
        	bcftools view -S "$non_missing_samples_file" -Oz -o "$genotyped_output_vcf" "$output_vcf"
        	tabix -p vcf "$genotyped_output_vcf"
    	else
        	cp "$output_vcf" "$genotyped_output_vcf"
        	tabix -p vcf "$genotyped_output_vcf"
    	fi

    	echo "Biallelic SNPs have been filtered, compressed, indexed, and cleaned for $filename successfully."
done

echo "All VCF files have been processed."
