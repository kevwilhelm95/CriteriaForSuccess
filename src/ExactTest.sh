#!/usr/bin/bash

# all file paths relative to taco cluster
set -e
snpsift=/storage/lichtarge/shared/snpEff/scripts/snpSift

# Set input variables
in_vcf=$1
gene_file=$2
sample_only_file=$3
sample_case_file=$4
sample_fam_file=$5
out_path=$6
cd $out_path

# Subset for defined samples, re-calculate fields, and perform exactTest #
echo "Performing ExactTest"
bcftools view -Ou -c 1 - -S ${sample_only_file} -R ${gene_file} ${in_vcf} | \
	bcftools +fill-tags -Ou -- -S ${sample_case_file} -t AC,AN,AF,AC_Het,AC_Hom | \
		/storage/lichtarge/shared/snpEff/scripts/snpSift CaseControl -tfam ${sample_fam_file} | \
				bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/SYMBOL\t%INFO/ENSP\t%INFO/Consequence\t%INFO/HGVSp\t%INFO/EA\t%INFO/AN_0\t%INFO/AN_1\t%INFO/Cases\t%INFO/Controls\t%INFO/AC_1\t%INFO/AC_Het_1\t%INFO/AC_Hom_1\t%INFO/AC_0\t%INFO/AC_Het_0\t%INFO/AC_Hom_0\t%INFO/CC_ALL\t%INFO/CC_DOM\t%INFO/CC_REC\t%INFO/AF\t%INFO/AC\n' \
					>CaseControl.Variants.OR.txt

echo "Done with grabbing variants"



# Need functions to create tables for samples only and fam file -- save in intermediate file folder in output path
# Chose not to limit to genes we are only interested in -- Can ammend and probably should
# Is there a way that I can make this automatically executable????


### QUESTIONS FOR DANNY
# Can I filter the vcf by gene names in a file
# Can I pipe from bcftools to snpSift?