# wrapper script for generating syntenyplots

# path to faidx file for reference genome (everything should be aligned to this reference)
REF_INDEX=path/to/reference.fa.fai

## Make a big syntenymap with all genomes
name=output_name
mkdir -p ${name}

# run the progressive_syntenymap.py script
python3 progressive_syntenymap.py \
	--reference ${REF_INDEX} \
	--index genome_1_index.fai \
	genome_2_index.fai \
	genome_3_index.fai \
	genome_4_index.fai \
	--alignment genome_1_alignment.paf \
	genome_2_alignment.paf \
	genome_3_alignment.paf \
	genome_4_alignment.paf \
	--output ${name}/${name}

# the above script creates two output files, suffixed with .chroms.txt and .syntenyMap.txt. These can be plotted with the R script below.
Rscript ../minimap2synteny/plotSyntenyMaps.R ${name}/${name}.chroms.txt ${name}/${name}.syntenyMap.txt \
 	"genome_1_label,genome_2_label,genome_3_label,genome_4_label"
