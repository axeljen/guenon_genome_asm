## Example workflow for identifying preliminary fissions and fusions


# first identify fissions/fusions
python3 find_fissions_and_fusions.py \
	--paf path/to/alignment.paf \
	--reference_index path/to/reference_index.fai \
	--query_index path/to/query_index.fai \
	--out-prefix output/prefix

# then plot them, prior to manual curation
Rscript plot_fissions_fusions.R \
	output/prefix.synmap.tsv \
	output/prefix.genomes.tsv \
	output/prefix.prel_fissions_and_fusions.tsv \
	output/prefix.raw_calls