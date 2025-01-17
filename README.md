# Custom scripts used in the guenon genome assembly project

This repo contains custom analytical scripts that were used in the guenon genome assembly project.

The script progressive_syntenymap.py and plotSyntenymap.R was used to generate the syntenymap of all guenon genomes (Figure 1A). A usage example is provided in make_progressive_synplot.sh.

Preliminary fissions and fusions were identified and plotted with find_fissions_and_fusions.py and plot_fissions_fusions.R. A usage example is given in find_fissions_and_fusions.sh. Note that these might need some manual curation afterwards.

To call SNPs, I first used angsd to generate haploid genotype calls against the rhesus macaque reference genome, and these were then combined per sample into diploid genotypes in vcf format using the hapcalls2vcf.py script.