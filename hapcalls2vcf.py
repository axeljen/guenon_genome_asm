import pysam
import gzip
import sys

hapcalls_file = sys.argv[1]

# using a map to connect the correct sample with the correct haplotype in the hapcalls file from angsd
# three cols: sequence name in hapcalls, sample name in output vcf, haplotype (1 or 2, or for haploid samples, haploid)
indmap_file = "indmap.txt"

# also using a file with sex assignments (two cols: sample name in vcf file, sex assignment (male/female) )
sex_assignments = "sex_assignments.txt"

# open output vcf
out_vcf_path = sys.argv[2]

# parse indmap file
def parse_indmap(indmap_file):
	indmap = {}
	with open(indmap_file) as f:
		for line in f:
			ind,sample,haplotype = line.strip().split()
			if haplotype == "haploid":
				indmap[sample] = [ind,ind]
			elif haplotype == "1":
				if not sample in indmap:
					indmap[sample] = [ind]
				else:
					hap2 = indmap[sample][0]
					indmap[sample][0] = ind
					indmap[sample].append(hap2)
			elif haplotype == "2":
				if not sample in indmap:
					indmap[sample] = [ind]
				else:
					indmap[sample].append(ind)
	return indmap

def parse_sex_assignments(file):
	sex_assignments = {}
	with open(file) as f:
		for line in f:
			sample,sex = line.strip().split()
			sex_assignments[sample] = sex
	return sex_assignments

# fetch indmap
indmap = parse_indmap(indmap_file)

# parse sex assignments (for differential filtering of X-chromosome genotypes)
sex_assignment_dict = parse_sex_assignments(sex_assignments)

# open hapcalls file
hapcalls = gzip.open(hapcalls_file, 'rt')

header = pysam.VariantHeader()

header.add_line('##fileformat=VCFv4.2')
header.add_line("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
header.add_line("##FILTER=<ID=PASS,Description=\"All filters passed\">")


# parse header line
line = hapcalls.readline()
chrom = line.split()[0]
if not chrom in header.contigs:
	header.contigs.add(chrom)
pos = line.split()[1]
major = line.split()[2]
samples = line.split()[3:]
for sample in indmap.keys():
	if not sample in header.samples:
		header.samples.add(sample)

outvcf = pysam.VariantFile(out_vcf_path, 'w', header=header)

# index-sample dictionary
indmap_dict = {}
for i,sample in enumerate(samples):
	indmap_dict[i] = sample

# function to parse hapcall line
def parseCall(indmap_dict,line):
	line = line.strip().split()
	chrom = line[0]
	pos = line[1]
	genotypes = line[3:]
	genotype_dict = {'chrom': chrom, 'pos': pos, 'genotypes': {ind: genotypes[i] for i,ind in indmap_dict.items()}}
	return genotype_dict

def callToVCFRecord(genotype_dict, vcf, indmap, sex_assignment_dict):
	if not genotype_dict['chrom'] in vcf.header.contigs:
		vcf.header.contigs.add(genotype_dict['chrom'])
	rec = outvcf.header.new_record(
		contig=genotype_dict['chrom'],
		start=int(genotype_dict['pos']) - 1,
		id='.',
		filter='PASS',
	)
	#print(rec)
	# get ref allele
	rec.ref = genotype_dict['genotypes'][indmap['reference'][0]]
	# get alt alleles
	alts = {}
	# and add binary genotypes
	genotype_dict['binary_genotypes'] = {}
	for sample in genotype_dict['genotypes']:
		if not genotype_dict['genotypes'][sample] == rec.ref and not genotype_dict['genotypes'][sample] == 'N':
			if not genotype_dict['genotypes'][sample] in alts:
				alts[genotype_dict['genotypes'][sample]] = len(alts) + 1
			genotype_dict['binary_genotypes'][sample] = alts[genotype_dict['genotypes'][sample]]
		elif genotype_dict['genotypes'][sample] == 'N':
			genotype_dict['binary_genotypes'][sample] = None
		elif genotype_dict['genotypes'][sample] == rec.ref:
			genotype_dict['binary_genotypes'][sample] = 0
	if len(alts) > 0:
		rec.alts = alts
	else:
		rec.alts = '.'
	# add genotypes to record
	for sample,inds in indmap.items():
		if sample == "reference":
			continue
		if not None in (genotype_dict['binary_genotypes'][inds[0]], genotype_dict['binary_genotypes'][inds[1]]):
			rec.samples[sample]['GT'] = (genotype_dict['binary_genotypes'][inds[0]], genotype_dict['binary_genotypes'][inds[1]])
		elif rec.contig == "NC_041774.1":
			# if on the x chromosome and the sample is a male, a missing genotype on one haplotype is expected
			if sex_assignment_dict[sample] == 'male':
				# check if either of the genotypes are not none
				if not genotype_dict['binary_genotypes'][inds[0]] == None:
					rec.samples[sample]['GT'] = (genotype_dict['binary_genotypes'][inds[0]], genotype_dict['binary_genotypes'][inds[0]])
				elif not genotype_dict['binary_genotypes'][inds[1]] == None:
					#print("Setting genotype for " + sample + " to " + str(genotype_dict['binary_genotypes'][inds[0]]))
					rec.samples[sample]['GT'] = (genotype_dict['binary_genotypes'][inds[1]], genotype_dict['binary_genotypes'][inds[1]])
				else:
					rec.samples[sample]['GT'] = (None, None)
			else:
				rec.samples[sample]['GT'] = (None, None)
		else:
			# otherwise set all genotypes with one missing to missing
			rec.samples[sample]['GT'] = (None, None)
		if not inds[0] == inds[1]:
			# phase genotypes from those with phased assemblies
			rec.samples[sample].phased = True
	return rec

line = hapcalls.readline()

genotype_dict = parseCall(indmap_dict,line)


for line in hapcalls:
	genotype_dict = parseCall(indmap_dict,line)
	vcfrec = callToVCFRecord(genotype_dict, outvcf, indmap, sex_assignment_dict)
	outvcf.write(vcfrec)

outvcf.close()