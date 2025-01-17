import sys
import numpy as np
import argparse
#import paf2synmap as ps
import importlib

""" This script attempts to identify fissions and fusions in a query genome based on an alignment file in paf format. It creates a preliminary list of fissions and fusions, and plots them on a syntenymap.
These then need to be manually curated as it sometimes ids false postives (specifically single events are sometimes identified as multiple events). Also contains the classes Genome(), Alignment() and PAF(), which were written for general parsing of paf alignments."""

class Genome():
	# class for storing simple genome info as chrom-length pairs
    def __init__(self):
        self.scaffolds = {}
        self.cumulative_startpos = {}
    # add chromosome to genome
    def add_chromosome(self, scaffold, length):
        self.scaffolds[scaffold] = length
    def add_cumulative_startpos(self,scaffold_gap):
        cumulative_pos = 0
        for scaffold, length in self.scaffolds.items():
            self.cumulative_startpos[scaffold] = cumulative_pos
            cumulative_pos += length + scaffold_gap
    def parse_genome(self, genome_file):
        # parse genome file and return genome object
        with open(genome_file) as f:
            for line in f.readlines():
                scaffold = line.strip().split()[0]
                length = line.strip().split()[1]
                self.add_chromosome(scaffold, int(length))
    def __str__(self):
        return "\n".join(["{}: {}".format(scaffold, length) for scaffold, length in self.scaffolds.items()])

class Alignment():
    def __init__(self):
        self.query_scaffold = None
        self.query_sequence_length = None
        self.query_start = None
        self.query_end = None
        self.query_aln_length = None
        self.strand = None
        self.target_scaffold = None
        self.target_length = None
        self.target_start = None
        self.target_end = None
        self.target_aln_length = None
        self.n_base_matches = None
        self.n_base_total = None
        self.mapping_quality = None
    def __repr__(self):
        return "{}:{}-{} -> {}:{}-{} ({} bp, strand: {})".format(self.query_scaffold, self.query_start, self.query_end, self.target_scaffold, self.target_start, self.target_end, self.query_aln_length, self.strand)
    def parse_paf_alignment(self, pafline):
        pafitems = pafline.strip().split()
        self.query_scaffold = pafitems[0]
        self.query_sequence_length = int(pafitems[1])
        self.query_start = int(pafitems[2])
        self.query_end = int(pafitems[3])
        self.query_aln_length = self.query_end - self.query_start
        self.strand = pafitems[4]
        self.target_scaffold = pafitems[5]
        self.target_sequence_length = int(pafitems[6])
        self.target_start = int(pafitems[7])
        self.target_end = int(pafitems[8])
        self.target_aln_length = self.target_end - self.target_start
        self.n_base_matches = int(pafitems[9])
        self.n_base_total = int(pafitems[10])
        self.mapping_quality = int(pafitems[11])
    def is_consecutive(self, previous_alignment, max_gap=10000):
        # 1. Check if they are on the same strand
        if previous_alignment.strand != self.strand:
            return False
        # 2. Query coordinates check based on the strand
        if previous_alignment.strand == '+':
            # Query start must be larger in the second alignment
            if self.query_start <= previous_alignment.query_start:
                return False
            # Check if query coordinates are within max distance
            if self.query_start - previous_alignment.query_end > max_gap:
                return False
        else:  # For minus strand
            # Query start must be smaller in the second alignment
            if self.query_start >= previous_alignment.query_start:
                return False
            # Check if query coordinates are within max distance
            if previous_alignment.query_start - self.query_end > max_gap:
                return False
        # 3. Target coordinates must be larger in the second alignment and within max_gap
        if self.target_start <= previous_alignment.target_start:
            return False
        if self.target_start - previous_alignment.target_end > max_gap:
            return False
        # 4. query and target scaffolds must be the same
        if not self.query_scaffold == previous_alignment.query_scaffold or not self.target_scaffold == previous_alignment.target_scaffold:
            return False
        return True
    def merge(self, aln2):
        # merge two alignments and return a new alignment
        merged_aln = Alignment()
        merged_aln.query_scaffold = self.query_scaffold
        merged_aln.query_sequence_length = self.query_sequence_length
        merged_aln.query_start = min(self.query_start, aln2.query_start)
        merged_aln.query_end = max(self.query_end, aln2.query_end)
        merged_aln.query_aln_length = merged_aln.query_end - merged_aln.query_start
        merged_aln.strand = self.strand
        merged_aln.target_scaffold = self.target_scaffold
        merged_aln.target_sequence_length = self.target_sequence_length
        merged_aln.target_start = min(self.target_start, aln2.target_start)
        merged_aln.target_end = max(self.target_end, aln2.target_end)
        merged_aln.target_aln_length = merged_aln.target_end - merged_aln.target_start
        merged_aln.n_base_matches = self.n_base_matches + aln2.n_base_matches
        merged_aln.n_base_total = self.n_base_total + aln2.n_base_total
        merged_aln.mapping_quality = (self.mapping_quality + aln2.mapping_quality) / 2
        return merged_aln
    def flip_query_alignment(self):
        self.query_start, self.query_end = self.query_sequence_length - self.query_end, self.query_sequence_length - self.query_start
        if self.strand == "-":
            self.strand = "+"
        else:
            self.strand = "-"
	
class PAF():
    def __init__(self):
        self.alignments = []
        self.target_genome = None
        self.query_genome = None
        self.sorted = None
    def __repr__(self):
        # print the first 10 alignments
        string_repr = "List of {} alignments:\n".format(len(self.alignments)) + "\n".join([str(a) for a in self.alignments[:10]]) + "\n...{}...".format(len(self.alignments)-10)
        return string_repr
    def parse_paf(self, paf_file, target_genome = None, query_genome = None, overwrite=False):
        # if paf is not empty, require overwrite
        if len(self.alignments) > 0 and not overwrite:
            raise ValueError("Alignments already present. Set overwrite=True to overwrite.")
        self.alignments = []
        with open(paf_file) as f:
            for line in f.readlines():
                alignment = Alignment()
                alignment.parse_paf_alignment(line)
                self.alignments.append(alignment)
        if target_genome:
            self.target_genome = Genome()
            self.target_genome.parse_genome(target_genome)
        if query_genome:
            self.query_genome = Genome()
            self.query_genome.parse_genome(query_genome)
    def filter_alignments(self, minlen=10000, minmapq=40, minpid=0.9, only_ref_chroms=False, only_query_chroms=False):
        # filter alignments based on user-defined criteria
        if only_ref_chroms:
            self.alignments = [a for a in self.alignments if a.target_scaffold in self.target_genome.scaffolds]
        if only_query_chroms:
            self.alignments = [a for a in self.alignments if a.query_scaffold in self.query_genome.scaffolds]
        # then filter based on the remaining criteria
        self.alignments = [a for a in self.alignments if a.query_aln_length >= minlen and a.mapping_quality >= minmapq and a.n_base_matches/a.n_base_total >= minpid]
    def sort_on_target(self):
        self.alignments = sorted(self.alignments, key=lambda x: (x.target_scaffold, x.target_start))
        self.sorted = "target"
    def sort_on_query(self):
        self.alignments = sorted(self.alignments, key=lambda x: (x.query_scaffold, x.query_start))
        self.sorted = "query"
    def merge_consecutive(self, max_gap=10000):
        if not self.sorted:
            raise ValueError("Alignments must be sorted before merging.")
        merged_alignments = [self.alignments[0]]
        for aln in self.alignments[1:]:
            if aln.is_consecutive(merged_alignments[-1], max_gap=max_gap):
                merged_alignments[-1] = aln.merge(merged_alignments[-1])
            else:
                merged_alignments.append(aln)
        self.alignments = merged_alignments
    def get_dominating_strand(self, query_chrom):
        # get the dominating strand for a query chromosome
        query_alignments = [a for a in self.alignments if a.query_scaffold == query_chrom]
        plus_strand = sum([a.query_aln_length for a in query_alignments if a.strand == "+"])
        minus_strand = sum([a.query_aln_length for a in query_alignments if a.strand == "-"])
        if plus_strand > minus_strand:
            return "+"
        else:
            return "-"
    def flip_query_alignments(self, query_chrom):
        # print("Flipping alignments for chromosome {}.".format(query_chrom))
        # flip the query alignments for a query chromosome
        flipped_alignments = []
        for aln in self.alignments:
            if aln.query_scaffold == query_chrom:
                aln.flip_query_alignment()
            flipped_alignments.append(aln)
        self.alignments = flipped_alignments
    def get_min_target_position(self, query_chrom):
        previous_position = 0
        previous_chrom = None
        for aln in self.alignments:
            if not previous_chrom == aln.target_scaffold:
                previous_position = 0
                previous_chrom = aln.target_scaffold
            if aln.query_scaffold == query_chrom:
                return (aln.target_scaffold, previous_position, aln.target_start)
            previous_position = aln.target_end
        return (None,None, None)
    def get_max_target_position(self, query_chrom):
        maxpos = 0
        for i,aln in enumerate(self.alignments):
            if aln.query_scaffold == query_chrom:
                maxpos = aln.target_end
                chrom = aln.target_scaffold
                try:
                    next_aln = self.alignments[i+1]
                    if next_aln.query_scaffold == query_chrom:
                        maxpos_next = next_aln.target_end
                    else:
                        maxpos_next = self.target_genome.scaffolds[aln.target_scaffold] 
                except IndexError:
                    maxpos_next = self.target_genome.scaffolds[aln.target_scaffold]
        return (chrom, maxpos, maxpos_next)
    def total_aln_length(self):
        return sum([a.query_aln_length for a in self.alignments])
    def get_query_position_from_target(self, target_chrom, target_pos):
        # see if any alignment overlaps the target position
        for aln in self.alignments:
            if aln.target_scaffold == target_chrom and aln.target_start <= target_pos <= aln.target_end:
                # check how far from the start of the alignment our target position is
                relative_pos = (target_pos - aln.target_start) / aln.target_aln_length
                # calculate the query position
                if aln.strand == "+":
                    return (aln.query_scaffold, round(aln.query_start + relative_pos * aln.query_aln_length, 0))
                else:
                    return (aln.query_scaffold, round(aln.query_end - relative_pos * aln.query_aln_length, 0))
        return None
    def switch_strand_sign(self,chrom):
        for aln in self.alignments:
            if aln.query_scaffold == chrom:
                if aln.strand == "+":
                    aln.strand = "-"
                else:
                    aln.strand = "+"
    def get_alignments(self,query_chrom = None, target_chrom = None):
        alignments = []
        if query_chrom:
            for a in self.alignments:
                if a.query_scaffold == query_chrom:
                    alignments.append(a)
        if target_chrom:
            for a in self.alignments:
                if a.target_scaffold == target_chrom:
                    alignments.append(a)
        return alignments
    def fetch(self, target_chrom = None, target_start = None, target_end = None):
        if target_chrom == None and not target_start == None:
            raise ValueError("If target_start is provided, target_chrom must also be provided.")
        if target_start == None and not target_end == None:
            raise ValueError("Target end can only be provided if target_start is also provided.")
        if target_start != None and target_end == None:
            raise ValueError("If target_start is provided, target_end must also be provided.")
        # fetch alignments that overlap with the target region
        overlapping_alignments = []
        for aln in self.alignments:
            if target_chrom == None:
                overlapping_alignments.append(aln)
            else:
                if target_start == None:
                    if aln.target_scaffold == target_chrom:
                        overlapping_alignments.append(aln)
                else:
                    if aln.target_scaffold == target_chrom and aln.target_start <= target_end and aln.target_end >= target_start:
                        overlapping_alignments.append(aln)
        return overlapping_alignments
    def fetch_query(self, query_chrom = None, query_start = None, query_end = None):
        if query_chrom == None and not query_start == None:
            raise ValueError("If query_start is provided, query_chrom must also be provided.")
        if query_start == None and not query_end == None:
            raise ValueError("Query end can only be provided if query_start is also provided.")
        if query_start != None and query_end == None:
            raise ValueError("If query_start is provided, query_end must also be provided.")
        # fetch alignments that overlap with the query region
        overlapping_alignments = []
        for aln in self.alignments:
            if query_chrom == None:
                overlapping_alignments.append(aln)
            else:
                if query_start == None:
                    if aln.query_scaffold == query_chrom:
                        overlapping_alignments.append(aln)
                else:
                    if aln.query_scaffold == query_chrom and aln.query_start <= query_end and aln.query_end >= query_start:
                        overlapping_alignments.append(aln)
        return overlapping_alignments
    def get_median_target_position(self, query_chrom):
        target_positions = []
        for aln in self.alignments:
            if aln.query_scaffold == query_chrom:
                try:
                    target_positions.append(aln.target_start + self.target_genome.cumulative_startpos[aln.target_scaffold])
                except KeyError:
                    print("No cumulative start position found for scaffold {}.".format(aln.target_scaffold))
                    sys.exit(1)
        return np.median(target_positions)
    def add_cumulative_query_start(self, scaffold_gap=10000000):
        # sort on target
        self.sort_on_target()
        # make sure cumulative target positions are set
        if self.target_genome.cumulative_startpos is None:
            self.target_genome.add_cumulative_startpos(scaffold_gap)
        # get median startpos of query scaffolds
        query_positions = []
        for scaffold in self.query_genome.scaffolds:
            query_positions.append({'scaffold':scaffold, 'median_start': self.get_median_target_position(scaffold)})
        # sort on median start position
        query_positions = sorted(query_positions, key=lambda x: x['median_start'])
        # add cumulative query start positions
        cumulative_pos = 0
        for q in query_positions:
            self.query_genome.cumulative_startpos[q['scaffold']] = cumulative_pos
            cumulative_pos += self.query_genome.scaffolds[q['scaffold']] + scaffold_gap
    def paf2synmap(self, scaffold_gap=10000000):
        # add cumulative startpositions od chromosomes
        self.add_cumulative_query_start(scaffold_gap=scaffold_gap)
        # make a big string
        synmap = ""
        synmap += "{tchrom}\t{tstart_original}\t{tend_original}\t{qchrom}\t{qstart_original}\t{qend_original}\t{tstart_cumulative}\t{tend_cumulative}\t{qstart_cumulative}\t{qend_cumulative}\t{strand}\n".format(tchrom="tchrom", tstart_original="tstart_original", tend_original="tend_original", qchrom="qchrom", qstart_original="qstart_original", qend_original="qend_original", tstart_cumulative="tstart_cumulative", tend_cumulative="tend_cumulative", qstart_cumulative="qstart_cumulative", qend_cumulative="qend_cumulative", strand="strand")
        # print alignments
        for aln in self.alignments:
            synmap += "{tchrom}\t{tstart_original}\t{tend_original}\t{qchrom}\t{qstart_original}\t{qend_original}\t{tstart_cumulative}\t{tend_cumulative}\t{qstart_cumulative}\t{qend_cumulative}\t{strand}\n".format(tchrom=aln.target_scaffold, tstart_original=aln.target_start, tend_original=aln.target_end, qchrom=aln.query_scaffold, qstart_original=aln.query_start, qend_original=aln.query_end, tstart_cumulative=aln.target_start + self.target_genome.cumulative_startpos[aln.target_scaffold], tend_cumulative=aln.target_end + self.target_genome.cumulative_startpos[aln.target_scaffold], qstart_cumulative=aln.query_start + self.query_genome.cumulative_startpos[aln.query_scaffold], qend_cumulative=aln.query_end + self.query_genome.cumulative_startpos[aln.query_scaffold], strand=aln.strand)
        return synmap


# Parse arguments
parser = argparse.ArgumentParser(description='Call preliminary fissions and fusions from a PAF file')

parser.add_argument('--paf', type=str, help='Path to the PAF file', required=True)
parser.add_argument('--reference_index', type=str, help='Path to the reference genome index (.fai or similar, first two columns should be chrom\tlength)', required=True)
parser.add_argument('--query_index', type=str, help='Path to the query genome index (.fai or similar, first two columns should be chrom\tlength)', required=True)
parser.add_argument('--minlen', type=int, help='Minimum alignment length', default=50000)
parser.add_argument('--minq', type=int, help='Minimum mapping quality', default=50)
parser.add_argument('--minpid', type=float, help='Minimum percent identity', default=.85)
parser.add_argument('--only_ref_chroms', type=bool, help='Only consider alignments to reference chromosomes', default=True)
parser.add_argument('--only_query_chroms', type=bool, help='Only consider alignments to query chromosomes', default=True)
parser.add_argument('--out-prefix', type=str, help='Prefix for output files', default="fission_fusion_calls")
parser.add_argument('--blocksize', type=int, help='Size of flanking regions to check for synteny', default=3000000)

args = parser.parse_args()



aln = args.paf
refindex = args.reference_index
queryindex = args.query_index

calls_out = args.out_prefix + ".prel_fissions_and_fusions.tsv"
synmap_out = args.out_prefix + ".synmap.tsv"
genomes_out = args.out_prefix + ".genomes.tsv"

minlen = args.minlen
minq = args.minq
minpid = args.minpid
only_ref_chroms = args.only_ref_chroms
only_query_chroms = args.only_query_chroms

def call_fissions(paf):
	fissions = []
	processed_alns = [paf.alignments[0]]
	# sort on target
	paf.sort_on_target()
	for a in paf.alignments[1:]:
		if a.query_scaffold != processed_alns[-1].query_scaffold:
			if a.target_scaffold == processed_alns[-1].target_scaffold:
				# this would be a fission
				fissions.append({
					'type': "fission",
					'target_chrom': a.target_scaffold,
					'min_pos': processed_alns[-1].target_end,
					'max_pos': a.target_start,
					'query_pos_1': (processed_alns[-1].query_scaffold,processed_alns[-1].query_end,processed_alns[-1].strand),
					'query_pos_2': (a.query_scaffold,a.query_end,a.strand),
				})
		processed_alns.append(a)
	return fissions

def call_fusions(paf):
	fusions = []
	processed_alns = [paf.alignments[0]]
	# sort on target
	paf.sort_on_query()
	for a in paf.alignments[1:]:
		if a.target_scaffold != processed_alns[-1].target_scaffold:
			if a.query_scaffold == processed_alns[-1].query_scaffold:
				# this would be a fusion
				fusions.append({
					'type': "fusion",
					'query_chrom': a.query_scaffold,
					'min_pos': processed_alns[-1].query_end,
					'max_pos': a.query_start,
					'target_pos_1': (processed_alns[-1].target_scaffold,processed_alns[-1].target_end,processed_alns[-1].strand),
					'target_pos_2': (a.target_scaffold,a.target_end,a.strand),
				})
		processed_alns.append(a)
	return fusions

def calculate_synteny(alignments):
	target_length = 0
	query_synteny = {}
	for a in alignments:
		if a.query_scaffold not in query_synteny:
			query_synteny[a.query_scaffold] = {'total_length': a.query_end - a.query_start}
		else:
			query_synteny[a.query_scaffold]['total_length'] += a.query_end - a.query_start
		target_length += a.query_end - a.query_start
	for q in query_synteny:
		query_synteny[q]['fraction'] = query_synteny[q]['total_length'] / target_length
	return query_synteny

def calculate_synteny_target(alignments):
	query_length = 0
	target_synteny = {}
	for a in alignments:
		if a.target_scaffold not in target_synteny:
			target_synteny[a.target_scaffold] = {'total_length': a.target_end - a.target_start}
		else:
			target_synteny[a.target_scaffold]['total_length'] += a.target_end - a.target_start
		query_length += a.target_end - a.target_start
	for t in target_synteny:
		target_synteny[t]['fraction'] = target_synteny[t]['total_length'] / query_length
	return target_synteny

def is_syntenic(paf,target_chrom,target_start,target_end,query_chrom,min_synteny = .8, min_totlen = 500000):
	if target_start < 0:
		return False
	if target_end > paf.target_genome.scaffolds[target_chrom]:
		return False
	alignments = paf.fetch(target_chrom,target_start,target_end)
	query_synteny = calculate_synteny(alignments)
	#print(query_synteny)
	if query_chrom in query_synteny:
		if query_synteny[query_chrom]['fraction'] >= min_synteny:
			if query_synteny[query_chrom]['total_length'] >= min_totlen:
				return True
	return False

def is_syntenic_query(paf,query_chrom,query_start,query_end,target_chrom,min_synteny = .8, min_totlen = 500000):
	if query_start < 0:
		return False
	if query_end > paf.query_genome.scaffolds[query_chrom]:
		return False
	alignments = paf.fetch_query(query_chrom,query_start,query_end)
	target_synteny = calculate_synteny_target(alignments)
	if target_chrom in target_synteny:
		if target_synteny[target_chrom]['fraction'] >= min_synteny:
			if target_synteny[target_chrom]['total_length'] >= min_totlen:
				return True
	return False


def filter_fission(paf,fission,blocksize = 3000000):
	if fission['min_pos'] - blocksize < 0:
		return False
	downstream_region = (fission['target_chrom'],fission['min_pos'] - blocksize,fission['min_pos'])
	upstream_region = (fission['target_chrom'],fission['max_pos'],fission['max_pos'] + blocksize)
	if is_syntenic(paf,downstream_region[0],downstream_region[1],downstream_region[2],fission['query_pos_1'][0], min_synteny=.5):
		if is_syntenic(paf,upstream_region[0],upstream_region[1],upstream_region[2],fission['query_pos_2'][0], min_synteny=.5):
			return True
	return False

def filter_fusions(paf,fusion,blocksize = 3000000):
	# first check upstream
	qchrom = fusion['query_chrom']
	qstart = fusion['min_pos'] - blocksize
	qend = fusion['min_pos']
	tchrom = fusion['target_pos_1'][0]
	if qstart < 0:
		return False
	if is_syntenic_query(paf,qchrom,qstart,qend,tchrom):
		qchrom = fusion['query_chrom']
		qstart = fusion['max_pos']
		qend = fusion['max_pos'] + blocksize
		tchrom = fusion['target_pos_2'][0]
		if is_syntenic_query(paf,qchrom,qstart,qend,tchrom):
			return True
	return False

# Parse the PAF file
paf = PAF()
paf.parse_paf(aln,target_genome=refindex,query_genome=queryindex)

# filter paf
paf.filter_alignments(minlen = minlen, minmapq=minq,minpid = minpid,only_ref_chroms=only_ref_chroms,only_query_chroms=only_query_chroms)

# sort on target 
paf.sort_on_target()

# flip query chroms that predominantly align to the minus strand
for chrom in paf.query_genome.scaffolds:
	if paf.get_dominating_strand(chrom) == "-":
		paf.flip_query_alignments(chrom)
# resort (should not be necessary)
paf.sort_on_target()

# call fissions
fissions = call_fissions(paf)
# sort on query
paf.sort_on_query()
# call fusions
fusions = call_fusions(paf)


# filter fissions, by requiring that up- and downstream regions are consistently syntenic with the same query/reference chromosomes
results = []
for f in fissions:
	if filter_fission(paf,f, blocksize=args.blocksize):
		results.append(f)

# filter fusions
for f in fusions:
	if filter_fusions(paf,f,blocksize=args.blocksize):
		results.append(f)

# print results
with open(calls_out,'w') as out:
	out.write("type\ttarget_chrom_1\ttarget_pos_1_min\ttarget_pos_1_max\ttarget_chrom_2\ttarget_pos_2_min\ttarget_pos_2_max\tquery_chrom_1\tquery_pos_1_min\tquery_pos_1_max\tquery_chrom_2\tquery_pos_2_min\tquery_pos_2_max\n")
	for r in results:
		if r['type'] == "fission":
			out.write("fission\t{target_chrom_1}\t{min_pos}\t{max_pos}\tNA\tNA\tNA\t{query_chrom_1}\t{query_pos_1}\tNA\t{query_chrom_2}\t{query_pos_2}\tNA\n".format(target_chrom_1=r['target_chrom'],min_pos=r['min_pos'],max_pos=r['max_pos'],query_chrom_1=r['query_pos_1'][0],query_pos_1=r['query_pos_1'][1],query_chrom_2=r['query_pos_2'][0],query_pos_2=r['query_pos_2'][1]))
		elif r['type'] == "fusion":
			out.write("fusion\t{target_chrom_1}\t{target_pos_1}\tNA\t{target_chrom_2}\t{target_pos_2}\tNA\t{query_chrom_1}\t{min_pos}\t{max_pos}\tNA\tNA\tNA\n".format(target_chrom_1=r['target_pos_1'][0],target_pos_1=r['target_pos_1'][1],target_chrom_2=r['target_pos_2'][0],target_pos_2=r['target_pos_2'][1],min_pos=r['min_pos'],max_pos=r['max_pos'],query_chrom_1=r['query_chrom']))
		
# prep and write synmap
paf.target_genome.add_cumulative_startpos(scaffold_gap=10000000)

# merge consecutive alignments
paf.merge_consecutive(max_gap=100000)

with open(synmap_out, 'w') as out:
	out.write(paf.paf2synmap(scaffold_gap=10000000))

# write genomes
with open(genomes_out,'w') as out:
	out.write("genome\tchrom\tlength\tcumulative_start\n")
	for chrom in paf.target_genome.scaffolds:
		out.write("target\t{chrom}\t{length}\t{cum_start}\n".format(chrom=chrom,length=paf.target_genome.scaffolds[chrom],cum_start=paf.target_genome.cumulative_startpos[chrom]))
	for chrom in paf.query_genome.scaffolds:
		out.write("query\t{chrom}\t{length}\t{cum_start}\n".format(chrom=chrom,length=paf.query_genome.scaffolds[chrom],cum_start=paf.query_genome.cumulative_startpos[chrom]))
