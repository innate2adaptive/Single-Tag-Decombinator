from IPython import embed

import os
import sys
import argparse
import operator
from functools import partial
import multiprocessing as mp

import re
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

import time


def args():
	parser = argparse.ArgumentParser( description='** script to find overlaps between fragments of TCR sequence and rebuild complete sequences. **')
	parser.add_argument('-f', '--filename', type=str, help='File of sequences to be analysed', required=False)
	parser.add_argument('-np', '--nproc', type=int, help='Number of cores for multprocessing alignment', required=False, default=None)

	return parser


class TCR:
	def __init__(self,vread = None):
		self.vread = vread
		self.jread = None
		self.v_id = vread[3]
		self.j_id = None
		self.chain = vread[0]
		self.alignments = []
		self.ranked_alignments = []
		self.chosen_alignment = None
		self.longest_overlap = None
		self.sequence = None

	def determineAlignments(self, jreads,min_o):
		s1 = getSequence(self.vread)
		for j in jreads:
			if j[0] != self.chain:
				continue
			s2 = getSequence(j)[:-20]					
			alignments = align(s1,s2,min_o)
			
			for a in alignments:
				self.alignments.append(Alignment(a,self.v_id, j))
		return self.alignments

	def rankAlignmentLengths(self):
		if self.alignments == []:
			self.longest_overlap = 0
			return 0
		lengths = sorted( union( map(lambda x: x.length,self.alignments), [] ), reverse=True )
		# prioritise based on upon length of overlap
		for a in self.alignments:
			a.lrank = lengths.index(a.length) + 1
		self.longest_overlap = max(lengths)
		return lengths

	def rankAlignmentPurities(self):
		# prioritise based on upon purity of overlap
		if self.alignments == []:
			return 0
		purities = sorted( union( map(lambda x: x.purity,self.alignments), [] ), reverse=True )
		for a in self.alignments:
			a.prank = purities.index(a.purity) + 1
		return purities

	def setPriorities(self):
		if self.alignments == []:
			return 0
		order = []
		#first rank by length
		lranked = sorted(self.alignments,key=operator.attrgetter('lrank'))
		
		#divide rank lengths into partitions for sub-ranking
		partition = []
		for l in range(len(lranked)):
			if partition == []:
				partition.append(lranked[l])
			elif lranked[l].length == partition[len(partition)-1].length:
				partition.append(lranked[l])
			else:
				order.append(partition)
				partition = [lranked[l]]
		order.append(partition)
		
		#then rank parititons by purity
		for j in order:
			j.sort(key=operator.attrgetter('prank'),reverse=True)

		self.ranked_alignments = flatten(order)

		return self.ranked_alignments

	def setSequence(self):

		overlap = list(self.chosen_alignment.alignment[0])
		for i in range(len(overlap)):
			if overlap[i] == "-":
				overlap[i] = self.chosen_alignment.alignment[1][i]
		overlap = "".join(overlap)
		self.sequence = self.vread[4][:-len(overlap)] + overlap + self.jread[4][len(overlap):]
		return self.sequence


class Alignment:
	def __init__(self, alignment, v_id, j):
		self.alignment = alignment
		self.length = alignment[4]
		self.score = alignment[2]
		self.purity = self.length - self.score
		self.v_id  = v_id
		self.jread = j
		self.j_id = j[3]
		self.lrank = None
		self.prank = None


def getSequence(read):
	return read[4]

def union(a, b):
	return list(set(a) | set(b))

def flatten(l):
	return [item for sublist in l for item in sublist]

def align(s1,s2,min_o):

	half1  = s1[-min_o:-min_o/2]
	half2 = s1[-min_o/2:]

	rel_overlaps = []
	good_alignments = []

	half1_matches = [m.start() for m in re.finditer('(?='+half1+')', s2)]
	half2_matches = [m.start() - min_o/2 for m in re.finditer('(?='+half2+')', s2)]

	half_matches = union(half1_matches,half2_matches)

	for i in half_matches:
		aligns = pairwise2.align.globalms(half1+half2,s2[i:i+min_o],1,0,-.5,-0.1)

		for k in aligns: 
			if ( k[4] - k[2] < 3 ) or ( k[4] - k[2] == 3 and "-" in k[1] ):
				#rel_overlaps.append([i,k])
				s2start =  s2[:i + min_o]
				s1end = s1[-(i + min_o):]
				aligns2 = pairwise2.align.globalms(s2start,s1end,1,0,-.5,-0.1)

				for j in aligns2:
					if ( j[4] - j[2] < 3 ) or ( j[4] - j[2] == 3 and "-" in j[1] ):
						good_alignments.append(j)
				

	# for j in rel_overlaps:
	# 	s2start =  s2[:j[0] + min_o]
	# 	s1end = s1[-(j[0] + min_o):]
	# 	aligns = pairwise2.align.globalms(s2start,s1end,1,0,-.5,-0.1)


	# 	for k in aligns: 
	# 		if ( k[4] - k[2] < 3 ) or ( k[4] - k[2] == 3 and "-" in k[1] ):
	# 			good_alignments.append(k)

	return good_alignments


def reconstruct(tcr,jreads):
	tcr.determineAlignments(jreads,8)
	tcr.rankAlignmentLengths()
	tcr.rankAlignmentPurities()
	tcr.setPriorities()
	return tcr




def main(args):
	total_time = time.time()
	file = args.filename	
	cores = args.nproc
	if not cores: cores = mp.cpu_count()

	lines = []
	with open(file) as f:
		for line in f:
			lines.append(line)	
	

	reads = [l.rstrip().split(", ") for l in lines]	
	
	vreads = []
	jreads = []	

	for r in reads:
		if r[1] == 'n/a':
			jreads.append(r)
		if r[2] == 'n/a':
			vreads.append(r)

	tcrs = []


	for v in vreads:
		tcrs.append(TCR(vread = v))

	usedjs = []

	reads_count = 0
	print "Aligning Reads..."
	print "pooling with " + str(cores)+" cores:"
	aligned_tcrs = []

	part_tcrs = [tcrs[i:i + 100] for i in xrange(0, len(tcrs), 100)]

	for t in part_tcrs:
		start = time.time()
		pool = mp.Pool(processes=cores)
		#results = [pool.apply(reconstruct, args=(tcr,jreads)) for tcr in tcrs[0:200]]
		results = pool.imap(partial(reconstruct,jreads=jreads),t)

		for x in results:
			aligned_tcrs.append(x)
		print str(len(aligned_tcrs)), "aligned"
		print time.time() - start
		pool.close()


	tcrs = aligned_tcrs

	
	# tcrs with longest overlap alignments get priority for matching
	print "Reconstructing TCRs..."

	for tcr in sorted(tcrs,key=operator.attrgetter('longest_overlap'),reverse=True):

		for a in tcr.ranked_alignments:
			if a.j_id in usedjs:
				continue
			else:
				tcr.j_id = a.j_id
				tcr.jread = a.jread
				tcr.chosen_alignment = a
				usedjs.append(a.j_id)
				break

		# tcr.chosen_alignment = tcr.ranked_alignments[0]
		# tcr.j_id = tcr.ranked_alignments[0].j_id
		# tcr.jread = tcr.ranked_alignments[0].jread
		# jreads.remove(tcr.jread)
		if not tcr.chosen_alignment:
			tcrs.remove(tcr)
		else:
			tcr.setSequence()

	outfile = 'bfd'+os.path.splitext(os.path.basename(args.filename))[0] + ".fastq"
	print "writing to", outfile
	with open(outfile, "w") as f:	
		for tcr in tcrs:
			new_read = "@"+tcr.v_id+"\n"
			new_read += tcr.sequence+"\n"
			new_read += "+\n"
			new_read += "~"*len(tcr.sequence)+"\n"	
			f.write(new_read)

	print "TOTAL TIME: "+str(time.time() - total_time)
	return outfile


if __name__ == '__main__':
	parser = args()
	args = parser.parse_args()
	main(args)

