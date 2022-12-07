#!/usr/bin/env python

# Caleb Lareau, Stanford University
# Implemented: 14 November 2020
# This program will error correct barcodes
# From 10x sequencing data from scATAC

##### IMPORT MODULES #####
import os
import re
import regex
import sys
import gzip
import itertools

from functools import partial
from optparse import OptionParser
from multiprocessing import Pool, freeze_support
from itertools import repeat
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from fuzzywuzzy import fuzz
from fuzzywuzzy import process
from fuzzysearch import find_near_matches

#### OPTIONS ####
opts = OptionParser()
usage = "usage: %prog [options] [inputs] Software to process raw .fastq reads and make data suitable for downstream processes"

opts.add_option("-a", "--fastq1", help="<Read1> Accepts fastq.gz")
opts.add_option("-b", "--fastq2", help="<Read2> Accepts fastq.gz")
opts.add_option("-c", "--fastq3", help="<Read3> Accepts fastq.gz")
opts.add_option("-f", "--iseven", help="i5 barcode to parse")

opts.add_option("-n", "--nreads", default = 10000000, help="Number of reads to process in a given chunk")
opts.add_option("-r", "--ncores", default = 8, help="Number of cores for parallel processing.")

opts.add_option("-o", "--output", help="Output sample convention")

options, arguments = opts.parse_args()


# return usage information if no argvs given
if len(sys.argv)==1:
	os.system(sys.argv[0]+" --help")
	sys.exit()

def batch_iterator(iterator, batch_size):
	"""
	Returns lists of tuples of length batch_size.
	"""
	entry = True  # Make sure we loop once
	while entry:
		batch = []
		while len(batch) < batch_size:
			try:
				entry = iterator.__next__()
			except StopIteration:
				entry = None
			if entry is None:
				# End of file
				break
			batch.append(entry)
		if batch:
			yield batch


def chunk_writer_gzip(filename, what):
	'''
	Basic function to write a chunk of a fastq file
	to a gzipped file
	'''
	with gzip.open(filename, 'wt') as out_write:
				out_write.writelines(what)
	return(filename)			

def prove_barcode_simple(bc, valid_set):
	'''
	Function that takes a putative barcode and returns the nearest valid one
	'''
		
	if(bc in valid_set):
		return(bc)
	else:
		return("NA")

			
def formatRead(title, sequence, quality):
	"""
	Takes three components of fastq file and stiches them together in a string
	"""
	return("@%s\n%s\n+\n%s\n" % (title, sequence, quality))

					
#------ 
	
def debarcode_trio(trio, i7seq):
	"""
	Function that is called in parallel
	"""
	# Parse out inputs
	listRead1 = trio[0]; listRead2 = trio[1]; listRead3 = trio[2]
	
	# parameters to return
	fq1 = ""
	fq2 = ""
	
	'''
	Function that takes a putative barcode and returns the nearest valid one
	'''

	# Grab attributes
	title1 = listRead1[0]; sequence1 = listRead1[1]; quality1 = listRead1[2]
	title2 = listRead2[0]; sequence2 = listRead2[1]; quality2 = listRead2[2]
	title3 = listRead3[0]; sequence3 = listRead3[1]; quality3 = listRead3[2]

	bc = sequence2[0:8]
	i7seq_revcomp = Seq(i7seq).reverse_complement()
	match_iseven = bc == i7seq_revcomp
	#print(bc)
	
	# Return the barcode with underscores + the biological sequence learned
	ofq1 =""
	ofq2 = ""

	if(match_iseven):
		ofq1 = formatRead(title1, sequence1, quality1)
		ofq2 = formatRead(title3, sequence3, quality3)
	return(ofq1, ofq2)


if __name__ == "__main__":


	##### INPUTS #####
	af = options.fastq1
	bf = options.fastq2
	cf = options.fastq3
	
	i7 = options.iseven

	outname = options.output
	o = options.output
	cpu = int(options.ncores)
	n = int(options.nreads)

	# Parse input files
	extension = af.split('.')[-1]
	if extension == "fastq" or extension == "fq":
		sys.exist("Quitting... GZIP your .fastq files!")
	elif extension == "gz":
		print("Found supplied .fastq.gz files")
	else:
		sys.exit("ERROR! The input files (-a , -b, -c) a *.fastq.gz")
	print(options)
	with gzip.open(af, "rt") as f1:
		with gzip.open(bf, "rt") as f2:
				with gzip.open(cf, "rt") as f3:

					# Establish iterators
					it1 = batch_iterator(FastqGeneralIterator(f1), n)
					it2 = batch_iterator(FastqGeneralIterator(f2), n)
					it3 = batch_iterator(FastqGeneralIterator(f3), n)

					# iterate over batches of length n
				
					for i, batch1 in enumerate(it1):
						batch2 = it2.__next__()
						batch3 = it3.__next__()
						output = o 
			
						# parallel process the barcode processing and accounting of failures.
						pool = Pool(processes=cpu)
						pm = pool.map(partial(debarcode_trio, i7seq=i7), zip(batch1, batch2, batch3))
						pool.close()
			
						# Aggregate output
						fastq1 = [item[0] for item in pm]
						fastq2 = [item[1] for item in pm]

						# Export one chunk in parallel
						filename1 = output +'_R1.fastq.gz'
						filename2 = output +'_R2.fastq.gz'
			
						pool = Pool(processes=2)
						toke = pool.starmap(chunk_writer_gzip, [(filename1, fastq1), (filename2, fastq2)])
						pool.close()
			
	