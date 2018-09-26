# Designed to demultiplex overlaped amplicon sequencing data according to primer seuqence. It checks if the forward and revers primer sequences are in the correct position for each read. In case they are it saves sequences containing the same primers into the same file.
# usage: python primer_demultiplex.py Directoty_containing_not_demultiplex_fastq_files primer_combination_file.csv maximum_mismatch_value(integer) minimum_read_length(integer) output_directory
# Primers combination file is a tab-seperated file with three collumns: 1st locus name, 2nd primer forward sequence, 3rd primer reverse sequence
# Files containign sequences from the same marker can be saved in different folders. To do that, uncoment the lines in the last block and use the command "output=open(outfile_dir+"/"+fastq_pref[0]+"_"+locus+".fastq", "w")" instead the one in use.

import sys
import os
from Bio import SeqIO


def rev_comp(seq):
	code={'A':'T','C':'G','T':'A','G':'C'}
	result_sequence=""
	for i in seq:
		complement=code.get(i.upper(),i)
		result_sequence=complement+result_sequence
	return result_sequence

def mismatch(seq_a, seq_b):
	len1= len(seq_a)
	len2= len(seq_b)
	mismatches = 0
	for pos in range (0,min(len1,len2)) :
		if seq_a[pos] != seq_b[pos]:
			mismatches+=1
	return mismatches

def extract_primer_comb(fastq, primerF, primerR, mm, l):
	for rec in SeqIO.parse(fastq, 'fastq'):
		read=rec.seq
		motifF=str(read[:len(primerF)])
		motifR=str(read[-len(primerR):])
		if mismatch(motifF, primerF) <= mm and mismatch(motifR, primerR) <= mm:
			if len(read) >= l:
				yield rec


primer_dict={}
primer_file=open(sys.argv[2])
for line in primer_file:
	line_list=line.rstrip('\n').split('\t')
	primers = line_list[1:]
	primer_dict[line_list[0]]=primers


indir=sys.argv[1]
mm=int(sys.argv[3])
lenght=int(sys.argv[4])
outdir=sys.argv[5]



for fastq in os.listdir(indir):
	print("processing sample "+fastq+" ...")
	fastq_pref=fastq.split(".")
#	outfile_dir=outdir+fastq_pref[0]+"/"
#	os.mkdir(outfile_dir)
	for locus, primers in primer_dict.items():
		primerF = primers[0]
		primerR = rev_comp(primers[1])
		output=open(outdir+fastq_pref[0]+"_"+locus+".fastq", "w") #output=open(outfile_dir+"/"+fastq_pref[0]+"_"+locus+".fastq", "w")
		print("processing primer "+locus)
		SeqIO.write(extract_primer_comb(indir+fastq, primerF, primerR, mm, lenght), output, "fastq")
		output.close()


