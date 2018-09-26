# Takes a fasta file containing multiple sequences and it creates a consensus sequence and a table containing the nucleotide frequency per position. The consensus threashold is user defined.
# Usage: python get_consensus_and_freq.py Directory_containing_input_fasta_files Directory_containing_output_files Consensus_threashold(float)

from __future__ import division
import sys
import os


def parse_fasta(fasta):
	f = open(fasta)
	data = {}
	for l in f:
		if l.startswith('>'):
			head = l.rstrip('\n')
			data[head] = ''
		else:
			oldSeq = data[head]
			data[head] = oldSeq + l.rstrip('\n')
	f.close()
	return data


def process_file(fasta, consThreashold, outDir):
	fasta_dir = parse_fasta(fasta)
	nr_seqs = len(fasta_dir)

	# get maximum length
	max_len = 0
	for rec, seq in fasta_dir.items():
		if len(seq) > max_len:
			max_len = len(seq)

	# complete the ends of smaller sequences with '-'
	fasta_new = {}
	for rec, seq in fasta_dir.items():
		lenSeq = len(seq)
		newSeq = seq + '-'*(max_len - lenSeq)
		fasta_new[rec] = newSeq

	# creates a dictionary with nucleotide composition per alignemt collumn
	data = {}
	for i in range(0,max_len):
		nuc = []
		for rec, seq in fasta_new.items():
			nuc.append(seq[i])
		data[i] = nuc

	outName = fasta.split('.')[0].split('/')[-1] + '_C_' + str(nr_seqs)
	out_freqs = open(outDir + outName + '.freq', 'w')
	out_freqs.write('pos.\tA\tC\tG\tT\n')
	OutConSeq = open(outDir + outName + str(int(consThreashold*100)) + '.fasta', 'w')
	conSeq = '>' + outName + '\n'

	# saves the frequency of each base in a freq file and the consensus sequences in another. A threashold should be defined for the consensus
	for rec, nuc in data.items():
		temp = str(rec)
		base = 'N'
		for n in ['A', 'C', 'G', 'T', '-']:
			freq = nuc.count(n)/nr_seqs
			temp += '\t' + str(freq)
			if freq >= float(consThreashold):
				base = n
		out_freqs.write(temp + '\n')
		conSeq += base
	OutConSeq.write(conSeq + '\n')

inDir = sys.argv[1]
outDir = sys.argv[2]
consThreashold = float(sys.argv[3])

for f in os.listdir(inDir):
	process_file(inDir + f, consThreashold, outDir)
