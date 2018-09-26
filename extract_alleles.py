# It uses a codominant matrix as an input containing alleles defined according to length. Then it sepeartes the sequences from each allele in a seperated file per sample.
# Usage: python extract_alleles.py codominat_matrix.txt(tab-delimited text) Directory_containing_sequence_per_marker_annd_per_sample(fasta or fastq) file_type(fasta/fastq) Output_Directory 

from Bio import SeqIO
import sys


# only keeps sequences with a defined length
def extractLengthFasta(fasta, length):
    for rec in SeqIO.parse(fasta, 'fasta'):
        if length == len(rec):
            yield rec

def extractLengthFastq(fastq, length):
    for rec in SeqIO.parse(fastq, 'fastq'):
        if length == len(rec):
            yield rec

# three inputs: 1) matrix containing lengths; 2) directory containing fasta files; 3) directory to save results
allelesFile = open(sys.argv[1])
FastaDir = sys.argv[2]
fileType = sys.argv[3]
outDir = sys.argv[4]


# this part parses the matrix and uses the function from before to extract and then save the correct sequences
for l in allelesFile:
    led = l.rstrip('\n').split(',')
    Sample = led[0]
    print('processing sample ' + Sample)
    if Sample == 'samplename': # if line starts with Mix it means it is the header
        loci = led[1:] # saves the loci information
    else:
        alleles = led[1:] # saves the alleles lengths for a certain sample
        for i in range(0, len(loci), 2): # goes locus by locus and ...
            locus = loci[i]
            al1 = alleles[i]
            al2 = alleles[i + 1] # extarcts information for each locus / sample combination from the matrix
            print('...processing locus ' + locus)
            if al1 != '0' and al2 != '0':
                print('......genotype: ' + al1 + ',' + al2)
                if fileType == 'fasta':
                    f = FastaDir + Sample + '_' + locus + '.fasta'
                elif fileType == 'fastq':
                    f = FastaDir + Sample + '_' + locus + '.fastq'
                else:
                    print('please provide a valid file format')
                for al in list(set([al1, al2])):
                    outFasta = open(outDir + locus + '_' + Sample +  '_Al_' + al + '.fasta', 'w') # creates an output file
                    if fileType == 'fasta':
                        SeqIO.write(extractLengthFasta(f, int(al)), outFasta, 'fasta') # sequences for each sample/allele combination
                    elif fileType == 'fastq':
                        SeqIO.write(extractLengthFastq(f, int(al)), outFasta, 'fasta') # sequences for each sample/allele combination
                    else:
                        print('please provide a valid file format')
                    outFasta.close()
            else:
                print('......genotype: missing')
