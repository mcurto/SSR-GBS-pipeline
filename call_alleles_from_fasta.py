# Converts corrected sequences into allelic information. To each unique sequence one number is atributed. All equal sequences will have the same number. This considers Ns. Because of Ns some sequences may be atributed to multiple numbers. In that case it will be saved as missing data.
# Tow outputs are saved. Fist, a codominante matrix containing the alleles numbers (tab delimited). In case of missing data the saved genotype is "0 0". Second, a list of alleles per marker.
# Usage: python call_alleles_from_fasta Directory_containing_alleles_sequences_in_fasta_format prefix_in_common_for_all_markers List_of_sample_name prefix_to_save_output
# Fasta files should contain all sequences from a sepecific marker. Sample list shoudl contain one sample name per line witout spaces or tabs

from __future__ import division
import sys
import os


# functions:

## parse fasta

def parse_fasta(fasta):
    f = open(fasta)
    data = {}
    count = 0
    for l in f:
        if l.startswith('>'):
            head = l.rstrip('\r\n')[1:]
            if head in data:
                head += '_2'
                data[head] = ''
            else:
                data[head] = ''
        else:
            oldSeq = data[head]
            data[head] = oldSeq + l.rstrip('\r\n').upper()
    f.close()
    return data

# msimatch funtion including ambiguous bases
def mismatch(seq_a, seq_b):
    IUPAC = {'A':['A', 'N'], 'C':['C', 'N'], 'G':['G', 'N'], 'T':['T', 'N'], '-':['-'], 'R':['A','G', 'N'], 'Y':['C','T', 'N'], 'S':['G','C', 'N'], 'W':['A','T', 'N'], 'K':['G','T', 'N'], 'M':['A','C', 'N'], 'B':['C','G','T', 'N'], 'D':['A','G','T', 'N'], 'H':['A','C','T', 'N'], 'V':['A','C','G', 'N'], 'N':['A','T','G','C', 'N'], 'I':['A','T','G','C', 'N']}
    len1= len(seq_a)
    len2= len(seq_b)
    mismatches = 0
    for pos in range (0,min(len1,len2)):
        base_b = seq_b[pos]
        base_b_IUPAC = IUPAC[base_b]
        if not seq_a[pos] in base_b_IUPAC:
            mismatches+=1
    return mismatches




# this fucntion compares the sequences from one fasta file and if they match completely to each other it atributes them to the same allele
def get_allele_dict(fasta_parsed):

# get unique sequences from fasta
    seq_list = []
    for name, seq in fasta_parsed.items():
        seq_list.append(seq)
    seq_list = list(set(seq_list))




# sort uniq sequences according to N content
    Ns_per_seq = {} # {N: [seq1,seq2]}
    sorted_seq_list = []
    Ns_per_seq_list = []
    for seq in seq_list:
        countN = seq.count('N')
        Ns_per_seq_list.append(countN)
        count_to_dict = Ns_per_seq.get(countN, []) + [seq]
        Ns_per_seq[countN] = count_to_dict
    Ns_per_seq_list = list(set(Ns_per_seq_list))
    Ns_per_seq_list.sort()
    for N in Ns_per_seq_list:
        sorted_seq_list += Ns_per_seq[N]



# makes a dictionary like folows : {allele: [sequences that match each other], ...}
## in this dict each sequence list correspond to all possible sequenes varianths for the same allele

    alleles_dict = {}
    for seqa in sorted_seq_list:
        match = 0
        current_allele = len(alleles_dict) + 1
        if len(alleles_dict) == 0:
            alleles_dict['1'] = [seqa]
        else:
            for alelle, seqs in alleles_dict.items():
                for seqb in seqs:
                    if mismatch(seqa, seqb) == 0:
                        current_allele = alelle
                        match += 1

            if match == 0:
                new_list = [seqa]
                alleles_dict[str(current_allele)] = new_list
            elif match == 1:
                new_list = alleles_dict[str(current_allele)] + [seqa]
                alleles_dict[current_allele] = new_list
            else:
                alleles_dict['MC'] = alleles_dict.get('MC', []) + [seqa]


    return alleles_dict

# saves a text file with which alleles exist per locus and which sequences correspond to those alleles
#locus1
#allele1{tab}seq1
#allele2{tab}seq2
#allele2{tab}seq3
#
#locus2
#...
def save_allele_database(alleles_dict, locus, out):
    out.write(locus + '\n')
    for allele, allele_seq in alleles_dict.items():
        for seq in allele_seq:
            out.write(allele + ':\t' + seq + '\n')
    out.write('\n')



# compare sequences from fasta to alleles and call them:
def compare_to_alleles(alleles_dict, fasta_parsed, locus):
    result = {}
    for head, seq in fasta_parsed.items():
        alleles_seq = []
        name = head.split('_')[2]
        for allele, allele_seq in alleles_dict.items():
            if seq in allele_seq:
                alleles_to_dict = result.get(name, []) + [allele]
                result[name] = alleles_to_dict

    return result



# directory containing fasta alignemts per locus
inDir = sys.argv[1]
# common prefix across all alognemt files
prefix = sys.argv[2]
# list of sample names:
samples_names = open(sys.argv[3])
# saves samples list into a list
samples_list = []
for l in samples_names:
    samples_list.append(l.rstrip('\n'))

#prefix for output files
out_pre  = sys.argv[4]

#opens dictionary to store genotypes
#it already saves the header which willcontain loci info
out_dict = {'samples' : '\t'}
for sample in samples_list:
    out_dict[sample] = ''

#it opens an output file to save the alleles list
out_list = open(out_pre + 'allelle_list.txt', 'w')

# iterates through all alignemet files and it saves the genotypes in a dictionary like this:
# {sample: genotype_locus1\tgenotype_locus2\t ...}
for f in os.listdir(inDir):
    if f.startswith(prefix):
        locus = f.split('_')[0] + '_' + f.split('_')[1] # get locus

        out_dict['samples'] += (locus + '\t')*2 #saves the current locus in the dictionary
        parsed_fasta = parse_fasta(inDir + f) #parse fasta

        alleles_dict = get_allele_dict(parsed_fasta) #get lsit of alleles
        names_alleles = compare_to_alleles(alleles_dict, parsed_fasta, locus) # get alleles for each sample
        save_allele_database(alleles_dict, locus, out_list) # save allele list
        # create genotypes: if two alleles -> allele1{tab}allele2; if one allele -> allele1{tab}allele1; if sample name not found in allele per sample list (missing data) -> 0{tab}0
        for sample in samples_list:
            if sample in names_alleles:
                alleles = names_alleles[sample]
                alleles.sort()
                if len(alleles) == 2:
                    genotype = '\t'.join(alleles)
                elif len(alleles) > 2:
                    print('locus ' + locus + ' for sample ' +  sample + '\thas ' + str(len(alleles)) + ' alleles: ' + str(alleles))
                    genotype = 'MC\tMC'
                else:
                    genotype = '\t'.join(alleles*2)
            else:
                genotype = '0\t0'
            out_dict[sample] += '\t' + genotype #it saves the genotypes in the dictionary
out_list.close()


out = open(out_pre + 'matrix.txt', 'w') #it opens output file to save the matrix
out.write('samples' + out_dict['samples'] + '\n') #it saves the matrix header with loci information

#saves the matrix the genotypes
for sample in samples_list:
    out.write(sample + out_dict[sample] + '\n')

out.close()
