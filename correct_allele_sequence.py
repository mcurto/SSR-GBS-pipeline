# It first records positions with Ns for each sequence of a fasta file containing all consensus sequences from one marker. Then it checks the frequency of the nucleotide combinations taken from the postions in the sequences extracted. The two most frequent combinations are used to devide and correct the sequences. In case a sample is heterozigote for length and it contains consensus sequences with Ns only the most frequent nucleotide is taken.
# Usage: python correct_allele_sequence.py fasta_file_containing_all_consensus_sequences_per_marker.fasta Directory_containing_all_seuqences_extracted_per_length_allele minimum number_of_reads_considered per_consensus output_fasta_file.fasta


from collections import Counter
import sys

## parse fasta and save it in the dict.: {header : Sequence, ...}

def parse_fasta(fasta):
    f = open(fasta)
    data = {}
    for l in f:
        if l.startswith('>'):
            head = l.rstrip('\r\n')[1:]
            data[head] = ''
        else:
            oldSeq = data[head]
            data[head] = oldSeq + l.rstrip('\r\n').replace('~', '-').upper()
    f.close()
    return data


def get_Ns_records(st):
    result = []
    c = 0
    for v in st:
        if v.upper() == 'N':
            result.append(c)
        c += 1
    result.sort()
    return result

def filter_Ns_records(positions_dict, fasta_parse):
    max_len = 0
    for r, s in fasta_parse.items():
        max_len =  len(s)
        break

    possible_positions = []
    for i in range(0,max_len):
        nuc = []
        for rec, seq in fasta_parse.items():
            nuc.append(seq[i].upper())
        nuc_un = list(set(nuc))
        if 'N' in nuc_un:
            nuc_un.remove('N')
            if '-' in nuc_un:
                nuc_un.remove('-')
            if len(nuc_un) > 1:
                possible_positions.append(i)
    result = {}
    for name, rec in positions_dict.items():
        new_recNoGap = []
        new_recGap = []
        rec_posNoGap = rec[0]
        rec_posGap = rec[1]
        for p in rec_posGap:
            if p in possible_positions:
                new_recNoGap.append(rec_posNoGap[rec_posGap.index(p)])
                new_recGap.append(p)
        new_recNoGap.sort()
        new_recGap.sort()
        result[name] = (rec_posNoGap, rec_posGap)
    return result



def get_seq_freq(fasta, pos_pair):
    haps = []
    for name, seq in parse_fasta(fasta).items():
        hap = ''
        for i in pos_pair:
            hap += seq[i]
        haps.append(hap)
    return dict(Counter(haps))

def get_to_el_from_dict(dict, nrElments):
    temp_count = []
    for hap, c in dict.items():
        temp_count.append(c)
    temp_count.sort()
    temp_count_new = temp_count[-nrElments:]
    result = []
    for hap, c in dict.items():
        if c in temp_count_new:
            result.append(hap)
    return result


fasta = sys.argv[1] # alignement
inDir = sys.argv[2] # directory containing alleles seperated
minCount = int(sys.argv[3]) # minimum number of sequences per fasta


positions = {}
names_info = {}
for name, seq in parse_fasta(fasta).items():
    positions[(name, seq)] = (get_Ns_records(seq.replace('-', '')), get_Ns_records(seq))
    positions_filt = filter_Ns_records(positions, parse_fasta(fasta))
    sample = name.split('_')[2]
    if sample in names_info:
        names_info[sample] = 1
    else:
        names_info[sample] = 2

out = open(sys.argv[4], 'w') # out fasta with corrected and filtered


for name_seq, positionS in positions_filt.items():
    name = name_seq[0]
    seq = name_seq[1]
    position = positionS[0]
    positionAl = positionS[1]
    nr_seqs = int(name.split('_')[-1])


    fileName = inDir + '_'.join(name.split('_')[:5]) + '.fasta'
    sample = name.split('_')[2]
    het_info = names_info[sample]
    counter_dict = get_seq_freq(fileName, position)
    haps = get_to_el_from_dict(counter_dict, het_info)

    if len(haps) > het_info:
        haps = ['N'*len(haps[0])]

    out_temp = []

    if nr_seqs > minCount:
        for hap in haps:
            temp = ''
            for p, i in enumerate(seq, 0):
                if i == 'N' and p in positionAl:
                    temp += hap[positionAl.index(p)]
                else:
                    temp += i
            out_temp.append(temp)



    if len(out_temp) <= 2:
        for c, seq in enumerate(out_temp, 0):
            out.write('>' + name + '_' + str(c) + '\n' + seq + '\n')
