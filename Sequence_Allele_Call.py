# Second part of the SSR-GBS pipeline and defines alleles based on seuqence information

# Version 0.2, from 24.06.2020

# It uses the corrected length allele matrix and the demultiplexed overlaped sequences as input (both fasta and fastq formats are suported)
# The script runs four steps:
## 1) Extract reads corresponding to allele length
## 2) Make consensus sequences with user defined threashold
## 3) Join seuqences per locus
## 4) Detect SNPs, and in case they exist it outputs two sequences (one per variant)
## 5) Call alleles based on sequence information. Two modes are possible:
##### 5.1) DeNovo, when it is the first time that the markers are being used
##### 5.2) List Based, when there is already data in the project and thus an allele list is available

# It does not require any third-party program

# To run it, it is necessary to provide:
### 1) matrix containing lengths;
### 2) the directory where to save the results;
### 3) consensus threashold;
### 4) minimum number of reads for a consensus sequence to be considered for SNP detection and allele call




from __future__ import division
from collections import Counter
import sys, os, re


# make directory
def makeDir(Directory):
    # if cannot make directory puts out error message
    try:
        os.mkdir(Directory)
    except OSError:
        print ("Creation of the directory %s failed" % Directory)
    else:
        print ("Successfully created the directory %s " % Directory)

# Parse fasta files by producing dictionary {heads: sequence}. Repeated heads will get the sufix _2
def parse_fasta(fasta):
    f = open(fasta)
    data = {}
    for l in f:
        if l.startswith('>'): # head starts with '>'
            head = l.rstrip('\r\n')[1:] # exclude new line and starting character '>'
            if head in data:
                head += '_2' # if head already existing in dictionary add _2 sufix so seuqneces with the same head are considered
                data[head] = '' # atribute empty sequnce to head
            else:
                data[head] = ''
        else:
            data[head] = data[head] + l.rstrip('\r\n').upper() # retrieve head and atribute seuqnece information. I did this to process interlieve formats
    f.close()
    return data

# parse fastq file. It makes the dictionary: {header: 'sequence'}
def parse_fastq(fastq):
    f = open(fastq)
    data = {} # dictionary to store result
    for c, l in enumerate(f): # iterate through lines enumerating them starting with 0
        if c % 4 == 0: # headers line numebers when devide by four should always have a remainder of 0
            name = l.rstrip('\r\n')[1:] # exclude new line and starting character '@'
        elif c % 4 == 1: # sequences line numeber when devide by four should always have a remainder of 1
            data[name] = l.rstrip('\r\n')
    f.close()
    return data

# function to write fasta file out of dictionary {heads: sequence}. head information should not have the starting charcter
# output file should already be open
def writeFasta(SeqDict, out):
    for rec, seq in SeqDict.items():
        out.write('>' + rec[1:] + '\n' + seq + '\n')

# detect sequence file format based on the extension. Input is the seperate merged seuqences. Both fasta and fastq are ok
def getFileType(Dir):
    for f in os.listdir(Dir):
        type = f.split('.')[-1] # get extension
        if type in ['fasta', 'fastq']: # extension should be either fasta or fastq
            return type
            break

# get marker and sample list from allele length codominat matrix
def getMarkerAndSampleList(allelesFile):
    allellesInfo =  open(allelesFile)
    sampleList = [] # variable list to save samples names

    for l in allellesInfo:
        led = re.split('\t|,|;', l.rstrip('\r\n').replace('"', '')) # parse matrix lines
        Sample = led[0]
        # if sample name is samplename it means that the line is the header and tus the locus information
        if Sample == 'samplename':
            loci = led[1:] # saves the loci information
            # exclude _1 and _2 from loci
            lociNew = ['']*len(loci)
            for i in range(0, len(loci)):
                lociNew[i] = '_'.join(loci[i].split('_')[:2])
            # get locus and make them unique
            lociNew = list(set(lociNew))
            # sort loci list
            lociNew.sort()
        else:
            # add sample name to sample list
            sampleList.append(Sample)

    return lociNew, sampleList
    allellesInfo.close()


# Fuctions to extract alleles from codominat matrix. The matrix should be coma separated
# Extracts reads from a fasta file with with a defined length.
def extractLengthFasta(fasta, length):
    result = {}
    for rec, seq in parse_fasta(fasta).items():
        if length == len(seq):
            result[rec] = seq
    return result

# Extracts reads from a fastq file with with a defined length.
def extractLengthFastq(fastq, length):
    result = {}
    for rec, seq in parse_fastq(fastq).items():
        if length == len(seq):
            result[rec] = seq
    return result


# this part parses the matrix and uses the function to extract reads with a specific length into a new file
def extract_alleles(allelesFile, FastaDir, fileType, outDir):
    allellesInfo =  open(allelesFile)

    for l in allellesInfo:
        led = re.split('\t|,|;', l.rstrip('\r\n').replace('"', '')) # parse matrix lines
        Sample = led[0].strip('"')
        print('processing sample ' + Sample)
        if Sample == 'samplename': # if line starts with Mix it means it is the header
            loci = led[1:] # saves the loci information
            # exclude _1 and _2 from loci
            lociNew = ['']*len(loci)
            for i in range(0, len(loci)):
                lociNew[i] = '_'.join(loci[i].split('_')[:2])
        else:
            alleles = led[1:] # saves the alleles lengths for a certain sample for all loci
            # correct alleles names
            allelesNew = ['']*len(alleles)
            for i in range(0, len(alleles)):
                al = alleles[i].strip('"')
                if al in ['empty', 'too little reads']:
                    allelesNew[i] = '0'
                else:
                    allelesNew[i] = al.split('_')[0]

            for i in range(0, len(loci), 2): # goes locus by locus and ...
                locus = lociNew[i]
                al1 = allelesNew[i] # get allele 1 (most reads)
                al2 = allelesNew[i + 1] # get allele 2 (less reads)
                print('...processing locus ' + locus)
                if al1 != '0' and al2 != '0': # missing data (0) is not considered
                    print('......genotype: ' + al1 + ',' + al2)
                    # get file type and input sequence file name
                    if fileType == 'fasta':
                        f = FastaDir + Sample + '_' + locus + '.fasta'
                    elif fileType == 'fastq':
                        f = FastaDir + Sample + '_' + locus + '.fastq'
                    else:
                        print('please provide a valid file format')
                    #extract sequences based on allele length
                    for al in list(set([al1, al2])):
                        outFasta = open(outDir + locus + '_' + Sample +  '_Al_' + al + '.fasta', 'w') # creates an output file. al is the allele length
                        if fileType == 'fasta':
                            ExtractedSeqs = extractLengthFasta(f, int(al)) # filter sequences (fasta input)
                            if len(ExtractedSeqs) == 0:
                                print('No sequence found for allele ' + al)
                            else:
                                writeFasta(ExtractedSeqs, outFasta) #  write to file
                        elif fileType == 'fastq':
                            ExtractedSeqs = extractLengthFastq(f, int(al)) # filter sequences (fastq input)
                            if len(ExtractedSeqs) == 0:
                                print('No sequence found for allele ' + al)
                            else:
                                writeFasta(ExtractedSeqs, outFasta) #  write to file
                        else:
                            print('please provide a valid file format')
                        outFasta.close()
                else:
                    print('......genotype: missing')
    allellesInfo.close()


# Functions to make consensus sequences
# make consensus per file. Input is the extracted sequences per allele length, a consensus threashold, and the output directory
def MakeConsensusPerFile(fasta, consThreashold, outDir):
    fasta_new = parse_fasta(fasta)
    nr_seqs = len(fasta_new) # number of sequences per file

    # get sequence length. This is necessary to iterate through seuqence positions
    SeqLength = 0
    for rec, seq in fasta_new.items():
        SeqLength = len(seq)

    # creates a dictionary with nucleotide composition per alignemt collumn {position : [nucleotides for all sequences]}
    data = {}
    for i in range(0,SeqLength): # iterate through positions
        nuc = []
        for rec, seq in fasta_new.items():
            nuc.append(seq[i]) # save nucleotides for all sequences into a list
        data[i] = nuc

    # prepare output files
    outName = fasta.split('/')[-1].split('.')[0] + '_C_' + str(nr_seqs) # name for all files
    out_freqs = open(outDir + outName + '.freq', 'w') # output to save freq file. Maybe don't do that anymore
    out_freqs.write('pos.\tA\tC\tG\tT\tN\n') # head for the freq file
    OutConSeq = open(outDir + outName + '_' + str(int(consThreashold*100)) + '.fasta', 'w') # output to save the consensu sequence

    # saves the frequency of each base in a freq file and the consensus sequences in another. A threashold should be defined for the consensus
    conSeq = '>' + outName + '\n' # consensus sequence head
    for rec, nuc in data.items():
        temp = str(rec + 1) # position in sequence
        # it assumes N as the defoult nucleotide and only replaces it if a nucleotyde has a frequency above the threashold
        base = 'N' # define nucleotide as N
        for n in ['A', 'C', 'G', 'T', 'N']:
            freq = nuc.count(n)/nr_seqs # get frequency from each nucleotide
            temp += '\t' + str(freq) # save fequency in freq file
            #check of nucleotide frequncy is above the threashold
            if freq >= float(consThreashold):
                base = n
        out_freqs.write(temp + '\n')
        conSeq += base # save nucleotide in consensus sequence
    OutConSeq.write(conSeq + '\n')

# Make consensus for all fasta files
# Input: Directory - Directory containing demultiplexed merged sequences, consThreashold - consenus threasgold in a float format, outDir - directory to save the consensus sequences
def RunConsensusAll(Directory, consThreashold, outDir):
    # just run the consensus function for all files
	for f in os.listdir(Directory):
		print('making consensus of ' + f)
		MakeConsensusPerFile(Directory + f, consThreashold, outDir)

# Join all consensus sequences from the same locus into the a single file
# Input: Directory - directory containing consensus sequences, markerList - list with all markers, OutDir - directory to save the results
def joinSamplesSameMarker(Directory, markerList, OutDir):
	for locus in markerList:
		out = open(OutDir + locus + '_together.fasta', 'w') # it makes one file per marker with the sufix together
		for file in os.listdir(Directory):
			if file.startswith(locus) and file.endswith('.fasta'): # file should start with marker name and end with the extension fasta. This way freq files are not considered
				f = open(Directory + file)
                # just write the content of the input fasta file into the output
				for l in f:
					out.write(l)
				f.close()
		out.close()

# Funtions to correct sequences and define seperate sequences with different SNPs
# get positions with Ns in a seuqence in a list format
def get_Ns_records(st):
    result = []
    c = 0 # get position
    #iterate through sequence
    for v in st:
        if v.upper() == 'N':
            result.append(c) # add position with N to list
        c += 1
    result.sort()
    return result

# iterate through the demultiplexed sequences and get frequences of linked nucleotides of the positions with Ns
def get_seq_freq(fasta, pos_pair):
    haps = [] # list to save all linked nucleotides from the Ns positions
    for name, seq in parse_fasta(fasta).items():
        # make linked nucleotide sequence (haplotype) and add it to list
        hap = ''
        for i in pos_pair:
            hap += seq[i]
        haps.append(hap)
        # use Counter() function to make a dictionary of linked nucleotides frequency: {linked nucleotides: frequncy}
    return dict(Counter(haps))

# get haplotypes (linked nucleotides) to replace the Ns with. Depending on the input it retrives one or two haplotypes.
# It chooses the two or the most frequent haplotypes. The sum frequency of those should be more the the sum frequency of the remaining (maybe use sylvia principla here or threashold)
# samples heterozygote for length can only have one haplotype, while others can have two
# in case more than the defined number of haplotypes are possible, retrive a sequence of Ns has haplotypes
# Input: dict - linked nucleotides frequency dictionary: {linked nucleotides: frequncy}, nrElments - maximum number of nucleotides that can be retrived
def get_to_el_from_dict(dict, nrElments):
    # save the haplotype counts into a list
    temp_count = []
    for hap, c in dict.items():
        temp_count.append(c)
    # sort count list
    temp_count.sort()
    # get most frequent counts
    temp_count_new = temp_count[-nrElments:]
    # get corresponding haplotypes and their relative frequency compared to the others
    result = [] # list to save included haplotypes
    incHapCount = 0 # count summ of the included haplotypes
    excHapCount = 0 # count summ of the excluded haplotypes
    for hap, c in dict.items():
        # get haplotypes corresponding to counts and save them in list
        if c in temp_count_new:
            result.append(hap)
            incHapCount += c
        else:
            excHapCount += c
    # Internal control #print('included: ' + str(incHapCount) + ' excluded: ' + str(excHapCount) + ' nr. haps: ' + str(len(result)))
    # decide if haplotypes should be considered. Include them if number of haplotypes are equal to number of elements and if their frequency is above than the remaining
    # sequences without Ns will retrieve empty result list and will go to the else condiction
    if len(result) == nrElments and incHapCount > excHapCount:
        return result
    else:
        return ['N' * len(result[0])] # save Ns has haplotypes in case the conditions are not matached
        # for sequences without Ns. since len(result[0]) == 0, the returned result will be ['']

# correct all sequences from on fasta file (marker). It also filters the data based on minimum number of sequences used to make the consensus. In this case, the complete sample will be excluded
# Input: fasta - file containing the all consensus sequences for one locus, AlleleDir - Directory containig sequences extracted based on allele length, minimum - number of reads used for consensus so that a sequence can be considered, outFasta - fasta file directory to save corrected sequences
def correctSequences(fasta, AlleleDir, minCount, outFasta):
    # get positions with Ns and how many haplotypes should be consider
    positions = {} # dictionary to save positions per sequence
    names_info = {} # dictionary to save number of alleles to be considered per sample
    for name, seq in parse_fasta(fasta).items():
        if seq != '':
            # key is a tuple with the name of the sample and the sequence, while the item the list with the postions conataing Ns
            positions[(name, seq)] = get_Ns_records(seq)
            # Check if a sample name was already process it by comparing it wiht the names_info dictionary.
            sample = name.split('_')[2]
            if sample in names_info:
                # If it alrady exists it can only have a maximum of one haplotype because it means it is heterozygote for length
                names_info[sample] = 1
            else:
                # Others are homozygote for length and can have more than one allele
                names_info[sample] = 2

    # File to save output
    out = open(outFasta, 'w')
    # List of samples marked to be excluded
    Samples_to_exclude = []
    # dictionary to save sequences to be saved
    DataOut = {}

    # Correct and save sequences
    for name_seq, position in positions.items():
        # get inpute data:
        name = name_seq[0] # sequence name
        seq = name_seq[1] # sequence
        nr_seqs = int(name.split('_')[-1]) # nr of sequences used in consensus
        sample = name.split('_')[2] # sample name
        fileName = AlleleDir + '_'.join(name.split('_')[:5]) + '.fasta' # length extracted sequences file
        het_info = names_info[sample] # number of haplotypes that can be recovered
        # get frequences of linked nucleotides (haplotypes)
        counter_dict = get_seq_freq(fileName, position)
        # get haplotypes to be used
        haps = get_to_el_from_dict(counter_dict, het_info)

        # for internal control #print(counter_dict)

        # for internal control #print(name + ' nrHapsAllowed: ' + str(het_info) + ' N positions: ' + str(position) + ' allelesFile: ' + fileName)

        # Check if haplotypes conatin Ns. If they do, marke the sample to be exluded (add to Samples_to_exclude)
        if haps[0].count('N') > 0:
            Samples_to_exclude.append(sample)
            print(sample + ' excluded for ambiguous SNPs')

        # replace the sequences with haplotype information and save them in data out. In case of two haplotypes save two sequences
        out_temp = [] # List to save corrected sequences
        # only consider sequences contructed with a read number above the defined threashold
        if nr_seqs > minCount:
            # iterate through haplotypes
            # for sequences without Ns haps == [''] and thus they are iterable.
            for hap in haps:
                temp = '' # varable to save the corect sequence
                # p is position and i is the original sequence nucleotide
                for p, i in enumerate(seq, 0):
                    # if sequence has an N repalce with corresponding nucleotide in the haplotype
                    if i == 'N' and p in position:
                        temp += hap[position.index(p)]
                    # otherwise just add original nucleotide
                    else:
                        temp += i
                # add sequence to list
                out_temp.append(temp)
            # save list of corrected sequences in dictionary
            DataOut[name] = out_temp
        # if read number used for consenus is too low mark sample to be exluded
        else:
            Samples_to_exclude.append(sample)
            print(sample + ' excluded for too little reads')

    # Save all corrected sequences for samples not marked to be excluded
    for name, outSeq in DataOut.items():
        sampleName = name.split('_')[2]
        if sampleName not in Samples_to_exclude:
            # All sequences will be save with '_0' at the end of the name. In case of heterozygotes SNPs the second sequence will be saved with '_1'
            for c, seq in enumerate(outSeq, 0):
                out.write('>' + name + '_' + str(c) + '\n' + seq + '\n')
    out.close()

# Correct consensus sequences for all fasta files
# Input: Directory - Directory containg consensus sequences file per marker, AlleleDir -  directory containing sequences extarcted based on allele length, minCount -  minimum number of reads used for consensus so that a sequence can be considered, OutDir - Directory to save output files
def correctAllSeqs(Directory, AlleleDir, minCount, OutDir):
    for file in os.listdir(Directory):
        print('correcting file ' + file)
        inFasta = Directory + file
        outFasta = OutDir + file.split('.')[0] + '_Corr.fasta'
        correctSequences(inFasta, AlleleDir, minCount, outFasta)


# Functions to call alleles based on sequence content
# get unique sequences for a given fasta file and atributes an allele number to them
def get_allele_dict(fasta_parsed):
    # data to save results
    alleles_dict = {}

    # get unique sequences from fasta
    seq_list = []
    # add all sequences to list
    for name, seq in fasta_parsed.items():
        seq_list.append(seq)
    # get unque sequences and sort them
    seq_list = list(set(seq_list))
    seq_list.sort()

    # make dict {allele: [seq]}. Allele start on 1
    for allele, seq in enumerate(seq_list, 1):
        alleles_dict[str(allele)] = [seq]
    return alleles_dict

# get alleles for all markers (deNovo)
def make_locusAlleleList_DeNovo(inDir):
    # save in dictionary: {locus: Dictionary containing alleles}
    locusAlleleList = {}
    for f in os.listdir(inDir):
        locus = f.split('_')[0] + '_' + f.split('_')[1] # get locus
        parsed_fasta = parse_fasta(inDir + f) #parse fasta
        alleles_dict = get_allele_dict(parsed_fasta) # get alleles dictionary
        locusAlleleList[locus] = alleles_dict
    return locusAlleleList

# get alleles from allele list and saves them into the alleles dictionary for all markes
def parseAlleleList(AlleleFile):
    result = {}
    f = open(AlleleFile)
    for l in f:
        led = l.rstrip('\r\n')
        # get marker information
        if len(led.split('_')) == 2: # This identify the marker name. All markers names have an underscore, thus if you slit it it will return a list of two ellements
            locus = '_'.join(led.split('_')[:2])
            alleles = {}
            seq_allele_List = []
        # get end of locus and save data
        elif led == '/' or led == '': # the end of all loci is defined by a line with \
            result[locus] = (alleles, seq_allele_List)
        # get allele and sequence information
        else:
            led = led.split('\t')
            allele = led[0].rstrip(':') # allele is saved with the format "name:"
            seq = led[1].replace('-', '') # in case loci were made from aligned data we strip gaps out
            # do not consier alleles marker for manual control. Tis means that they have an N and match to multiple alleles. These are marked with MC. This only aplies for the first versions of the script where
            if allele != 'MC':
                alleles[allele] = alleles.get(allele, []) + [seq]
                seq_allele_List.append(seq)
    f.close()
    return result

# compare sequences from new data with already existing alleles and in case new alleles are found add them to alleles dictionary
# Input: alleles_dict_per_locus - dictionary containing alleles and sequences for one marker, seq_allele_List_locus - list of sequences used to call previous alleles, fasta_parsed - new fasta data
def Complete_AlleleList_PerLocus(alleles_dict_per_locus, seq_allele_List_locus, fasta_parsed):
    # define new lists and dictionraies to save the data with preexisting data
    alleles_dict_per_locus_New = alleles_dict_per_locus
    seq_allele_List_locus_New = seq_allele_List_locus
    # compare new sequences with existing sequence list
    for head, seq in fasta_parsed.items():
        allele_Number = len(alleles_dict_per_locus) # get number of alleles already existing. This is necessary to define new alleles, which will be allele_Number + 1
        seqDeGap = seq.replace('-', '') # degap sequence in case of alignement
        # check if sequence exists
        if not seqDeGap in seq_allele_List_locus_New:
            allele_Number += 1 # make new allele name
            alleles_dict_per_locus_New[str(allele_Number)] = [seqDeGap] # add new allele and corresponding sequence to dictionary
            seq_allele_List_locus_New.append(seqDeGap) # add new sequence to sequence list
    return alleles_dict_per_locus_New

# Make dictionary with all allele information for all markers when allele list exist already
# Input: inDir - Directory containing new fasta files, AllelesInfo - parsed allele list
def Complete_AlleleList_All(inDir, AllelesInfo):
    locusAlleleList = {}
    for f in os.listdir(inDir):
        locus = f.split('_')[0] + '_' + f.split('_')[1] # get locus name
        parsed_fasta = parse_fasta(inDir + f) #parse fasta
        # go through files and and check if marker exists or not
        # non existing loci are treated as denovo data
        if locus not in AllelesInfo:
            print('locus ' + locus + ' not present in allele list') # report non existing loci
            alleles_dict = get_allele_dict(parsed_fasta) # if it does not exist alleles are defined denovo
            locusAlleleList[locus] = alleles_dict
        # if it exists comple allele inforamtion
        else:
            alleles_info_locus = AllelesInfo[locus] # get infomation for a partivular locus
            alleles_dict_locus = alleles_info_locus[0] # get alleles dictionary
            seq_allele_locus = alleles_info_locus[1] # get list of already existing seuqences for that locus
            new_alleles_dict_locus = Complete_AlleleList_PerLocus(alleles_dict_locus, seq_allele_locus, parsed_fasta) # complete allele information
            locusAlleleList[locus] = new_alleles_dict_locus
    return locusAlleleList

# saves a text file with which alleles exist per locus and which sequences correspond to those alleles.
#Format:
#locus1
#1:{tab}seq1
#2:{tab}seq2
#3{tab}seq3
#\
#locus2
#...
# Input: alleles_dict - dictionary containing alleles and sequences for one marker, locus - marker name, out - output file name
def save_allele_database(alleles_dict, locus, out):
    out.write(locus + '\n')
    for allele, allele_seq in alleles_dict.items():
        for seq in allele_seq:
            out.write(allele + ':\t' + seq + '\n')
    out.write('/\n')

# compare sequences from fasta to alleles and call them
# Input: alleles_dict_per_locus - alleles dictionary for one locus, fasta_parsed - new seuqnece data
def Define_genotypes(alleles_dict_per_locus, fasta_parsed):
    result = {}
    for head, seq in fasta_parsed.items():
        name = head.split('_')[2] # get sample name
        seqDeGap = seq.replace('-', '') # degap sequences in case of aligned files
        for allele, allele_seq in alleles_dict_per_locus.items():
            # compare with sequences
            if seqDeGap in allele_seq:
                # save information per sample
                alleles_to_dict = result.get(name, []) + [allele]
                result[name] = alleles_to_dict
    return result

# call alleles and save codominant matrix
# Input: inDir - Directory conatining corrected consensus sequences, samples_list - list of all samples, locusAlleleList - dictionary with allele information for all markers, out_pre - prefix to save output file
def CallAlleles(inDir, samples_list, locusAlleleList, out_pre):
    #opens dictionary to store genotypes
    #it already saves the header which willcontain loci info
    out_dict = {'samples' : '\t'}
    for sample in samples_list:
        out_dict[sample] = ''

    #it opens an output file to save the alleles list
    out_list = open(out_pre + '_allelle_list.txt', 'w')

    # iterates through all alignemet files and it saves the genotypes in a dictionary like this:
    # {sample: genotype_locus1\tgenotype_locus2\t ...}
    for f in os.listdir(inDir):
        locus = f.split('_')[0] + '_' + f.split('_')[1] # get locus
        out_dict['samples'] += (locus + '\t')*2 #saves the current locus in the dictionary
        alleles_dict = locusAlleleList[locus]
        parsed_fasta = parse_fasta(inDir + f) #parse fasta
        names_alleles = Define_genotypes(alleles_dict, parsed_fasta) # get alleles for each sample
        save_allele_database(alleles_dict, locus, out_list) # save new allele list
        # create genotypes: if two alleles -> allele1{tab}allele2; if one allele -> allele1{tab}allele1; if sample name not found in allele per sample list (missing data) -> 0{tab}0
        for sample in samples_list:
            # call only if sample exists. Otherwise the genotype will be 0/t0 for missing data
            if sample in names_alleles:
                alleles = names_alleles[sample]
                alleles.sort()
                # if two alleles saves it as heterozygote
                if len(alleles) == 2:
                    genotype = '\t'.join(alleles)
                # if more report the genotype as ambiguous and mark genotype for manual control. It does not aply anymore
                elif len(alleles) > 2:
                    print('locus ' + locus + ' for sample ' +  sample + '\thas ' + str(len(alleles)) + ' alleles: ' + str(alleles))
                    genotype = 'MC\tMC'
                # if one saves as homozygote
                else:
                    genotype = '\t'.join(alleles*2)
            else:
                genotype = '0\t0'
            out_dict[sample] += '\t' + genotype #it saves the genotypes in the dictionary
    out_list.close()

    #Saving outputs
    out = open(out_pre + '_matrix.txt', 'w') #it opens output file to save the matrix
     #it saves the matrix header with loci information
    out.write('samples' + out_dict['samples'] + '\n')
    #saves the matrix the genotypes
    for sample in samples_list:
        out.write(sample + out_dict[sample] + '\n')

    out.close()

# Call alleles de novo
# Input: FastaDir - Directory conatining corrected consensus sequences, samples_list - list of all samples, out_pre - prefix to save output file
def deNovoCall(FastaDir, samples_list, out_pre):
    locusAlleleList = make_locusAlleleList_DeNovo(FastaDir)
    CallAlleles(FastaDir, samples_list, locusAlleleList, out_pre)

# Call alleles with existing allele list
# Input: FastaDir - Directory conatining corrected consensus sequences, samples_list - list of all samples, out_pre - prefix to save output file, allelesFile - file with allele list
def ListBasedCall(FastaDir, samples_list, out_pre, allelesFile):
    AllelesInfo = parseAlleleList(allelesFile)
    NewLocusAlleleList = Complete_AlleleList_All(FastaDir, AllelesInfo)
    CallAlleles(FastaDir, samples_list, NewLocusAlleleList, out_pre)

# getting information from parameter file
def parseParameterFile(param):
    expected_params = ["LengthAllelels", "WorkingDir", "SeparatOut", "consThreashold", "minCount", "AlleleList", "out_pre"]
    input = {}
    with open(param) as p:
        for l in p:
            if not l.startswith("#") and l not in ["\n", "\r\n"]:
                led = l.rstrip("\r\n").split(" = ") # parameters names should be followed by a " = " to parse them
                if led[0] in expected_params:
                    input[led[0].strip(" ")] = led[-1].strip(" ")
                else:
                    print("Parameter " + led[0] + "not in conformaty")
    return input



# Main Pipeline
# reading parameter file

intputs = parseParameterFile(sys.argv[1])


# input information
allelesFile = intputs["LengthAllelels"] # matrix containing lengths
WorkingDir = intputs["WorkingDir"] + '/' # Directory containing all outputs
SeparatOut = intputs["SeparatOut"] + '/'
consThreashold = float(intputs["consThreashold"]) # consensus threashold
minCount = int(intputs["minCount"]) # minimum number of reads per consensus sequence

# allele list information already exists
AlleleList = intputs["AlleleList"]

out_pre = intputs["out_pre"] # name in common across both output files

# get file type
fileType = getFileType(SeparatOut)
# get marker and sample list
LociAndSamplesUsed = getMarkerAndSampleList(allelesFile)
Loci = LociAndSamplesUsed[0] # marker list
Samples = LociAndSamplesUsed[1] # sample list

# make output directories
for dir in ['AllelesOut', 'ConsensusOut', 'ConsensusTogether', 'Corrected', 'AlleleCall']:
    makeDir(WorkingDir + dir)

# run pipleine
# extract sequence corresponding to alleles based on length
extract_alleles(allelesFile, SeparatOut, fileType, WorkingDir + 'AllelesOut/')
# make consensus sequences
RunConsensusAll(WorkingDir + 'AllelesOut/', consThreashold, WorkingDir + 'ConsensusOut/')

# put consensus sequences together into fasta files per marker
joinSamplesSameMarker(WorkingDir + 'ConsensusOut/', Loci, WorkingDir + 'ConsensusTogether/')
# correct consensus seuences and define SNPs
correctAllSeqs(WorkingDir + 'ConsensusTogether/', WorkingDir + 'AllelesOut/', minCount, WorkingDir + 'Corrected/')
# call alleles based on sequence information
# run denovo or list based call depending on AlleleFile variable
if AlleleList == '':
    deNovoCall(WorkingDir + 'Corrected/', Samples, WorkingDir + 'AlleleCall/' + out_pre)
else:
    ListBasedCall(WorkingDir + 'Corrected/', Samples, WorkingDir + 'AlleleCall/' + out_pre, AlleleList)
