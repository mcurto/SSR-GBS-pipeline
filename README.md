Script number 1
File name:  primer_demultiplex.py
Requirements: python 2 or 3; Biopython
Description: This script demultiplexes merged fastq files according to primer content. The outputs are one fastq file per sample and locus. This script allows for a user specified maximum number of mismatches between the primer and the reads. Only reads with a mismatch to both primers below to the defined are kept. In this case they are saved in a separate file. Moreover, sequences below a certain length can be excluded. 
How to run: python extract_reads_correct_primer_merged.py [1] [2] [3] [4] [5]
[1] Directory containing input fastq files. Files should be named without the underscored character (ex: Sample1.fastq)
[2] File containing primer information: this should be a tab separated text file containing one locus per line with the following information:
maker name[TAB]sequence primer forward[TAB] sequence reverse primer
Markers should be named in the following way: MarkerName_RepetitionMotif (ex: HH1_AT)
[3] Maximum number of mismatches
[4] Minimum sequence length
[5] Directory to save output files. Files names are save in the following format: RepetitionMotif_SampleName_MarkerName.fasq (ex: Sample1_HH1_AT.fastq) 


Script number 2
File name: CountLengths.sh
Requirements: Unix system
Description: Per fastq file it counts the number of occurrences of each sequence length present. It outputs this information in a space separated text file being the first column the length and the second the number of occurrences. Example:
417	14
422	18
418	276
423	282
How to run: sh CountLengths.sh [1] [2]
[1] Directory containing input demultiplexed fastq files. Files should be as described in the output from script 1: RepetitionMotif_SampleName_MarkerName.fasq (ex: Sample1_HH1_AT. fastq)
[2] Directory to save output files. Output file is saved in the following format: MarkerName_RepetitionMotif_SampleName statistics (ex: HH1_AT_Sample1_.statistics)

Script number 3                                                                                                                                            
File name: Rscript_Markerlength_develop_Color.R
Requirements: R version 3.4.1 or newer, packages reshape and reshape2 (Wickham 2007) will be installed automatically by the R script.
Description: This script uses sequence length counts to call alleles and define genotypes. A homozygote genotype is considered if the relative read count of the most abundant allele is larger than the user defined alpha value (e.g. 0.7). A heterozygote genotype is considered when the sum of the relative abundance of two most abundant read lengths exceeds the user defined alpha value 
(i) and if the difference in length of the potential alleles is larger than one time the repeat motif length. 
(ii) and if the two most abundant potential alleles only differ in one repeat size length, the less abundant allele is considered true if it is longer than the most abundant allele. 
(iii) and if the second largest potential allele is just one time smaller than the repeat motif length, it is only considered if its relative abundance equals 0.75 of the most abundant allele. 
(iv) in addition, we also considered point mutations which lead to non-integer multiples of the repeat length, if their relative abundance was at least 0.6 of the most abundant allele. 
If these criteria are not met the sample is considered for manual control. Genotypes are saved in a coma separate file and length histograms are plotted in a pdf file highlighting the called alleles with different colors: blue for the homozygous allele, green for the alleles in heterozygous samples and black for the other sequence lengths.
How to run: Rscript --vanilla Rscript_Markerlength_develop_Color.R [1] [2] [3] [4] [5] [6] [7] [8]
[1] Directory containing length counts. This is the output from script 2.
[2] Minimum count of the total reads for a genotype to be considered.
[3] Path to directory where pdf file with length counts histograms will be saved.
[4] Alpha threshold for defining homo- and heterozygotes genotypes.
[5] Path to Csv matrix file containing genotypes of markers of the respective samples
[6] Path to Csv file containing three columns with sample and marker name and respective allele
[7] Minimum allele size length for plotting the results
[8] Maximum allele size length for plotting the results

References
Wickham, H. (2007). Reshaping Data with the reshape Package. Journal of Statistical Software, 21(12), 1-20. URL http://www.jstatsoft.org/v21/i12/.


Script number 4
File name: extract_alleles.py
Requirements: python 2 or 3; Biopython.
Description: It extracts all sequences with the same length of the alleles saved in the csv file produced by script number 2 and saves them in a fasta file per sample and allele.
How to run: python extract_alleles_of_a_certain_length_v2.py [1] [2] [3]
[1] Coma separated text file containing genotypes. This file should contain a header containing the marker information. The next lines it contains the genotype information. Th first column should contain the sample names while the remaining the genotype information. Two columns per markers should be added allowing for heterozygote genotypes. Missing data should be coded as 0. Example:
samplename	HH1_TA	HH1_TA	HH2_TGT	HH2_TGT	HH3_AAAC	HH4_AAAC
Sample1	425	425	429	426	418	418
Sample2	423	425	429	429	0	0
Sample3	0	0	429	429	414	414
Sample4	427	427	426	426	418	418
[2] Directory containing input fastq files. These are the output from script 1.
[3] Directory to save output fasta files. These files are named in the following format: MarkerName_SampleName_Al_SequenceLength.fasta (ex: HH1_AT_Sample1_Al_425.fasta)


Script number 5
File name: get_consensus_and_freq.py
Requirements: python 2 or 3.
Description: It produces a consensus sequence per file keeping bases above a certain similarity threshold. For positions where this is not meet the script outputs a “N”.
How to run: python get_consensus_and_freq.py [1] [2] [3]
[1] Directory containing input fasta extracted based on length genotype information (output from script 4)
[2] Directory to save the consensus files. Files will be named in the following format: Marker Name_Sample Name_Al_Sequence Length_C Numer Of Sequences Used_Consensus Threashold.fasta (ex: HH1_AT_Sample1_Al_425_C881_70. fasta)
[3] Similarity of frequency threshold in an integer form. The value of 0.7 will do a 70% consensus.

Script 6
File name: correct_allele_sequence.py
Requirements: python 2 or 3
Description: In case a sequence has an ambiguous base (“N”) after the consensus and it is homozygote based on sequence length (SL), it divides the sequence into two new ones. The Ns are corrected based on the frequency that the N position shows up in the reads extracted from script 4 taking the two most frequent nucleotide combinations. In case of heterozygote genotype based on SL the sequence is not divided but only corrected with the most frequent nucleotide information.
How to run: python correct_allele_sequence.py [1] [2] [3] [4]
[1] fasta file containing all consensus sequences from one marker. File name should be MarkerName.fasta (ex: HH1_AT. fasta). Sequences should be named in the following format: MarkerName_SampleName_Al_SequenceLength_CNumerOfSequencesUsed_ConsensusThreashold (ex: >HH1_AT_Sample1_Al_425_C881_70)
[2] Directory containing sequences extracted based on length genotypes (output from script 4)
[3] Minimum number of sequence counts required for an allele to be condidered
[4] Name of the output fasta file


Script number 7
File name: call_alleles_from_fasta.py
Requirements: python 2 or 3.
Description: Uses the haplotypes obtained from the SNP correction process and converts them into allele’s numbers. If a haplotype can be assigned to more than one allele it is saved as missing data. The results are saved in tab separated text file in the format of a codominant matrix (*matrix.txt). Allele’s numbers and which haplotypes they correspond to are saved in file ending with *allelle_list.txt with the following format:
Marker 1
Allele 1:	Haplotype
Allele 2:	Haplotype
…
Marker 2
…
How to run: python call_alleles_from_fasta.py [1] [2] [3] [4] [5]
[1] Directory containing haplotypes per locus in fasta format (output from script 7)
[2] Prefix common to all input files (ex: HH in marker HH1)
[3] List of samples names to be considered. This should be a text file with one sample name per line)
[4] Prefix that should be used to save output files
[5] Minimum number of sequences required to for a haplotype to be consider