# SSR-GBS pipeline

﻿The SSR-GBS pipeline contains a set of scripts that can be used to analyze SSR-GBS sequence files as described in Tibihika et al (2019) and Curto et al. (2019).
Impute files should be in fastq format where paired reads have been merged (see Curto et al. 2019 for more details). Each file should correspond to a different sample.
The way the pipeline is made it will only work in unix based operating systems such as Linux and MacOS. However, a windows compatible pipeline is being produced.

The pipeline is devided in two parts. A first one where a codominat matrix is produced with genotypes based on allele length is produced and a part where the gentypes are based on sequence composition. It is strognly advised that the length based matrix is manually quality controled using the length profiles histograms produced in the end of the first step. For more details please see Curto et al. (2019). The first part is done by running the scripts primer_demultiplex.py, CountLengths.sh, and Rscript_Markerlength_develop_Color.R. The senconf part by running Sequence_Allele_Call.py.


## References
Tibihika, P. D., Curto, M., Dornstauder-Schrammel, E., Winter, S., Alemayehu, E., Waidbacher, H., & Meimberg, H. (2019). Application of microsatellite genotyping by sequencing (SSR-GBS) to measure genetic diversity of the East African Oreochromis niloticus. Conservation Genetics, 20(2), 357-372.
Curto, M., Winter, S., Seiter, A., Schmid, L., Scheicher, K., Barthel, L. M., ... & Meimberg, H. (2019). Application of a SSR‐GBS marker system on investigation of European Hedgehog species and their hybrid zone dynamics. Ecology and evolution, 9(5), 2814-2832.


## Requirements
Unix based operating system
python 2 or 3
Biopython
R version 3.4.1 or newer, packages reshape and reshape2 (Wickham 2007) will be installed automatically by the R script number 3.


 
## primer_demultiplex.py

### Description
This script demultiplexes merged fastq files according to primer content. The outputs are one fastq file per sample and locus. This script allows for a user specified threshold of maximum number of mismatches between the primer and the reads. Only reads with a mismatch between sequence and both primers below to threshold are kept. A minimum length threshold is also necessary to be defined. 


### How to run:

python primer_demultiplex.py [1] [2] [3] [4] [5]

[1] Directory containing input fastq files. Files should be named without the underscored character (ex: Sample1.fastq)

[2] File containing primer information: this should be a tab separated text file containing one locus per line with the following information:
maker name[TAB]sequence primer forward[TAB] sequence reverse primer
Markers should be named in the following way: MarkerName_RepetitionMotif (ex: HH1_AT)

[3] Maximum number of mismatches

[4] Minimum sequence length

[5] Directory to save output files. Files names are save in the following format: RepetitionMotif_SampleName_MarkerName.fastq (ex: Sample1_HH1_AT.fastq) 




## CountLengths.sh

### Description

Per fastq file, it counts the number of occurrences of each sequence length present. It outputs this information in a space separated text file being the first column the length and the second the number of occurrences. Example:

417	14
422	18
418	276
423	282

### How to run

sh CountLengths.sh [1] [2]

[1] Directory containing input demultiplexed fastq files. Files should be as described in the output from primer_demultiplex.py: RepetitionMotif_SampleName_MarkerName.fasq (ex: Sample1_HH1_AT. fastq)

[2] Directory to save output files. Output file is saved in the following format: MarkerName_RepetitionMotif_SampleName statistics (ex: HH1_AT_Sample1_.statistics)



## Rscript_Markerlength_develop_Color.R

Requirements: R version 3.4.1 or newer, packages reshape and reshape2 (Wickham 2007) will be installed automatically by the R script.


### Description
This script uses sequence length counts to call alleles and define genotypes. A homozygote genotype is considered if the relative read count of the most abundant allele is larger than the user defined alpha value (e.g. 0.7). A heterozygote genotype is considered when the sum of the relative abundance of two most abundant read lengths exceeds the user defined alpha value 
(i) and if the difference in length of the potential alleles is larger than one time the repeat motif length. 
(ii) and if the two most abundant potential alleles only differ in one repeat size length, the less abundant allele is considered true if it is longer than the most abundant allele. 
(iii) and if the second largest potential allele is just one time smaller than the repeat motif length, it is only considered if its relative abundance equals 0.75 of the most abundant allele. 
(iv) in addition, we also considered point mutations which lead to non-integer multiples of the repeat length, if their relative abundance was at least 0.6 of the most abundant allele. 
If these criteria are not met the sample is considered for manual control. Genotypes are saved in a coma separate file and length histograms are plotted in a pdf file highlighting the called alleles with different colors: blue for the homozygous allele, green for the alleles in heterozygous samples and black for the other sequence lengths.

### How to run:

Rscript --vanilla Rscript_Markerlength_develop_Color.R [1] [2] [3] [4] [5] [6] [7] [8]


[1] Directory containing length counts. This is the output from CountLengths.sh.

[2] Minimum count of the total reads for a genotype to be considered.

[3] Path to pdf file where the length counts histograms will be saved.

[4] Alpha threshold for defining homo- and heterozygotes genotypes.

[5] Path to csv file where the length gentotype matrix file will be saved

[6] Path to csv file where sample and marker and allele information will be saved

[7] Minimum allele size length for plotting the results

[8] Maximum allele size length for plotting the results

### References
Wickham, H. (2007). Reshaping Data with the reshape Package. Journal of Statistical Software, 21(12), 1-20. URL http://www.jstatsoft.org/v21/i12/.

### Outputs
Two output files are produced. A csv file contaning containg the length alleles. Alleles correspond amplicon length. If there no reads a genotype with "empty" will be save, if the read count for that marker sample combination is too low "too little reads" will be saved. In case the gentype does not follow the the critera defined in Curto et al. 2019 is marked for manual control by adding "_man check". 




## Sequence_Allele_Call.py

### Description
This scripts replaces the script 4 to 7 from Tibihika et al. 2028. It takes alleles saved in the csv file produced by Rscript_Markerlength_develop_Color.R and extract the corresponding reads, makes consensus sequence for each length-based allele, detects possible SNP variation and calls alleles based on sequence information. For each of these steps a new directory is created to save output files. 

### How to run

python extract_alleles_of_a_certain_length_v2.py [1]


[1] parameter text file in the following format with the opetions “LengthAllelels = “, “ WorkingDir = “, “SeparatOut = “, “consThreashold = “, “minCount = “, “AlleleList = “, “out_pre = “, followed by used defined input.

Example:


LengthAllelels = OutDir_20200520/MarkerPlots/markermatrix.csv

WorkingDir = OutDir_20200520/

SeparatOut = OutDir_20200520/SeparatOut/

consThreashold = 0.7

minCount = 10

AlleleList = 

out_pre = Test22-06-2020


##### Description of parameter file options:

The SSR-GBS_pipeline_paramter_file.txt is an example of a parameter file

“LengthAllelels” is a coma separated text file containing length genotypes. This file should contain a header containing the marker information. The next lines it contains the genotype information. Th first column should contain the sample names while the remaining the genotype information. Two columns per markers should be added allowing for heterozygote genotypes. Missing data should be coded as 0.

Example:
samplename	HH1_TA	HH1_TA	HH2_TGT	HH2_TGT	HH3_AAAC	HH4_AAAC
Sample1	425	425	429	426	418	418
Sample2	423	425	429	429	0	0
Sample3	0	0	429	429	414	414
Sample4	427	427	426	426	418	418

If the matrix is produced by Rscript_Markerlength_develop_Color.R other characters besides allele calls such as “empty” when not read is found or “_man check” when the genotype does not follow the expectation. The script can parse the file even with these characters are read the genotypes.


“WorkingDir” is the directory where the output directories and consequently all output files will be saved.

“SeparatOut” is the directory containing the demultiplexed fastq files produced by the script primer_demultiplex.py.

“consThreashold” threshold for consensus sequence. Nucleotides with a frequency per position above this will be saved otherwise an “N” will be outputted. This should be float. For example, 0.7.

“minCount” is the minimum number of reads necessary for an allele to be called. If a consensus seuqnece was made by a number of reads below this value the genotype will be considered as missing.

“AlleleList” is a file containing alleles and corresponding sequences to be used as reference for the allele call. If there are already genotypes defined based on sequence information for other samples and allele list can be used so that the new samples are directly comparable with the older ones. If that is not the case, just leave this option empty. The allele list should be in the following format:

Marker name
1: Sequence1
2: Sequence2
\
Where 1 and 2 are allele corresponding to Sequence1 and Sequence2, respectively. “\” is used to separate markers.


For example

AF10_TGCC

1:	CTCCCCATCGACGGTAACGCTCTCTCCTGCCCCT…

2:	CTCCCCATCGACGGTAACGCTCTCTCCTGCCCCT…

3:	CTCCCCATCGACGGTAACGCTCTCTCCTGCCCCT…

\


AF12_TTTC

1:	GGCCTATTGTGTTCGAAATTATGCAGGCTCACCGAAGCTTCTCGCTTATCGTTGT…

2:	GGCCTATTGTGTTCGAAATTATGCAGGCTCACCGAAGCTTCTCGCTTATCGTTGTG…

\


“out_pre” is a prefix to save the final output files.


### Output files:

There are four directories created to save output files:

“AllelesOut”: fasta files with sequences corresponding to each allele defined in the amplicon length call. The allele information is added to the file name. For example, AF1_AGTG_EO-415_Al_356.fasta, correspond to all sequences in the sample EO-415 and marker AF1_AGTG that have a length of 356 bp.

“ConsensusOut”: Files containing the consensus sequence for each length allele. Two files per consensus sequence are produced. A fasta file with sequence, and a .freq file with the nucleotide frequency per position. The latter is a tab separated file with six columns: First the sequence position; Second to fifth columns the relative frequency of each nucleotide; Sixth column the relative frequency of uncalled bases (“N”). 

Here it is an example

pos.	A	C	G	T	N

1	1.0	0.0	0.0	0.0	0.0

2	0.0	0.0	1.0	0.0	0.0

3	0.0	0.0	0.0	1.0	0.0

4	0.0	1.0	0.0	0.0	0.0

5	0.0	0.0	0.0	1.0	0.0

6	1.0	0.0	0.0	0.0	0.0

7	0.0	0.0	1.0	0.0	0.0

8	0.0	0.0	0.0	1.0	0.0

9	0.0	0.0	0.0	1.0	0.0

10	0.0	0.0	0.0	1.0	0.0


The number of reads used to make the consensus sequence and the consensus threshold is added in the file name. For example: AF2_GGAA_EO-616_Al_352_C_1917_70.fasta corresponds to a 70% consensus sequence produced from 1917 reads.

“ConsensusTogether”: Here one fasta file per marker is saved containing all consensus sequences for all sample for a particular marker. 

“Corrected”: Here the corrected files with the final sequence information are saved. In case of SNPs two sequences per sample are found with the corrected variant. 


“AlleleCall”: For the allele call two files are saved. 1) The *_matrix.txt, which is a tab separated text file corresponding to the codominant matrix. The header has the marker information while the remaining lines the genotypes. The first column corresponds to the sample information. Two columns per marker are saved. 2)The *_allelle_list.txt file, which contains the information of which sequence each allele corresponds to. It is saved in the following format 

