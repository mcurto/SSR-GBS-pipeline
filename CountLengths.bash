
Inputdir=$1;
Outputdir=$2

for filename in $(ls $Inputdir); do
    awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' $Inputdir/$filename | sort -k2 -n  > $Outputdir/$filename.statistics
    newname=$(echo $filename | awk 'BEGIN { FS = "_" }; { print $2"_"$3"_"$1"_" }' | sed 's/\.fastq//g' )      #FASTQ files need here fa
    echo $newname
    mv $Outputdir/$filename.statistics $Outputdir/$newname.statistics
done
