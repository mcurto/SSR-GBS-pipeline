
Inputdir=$1;
Outputdir=$2

for filename in $(ls $Outputdir/SeparatOut); do
    awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' $Outputdir/SeparatOut/$filename | sort -k2 -n  > $Outputdir/MarkerStatistics/$filename.statistics
    newname=$(echo $filename | awk 'BEGIN { FS = "_" }; { print $2"_"$3"_"$1"_" }' | sed 's/\.fastq//g' )      #FASTQ files need here fa
    echo $newname
    mv $Outputdir/MarkerStatistics/$filename.statistics $Outputdir/MarkerStatistics/$newname.statistics
done
