#!/bin/bash

runLocal=1

outdir=/afs/cern.ch/work/j/jlawhorn/public/wwz-samples/

while IFS=" " read -r f1 f2
do
    echo $f1
    for file in `ls "$f2"*root`
    do
	outfile=${outdir}${f1}_${file##*_}
	if [ $runLocal = 1 ]; then 
	    echo root -l -b -q Selection.C+\(\"${file}\",\"${outfile}\"\)
	    root -l -b -q Selection.C+\(\"${file}\",\"${outfile}\"\)
	fi
    done
done < wwz.conf
