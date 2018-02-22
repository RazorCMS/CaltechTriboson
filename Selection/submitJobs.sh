#!/bin/bash

#ARGS=("$@")
# scramdir=${ARGS[0]}
#outputdir=${ARGS[1]}
# runmacro=${ARGS[2]}
#   sample=${ARGS[3]}
#    jobid=${ARGS[4]}
#  infiles=${ARGS[@]:5}

runLocal=1

outdir=/afs/cern.ch/work/j/jlawhorn/public/wwz-samples/
runmacro=SelectionFromNanoReader.C

filesPerJob=15

#if [[ ! -e ${runmacro%.*}_C.so ]] || [[ ! -e NanoReader_cc.so ]]; then
#    echo "Warning! You don't have a .so file for your macro!"
#    exit 1
#fi

while IFS=" " read -r f1 f2
do
    masterfiles=( `cat ${f2}` )

    i=1
    while [ ${#masterfiles[@]} -gt 0 ]; do
    	files=( "${masterfiles[@]:0:${filesPerJob}}" )
    	
    	echo `pwd`/runJob.sh ${CMSSW_BASE} ${outdir} ${runmacro} ${f1} ${i} ${files[@]}
	bsub -q 2nd `pwd`/runJob.sh ${CMSSW_BASE} ${outdir} ${runmacro} ${f1} ${i} ${files[@]}
    	masterfiles=( "${masterfiles[@]:${filesPerJob}}")
    	let "i=i + 1"
    	break
    done
    
done < wwzNano.conf
