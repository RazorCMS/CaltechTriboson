#!/bin/bash
ARGS=("$@")
 scramdir=${ARGS[0]}
outputdir=${ARGS[1]}
 runmacro=${ARGS[2]}
   sample=${ARGS[3]}
    jobid=${ARGS[4]}
  infiles=${ARGS[@]:5}

workdir=`pwd`
echo `hostname`
echo "args: " ARGS

cd ${scramdir}
eval `scramv1 runtime -sh`
cd $workdir

export XRD_NETWORKSTACK=IPv4
export X509_USER_PROXY=/afs/cern.ch/user/j/jlawhorn/x509up_u69636

cp ${scramdir}/src/CaltechTriboson/Selection/rootlogon.C .
cp ${scramdir}/src/CaltechTriboson/Selection/NanoReader.* .
#cp ${scramdir}/src/CaltechTriboson/Selection/*.so .
#cp ${scramdir}/src/CaltechTriboson/Selection/*.pcm .
cp ${scramdir}/src/CaltechTriboson/Selection/${runmacro} .

for file in ${infiles[@]}
do
    echo root -l -b -q ${runmacro}+\(\"${file}\",\"${file##*/}\"\)
    root -l -b -q ${runmacro}+\(\"${file}\",\"${file##*/}\"\)
done

hadd ${sample}_${jobid}.root *.root

mv ${sample}_${jobid}.root ${outputdir}

status=`echo $?`
echo "Status - $status"

exit $status