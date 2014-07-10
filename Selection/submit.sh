#!/bin/bash
#------------------------------------------------------------
# Submit all seletion jobs to lxbatch.
#
# To run use the following command: ./submit.sh xsec.txt
#------------------------------------------------------------

root_script=$1
  conf_file=$2

rm *_14TEV.txt
rm -r /afs/cern.ch/work/a/ariostas/private/bbgg_temp/*
rm -r LSFJOB_*

mkdir -v /afs/cern.ch/work/a/ariostas/private/bbgg_temp/HHToGGBB_14TeV
bsub -q 1nd -W 1440 -J Signal run.sh ${root_script} HHToGGBB_14TeV 0.000089

while read line #loop over lines in ${conf_file}
do
  array=($line)
    if [ "${array[0]}" != "#" ]; then 
	# get list of files in eos for that sample+configuration combination

	      mkdir -v /afs/cern.ch/work/a/ariostas/private/bbgg_temp/${array[0]}
	      /afs/cern.ch/project/eos/installation/0.2.41/bin/eos.select ls ${conf}/${array[0]}/ | grep root > "${array[0]}.txt"
	      bsub -q 1nd -W 1440 -J ${array[0]} run.sh ${root_script} ${array[0]} ${array[1]} # submit to lxbtch

    fi
done < ${conf_file}

