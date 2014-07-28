#!/bin/bash
#------------------------------------------------------------
# Submit a batch of jobs to lxbtch
#
# example command:
#             root_script     conf_file output_location
# ./submit.sh selectDelphes.C xsec.txt  /afs/cern.ch/work/a/ariostas/private/HHToGGBB/
# 
# conf_file has format (no leading "#")
# sample_type cross_section
#
# Note: conf file needs empty last line... clunky I know.
# Jay Lawhorn 11/4/13
#------------------------------------------------------------

  conf_file=$1

./mergeSelFiles.sh HHToGGBB_14TeV /afs/cern.ch/work/a/ariostas/private/bbgg_temp/HHToGGBB_14TeV/ /afs/cern.ch/work/a/ariostas/private/bbgg/ 

while read line #loop over lines in ${conf_file}
do
  array=($line)
  #for conf in PhaseI/Configuration0 PhaseII/Configuration3 PhaseII/Configuration4v2 # to process all three configurations
  
   if [ "${array[0]}" != "#" ]; then 

	   ./mergeSelFiles.sh ${array[0]} /afs/cern.ch/work/a/ariostas/private/bbgg_temp/${array[0]}/ /afs/cern.ch/work/a/ariostas/private/bbgg/   

    fi

done < ${conf_file}
