#----------------------------------------------------
# script to merge and clean up output from lxbtch submissions
#
# example command
# ./mergeSelFiles.sh HHToGGBB_14TeV /afs/cern.ch/work/a/ariostas/private/HHToGGBB/ /afs/cern.ch/work/a/ariostas/private/HHToGGBB/
#
# Jay Lawhorn 11/4/13
#----------------------------------------------------
#!/bin/bash

  sample_type=$1
file_location=$2
 out_location=$3

count=0

for file in `ls ${file_location}/${sample_type}*`
do
  file_array[$count]=$file
  let "count+=1"
done

hadd -f ${out_location}/${sample_type}_temp.root ${file_array[@]}

root -l -b -q "cleanUpMergedFiles.C+(\"${out_location}/${sample_type}_temp.root\",\"${out_location}/${sample_type}.root\")"

#rm ${out_location}/${sample_type}_temp.root
