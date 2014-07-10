#!/bin/bash
#------------------------------------------------------------
#
# Cancel jobs submitted to lxbatch
#
# To use this run "bjobs" on lxplus. Then copy the output and paste
# it on cancel.txt . Finally, run this with ./cancel.sh cancel.txt
#
#------------------------------------------------------------

conf_file=$1

while read line #loop over lines in ${conf_file}
do
  array=($line)
    if [ "${array[0]}" != "#" ]; then 
	  
	      bkill ${array[0]} # kill job

    fi
done < ${conf_file}
