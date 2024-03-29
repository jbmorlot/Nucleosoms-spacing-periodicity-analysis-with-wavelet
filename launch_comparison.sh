#!/bin/bash


#Parameters if Not Interactive
data_file_1=$2
data_file_2=$3
name_data_1=$4
name_data_2=$5
data_min=0
data_max=0
save_txt_var=1
local_analysis=0
histogram_analysis=1
reverse=0 

#We set the output file in the repertory PRINT
mkdir -p PRINT
DIR=$(pwd)
DIR_PRINT="$DIR/PRINT"

#File in which the output will be redirected if Non Interactive OR Restart
#(Trick use:  tail -f print_no_interactive.txt  in order to follow the evolution of the program)
OUTFILE="$DIR_PRINT/print_no_interactive$4-$5.txt"

#Check the option Interactive (-I) or Not Interactive (-N)
while getopts ":I:N:" opt; do
  case $opt in
    I)
      python Nucleosom_local_comparison.py
      ;;
      
    N)
      python Nucleosom_local_comparison_no_interactive.py $data_file_1 $data_file_2 $name_data_1 $name_data_2 $data_min $data_max $save_txt_var $local_analysis $histogram_analysis $reverse > $OUTFILE 2>&1
      ;;
      
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

#rm $OUTFILE

exit 0

