#!/bin/bash


#Parameters if Not Interactive
launch_all=0
data_min=109500
data_max=110500
period_min=50
period_max=300
width_min=5
width_max=200
bin_period=5
plot=0
save=1
save_txt=0
local_var=1
global_var=0

#We set the output file in the repertory PRINT
mkdir -p PRINT
DIR=$(pwd)
DIR_PRINT="$DIR/PRINT"

#File in which the output will be redirected if Non Interactive OR Restart
#(Trick use:  tail -f print_no_interactive.txt  in order to follow the evolution of the program)
OUTFILE="$DIR_PRINT/print_no_interactive$3.txt"

#Check the option Interactive (-I) or Not Interactive (-N)
while getopts ":I:N:R:" opt; do
  case $opt in
    I)
      python Nucleosom_position_analysis.py
      ;;
      
    N)
      python Nucleosom_position_analysis_no_interactive.py $2 $launch_all $data_min $data_max $period_min $period_max $width_min $width_max $bin_period $plot $save $save_txt $local_var $global_var 0 > $OUTFILE 2>&1
      ;;
      
    R)
      python Nucleosom_position_analysis_no_interactive.py $2 $launch_all $data_min $data_max $period_min $period_max $width_min $width_max $bin_period $plot $save $save_txt $local_var $global_var 1 > $OUTFILE 2>&1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

#rm $OUTFILE

exit 0

