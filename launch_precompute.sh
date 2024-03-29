#!/bin/bash

launch_all=0


#We set the output file in the repertory PRINT
mkdir -p PRINT
DIR=$(pwd)
DIR_PRINT="$DIR/PRINT"

#File in which the output will be redirected if Non Interactive OR Restart
#(Trick use:  tail -f print_no_interactive.txt  in order to follow the evolution of the program)
OUTFILE="$DIR_PRINT/print_no_interactive$3.txt"


python precompute.py $1 $launch_all $2 #> $OUTFILE 2>&1

#rm $OUTFILE

exit 0

