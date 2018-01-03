# Nucleosoms_spacing_periodicity_analysis_with_wavelet
Nucleosomes are proteins, the histones, around which the DNA is enrolled. This complex is present every 200 bp in human cells
in average and it plays an important role in the control of gene expression [1].
The local period of the nucleosomes spacing can gives access to different layers of chromatine organisation and therefore genic regulation.

## Purpose
This code aim to produce study the local period of nucleosome spacing and the width of the peaks.

## How to use it
### For interactive Start ( with a window to enter data ): --> Recommended version
	./launch.sh -I ""
Then fill the fields

### For Non Interactive Start (for clusters for exemple):
Open launch.sh and change directly the parameters on the file, then type in a terminal
	./launch.sh -N "/path_to_the_file/file.wig" Number_or_string
with Number_or_string a number or a string which will be inserted in the output file name print_no_interaction(Number_or_string).txt
(allow you to follow multiple parallel launch of the program) 
if you want to check the progress of the work, type in another terminal (in repertory/PRINT): 
	tail -f print_no_interactive.txt
  
 ## Restore previous calculation:
 A restart implementation is provided if the program crashes. You can restart the program at the last 
 automatic checkpoint with:

	./launch.sh -R "/path_to_the_file/file.wig"
OR
	./launch.sh -I ""
Then select your previous data file and Click on "Restart"

In both case, the parameters of the previous simulation will be loaded.


## Exemple 1

Launch the interactive mode: 
  ./launch.sh -I "" 
and select the exemple file in ./exemple/GSM724105_AG222_min100_max200_PE_5milFMR.wig

with

launch_all=0\
data_min=0\
data_max=0\
period_min=50\
period_max=300\
width_min=5\
width_max=200\
plot=0\
save=1\
save_txt=0\
local_var=0\
global_var=1\

## Exemple 2

Launch the interactive mode: ./launch.sh -I "" and select the exemple file in ./exemple/GSM724105_AG222_min100_max200_PE_5milFMR.wig

with

  launch_all=0\
  data_min=109500\
  data_max=110500\
  period_min=50\
  period_max=300\
  width_min=5\
  width_max=200\
  bin_period=5\
  plot=0\
  save=1\
  save_txt=0\
  local_var=1\
  global_var=0\



![alt text](http://url/to/img.png)

[1] Audit, Benjamin, et al. "Long-range correlations in genomic DNA: a signature of the nucleosomal structure." Physical Review Letters 86.11 (2001): 2471.





