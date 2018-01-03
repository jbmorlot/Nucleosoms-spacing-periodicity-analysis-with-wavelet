#!/bin/bash

data_min=10000
data_max=1000000
gauss=0

python plot_histogram.py $1 $data_min $data_max $gauss

exit 0
