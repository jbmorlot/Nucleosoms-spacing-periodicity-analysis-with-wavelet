#!/bin/python

import numpy as np
from math import *
from sys import stdout
import sys

import os
import re
import gzip
import linecache


def input_file_analysis(data_file,data_min,data_max):
  
  regexp = re.compile(r'[A-Za-z]')
  regexpscientific = re.compile(r'e\+')
  regexpgz = re.compile(r'.gz')
  regexpheader = re.compile(r'(data_length,data_min,data_max)')
  
  
  with open(data_file) as f:
    first_line = f.readline()
  
  if regexpheader.search(first_line) is not None:
    print('Length of file --> In the header')
   
    data_input = (first_line.split(" ")[2]).split(",")
    size = int(float(data_input[0]))
      
  else:
    if regexpgz.search(data_file) is not None:
      i = 0
      with gzip.open(data_file) as f:
	for i, l in enumerate(f):
	  pass
	
	size = i +1
	  
    else:
      with open(data_file) as f:
	for i, l in enumerate(f):
	  pass
	
	size = i + 1
      
  #We set the parameters from data_file
  data_max_default = size
  data_min_default = 0
  
  #Default value
  if data_max == 0:
    data_max  = data_max_default
    
  if data_min == 0:
    data_min = data_min_default
  
  #Verifying if the input data_min/data_max from user are out of range
  if data_max > size or data_min < 0:
    sys.exit('WARNING: Invalid input limits: data_max > {0} or data_min < {1}'.format(int(size),0))
    
  if size == 0:
    sys.exit('WARNING: Invalid Input File: Size of the histogram = 0')
  
  data_length = data_max - data_min
    
  return data_length, data_min, data_max


def input_file_comparison(data_file_1,data_file_2,data_min,data_max):
  
  regexp = re.compile(r'[A-Za-z]')
  regexpscientific = re.compile(r'e\+')
  regexpgz = re.compile(r'.gz')
  regexpheader = re.compile(r'(data_length,data_min,data_max)')
  
  #file 1
  min_data_file_1 = 0
  max_data_file_1 = 0
  data_min_1 = 0
  data_max_1 = 0
  
  with open(data_file_1) as f:
    first_line = f.readline()
  
  if regexpheader.search(first_line) is not None:
    print('Length of file 1 --> In the header')
   
    data_input_1 = (first_line.split(" ")[2]).split(",")
    size_1 = int(float(data_input_1[0]))
    min_data_file_1 = int(float(data_input_1[1]))
    max_data_file_1 = int(float(data_input_1[2]))
      
  else:
    if regexpgz.search(data_file_1) is not None:
      i = 0
      with gzip.open(data_file_1) as f:
	for i, l in enumerate(f):
	  pass
	
	size_1 = i +1
	  
    else:
      with open(data_file_1) as f:
	for i, l in enumerate(f):
	  pass
	
	size_1 = i + 1
  
  #We set the parameters from data_file
  data_max_1_default = size_1
  data_min_1_default = 0
  
  
  #Default value
  if data_max == 0:
    data_max_1 = data_max_1_default
  else: 
    data_max_1 = data_max
    
  if data_min == 0:
    data_min_1 = data_min_1_default
  else:
    data_min_1 = data_min
    
    
  #Verifying if the input data_min/data_max from user are out of range
  if data_max_1 > size_1 or data_min_1 < 0:
    sys.exit('WARNING: Invalid input limits: data_max > {0} or data_min < {1}'.format(int(size_1),0))
    
  if size_1 == 0:
    sys.exit('WARNING: Invalid Input File: Size of the histogram = 0')
  
  #file 2
  min_data_file_2 = 0
  max_data_file_2 = 0
  data_min_2 = 0
  data_max_2 = 0
  
  with open(data_file_2) as f:
    first_line = f.readline()
    
  if regexpheader.search(first_line) is not None:
    print('Length of file 2 --> In the header')
    data_input_2 = (first_line.split(" ")[2]).split(",")
    size_2 = int(float(data_input_2[0]))
    min_data_file_2 = int(float(data_input_2[1]))
    max_data_file_2 = int(float(data_input_2[2]))
    
  else:
    if regexpgz.search(data_file_2) is not None:
      i = 0
      with gzip.open(data_file_2) as f:
	for i, l in enumerate(f):
	  pass
	
	size_2 = i +1
	  
    else:
      with open(data_file_2) as f:
	for i, l in enumerate(f):
	  pass
	
	size_2 = i + 1
    
  #We set the parameters from data_file
  data_max_2_default = size_2
  data_min_2_default = 0
  
  #Default value
  if data_max == 0:
    data_max_2 = data_max_2_default
  else: 
    data_max_2 = data_max
  
  if data_min == 0:
    data_min_2 = data_min_2_default
  else:
    data_min_2 = data_min
    
  #Verifying if the input data_min/data_max from user are out of range
  if data_max_2 > size_2 or data_min_2 < 0:
    sys.exit('WARNING: Invalid input limits: data_max > {0} or data_min < {1}'.format(int(size_2),0))
    
  if size_2 == 0:
    sys.exit('WARNING: Invalid Input File: Size of the histogram = 0')
    
  
  #We adjust the min and the max of the 2 histograms in order to match
  print data_min_1,data_min_2
  print data_max_1,data_max_2
  print(min_data_file_1,min_data_file_2)
  print(max_data_file_1,max_data_file_2)
  
  min_data_file_diff = min_data_file_1 - min_data_file_2
  max_data_file_diff = max_data_file_1 - max_data_file_2
  
  print "Difference between the 2 data: ( delta min , delta_max ) = ( {0} , {1} )".format(min_data_file_diff,max_data_file_diff)
  
  if min_data_file_diff < 0:
    data_min_1 +=  abs(min_data_file_diff)
    data_min_2 += 0
    
    #We shift all the histogram in order to keep the same length
    if data_max_1 != data_max_1_default :
      data_max_1 += abs(min_data_file_diff)
      
      #we make sure that data_max_1 > data_max_1_default
      if data_max_1 > data_max_1_default:
	data_max_1 = data_max_1_default
	
  
  if min_data_file_diff > 0:
    data_min_1 += 0
    data_min_2 += abs(min_data_file_diff)
    
    #We shift all the histogram in order to keep the same length
    if data_max_2 != data_max_2_default :
      data_max_2 += abs(min_data_file_diff)
      
      #we make sure that data_max_2 > data_max_2_default
      if data_max_2 > data_max_2_default:
	data_max_2 = data_max_2_default
  
  if max_data_file_diff < 0 and data_max_1 == data_max_1_default :
    data_max_1 += 0 
    data_max_2 -= abs(max_data_file_diff)
  
  if max_data_file_diff > 0 and data_max_2 == data_max_2_default :
    data_max_1 -= abs(max_data_file_diff)
    data_max_2 += 0
   
  
  #data_min = max(data_min_1,data_min_2)
  #data_max = min(data_max_1,data_max_2)
  data_length_1 = data_max_1 - data_min_1
  data_length_2 = data_max_2 - data_min_2
  print (data_length_1,data_length_2)
    
  return data_length_1, data_min_1, data_min_2, data_max_1, data_max_2


def load_data(data_file,data_file_min,data_file_max):
 
  #Loading data_file
  print("Loading data file...")

  regexp = re.compile(r'[A-Za-z]')
  regexpscientific = re.compile(r'e\+')

  regexpgz = re.compile(r'.gz')
  
  data_input_1 = []

  #Au cas ou le fichier est compresse
  if regexpgz.search(data_file) is not None:
    with gzip.open(data_file) as f:
      for i,line in enumerate(f):
	if i < data_file_min:
	  pass
	if data_file_min <= i and i < data_file_max:
	  if regexp.search(line) is None:
	    data_input_1.append(float(line.replace("\n", "")))
	  if regexpscientific.search(line) is not None:
	    data_input_1.append(float(line.replace("\n", "")))
	if i > data_file_max:
	  break	
  else:
    with open(data_file) as f:
      for i,line in enumerate(f):
	if i < data_file_min:
	  pass
	if data_file_min <= i and i < data_file_max:
	  if regexp.search(line) is None:
	    data_input_1.append(float(line.replace("\n", "")))
	  if regexpscientific.search(line) is not None:
	    data_input_1.append(float(line.replace("\n", "")))
	if i > data_file_max:
	  break
  
  data_input = np.asarray(data_input_1,np.float32)
  print("DONE!\n")
  
  histogram = data_input
  histogram_length = histogram.shape[0]
   
  return histogram, histogram_length