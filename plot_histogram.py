#!/bin/python

import scipy
from scipy import signal
import numpy as np
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from input_file_analysis import *


def plot_histogram(data_file,data_min,data_max,gauss):
  
  #Repertoire dans lequel sont les donnes
  directory = os.path.dirname(data_file)
  #Creating new directory and executing the program in it if saving
  filename = os.path.split(data_file)[1]
  filename_split = filename.split('.')[0]
  sub_directory = directory +'/'+ filename_split
  
  if not os.path.exists(sub_directory):
    os.mkdir(sub_directory)
  #Changing Current directory
  os.chdir(sub_directory)
  
  histogram, histogram_length = load_data(data_file,data_min,data_max)
  
  #histogram_sum = np.sum(histogram)
  #histogram = histogram/histogram_sum
  
  os.chdir(sub_directory)
  print ("Ploting Histogram...")
  fig1 = plt.figure()
  ax11 = fig1.add_subplot(111)
  ax11.plot(np.linspace(data_min,data_max, histogram_length), histogram)
  
  if gauss == 1:
    #Signal convolved with gaussian
    sigma = 30
    
    histogram_gaussian = scipy.signal.fftconvolve(histogram,signal.gaussian(histogram_length, std = sigma),"same")
    
    histogram_gaussian_sum = np.sum(histogram_gaussian)
    histogram_gaussian = histogram_gaussian/histogram_gaussian_sum
    ax11.plot(np.linspace(data_min,data_max, histogram_length), histogram_gaussian,'r')
    
  ax11.set_title('Histogram of the nucleosoms position', fontsize=10)
  ax11.set_xlabel('Position (bp)', fontsize=8)
  ax11.set_ylabel('Frequency', fontsize=8)
  ax11.tick_params(axis='x', labelsize=8)
  ax11.tick_params(axis='y', labelsize=8)
  ax11.grid() 
  ax11.axis('tight')
  ax11.set_xlim([data_min,data_max])
  max_histo = np.amax(histogram)
  #ax11.set_ylim([0,max_histo + 0.1* max_histo])
  
  fig1.set_tight_layout(fig1)
  fig1.savefig('Histogram_data_min_{0}_data_max_{1}.pdf'.format(data_min,data_max) )

  plt.clf()
  plt.close()
  
  print ('Histogram Ploted!')
  
  return 0


if __name__ =='__main__':
  plot_histogram(sys.argv[1],int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]))
  