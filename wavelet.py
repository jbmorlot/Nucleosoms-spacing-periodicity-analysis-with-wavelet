#!/bin/python

import scipy
from sys import stdout
from scipy import signal
from scipy import ndimage
import numpy as np
import multiprocessing as mp
#from multiprocessing.sharedctypes import Value
from math import *
import psutil
import time as tim
import os

from wavelet_plot import *
from input_file_analysis import *
from restart import *

def wavelet_convolution(histogram,w,period,histogram_length,data_range,period_length,mp_print,j):
  
    cwt = np.zeros(histogram_length)
    convol = scipy.signal.fftconvolve(histogram,scipy.signal.morlet(histogram_length,w,float(data_range)/(2.*period*w),False),"same")
    cwt = np.absolute(convol) 
    
    if int(100*float(j+1)/float(period_length))%(20) == 0 and mp_print.cwt < int(100*float(j+1)/float(period_length)):
      mp_print.cwt = int(100*float(j+1)/float(period_length))
      stdout.write("\r%d /100" % (100*float(j+1)/float(period_length)))
      stdout.flush()
      
    return ((float(data_range)/(2.*period*w)))*cwt

def wavelet_convolution_wrapper(args):
  return wavelet_convolution(*args)

#def maximum_map(cwt_array,histogram_length,period_max,print_max_map,position):
    
  #array_max_index = np.argmax(cwt_array)  
  
  #if int(100*float(position+1)/float(histogram_length))%(20) == 0 and print_max_map < int(100*float(position+1)/float(histogram_length)):
    #print_max_map = int(100*float(position+1)/float(histogram_length))
    #stdout.write("\r%d /100" % (100*float(position+1)/float(histogram_length)))
    #stdout.flush()
    
  #return period_max - array_max_index 

#def maximum_map_wrapper(args):
  #return maximum_map(*args)

def maximum_map_period(cwt_array,histogram_length,period_length,compare_value_map,mp_print,position):
    
  cwt_middle = np.zeros(period_length)    
  cwt_middle_index = scipy.signal.argrelmax(cwt_array, order = 10)
  cwt_middle_index = np.asarray(cwt_middle_index,np.int32)
  
  if compare_value_map == 0:
    cwt_middle[cwt_middle_index] = cwt_array[cwt_middle_index]
  else:
    cwt_middle[cwt_middle_index] = compare_value_map

  if int(100*float(position+1)/float(histogram_length))%(20) == 0 and mp_print.max_map_period < int(100*float(position+1)/float(histogram_length)):
    
    #print print_max_map_period
    #print int(100*float(position+1)/float(histogram_length))
    #print print_max_map_period < int(100*float(position+1)/float(histogram_length))
    mp_print.max_map_period = int(100*float(position+1)/float(histogram_length))
    stdout.write("\r%d /100" % (100*float(position+1)/float(histogram_length)))
    stdout.flush()
  
  return cwt_middle

def maximum_map_period_wrapper(args):
  return maximum_map_period(*args)  
  
def maximum_period_index(cwt_array,cwt_middle,histogram_length,mp_print,position):

  cwt_array_cpy = np.copy(cwt_array)
  index = np.where(cwt_middle > 0)[0]
  if index.shape[0] == 0:
    index = [0]
  
  max_1_index = index[np.argmax(cwt_array_cpy[index])]
  cwt_array_cpy[max_1_index] = 0
  max_2_index = index[np.argmax(cwt_array_cpy[index])]
  if cwt_array_cpy[max_2_index] == 0:
    max_2_index = 0
  cwt_array_cpy[max_2_index] = 0
  max_3_index = index[np.argmax(cwt_array_cpy[index])]
  if cwt_array_cpy[max_3_index] == 0:
    max_3_index = 0
    
  if int(100*float(position+1)/float(histogram_length))%(20) == 0 and mp_print.max_period_index < int(100*float(position+1)/float(histogram_length)):
    mp_print.max_period_index = int(100*float(position+1)/float(histogram_length))
    stdout.write("\r%d /100" % (100*float(position+1)/float(histogram_length)))
    stdout.flush()
  
  return max_1_index,max_2_index,max_3_index

def maximum_period_index_wrapper(args):
  return maximum_period_index(*args)

def maximum_map_period_convolution(cwt_middle,local_period_index,histogram_length,period_length,period_max,mp_print,position):
    
  cwt_out = np.zeros(period_length)
  cwt_out[int(local_period_index[0])] = 3
  cwt_out[int(local_period_index[1])] = 2
  cwt_out[int(local_period_index[2])] = 1
  step = np.zeros(period_length)
  step[ int(period_length/2) - 1 : int(period_length/2) + 1 ] = 1
  
  cwt_out = scipy.signal.fftconvolve(cwt_out,step,"same")
  
  if int(100*float(position+1)/float(histogram_length))%(20) == 0 and mp_print.max_map_conv < int(100*float(position+1)/float(histogram_length)) :
    mp_print.max_map_conv = int(100*float(position+1)/float(histogram_length))
    stdout.write("\r%d /100" % (100*float(position+1)/float(histogram_length)))
    stdout.flush()
      
  return cwt_out
  

def maximum_map_period_convolution_wrapper(args):
  return maximum_map_period_convolution(*args)  
  
def peak_hole_width(histogram,histogram_length, width_min, width_max, width_length, hole_width_min, hole_width_max, hole_width_length):
  #Signal convolved with gaussian
  sigma = 30
  
  histogram_min = np.amin(histogram)
  histogram = histogram - histogram_min
  histogram_max = np.amax(histogram)
  histogram_norm = np.zeros(histogram_length)
  list_peaks = np.where(histogram > histogram_max * 0.10)[0]
  histogram_norm[list_peaks] = histogram[list_peaks]
  
  histogram_gaussian = scipy.signal.fftconvolve(histogram_norm,signal.gaussian(histogram_length, std = sigma),"same")
  #histogram_gaussian_max = np.amax(histogram_gaussian)
  #histogram_gaussian_min = np.amin(histogram_gaussian)
  #histogram_gaussian = histogram_gaussian  - histogram_gaussian_min
  
  #deriv_gauss = np.diff(signal.gaussian(histogram_length + 1, std = sigma))
  histogram_gaussian_derivative = scipy.signal.fftconvolve(histogram_norm,np.diff(signal.gaussian(histogram_length+1, std = sigma)),"same")
  max_histogram_gaussian_derivative = np.amax(histogram_gaussian_derivative)

  maximum = 0
  imax = 0
  minimum = 0
  imin = 0
  width = 0
  prev_peak_position = 0
  
  histogram_width_peak = np.zeros(width_length)
  local_width_peak = np.zeros(histogram_length)
  
  histogram_width_hole = np.zeros(hole_width_length)
  local_width_hole = np.zeros(histogram_length)
  
  #Used to print only once each step
  print_local_width = 0
  
  for i in xrange(histogram_length-1):
    
    if histogram_gaussian_derivative[i] > 0 and histogram_gaussian_derivative[i] > maximum and minimum == 0:
      #We remove all small noise (< 10% max od derivative)
      if  abs(histogram_gaussian_derivative[i]) > 0.1*max_histogram_gaussian_derivative:
	maximum = histogram_gaussian_derivative[i]
	imax = i

    if histogram_gaussian_derivative[i] < 0 and histogram_gaussian_derivative[i] < minimum and maximum != 0:
      #We remove all small noise (< 10% max od derivative)
      if  abs(histogram_gaussian_derivative[i]) > 0.1*max_histogram_gaussian_derivative:
	minimum = histogram_gaussian_derivative[i]
	imin = i
	
    #We check that we are at the end of the gaussian with histogram_gaussian_derivative[i] * histogram_gaussian_derivative[i+1] < 0:
    if maximum != 0 and minimum != 0 and histogram_gaussian_derivative[i] * histogram_gaussian_derivative[i+1] < 0:
      width = abs(imax - imin)
      peak_position = int(abs(imin - imax)/2) + imax
      maximum = 0     
      minimum = 0
      #We build the width histogram and local function
      #We verify that the point is not out of range (due to the int conversion)
      if int(abs(imax - imin)/2) + imax < histogram_length : 
	local_width_peak[peak_position] = width
	if width < width_max and width >= width_min:
	  histogram_width_peak[int(width - width_min)] += 1
	
      #histogram of hole width
      hole_width = peak_position - prev_peak_position
      prev_peak_position = peak_position
      
      #prev_peak_position = 0
      if int(hole_width/2) + prev_peak_position < histogram_length : 
	local_width_hole[int(hole_width/2) + prev_peak_position ] = hole_width
	
      if hole_width < hole_width_max and hole_width >= hole_width_min:
	histogram_width_hole[int(hole_width - hole_width_min)] += 1

      
    if int(100*float(i+1)/float(histogram_length-1))%(10) == 0 and print_local_width < int(100*float(i+1)/float(histogram_length-1)):
      print_local_width = int(100*float(i+1)/float(histogram_length-1))
      stdout.write("\r%d /100" % (100*float(i+1)/float(histogram_length-1)))
      stdout.flush()

  return histogram_gaussian, histogram_width_peak, local_width_peak, histogram_width_hole, local_width_hole

  

def wavelet_analysis(histogram,period_min,period_max,period_length,save_var,plot_var,width_min,width_max,width_length,hole_width_min, hole_width_max,data_min,data_max,histogram_length,compare_value_map,bin_period,local_signal_var):
  
    
  print("Computing local width and width histogram ...") 
  
  histogram_gaussian, histogram_width, local_width, histogram_width_hole, local_width_hole = peak_hole_width(histogram,histogram_length, width_min, width_max, width_length, hole_width_min, hole_width_max, hole_width_max - hole_width_min)
      
  print ('\nDONE!\n')
  
  #Period
  w = 5		#Frequency of the wavelet
  data_range = data_max - data_min
  print("Computing Wavelet Convolution...")
  cwt = np.zeros((period_length,histogram_length))
  
  #We share a parameter throught all process in order to print only once each step (avoid too much I/O)
  mp_print = mp.Manager().Namespace()
    
  #Used to print only once each step
  mp_print.cwt = 0
  
  pool_cwt = mp.Pool()
  cwt_arg = [(histogram,w,(period_max-period*bin_period),histogram_length,data_range,period_length,mp_print,period) for period in range(period_length)]
  cwt_1 = pool_cwt.map(wavelet_convolution_wrapper, cwt_arg)
  pool_cwt.close()
  pool_cwt.join()
  cwt_2 = np.hstack(cwt_1)
  cwt = cwt_2.reshape((period_length,histogram_length))

  print("\nDONE!\n")
  
  print("Searching Local Maxima in Wavelet Map ...")
  map_local_period = np.zeros((period_length,histogram_length))
  
  #Used to print only once each step
  mp_print.max_map_period = 0
  
  pool_period_max_map = mp.Pool()
  period_max_map_arg = [(cwt[:,position],histogram_length,period_length,compare_value_map,mp_print,position) for position in range(histogram_length)]
  map_local_period_1 = pool_period_max_map.map(maximum_map_period_wrapper, period_max_map_arg)
  pool_period_max_map.close()
  pool_period_max_map.join()
  map_local_period_2 = np.hstack(map_local_period_1)
  map_local_period_transpose = map_local_period_2.reshape((histogram_length,period_length))
  map_local_period = np.transpose(map_local_period_transpose)

  print("\nDONE!\n")
  
  print("Identifying Local Period from Local Maxima ...")
  local_period = np.zeros((3,histogram_length))

  #Used to print only once each step
  mp_print.max_period_index = 0
  
  pool_period_index_max = mp.Pool()
  period_index_max_arg = [(cwt[:,position],map_local_period[:,position],histogram_length,mp_print,position) for position in range(histogram_length)]
  local_period_index_1 = pool_period_index_max.map(maximum_period_index_wrapper, period_index_max_arg)
  pool_period_index_max.close()
  pool_period_index_max.join()
  local_period_index_2 = np.hstack(local_period_index_1)
  local_period_index_transpose = local_period_index_2.reshape((histogram_length,3))
  local_period_index = np.transpose(local_period_index_transpose)
    
  print("\nDONE!\n")
  
  print("Computing period histogram ...")
  
  histogram_period_1 = np.zeros(period_length)
  histogram_period_2 = np.zeros(period_length)
  histogram_period_3 = np.zeros(period_length)
  
  #Used to print only once each step
  print_histo_period = 0
  
  for i in xrange(histogram_length):
    
    if local_period_index[0,i] != 0:
      histogram_period_1[period_length - local_period_index[0,i]] += 1
      
    if local_period_index[1,i] != 0:
      histogram_period_2[period_length - local_period_index[1,i]] += 1
      
    if local_period_index[2,i] != 0:
      histogram_period_3[period_length - local_period_index[2,i]] += 1
      
    if int(100*float(i+1)/float(histogram_length))%(20) == 0 and print_histo_period < int(100*float(i+1)/float(histogram_length)) :
      print_histo_period = int(100*float(i+1)/float(histogram_length))
      stdout.write("\r%d /100" % (100*float(i+1)/float(histogram_length)))
      stdout.flush()
      
  print("\nDONE!\n")
  
  #This part is in order to make the map of the maximum more readable
  map_local_period_convolution = np.zeros((period_length,histogram_length))
  if local_signal_var == 1:
    print("Computing Convolution of the Map with step for the display ...")

    #Used to print only once each step
    mp_print.max_map_conv  = 0
    
    pool_period_max_map_convolution = mp.Pool()
    period_max_map_convolution_arg = [(map_local_period[:,position],local_period_index[:,position],histogram_length,period_length,period_max,mp_print,position) for position in range(histogram_length)]
    map_local_period_convolution_1 = pool_period_max_map_convolution.map(maximum_map_period_convolution_wrapper, period_max_map_convolution_arg)
    pool_period_max_map_convolution.close()
    pool_period_max_map_convolution.join()
    map_local_period_convolution_2 = np.hstack(map_local_period_convolution_1)
    map_local_period_convolution_transpose = map_local_period_convolution_2.reshape((histogram_length,period_length))
    map_local_period_convolution = np.transpose(map_local_period_convolution_transpose)

    print("\nDONE!\n")
    
    #We free the memory
    del mp_print
  
  return histogram_period_1,histogram_period_2,histogram_period_3,histogram_width,cwt,map_local_period_convolution,local_width,histogram_gaussian, histogram_width_hole, local_width_hole

  
def wavelet_analysis_pieces(data_file,launch_all,data_min,data_max,period_min,period_max,width_min,width_max,hole_width_min,hole_width_max,bin_period,plot_var,save_var,save_txt_var,local_signal_var,global_signal_var,restart_var,histogram_period_1,histogram_period_2,histogram_period_3,histogram_width,sub_directory):
  
  if (data_file == ''):
    print("WARNING: No input file ")
    return 0

  if period_max < period_min:
    print("WARNING: period max < period min")
    return 0
    
  if width_max < width_min :
    print("WARNING: width max < width min")
    return 0
  
  if hole_width_max < hole_width_min :
    print("WARNING: hole width max < hole width min")
    return 0
  
  if  data_max < data_min :
    print("WARNING: data max < data min")
    return 0

  if  bin_period  <= 0:
    print("WARNING: bin period length <= 0")
    return 0
    
  if (period_max - period_min) % bin_period != 0:
    print ('WARNING: Select a bin_period which could divide (period_max - period_min) ')
    return 0
    
  #Analysing the size of the file
  print('\nDetermining the length of the Data ...')
  data_length, data_min, data_max = input_file_analysis(data_file,data_min,data_max)
  print('DONE!\n')
  
  #Parameters
  period_length = int((period_max - period_min)/bin_period)
  width_length = width_max - width_min
  hole_width_length = hole_width_max - hole_width_min
  
  print("\n--------------- Parameters ---------------------\n\nData Length = {0}\nData Min = {1}\nData Max = {2}\nPeriod Min = {3}\nPeriod Max = {4}\nPeriod Length = {5}\nBin Period = {6}\nPeak Width Min = {7}\nPeak Width Max = {8}\nPeak Width Length = {9}\nHole Width Min = {10}\nHole Width Max = {11}\nHole Width Length = {12}\n\n----------------------------------------------\n").format(data_length,data_min,data_max,period_min,period_max,period_length,bin_period,width_min,width_max,width_length,hole_width_min,hole_width_max,hole_width_max - hole_width_min)
  
  #Calculating the number of pieces (The number of split is calculated in order have memory of wavelet map = 2% of avalaible memory (Biggest object in memory) )
  memory_available = int(psutil.virtual_memory().available)
  split = int(4*period_length * data_length/(memory_available*0.5/100))
  
  if split <= 0:
    split = 1

  print('\n-------- Memory Information ------------\n')
  
  print ('Avalaible Memory = {0} GB'.format(int(memory_available/10e8)))
  print ('Number of Pieces = {0}'.format(split))
  print ('Avalaible Memory for Each Pieces = {0} MB'.format(int(memory_available/(split*10e5))))
  
  print('\n----------------------------------------\n')
  
  
  
  #Launching the analysis pieces by pieces
  
  print('\n-------- Analysis Launched! ------------\n')
  
  #Some needed variables
  piece_length = int(data_length/split)
  
  histogram_period_1 = np.zeros(period_length)
  histogram_period_2 = np.zeros(period_length)
  histogram_period_3 = np.zeros(period_length)
  histogram_width = np.zeros(width_length)
  histogram_width_hole = np.zeros(hole_width_length)
  
  #Used to save only once each step
  save_restart_control = 0
  
  time_start = tim.time()
  
  for i in xrange(split):
    
    #Timer
    if i != 0:
      time_stop = tim.time()
      print('\n--------------------------------------------------------------')
      print('/ Estimated Remaining Time: {0} day(s) {1} hour(s) {2} minute(s)  /'.format(int((split-i)*(time_stop - time_start)/(3600*24)),int((split-i)*(time_stop - time_start)/3600)%24,int((split-i)*(time_stop - time_start)/60)%60))
      print('----------------------------------------------------------------\n')
      time_start = tim.time()
    
    print('\n-------- Computing Piece {0} / {1} ------------\n'.format(i+1,split))
    
    data_file_min = data_min + i*piece_length
    data_file_max = data_min + (i+1)*piece_length
    
    #Loading data_file
    histogram, histogram_length = load_data(data_file,data_file_min,data_file_max)
    #histogram = np.cos(2*pi*200*np.linspace(0,10000,10000))
    #histogram_length = 10000
    
    
    #Wavelet Analysis
    histogram_period_1_current,histogram_period_2_current,histogram_period_3_current,histogram_width_current,cwt,map_local_period_convolution,local_width,histogram_gaussian, histogram_width_hole_current, local_width_hole = wavelet_analysis(histogram,period_min,period_max,period_length,save_var,plot_var,width_min,width_max,width_length,hole_width_min,hole_width_max,data_file_min,data_file_max,histogram_length,0,bin_period,local_signal_var)
    
    #Plotting local signal if asked
    if int(local_signal_var) == 1:
      wavelet_local_plot(data_file_min, data_file_max,histogram_length,period_min,period_max,width_min,width_max,hole_width_min,hole_width_max,plot_var,save_var,save_txt_var,local_width,local_width_hole,map_local_period_convolution,cwt,histogram,histogram_gaussian,sub_directory)
    
    #Adding computed histograms to older histograms
    histogram_period_1 = np.add(histogram_period_1,histogram_period_1_current)
    histogram_period_2 = np.add(histogram_period_2,histogram_period_2_current)
    histogram_period_3 = np.add(histogram_period_3,histogram_period_3_current)
    histogram_width = np.add(histogram_width,histogram_width_current)
    histogram_width_hole = np.add(histogram_width_hole,histogram_width_hole_current)
    
    #Saving in case of a crash
    if int(100*float(i+1)/float(split))%(10) == 0 and save_restart_control < int(100*float(i+1)/float(split)):
      save_restart_control = int(100*float(i+1)/float(split))
      input_data = (launch_all,data_file_min,data_max,period_min,period_max,width_min,width_max,hole_width_min,hole_width_max,bin_period,plot_var,save_var,save_txt_var,local_signal_var,global_signal_var,restart_var)
      save_restart(data_file,input_data,histogram_period_1,histogram_period_2,histogram_period_3,histogram_width,histogram_width_hole)
    
    #We free the memory
    del cwt
    del local_width
    del map_local_period_convolution
    del histogram
    del histogram_gaussian
  
  if data_length % split > 1000:
    print('\n-------- Computing Last Piece ------------\n')
    
    data_file_min = data_min + piece_length*split
    data_file_max = data_max
    
    #Loading data_file
    histogram, histogram_length = load_data(data_file,data_file_min,data_file_max)
    
    #Numeric Analysis
    histogram_period_1_current,histogram_period_2_current,histogram_period_3_current,histogram_width_current,cwt,map_local_period_convolution,local_width,histogram_gaussian = wavelet_analysis(histogram,period_min,period_max,period_length,save_var,plot_var,width_min,width_max,width_length,hole_width_min,hole_width_max,data_file_min,data_file_max,histogram_length,0,bin_period,local_signal_var)
    
    #Plotting local signal if asked
    if int(local_signal_var) == 1:
      wavelet_local_plot(data_file_min, data_file_max,histogram_length,period_min,period_max,width_min,width_max,hole_width_min,hole_width_max,plot_var,save_var,save_txt_var,local_width,local_width_hole,map_local_period_convolution,cwt,histogram,histogram_gaussian,sub_directory)
    
    
    #We free the memory
    del cwt
    del local_width
    del map_local_period_convolution
    del histogram
    del histogram_gaussian
    
    #Adding computed histograms to older histograms
    histogram_period_1 = np.add(histogram_period_1,histogram_period_1_current)
    histogram_period_2 = np.add(histogram_period_2,histogram_period_2_current)
    histogram_period_3 = np.add(histogram_period_3,histogram_period_3_current)
    histogram_width = np.add(histogram_width,histogram_width_current)
    histogram_width_hole = np.add(histogram_width_hole,histogram_width_hole_current)

    
  print('\n-------- Analysis DONE! ------------\n')
  
  print('Normalization of Histograms')
  histogram_period_1_sum = np.sum(histogram_period_1)
  histogram_period_2_sum = np.sum(histogram_period_2)
  histogram_period_3_sum = np.sum(histogram_period_3)
  histogram_width_sum = np.sum(histogram_width)
  histogram_width_hole_sum = np.sum(histogram_width_hole)
  if histogram_period_1_sum != 0:
    histogram_period_1 = histogram_period_1 / histogram_period_1_sum
  if histogram_period_2_sum != 0:
    histogram_period_2 = histogram_period_2 / histogram_period_2_sum
  if histogram_period_3_sum != 0:
    histogram_period_3 = histogram_period_3 / histogram_period_3_sum
  if histogram_width_sum != 0:
    histogram_width =  histogram_width / histogram_width_sum
  if histogram_width_hole_sum != 0:
    histogram_width_hole =  histogram_width_hole / histogram_width_hole_sum
  
  print('DONE!\n')
  
  if int(global_signal_var) == 1:
    wavelet_histogram_plot(data_min,data_max,period_min,period_max,period_length,width_min,width_max,width_length,hole_width_min,hole_width_max,hole_width_length,plot_var,save_var,save_txt_var,histogram_period_1,histogram_period_2,histogram_period_3,histogram_width,histogram_width_hole,sub_directory)
  
  #We delete the Checkpoints if the restart is not launched
  if restart_var == 0:
    delete_restart(data_file)
  

  return 0