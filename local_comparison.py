#!/bin/python

import scipy
from scipy import signal
from scipy import ndimage
import numpy as np
from math import *
import psutil
import time as tim
import os
import sys
from matplotlib.patches import Circle
import matplotlib.transforms as transforms

from input_file_analysis import *

def nucleosom_match_local_plot(histogram_gaussian_1,histogram_gaussian_2,nucl_state_1,nucl_state_2,nucl_state_1_reverse,nucl_state_2_reverse,peaks_1_position_length,peaks_2_position_length,name_data_1,name_data_2,data_min,data_max,reverse_var,sub_directory,save_txt_var):
  
  #We change the directory
  os.chdir(sub_directory)
  
  if save_txt_var == 1:
    np.savetxt('Nucleosom_comparaison_states_1_data_min{0}_data_max{1}.txt'.format(data_min,data_max),nucl_state_1)
    np.savetxt('Nucleosom_comparaison_states_2_data_min{0}_data_max{1}.txt'.format(data_min,data_max),nucl_state_2)
  
  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt

  fig1 = plt.figure()
  ax1 = fig1.add_subplot(111)
  
  
  #We normalize the histograms before plotting
  print("Normalizing the histograms before plotting...")
  histogram_gaussian_1_sum = np.sum(histogram_gaussian_1)
  histogram_gaussian_2_sum = np.sum(histogram_gaussian_2)
  histogram_gaussian_1 = histogram_gaussian_1/histogram_gaussian_1_sum
  histogram_gaussian_2 = histogram_gaussian_2/histogram_gaussian_2_sum
  print("Done!\n")
  
  #We adjuste the window to be symmetric on y
  print("Adjusting the size of the window...")
  histogram_gaussian_1_max = np.amax(histogram_gaussian_1)
  histogram_gaussian_2_max = np.amax(histogram_gaussian_2)
  histogram_gaussian_max = max(histogram_gaussian_1_max,histogram_gaussian_2_max)
  print("Done!\n")
  
  #print histogram_gaussian_1.shape
  #print data_max - data_min
  
  ax1.plot(np.linspace(data_min,data_max, data_max - data_min ),histogram_gaussian_1)
  ax1.plot(np.linspace(data_min,data_max, data_max - data_min ),-histogram_gaussian_2)
  ax1.set_title('Nucleosom Position (Histogram Up: '+ name_data_1 + ' ; Histogram Down: ' + name_data_2 + ' )', fontsize=10)
  ax1.set_xlabel('Position (bp)', fontsize=8)
  ax1.set_ylabel('Frequency of Histograms Normalized', fontsize=8)
  ax1.tick_params(axis='x', labelsize=8)
  ax1.tick_params(axis='y', labelsize=8)
  ax1.grid() 
  ax1.axis('tight')
  ax1.set_xlim([data_min,data_max])
  ax1.set_ylim([-histogram_gaussian_max - 0.1*histogram_gaussian_max,histogram_gaussian_max + 0.1*histogram_gaussian_max])
  
  trans_1 = ax1.transAxes
  
  data_length = data_max - data_min
  high = 0.01
  diameter = 0.002
  shift = 0.05
  
  for i in xrange(peaks_1_position_length):
    
    #el1 = Circle( (nucl_state_1[i][3]/data_length,0.5 + 2*high),diameter, color='w',ec = 'k',transform=trans_1)
    #ax1.add_artist(el1)
    
    if nucl_state_1[i][0] != 0:
      el1 = Circle((nucl_state_1[i][3]/data_length, 0.5 + high ),diameter ,color='r',transform=trans_1)
      ax1.add_artist(el1)
      el2 = Circle((nucl_state_1[i][3]/data_length, 0.5 - high ),diameter ,color='b',transform=trans_1)
      ax1.add_artist(el2)
    
    if nucl_state_1[i][1] != 0:
      el2 = Circle((nucl_state_1[i][3]/data_length, 0.5  + high),diameter,color='c',transform=trans_1)
      ax1.add_artist(el2)
    
    if nucl_state_1[i][2] != 0:
      el2 = Circle((nucl_state_1[i][3]/data_length,  0.5 + high),diameter,color='k',transform=trans_1)
      ax1.add_artist(el2)
    
    if reverse_var == 1:
      
      #el1 = Circle( (nucl_state_1_reverse[i][3]/data_length,0.5 + 2*high),diameter, color='w',ec = 'k',transform=trans_1)
      #ax1.add_artist(el1)
    
      if nucl_state_1_reverse[i][0] != 0:
	el1 = Circle((nucl_state_1_reverse[i][3]/data_length, 0.5 + 2*high ),diameter ,color='r',transform=trans_1)
	ax1.add_artist(el1)
	el2 = Circle((nucl_state_1_reverse[i][3]/data_length,  0.5 - 2*high ),diameter,color='b',transform=trans_1)
	ax1.add_artist(el2)
      
      if nucl_state_1_reverse[i][1] != 0:
	el2 = Circle((nucl_state_1_reverse[i][3]/data_length, 0.5  + 2*high),diameter,color='c',transform=trans_1)
	ax1.add_artist(el2)
      
      if nucl_state_1_reverse[i][2] != 0:
	el2 = Circle((nucl_state_1_reverse[i][3]/data_length,  0.5 + 2*high ),diameter,color='k',transform=trans_1)
	ax1.add_artist(el2)
    
  for i in xrange(peaks_2_position_length):
    
    if nucl_state_2[i][0] != 0:
      el2 = Circle((nucl_state_2[i][3]/data_length,  0.5 - high),diameter,color='r',transform=trans_1)
      ax1.add_artist(el2)
      el1 = Circle((nucl_state_2[i][3]/data_length, 0.5 + high ),diameter ,color='b',transform=trans_1)
      ax1.add_artist(el1)
      
    if nucl_state_2[i][1] != 0:
      el2 = Circle((nucl_state_2[i][3]/data_length, 0.5  - high),diameter,color='c',transform=trans_1)
      ax1.add_artist(el2)
    
    if nucl_state_2[i][2] != 0:
      el2 = Circle((nucl_state_2[i][3]/data_length,  0.5 - high),diameter,color='k',transform=trans_1)
      ax1.add_artist(el2)
      
    if reverse_var == 1:
      
      if nucl_state_2_reverse[i][0] != 0:
	el1 = Circle((nucl_state_2_reverse[i][3]/data_length, 0.5 + 2*high ),diameter ,color='b',transform=trans_1)
	ax1.add_artist(el1)
	el2 = Circle((nucl_state_2_reverse[i][3]/data_length,  0.5 - 2*high ),diameter,color='r',transform=trans_1)
	ax1.add_artist(el2)
	
      if nucl_state_2_reverse[i][1] != 0:
	el2 = Circle((nucl_state_2_reverse[i][3]/data_length, 0.5  - 2*high),diameter,color='c',transform=trans_1)
	ax1.add_artist(el2)
      
      if nucl_state_2_reverse[i][2] != 0:
	el2 = Circle((nucl_state_2_reverse[i][3]/data_length,  0.5 - 2*high ),diameter,color='k',transform=trans_1)
	ax1.add_artist(el2)
  
  
  fig1.set_tight_layout(fig1)
  fig1.savefig('Nucleosom Position'+ name_data_2 + '_data_min_{0}_data_max_{1}.pdf'.format(data_min,data_max) )
  
  plt.clf()
  plt.close()
  
  return 0
  
  
def nucleosom_match_global_plot(number_nucleosom_detected,gain,loss,occupancy,shift_tot,histogram_shift,name_data_2,data_min,data_max,sub_directory,save_txt_var):
  
  os.chdir(sub_directory)
  if save_txt_var == 1:
    np.savetxt('Nucleosom_comparaison_states_gain_loss_occupancy_data_min{0}_data_max{1}.txt'.format(data_min,data_max),(gain,loss,occupancy))
    np.savetxt('Nucleosom_comparaison_states_shift_histogram_data_min{0}_data_max{1}.txt'.format(data_min,data_max),histogram_shift)

  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt
  
  
  fig1 = plt.figure()
  ax11 = fig1.add_subplot(111)
  #ax11.hist(histogram_shift,bins=160,range=(-80,80),facecolor="r", histtype = 'step')
  ax11.plot(np.linspace(-80,80, 160), histogram_shift)
  ax11.set_title('Shift Histogram', fontsize=10)
  ax11.set_xlabel('Shift (bp)', fontsize=8)
  ax11.set_ylabel('Frequency', fontsize=8)
  ax11.tick_params(axis='x', labelsize=8)
  ax11.tick_params(axis='y', labelsize=8)
  ax11.grid() 
  ax11.axis('tight')
  ax11.set_xlim([-80,80])
  max_shift = np.amax(histogram_shift)
  if max_shift == 0:
    max_shift = 1
  ax11.set_ylim([-0.1 * max_shift,max_shift + 0.4 * max_shift])
  
  text = "Comparaison Summary: Nucleosom Detected = {0}\nGain = {1} ; Loss = {2} ; Occupancy = {3} ; Shift = {4}".format(number_nucleosom_detected,gain,loss,occupancy,shift_tot)
  
  print('\n----------------------------------------\n')
  print (text)
  print('\n----------------------------------------\n')
  
  ax11.text(0.05, 0.95, text, transform=ax11.transAxes, fontsize=8, fontweight='bold', va='top')
  
  fig1.set_tight_layout(fig1)
  fig1.savefig('Histogram_shift'+ name_data_2 + '_data_min_{0}_data_max_{1}.pdf'.format(data_min,data_max) )

  plt.clf()
  plt.close()
  
  return 0

 
def nucleosom_match(histogram_1,histogram_2,histogram_length,name_data_1,name_data_2,data_min,data_max,save_txt_var,local_signal_var,global_signal_var,reverse_var,sub_directory):

  #Signal convolved with gaussian
  sigma = 30
  
  print("Computing peaks position of histogram 1 ...")
  
  #histogram_min_1 = np.amin(histogram_1)
  #histogram_1 = histogram_1 - histogram_min_1
  #histogram_mean_1 = np.sum(histogram_1) / histogram_length
  #list_peaks_1 = np.where(histogram_1 > histogram_mean_1 * 1.9)[0]
  #histogram_1_norm = np.zeros(histogram_length)
  #histogram_1_norm[list_peaks_1] = histogram_1[list_peaks_1]
  
  histogram_min_1 = np.amin(histogram_1)
  histogram_1 = histogram_1 - histogram_min_1
  histogram_max_1 = np.amax(histogram_1)
  histogram_1_norm = np.zeros(histogram_length)
  list_peaks_1 = np.where(histogram_1 > histogram_max_1 * 0.1)[0]
  histogram_1_norm[list_peaks_1] = histogram_1[list_peaks_1]
  
  histogram_gaussian_1 = scipy.signal.fftconvolve(histogram_1_norm,signal.gaussian(histogram_length, std = sigma),"same")
  
  histogram_gaussian_derivative_1 = scipy.signal.fftconvolve(histogram_1_norm,np.diff(signal.gaussian(histogram_length+1, std = sigma)),"same")

  max_histogram_gaussian_derivative_1 = np.amax(histogram_gaussian_derivative_1)
  
  maximum = 0
  imax = 0
  minimum = 0
  imin = 0
  width = 0
  counter_width = 0
  
  #Used to print only once each step
  print_once = 0

  k=0
  peaks_1 = np.zeros(histogram_length)
  
  for i in xrange(histogram_length-1):
    
    #We want first a max and the a min in order to have a gaussian and not the space between 2 gaussian
    if histogram_gaussian_derivative_1[i] > 0 and histogram_gaussian_derivative_1[i] > maximum and minimum == 0:
      #We remove all small noise (< 10% max od derivative)
      if  abs(histogram_gaussian_derivative_1[i]) > 0.1*max_histogram_gaussian_derivative_1:
	maximum = histogram_gaussian_derivative_1[i]
	imax = i
    
    if histogram_gaussian_derivative_1[i] < 0 and histogram_gaussian_derivative_1[i] < minimum and maximum != 0:
      #We remove all small noise (< 10% max od derivative)
      if  abs(histogram_gaussian_derivative_1[i]) > 0.1*max_histogram_gaussian_derivative_1:
	minimum = histogram_gaussian_derivative_1[i]
	imin = i
      
    if maximum != 0 and minimum != 0 and histogram_gaussian_derivative_1[i] * histogram_gaussian_derivative_1[i+1] < 0:
      maximum = 0     
      minimum = 0
      
      if int(abs(imax - imin)/2) + imin < histogram_length : 
	peaks_1[int(abs(imax - imin)/2) + imax] = int(abs(imax - imin)/2) + imax

            
    if int(100*float(i+1)/float(histogram_length-1))%(10) == 0 and print_once < int(100*float(i+1)/float(histogram_length-1)):
      print_once = int(100*float(i+1)/float(histogram_length-1))
      stdout.write("\r%d /100" % (100*float(i+1)/float(histogram_length-1)))
      stdout.flush()
      
  print ('\nDONE!\n')
  
  print("Computing peaks position of histogram 2 ...")
  
  #histogram_min_2 = np.amin(histogram_2)
  #histogram_2 = histogram_2 - histogram_min_2
  #histogram_mean_2 = np.sum(histogram_2) / histogram_length
  #list_peaks_2 = np.where(histogram_2 > histogram_mean_2 * 1.9)[0]
  #histogram_2_norm = np.zeros(histogram_length)
  #histogram_2_norm[list_peaks_2] = histogram_2[list_peaks_2]
  
  histogram_min_2 = np.amin(histogram_2)
  histogram_2 = histogram_2 - histogram_min_2
  histogram_max_2 = np.amax(histogram_2)
  histogram_2_norm = np.zeros(histogram_length)
  list_peaks_2 = np.where(histogram_2 > histogram_max_2 * 0.1)[0]
  histogram_2_norm[list_peaks_2] = histogram_2[list_peaks_2]
  
  histogram_gaussian_2 = scipy.signal.fftconvolve(histogram_2_norm,signal.gaussian(histogram_length, std = sigma),"same")
  #histogram_gaussian_max_2 = np.amax(histogram_gaussian_2)

  
  histogram_gaussian_derivative_2 = scipy.signal.fftconvolve(histogram_2_norm,np.diff(signal.gaussian(histogram_length+1, std = sigma)),"same")

  max_histogram_gaussian_derivative_2 = np.amax(histogram_gaussian_derivative_2)
  
  maximum = 0
  imax = 0
  minimum = 0
  imin = 0
  width = 0
  counter_width = 0
  
  
  #Used to print only once each step
  print_once = 0
  
  peaks_2 = np.zeros(histogram_length)
  
  for i in xrange(histogram_length-1):
    
    #We want first a max and the a min in order to have a gaussian and not the space between 2 gaussian
    if histogram_gaussian_derivative_2[i] > 0 and histogram_gaussian_derivative_2[i] > maximum and minimum == 0:
      #We remove all small noise (< 10% max od derivative)
      if  abs(histogram_gaussian_derivative_2[i]) > 0.1*max_histogram_gaussian_derivative_2:
	maximum = histogram_gaussian_derivative_2[i]
	imax = i
    
    if histogram_gaussian_derivative_2[i] < 0 and histogram_gaussian_derivative_2[i] < minimum and maximum != 0:
      #We remove all small noise (< 10% max od derivative)
      if  abs(histogram_gaussian_derivative_2[i]) > 0.1*max_histogram_gaussian_derivative_2:
	minimum = histogram_gaussian_derivative_2[i]
	imin = i
    
    if maximum != 0 and minimum != 0 and histogram_gaussian_derivative_2[i] * histogram_gaussian_derivative_2[i+1] < 0:
      maximum = 0     
      minimum = 0
      
      if int(abs(imax - imin)/2) + imin < histogram_length : 
	peaks_2[int(abs(imax - imin)/2) + imax] = int(abs(imax - imin)/2) + imax 
	      
    if int(100*float(i+1)/float(histogram_length-1))%(10) == 0 and print_once < int(100*float(i+1)/float(histogram_length-1)):
      print_once = int(100*float(i+1)/float(histogram_length-1))
      stdout.write("\r%d /100" % (100*float(i+1)/float(histogram_length-1)))
      stdout.flush()
      
  print ('\nDONE!\n')
  
  
  print('Matching the 2 histograms Forward...')
  
  #We locate in each the postion of every peak
  peaks_1_position = np.where(peaks_1 > 0)[0]
  peaks_2_position = np.where(peaks_2 > 0)[0]
  
  peaks_1_position_length = peaks_1_position.shape[0]
  peaks_2_position_length = peaks_2_position.shape[0]
  
  if peaks_1_position_length <= 1 or peaks_2_position_length <= 1:
    print("\nWarning: No peak detected in this section!")
  
  #Used to print only once each step
  print_once = 0
  
  k=0
  j=0
  
  while j < peaks_1_position_length:
    
    if j == 0:
      nucl_state_1 = np.zeros((peaks_1_position_length,4))		#0:gain, 1:shift 2:occupancy 3:position
      nucl_state_2 = np.zeros((peaks_2_position_length,4))		#0:gain, 1:shift 2:occupancy 3:position
    
    position_1 = peaks_1_position[j]
    position_2 = peaks_2_position[k]
    
    nucl_state_1[j][3] = position_1
    nucl_state_2[k][3] = position_2
    
    #print (j , k)
    #print (position_1,position_2)
    
    if position_2 <= position_1 - 80:
      nucl_state_2[k][0] = +1							#gain for 2 or loss for 1 at position_2
      
      if k == peaks_2_position_length - 1:					#If it is the last peak_2 , all peak_1 beyond will be loss for 2
	nucl_state_1[j][0] = +1
      
      #Go to the next peak_2 and stay for peak_1
      if k == peaks_2_position_length - 1:					#If it is the last peak_2, only 1 moves
	j += 1
      if k < peaks_2_position_length - 1:					
	k += 1
     
      continue
     
    
    if position_2 >= position_1 + 80:
      nucl_state_1[j][0] = +1							#gain for 1 or loss for 2 at position_1
      
      #Go to the next peak_1 and stay for peak_2
      j += 1									
      continue
      
    
    if position_1 + 80 > position_2 and position_2 > position_1 - 80:
      
      if j < peaks_1_position_length - 1:
	position_1_next = peaks_1_position[j+1]
	
	if abs(position_1 - position_2) <= abs(position_1_next - position_2):	#The shortest distance = same nucleosom
	  if position_1 + 30 >= position_2 and position_2 >= position_1 - 30:
	    nucl_state_1[j][2] = 1						#Occupancy
	    nucl_state_2[k][2] = 1						#Occupancy
	
	  if position_1 + 30 < position_2 or position_2 < position_1 - 30:
	    nucl_state_1[j][1] = position_2 - position_1			#Shift 2 in comparaison of 1
	    nucl_state_2[k][1] = position_1 - position_2			#Shift 1 in comparaison of 2
	  #Go to the next for peak_1 and peak_2
	  j += 1 
	  if k < peaks_2_position_length - 1:
	    k += 1
	  continue
	
	if abs(position_1 - position_2) > abs(position_1_next - position_2):	#The shortest distance = same nucleosom
	  nucl_state_1[j][0] = 1						#Loss for 2 at position_1
	  
	  if position_1_next + 30 >= position_2 and position_2 >= position_1_next - 30:
	    nucl_state_1[j+1][2] = 1						#Occupancy
	    nucl_state_2[k][2] = 1						#Occupancy
	
	  if position_1_next + 30 < position_2 or position_2 < position_1_next - 30: 
	    nucl_state_1[j+1][1] = position_2 - position_1_next		#Shift 2 in comparaison of 1
	    nucl_state_2[k][1] = position_1_next - position_2			#Shift 1 in comparaison of 2
	    
	  #Go to the next for peak_1 and peak_2
	  j += 2 
	  if k < peaks_2_position_length - 1:
	    k += 1
	  continue
	
      if j == peaks_1_position_length - 1:

	if position_1 + 30 >= position_2 and position_2 >= position_1 - 30:
	  nucl_state_1[j][2] = 1						#Occupancy
	  nucl_state_2[k][2] = 1						#Occupancy
      
	if position_1 + 30 < position_2 or position_2 < position_1 - 30: 
	  nucl_state_1[j][1] = position_2 - position_1				#Shift 2 in comparaison of 1
	  nucl_state_2[k][1] = position_1 - position_2				#Shift 1 in comparaison of 2
	  
	#Go to the next for peak_1 and peak_2
	j += 1 
	if k < peaks_2_position_length - 1:
	  k += 1
	continue
	
    if int(100*float(j+1)/float(peaks_1_position_length-1))%(10) == 0 and print_once < int(100*float(j+1)/float(peaks_1_position_length-1)):
      print_once = int(100*float(j+1)/float(peaks_1_position_length-1))
      stdout.write("\r%d /100" % (100*float(j+1)/float(peaks_1_position_length-1)))
      stdout.flush()
    
  print ('\nDONE!\n')
  
  if reverse_var == 1:
    
    print('Matching the 2 histograms Backward...')
    
    
    if peaks_1_position_length <= 1 or peaks_2_position_length <= 1:
      print("\nWarning: No peak detected in this section!")

    #Used to print only once each step
    print_once = 0
    
    k=0
    j=0
    
    while j < peaks_1_position_length:
      
      if j == 0:
	#reverse
	nucl_state_1_reverse = np.zeros((peaks_1_position_length,4))		#0:gain, 1:shift 2:occupancy 3:position
	nucl_state_2_reverse = np.zeros((peaks_2_position_length,4))		#0:gain, 1:shift 2:occupancy 3:position
	
	peaks_1_position_reverse = peaks_1_position[::-1]
	peaks_2_position_reverse = peaks_2_position[::-1]
    
      
      position_1 = peaks_1_position_reverse[j]
      position_2 = peaks_2_position_reverse[k]
      
      nucl_state_1_reverse[j][3] = position_1
      nucl_state_2_reverse[k][3] = position_2
      
      if position_2 >= position_1 + 80:
	if k < peaks_2_position_length - 1:
	  nucl_state_2_reverse[k][0] = +1						#gain for 2 or loss for 1 at position_2
	
	if k == peaks_2_position_length - 1:
	  nucl_state_1_reverse[j][0] = +1						#If it is the last peak_2 , all peak_1 beyond will be loss for 2
	
	#Go to the next peak_2 and stay for peak_1
	if k == peaks_2_position_length - 1:						#If it is the last peak_2, only 1 moves
	  j += 1
	if k < peaks_2_position_length - 1:					
	  k += 1
	
	continue
      
      if position_2 <= position_1 - 80:
	nucl_state_1_reverse[j][0] = +1						#gain for 1 or loss for 2 at position_1
	
	#Go to the next peak_1 and stay for peak_2
	j += 1									
	continue
	
      
      if position_1 + 80 > position_2 and position_2 > position_1 - 80:
	
	if j < peaks_1_position_length - 1:
	  position_1_next = peaks_1_position_reverse[j+1]
	  
	  if abs(position_1 - position_2) <= abs(position_1_next - position_2):	#The shortest distance = same nucleosom
	    if position_1 + 30 >= position_2 and position_2 >= position_1 - 30:
	      nucl_state_1_reverse[j][2] = 1						#Occupancy
	      nucl_state_2_reverse[k][2] = 1						#Occupancy
	  
	    if position_1 + 30 < position_2 or position_2 < position_1 - 30:
	      nucl_state_1_reverse[j][1] = position_2 - position_1			#Shift 2 in comparaison of 1
	      nucl_state_2_reverse[k][1] = position_1 - position_2			#Shift 1 in comparaison of 2
	    #Go to the next for peak_1 and peak_2
	    j += 1 
	    if k < peaks_2_position_length - 1:
	      k += 1
	    continue
	  
	  if abs(position_1 - position_2) > abs(position_1_next - position_2):	#The shortest distance = same nucleosom
	    nucl_state_1_reverse[j][0] = 1						#Loss for 2 at position_1
	    
	    if position_1_next + 30 >= position_2 and position_2 >= position_1_next - 30:
	      nucl_state_1_reverse[j+1][2] = 1						#Occupancy
	      nucl_state_2_reverse[k][2] = 1						#Occupancy
	  
	    if position_1_next + 30 < position_2 or position_2 < position_1_next - 30: 
	      nucl_state_1_reverse[j+1][1] = position_2 - position_1_next		#Shift 2 in comparaison of 1
	      nucl_state_2_reverse[k][1] = position_1_next - position_2		#Shift 1 in comparaison of 2
	      
	    #Go to the next for peak_1 and peak_2
	    j += 2 
	    if k < peaks_2_position_length - 1:
	      k += 1
	    continue
	  
	if j == peaks_1_position_length - 1:

	  if position_1 + 30 >= position_2 and position_2 >= position_1 - 30:
	    nucl_state_1_reverse[j][2] = 1						#Occupancy
	    nucl_state_2_reverse[k][2] = 1						#Occupancy
	
	  if position_1 + 30 < position_2 or position_2 < position_1 - 30: 
	    nucl_state_1_reverse[j][1] = position_2 - position_1			#Shift 2 in comparaison of 1
	    nucl_state_2_reverse[k][1] = position_1 - position_2			#Shift 1 in comparaison of 2
	    
	  #Go to the next for peak_1 and peak_2
	  j += 1 
	  if k < peaks_2_position_length - 1:
	    k += 1
	  continue
    
      if int(100*float(j+1)/float(peaks_1_position_length-1))%(10) == 0 and print_once < int(100*float(j+1)/float(peaks_1_position_length-1)):
	print_once = int(100*float(j+1)/float(peaks_1_position_length-1))
	stdout.write("\r%d /100" % (100*float(j+1)/float(peaks_1_position_length-1)))
	stdout.flush()
  
  
      print ('\nDONE!\n')
  
  else:
    nucl_state_1_reverse = np.zeros((peaks_1_position_length,4))
    nucl_state_2_reverse = np.zeros((peaks_2_position_length,4))
  
  histogram_shift = np.zeros(160)
  global_gain_1 = 0
  global_gain_2 = 0
  global_occupancy = 0
  
  if int(global_signal_var) == 1:

    #We compute the global changes of nucleosom state    
    nucl_state_gain_1 = np.hsplit(nucl_state_1, 4)[0]
    nucl_state_gain_2 = np.hsplit(nucl_state_2, 4)[0]
    nucl_state_shift_1 = np.hsplit(nucl_state_1, 4)[1]
    nucl_state_occupancy_2 = np.hsplit(nucl_state_2, 4)[2]
    if reverse_var == 1:
      nucl_state_gain_1_reverse = np.hsplit(nucl_state_1_reverse, 4)[0]
      nucl_state_gain_2_reverse = np.hsplit(nucl_state_2_reverse, 4)[0]
      nucl_state_shift_1_reverse = np.hsplit(nucl_state_1_reverse, 4)[1]
      nucl_state_occupancy_2_reverse = np.hsplit(nucl_state_2_reverse, 4)[2]

      #Global gain
      global_gain_1 = (np.sum(nucl_state_gain_1) + np.sum(nucl_state_gain_1_reverse))/2
      global_gain_2 = (np.sum(nucl_state_gain_2) + np.sum(nucl_state_gain_2_reverse))/2
      
      #Global occupancy
      global_occupancy = (np.sum(nucl_state_occupancy_2) + np.sum(nucl_state_occupancy_2_reverse))/2
      
      #Shift histogram
      for k in xrange(peaks_1_position_length):
	if int(nucl_state_shift_1[k]) != 0:
	  histogram_shift[int(nucl_state_shift_1[k]) + 80] += 1
	  histogram_shift[int(nucl_state_shift_1_reverse[k]) + 80] += 1
    else:
      #Global gain
      global_gain_1 = np.sum(nucl_state_gain_1)
      global_gain_2 = np.sum(nucl_state_gain_2)
      
      #Global occupancy
      global_occupancy = np.sum(nucl_state_occupancy_2)
      
      #Shift histogram
      for k in xrange(peaks_1_position_length):
	if int(nucl_state_shift_1[k]) != 0:
	  histogram_shift[int(nucl_state_shift_1[k]) + 80] += 1
    
  if local_signal_var == 1:
    nucleosom_match_local_plot(histogram_gaussian_1,histogram_gaussian_2,nucl_state_1,nucl_state_2,nucl_state_1_reverse,nucl_state_2_reverse,peaks_1_position_length,peaks_2_position_length,name_data_1,name_data_2,data_min,data_max,reverse_var,sub_directory,save_txt_var)
  
  #We free manually the space (don't trust python!)
  del histogram_gaussian_1
  del histogram_gaussian_2
  
  return histogram_shift, global_gain_2, global_gain_1, global_occupancy

  
def nucleosom_match_pieces(data_file_1,data_file_2,name_data_1,name_data_2,data_min,data_max,save_txt_var,local_signal_var,global_signal_var,reverse_var,sub_directory):
  
  if (data_file_1 == ''):
    print("WARNING: No input file ")
    return 0
    
  if (data_file_2 == ''):
    print("WARNING: No input file ")
    return 0
    
  if  data_max < data_min :
    print("WARNING: data max < data min")
    return 0

  #Analysing the size of the file
  print('\nDetermining the length of the Data ...')
  data_length, data_min_1, data_min_2, data_max_1, data_max_2 = input_file_comparison(data_file_1,data_file_2,data_min,data_max)
  print('DONE!\n')
  
  #Parameters
  print("\n--------------- Parameters ---------------------\n\nData Length = {0}\nData Min 1 = {1}\nData Max 1 = {2}\nData Min 2 = {3}\nData Max 2 = {4}\n----------------------------------------------\n").format(data_length,data_min_1,data_max_1,data_min_2,data_max_2)
  
  #Calculating the number of pieces
  memory_available = int(psutil.virtual_memory().available)
  split = int(4 * data_length/(memory_available*0.5/100))*2
  #split = int(data_length / 50000)	#We divise in small blocs in order to adapt the peak selection to the nearest neighbour (5e7 for example)
  if split <= 0:
    split = 1

  #print('\n-------- Memory Information ------------\n')
  
  #print ('Avalaible Memory = {0} GB'.format(int(memory_available/10e8)))
  #print ('Number of Pieces = {0}'.format(split))
  #print ('Avalaible Memory for Each Pieces = {0} MB'.format(int(memory_available/(split*10e5))))
  
  #print('\n----------------------------------------\n')
  
  #Launching the analysis pieces by pieces
  
  print('\n-------- Analysis Launched! ------------\n')
  
  #Some needed variables
  piece_length = int(data_length/split)
  
  histogram_shift = np.zeros(160)
  gain = 0
  loss = 0
  occupancy = 0
  
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
    
    data_file_min_1 = data_min_1  + i*piece_length
    data_file_max_1 = data_min_1 + (i+1)*piece_length
    
    data_file_min_2 = data_min_2  + i*piece_length
    data_file_max_2 = data_min_2 + (i+1)*piece_length
    
    
    #Loading data_file
    histogram_1, histogram_length_1 = load_data(data_file_1,data_file_min_1,data_file_max_1)
    histogram_2, histogram_length_2 = load_data(data_file_2,data_file_min_2,data_file_max_2)
    
    #In case there is text in the middle of the data, it change the length of histograms
    #And they need to be the same in order to be compared
    if histogram_length_1 < histogram_length_2:
      histogram_length = histogram_length_1
      histogram_2 = histogram_2[0:histogram_length]
      data_file_max_1 = data_file_min_1 + histogram_length
      data_file_max_2 = data_file_min_2 + histogram_length
      
      if histogram_length == 0:
	print("The file 1 reach it's end: End of the analysis!")
    
    if histogram_length_1 > histogram_length_2:
      histogram_length = histogram_length_2
      histogram_1 = histogram_1[0:histogram_length]
      data_file_max_1 = data_file_min_1 + histogram_length
      data_file_max_2 = data_file_min_2 + histogram_length
      
      if histogram_length == 0:
	print("The file 2 reach it's end: End of the analysis!")
    
    if histogram_length_1 == histogram_length_2:
      histogram_length = histogram_length_2
    
    #Nucleosom matching
    histogram_shift_current, gain_current, loss_current, occupancy_current = nucleosom_match(histogram_1,histogram_2,histogram_length,name_data_1,name_data_2,data_file_min_1,data_file_max_1,save_txt_var,local_signal_var,global_signal_var,reverse_var,sub_directory)
    
    #Adding computed histograms to older histograms
    histogram_shift = np.add(histogram_shift,histogram_shift_current)
    gain += gain_current
    loss += loss_current
    occupancy += occupancy_current
    
  if data_length % split > 1000:
    print('\n-------- Computing Last Piece ------------\n')
    
    data_file_min_1 = data_min_1  + piece_length*split
    data_file_max_1 = data_max_1
    
    data_file_min_2 = data_min_2  + piece_length*split
    data_file_max_2 = data_max_2
    
     #Loading data_file
    histogram_1, histogram_length_1 = load_data(data_file_1,data_file_min_1,data_file_max_1)
    histogram_2, histogram_length_2 = load_data(data_file_2,data_file_min_2,data_file_max_2)
    
    #In case there is text in the middle of the data, it change the length of histograms
    #And they need to be the same in order to be compared
    if histogram_length_1 < histogram_length_2:
      histogram_length = histogram_length_1
      histogram_2 = histogram_2[0:histogram_length]
      data_file_max_1 = data_file_min_1 + histogram_length
      data_file_max_2 = data_file_min_2 + histogram_length
    
    if histogram_length_1 > histogram_length_2:
      histogram_length = histogram_length_2
      histogram_1 = histogram_1[0:histogram_length]
      data_file_max_1 = data_file_min_1 + histogram_length
      data_file_max_2 = data_file_min_2 + histogram_length
    
    if histogram_length_1 == histogram_length_2:
      histogram_length = histogram_length_2
    
    #Nucleosom matching
    histogram_shift_current, gain_current, loss_current, occupancy_current = nucleosom_match(histogram_1,histogram_2,histogram_length,name_data_1,name_data_2,data_file_min_1,data_file_max_1,save_txt_var,local_signal_var,global_signal_var,reverse_var,sub_directory)
    
    #Adding computed histograms to older histograms
    histogram_shift = np.add(histogram_shift,histogram_shift_current)
    gain += gain_current
    loss += loss_current
    occupancy += occupancy_current

  print('\n-------- Analysis DONE! ------------\n')
  
  print('Normalization of Histograms')
  histogram_shift_sum = np.sum(histogram_shift)
  if histogram_shift_sum != 0:
    histogram_shift = histogram_shift / histogram_shift_sum
  print('DONE!\n')
  
  if int(global_signal_var) == 1:
    number_nucleosom_detected = histogram_shift_sum + gain + occupancy
    nucleosom_match_global_plot(number_nucleosom_detected,gain,loss,occupancy,histogram_shift_sum,histogram_shift,name_data_2,data_min,data_max,sub_directory,save_txt_var)
  
  return 0
