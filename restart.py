#!/bin/python

import os
import shutil
import sys

import numpy as np

def load_restart_input(data_file):
  print ("\nLoading Input Data for Restart...")
  #Repertoire dans lequel sont les donnes
  directory = os.path.dirname(data_file)
  
  filename = os.path.split(data_file)[1]
  filename_split = filename.split('.')[0]
  sub_directory = directory +'/'+ filename_split + "/RESTART"
  
  if not os.path.exists(sub_directory):
    print("WARNING: No Restart File !")
    return 0
  
  restart_file_input = sub_directory + "/RESTART_input.txt"

  file_input = np.loadtxt(restart_file_input)
  
  print("DONE!\n")
  
  return file_input

def load_restart_histograms(data_file):
  
  print ("\nLoading Histograms for Restart...")
  
  #Repertoire dans lequel sont les donnes
  directory = os.path.dirname(data_file)
  
  filename = os.path.split(data_file)[1]
  filename_split = filename.split('.')[0]
  sub_directory = directory +'/'+ filename_split + "/RESTART"
  
  if not os.path.exists(sub_directory):
    sys.exit("WARNING: No Restart File !")
  
  restart_file_histogram_1 = sub_directory  + "/RESTART_histogram_1.txt"
  restart_file_histogram_2 = sub_directory  + "/RESTART_histogram_2.txt"
  restart_file_histogram_3 = sub_directory  + "/RESTART_histogram_3.txt"
  restart_file_width_peak = sub_directory + "/RESTART_width_peak.txt"
  restart_file_width_hole = sub_directory + "/RESTART_width_hole.txt"
  
  file_histogram_1 = np.loadtxt(restart_file_histogram_1)
  file_histogram_2 = np.loadtxt(restart_file_histogram_2)
  file_histogram_3 = np.loadtxt(restart_file_histogram_3)
  file_width_peak = np.loadtxt(restart_file_width_peak)
  file_width_hole = np.loadtxt(restart_file_width_hole)
  
  print("DONE!\n")
  
  return file_histogram_1,file_histogram_2,file_histogram_3,file_width_peak,file_width_hole
  
def save_restart(data_file,input_data,histogram_period_1,histogram_period_2,histogram_period_3,histogram_width,histogram_width_hole):
  
  print ("\nSaving Checkpoint in case of crash...")
  
  #Repertoire dans lequel sont les donnes
  directory = os.path.dirname(data_file)
  
  filename = os.path.split(data_file)[1]
  filename_split = filename.split('.')[0]
  sub_directory = directory +'/'+ filename_split + "/RESTART"
    
  if not os.path.exists(sub_directory):
    os.mkdir(sub_directory)
  
  restart_file_input = sub_directory + "/RESTART_input.txt"
  restart_file_histogram_1 = sub_directory + "/RESTART_histogram_1.txt"
  restart_file_histogram_2 = sub_directory  + "/RESTART_histogram_2.txt"
  restart_file_histogram_3 = sub_directory  + "/RESTART_histogram_3.txt"
  restart_file_width_peak = sub_directory + "/RESTART_width_peak.txt"
  restart_file_width_hole = sub_directory + "/RESTART_width_hole.txt"

  np.savetxt(restart_file_input,input_data)
  np.savetxt(restart_file_histogram_1,histogram_period_1)
  np.savetxt(restart_file_histogram_2,histogram_period_2)
  np.savetxt(restart_file_histogram_3,histogram_period_3)
  np.savetxt(restart_file_width_peak,histogram_width)
  np.savetxt(restart_file_width_hole,histogram_width_hole)
  
  print("DONE!\n")
  
def delete_restart(data_file):
    
  print ("\nDeleting Checkpoint files...")
  
  #Repertoire dans lequel sont les donnes
  directory = os.path.dirname(data_file)
  
  filename = os.path.split(data_file)[1]
  filename_split = filename.split('.')[0]
  sub_directory = directory +'/'+ filename_split + "/RESTART"
  
  #We delete the RESTART File
  shutil.rmtree(sub_directory)
  
  print("DONE!\n")
  
