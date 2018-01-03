#!/bin/python

import Tkinter as tk
import tkFileDialog
import numpy as np
import os
import re


from wavelet import *
from wavelet_plot import *
from input_file_analysis import *
from restart import *
from plot_histogram import *

  
class Window(tk.Frame):
  def __init__(self):
    tk.Frame.__init__(self)
    self.master.title("Nucleosom Analysis")
    self.master.columnconfigure(0, weight=1)
    self.master.rowconfigure(0, weight=1)
    self.columnconfigure(0, weight=1)
    self.rowconfigure(0, weight=1)
    self.grid(sticky="NSEW")
    self.createWidgets()
    

  def createWidgets(self):
    #creation des objets variables
    self.data_file_var = tk.StringVar()
    self.data_file_var.set("")

    self.data_min_var = tk.StringVar()
    self.data_min_var.set("0")
    
    self.data_max_var = tk.StringVar()
    self.data_max_var.set("0")
    
    self.period_min_var = tk.StringVar()
    self.period_min_var.set("")
    
    self.period_max_var = tk.StringVar()
    self.period_max_var.set("") 
    
    self.bin_period_var = tk.StringVar()
    self.bin_period_var.set("1")
        
    self.width_min_var = tk.StringVar()
    self.width_min_var.set("")
    
    self.width_max_var = tk.StringVar()
    self.width_max_var.set("")
    
    self.hole_width_min_var = tk.StringVar()
    self.hole_width_min_var.set("")
    
    self.hole_width_max_var = tk.StringVar()
    self.hole_width_max_var.set("")
    
    self.launch_all_var = tk.StringVar()
    self.save_var = tk.StringVar()
    self.save_txt_var = tk.StringVar()
    self.plot_var = tk.StringVar()
    self.local_signal_var = tk.StringVar()
    self.global_signal_var = tk.StringVar()
    
    self.restart_var = tk.StringVar()
    self.restart_var.set("0")
    
    
    #creation des widgets
    mainFrame = tk.Frame(self, borderwidth=2, relief="groove")
    
    launch_all_label = tk.Label(mainFrame, text="Compute all files in directory")
    data_min_label = tk.Label(mainFrame, text="Data Min (in bp)\n(if 0->min) ")
    data_max_label = tk.Label(mainFrame, text="Data Max (in bp)\n(if 0->max) ")
    period_min_label = tk.Label(mainFrame, text="Period Min (in bp) ")
    period_max_label = tk.Label(mainFrame, text="Period Max (in bp)")
    bin_period_label = tk.Label(mainFrame, text="Period Bin length") 
    width_min_label = tk.Label(mainFrame, text="Width Min (in bp) ")
    width_max_label = tk.Label(mainFrame, text="Width Max (in bp)")
    hole_width_min_label = tk.Label(mainFrame, text="Hole Width Min (in bp) ")
    hole_width_max_label = tk.Label(mainFrame, text="Hole Width Max (in bp)")
    plot_save_label = tk.Label(mainFrame, text="Display OR Save Plots:") 
    save_txt_label = tk.Label(mainFrame, text="Save Data into Text:") 
    local_global_label = tk.Label(mainFrame, text="Local Signal\nAND/OR Histograms") 
    

    data_file_entry = tk.Button(mainFrame, text='Data File ( histogram of nucleosom position only)', command=self.askopenfilename)
    data_min_entry = tk.Entry(mainFrame, textvariable=self.data_min_var)
    data_max_entry = tk.Entry(mainFrame, textvariable=self.data_max_var)
    period_min_entry = tk.Entry(mainFrame, textvariable=self.period_min_var)
    period_max_entry = tk.Entry(mainFrame, textvariable=self.period_max_var)
    bin_period_entry = tk.Entry(mainFrame, textvariable=self.bin_period_var)
    width_min_entry = tk.Entry(mainFrame, textvariable=self.width_min_var)
    width_max_entry = tk.Entry(mainFrame, textvariable=self.width_max_var)
    hole_width_min_entry = tk.Entry(mainFrame, textvariable=self.hole_width_min_var)
    hole_width_max_entry = tk.Entry(mainFrame, textvariable=self.hole_width_max_var)
    
    button = tk.Button(mainFrame, text="Launch", command=self.action)
    button_restart = tk.Button(mainFrame, text="Restart", command=self.restart)
    
    launch_all_checkbutton = tk.Checkbutton(mainFrame,   text="Activate       ",variable = self.launch_all_var)
    save_checkbutton = tk.Checkbutton(mainFrame,         text="Save Plots    ",variable = self.save_var)
    save_txt_checkbutton = tk.Checkbutton(mainFrame,     text="Save Text    ",variable = self.save_txt_var)
    plot_checkbutton = tk.Checkbutton(mainFrame,         text="Display Plot  ",variable = self.plot_var)
    local_signal_checkbutton = tk.Checkbutton(mainFrame, text="Local Signal ",variable = self.local_signal_var)
    global_signal_checkbutton = tk.Checkbutton(mainFrame,text="Histograms  ",variable = self.global_signal_var)
    
    
    # define options for opening or saving a file
    self.file_opt = options = {}
    options['defaultextension'] = '.wig'
    options['filetypes'] = [('all files', '.*'), ('wig files', '.wig')]    

    #position des widgets
    mainFrame.grid(column=0, columnspan=2, row=0, sticky="NSEW")
    
    data_file_entry.grid(column=0,columnspan=2, row=1, sticky="EW")
    
    launch_all_label.grid(column=0, row=2, sticky="EW")
    launch_all_checkbutton.grid(column=1, row=2, sticky="EW")
    
    data_min_label.grid(column=0, row=4, sticky="EW")
    data_min_entry.grid(column=1, row=4, sticky="EW")
    
    data_max_label.grid(column=0, row=5, sticky="EW")
    data_max_entry.grid(column=1, row=5, sticky="EW")
    
    period_min_label.grid(column=0, row=6, sticky="EW")
    period_min_entry.grid(column=1, row=6, sticky="EW")
    
    period_max_label.grid(column=0, row=7, sticky="EW")
    period_max_entry.grid(column=1, row=7, sticky="EW")
    
    bin_period_label.grid(column=0, row=8, sticky="EW")
    bin_period_entry.grid(column=1, row=8, sticky="EW")
    
    width_min_label.grid(column=0, row=9, sticky="EW")
    width_min_entry.grid(column=1, row=9, sticky="EW")
    
    width_max_label.grid(column=0, row=10, sticky="EW")
    width_max_entry.grid(column=1, row=10, sticky="EW")
    
    hole_width_min_label.grid(column=0, row=11, sticky="EW")
    hole_width_min_entry.grid(column=1, row=11, sticky="EW")
    
    hole_width_max_label.grid(column=0, row=12, sticky="EW")
    hole_width_max_entry.grid(column=1, row=12, sticky="EW")
    
    plot_save_label.grid(column=0,rowspan=2, row=13, sticky="EW")
    save_checkbutton.grid(column=1, row=13, sticky="EW")
    plot_checkbutton.grid(column=1, row=14, sticky="EW")
    
    
    save_txt_label.grid(column=0,rowspan=1, row=15, sticky="EW")
    save_txt_checkbutton.grid(column=1, row=15, sticky="EW")
    
    local_global_label.grid(column=0,rowspan=2, row=16, sticky="EW")
    local_signal_checkbutton.grid(column=1, row=16, sticky="EW")
    global_signal_checkbutton.grid(column=1, row=17, sticky="EW")
    
    button.grid(column=1, columnspan=1, row=19,rowspan=1, sticky="NSEW")
    button_restart.grid(column=0, columnspan=1, row=19, sticky="NSEW")
    

    
  def action(self):
     
    #In case we want to restart the program to the last checkpoint
    print self.restart_var.get()
    if int(self.restart_var.get()) == 1:
      file_input = load_restart_input(self.data_file_var.get())
      
      self.launch_all_var.set("{0}".format(int(file_input[0])))
      self.data_min_var.set("{0}".format(int(file_input[1])))
      self.data_max_var.set("{0}".format(int(file_input[2])))
      self.period_min_var.set("{0}".format(int(file_input[3])))
      self.period_max_var.set("{0}".format(int(file_input[4])))
      self.width_min_var.set("{0}".format(int(file_input[5])))
      self.width_max_var.set("{0}".format(int(file_input[6])))
      self.bin_period_var.set("{0}".format(int(file_input[7])))
      self.plot_var.set("{0}".format(int(file_input[8])))
      self.save_var.set("{0}".format(int(file_input[9])))
      self.save_txt_var.set("{0}".format(int(file_input[10])))
      self.local_signal_var.set("{0}".format(int(file_input[11])))
      self.global_signal_var.set("{0}".format(int(file_input[12])))
      
      #Load Histograms
      period_length = int((int(self.period_max_var.get()) - int(self.period_min_var.get()))/int(self.bin_period_var.get()))
      width_length = int(int(self.width_max_var.get()) - int(self.width_min_var.get()))
      
      histogram_period_1 = np.zeros(period_length)
      histogram_period_2 = np.zeros(period_length)
      histogram_period_3 = np.zeros(period_length)
      histogram_period_width = np.zeros(width_length)
      
      histogram_period_1,histogram_period_2,histogram_period_3,histogram_width = load_restart_histograms(self.data_file_var.get())
      
    if int(self.restart_var.get()) == 0:
      #We set the histogram to zero per default
      period_length = int((int(self.period_max_var.get()) - int(self.period_min_var.get()))/int(self.bin_period_var.get()))
      width_length = int(int(self.width_max_var.get()) - int(self.width_min_var.get()))
      
      histogram_period_1 = np.zeros(period_length)
      histogram_period_2 = np.zeros(period_length)
      histogram_period_3 = np.zeros(period_length)
      histogram_width = np.zeros(width_length)
      
    print('\n------------------ Data File -----------------------\n')
    
    
    #Repertoire dans lequel sont les donnes
    directory = os.path.dirname(self.data_file_var.get())
    

    if int(self.launch_all_var.get()) == 0:
     
      #Creating new directory and executing the program in it if saving
      filename = os.path.split(self.data_file_var.get())[1]
      print ('Data File: '+ filename)
      filename_split = filename.split('.')[0]
      sub_directory = directory +'/'+ filename_split
      
      if not os.path.exists(sub_directory) and int(self.save_var.get()) == 1:
	os.mkdir(sub_directory)
	#Changing Current directory
	os.chdir(sub_directory)
      
      print ( 'Directory Path: ' + directory )
      print('\n-----------------------------------------------------\n')
      
      wavelet_analysis_pieces(self.data_file_var.get(),int(self.launch_all_var.get()),int(self.data_min_var.get()),int(self.data_max_var.get()),int(self.period_min_var.get()),int(self.period_max_var.get()),int(self.width_min_var.get()),int(self.width_max_var.get()),int(self.hole_width_min_var.get()),int(self.hole_width_max_var.get()),int(self.bin_period_var.get()),int(self.plot_var.get()),int(self.save_var.get()),int(self.save_txt_var.get()),int(self.local_signal_var.get()),int(self.global_signal_var.get()),int(self.restart_var.get()),histogram_period_1,histogram_period_2,histogram_period_3,histogram_width,sub_directory)
    
    if int(self.launch_all_var.get()) == 1:
      
      for filename in os.listdir(directory):
	
	regexp = re.compile(r'.wig')
	if regexp.search(filename) is not None:
	  print ('Data File: '+ filename)
	  
	  #Creating new directory and executing the program in it if saving
	  filename_split = filename.split('.')[0]
	  sub_directory = directory +'/'+ filename_split
	  if not os.path.exists(sub_directory) and int(self.save_var.get()) == 1:
	    os.mkdir(sub_directory)
	    #Changing Current directory
	    os.chdir(sub_directory)
	    
	  print ( 'Directory Path: ' + directory )
	  print('\n-----------------------------------------------------\n')
	  
	  #The analysis need the absolute path in order to open data file
	  filename_absolute_path = directory + '/'+ filename
	
	  wavelet_analysis_pieces(filename_absolute_path,int(self.launch_all_var.get()),int(self.data_min_var.get()),int(self.data_max_var.get()),int(self.period_min_var.get()),int(self.period_max_var.get()),int(self.width_min_var.get()),int(self.width_max_var.get()),int(self.bin_period_var.get()),int(self.plot_var.get()),int(self.save_var.get()),int(self.save_txt_var.get()),int(self.local_signal_var.get()),int(self.global_signal_var.get()),int(self.restart_var.get()),histogram_period_1,histogram_period_2,histogram_period_3,histogram_width,sub_directory)
	      
	else:
	  print (filename + ' : Not taken into account')  
	  
  def restart(self):
    self.restart_var.set("1")
    self.action()

  #Browser for input data file
  def askopenfilename(self):
    filename = tkFileDialog.askopenfilename(**self.file_opt)
    self.data_file_var.set(filename)
  
if __name__ =='__main__':
  Window().mainloop()
   
