#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 10:28:13 2021

@author: dla


SCRIPT PRACTICA ALFA MASTER, PASADO DE MATLAB A PYTHON
y
"""

#reset to manually clear all the variables
#clear               #to clear the command windows
#%reset -f          #to clear all the variables without confirmation
#magic('reset -sf')


#######General packages useful##

import matplotlib.pyplot as plt  #for simplicity, to not write matplotlib.pyplot
        #everytime we want to plot something

import math 
#from scipy.stats import norm               ##norm.fit() fit to gaussian
import numpy as np
    #np contain linspaces as np.linspace(a,b,N)
import csv
        
from plotly.graph_objs import Bar, Layout
from plotly import offline

######3

#plt.close("all")

#%%
#########################1), Data loading #####################3
#The files to load are in txt. The best way to read is:


##Data of the natural signal, from the oscilloscope, in.csv. 

voltage_stored = []
time_stored = []


#################################
###########RAW##################3
################################


with open('TEK0001.CSV') as file_object:            #LYSO, raw
            
        reader = csv.reader(file_object) #reader object assoaciated with the 
        #filerow_count = sum(1 for row in reader)            #number of rows
        #header_row = next(reader) #next return the next line. Since we only call it
        #once, we only get the 1st line
        #print(header_row)
        #n_columns = len(header_row)     #number of columns
    
    #Storing of the voltage
        time_help = []
        voltage_help = []
        for row in reader:
            voltage_help.append(float(row[-1 -1]))     #voltage; -1 = last , so -1 -1 is the
                            #second to last
            time_help.append(float(row[-1-2]))     #time

        #Storing of the final values
        voltage_stored.append(voltage_help)
        time_stored.append(time_help)
              
# #Plot

plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
plt.plot([1e6 * x for x in time_stored[-1]],[1e3 * x for x in voltage_stored[-1]], 'bo-')    #-1 chooses last element, which is the
        #one that have been added to the lsit the lastest ;)    
        #widht so that each bar touches each other!
plt.title("Raw Waveform of Cs137 with LYSO", fontsize=22)           #title
plt.xlabel("time [us]", fontsize=14)                        #xlabel
plt.ylabel("voltage [mV]", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
#plt.xlim(0,max(ADC_channel))                       #limits of x axis     
plt.savefig('Raw_signal_LYSO.png', format='png')




with open('TEK0005.CSV') as file_object:            #CsI, raw
            
        reader = csv.reader(file_object) 
    #Storing of the voltage
        time_help = []
        voltage_help = []
        for row in reader:
            voltage_help.append(float(row[-1 -1]))     #voltage; -1 = last , so -1 -1 is the
                            #second to last
            time_help.append(float(row[-1-2]))     #time

        #Storing of the final values
        voltage_stored.append(voltage_help)
        time_stored.append(time_help)
              
# #Plot

plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
plt.plot([1e6 * x for x in time_stored[-1]],[1e3 * x for x in voltage_stored[-1]], 'bo-')    #-1 chooses last element, which is the
        #widht so that each bar touches each other!
plt.title("Raw waveform of Cs137 with CsI", fontsize=22)           #title
plt.xlabel("time [us]", fontsize=14)                        #xlabel
plt.ylabel("voltage [mV]", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
#plt.xlim(0,max(ADC_channel))                       #limits of x axis     
plt.savefig('Raw_signal_CsI.png', format='png')





with open('TEK0007.CSV') as file_object:            #BGo, raw
            
        reader = csv.reader(file_object)
    #Storing of the voltage
        time_help = []
        voltage_help = []
        for row in reader:
            voltage_help.append(float(row[-1 -1]))     #voltage; -1 = last , so -1 -1 is the
                            #second to last
            time_help.append(float(row[-1-2]))     #time

        #Storing of the final values
        voltage_stored.append(voltage_help)
        time_stored.append(time_help)
              
# #Plot

plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
plt.plot([1e6 * x for x in time_stored[-1]],[1e3 * x for x in voltage_stored[-1]], 'bo-')    #-1 chooses last element, which is the
        #widht so that each bar touches each other!
plt.title("Raw waveform of Cs137 with BGO", fontsize=22)           #title
plt.xlabel("time [us]", fontsize=14)                        #xlabel
plt.ylabel("voltage [mV]", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
#plt.xlim(0,max(ADC_channel))                       #limits of x axis     
plt.savefig('Raw_signal_BGO.png', format='png')


##########3
############Waves from preampl############################
##########################
#The new data comes from split, ensuring with the threshold that the data are mostly
#from the gamma peak.


#with open('TEK0002.CSV') as file_object:            #LYSO, pre
with open('TEK0002_LYSO_pre.CSV') as file_object:           #the new data         
 
        reader = csv.reader(file_object)
    
    #Storing of the voltage
        time_help = []
        voltage_help = []
        for row in reader:
            voltage_help.append(float(row[-1 -1]))     #voltage; -1 = last , so -1 -1 is the
                            #second to last
            time_help.append(float(row[-1-2]))     #time

        #Storing of the final values
        voltage_stored.append(voltage_help)
        time_stored.append(time_help)
              
# #Plot

plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
plt.plot([1e6 * x for x in time_stored[-1]],[1e3 * x for x in voltage_stored[-1]], 'bo-')    #-1 chooses last element, which is the
        #widht so that each bar touches each other!
plt.title("Waveform (from Pre) of Cs137 with LYSO", fontsize=22)           #title
plt.xlabel("time [us]", fontsize=14)                        #xlabel
plt.ylabel("voltage [mV]", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
#plt.xlim(0,max(ADC_channel))                       #limits of x axis     
plt.savefig('Signal_LYSO.png', format='png')


#with open('TEK0006.CSV') as file_object:            #Csi, pre
with open('TEK0001_CsI_pre.CSV') as file_object:           #the new data          
  
        reader = csv.reader(file_object) 
    #Storing of the voltage
        time_help = []
        voltage_help = []
        for row in reader:
            voltage_help.append(float(row[-1 -1]))     #voltage; -1 = last , so -1 -1 is the
                            #second to last
            time_help.append(float(row[-1-2]))     #time

        #Storing of the final values
        voltage_stored.append(voltage_help)
        time_stored.append(time_help)
              
# #Plot

plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
plt.plot([1e6 * x for x in time_stored[-1]],[1e3 * x for x in voltage_stored[-1]], 'bo-')    #-1 chooses last element, which is the
        #widht so that each bar touches each other!
plt.title("Waveform (pre) of Cs137 with CsI", fontsize=22)           #title
plt.xlabel("time [us]", fontsize=14)                        #xlabel
plt.ylabel("voltage [mV]", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
#plt.xlim(0,max(ADC_channel))                       #limits of x axis     
plt.savefig('Signal_CsI.png', format='png')

###########
###########

#with open('TEK0008.CSV') as file_object:            #BGo, pre
with open('TEK0000_BGO_pre.CSV') as file_object:           #the new data
          
        reader = csv.reader(file_object) 
    
    #Storing of the voltage
        time_help = []
        voltage_help = []
        for row in reader:
            voltage_help.append(float(row[-1 -1]))     #voltage; -1 = last , so -1 -1 is the
                            #second to last
            time_help.append(float(row[-1-2]))     #time

        #Storing of the final values
        voltage_stored.append(voltage_help)
        time_stored.append(time_help)
              
# #Plot

plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
plt.plot([1e6 * x for x in time_stored[-1]],[1e3 * x for x in voltage_stored[-1]], 'bo-')    #-1 chooses last element, which is the
        #widht so that each bar touches each other!
plt.title("Waveform (pre) of Cs137 with BGO", fontsize=22)           #title
plt.xlabel("time [us]", fontsize=14)                        #xlabel
plt.ylabel("voltage [mV]", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
#plt.xlim(0,max(ADC_channel))                       #limits of x axis     
plt.savefig('Signal_BGO.png', format='png')





#%% 2) Calcs of the light yield

#Since the signals from the Pre are way better, will use them. Remember they are
#the last stored (6 total stored), so from 3 to 5 (0 the first). LYSO, CsI, BGO the order


####LYSO###########3

#TO find the peak I could easily do:
peak = min(voltage_stored[3])                 #Peak value
index_peak = voltage_stored[3].index(peak)        #Index of the peak value

#using the index of the peak I see the interval by looking at the .csv file
#and the plot
index_min = 499         #index that starts the peak (by looking the graph)
index_max = 820#774         #index that ends the peak (by looking the graph)


voltage_peak_neg = voltage_stored[3][index_min-1:index_max-1]  
time_peak = time_stored[3][index_min-1:index_max-1]
voltage_peak = [abs(x) for x in voltage_peak_neg]     #abs() because this
        #peak contains both > and ,0 values, so to add them in order to count
        #them, I have to put everything in the positive value
len_peak =  len(voltage_stored[3][index_min-1:index_max-1])           #len of the peak   

delta_V = 1/5 *100e-3           #[V] error of the voltage emasurements, from the scale of
                            #the oscilloscope, which come from the photos
        
sum_voltage = sum(voltage_peak)         #[V] total voltage of the peak

#Storing of the current

voltage_peak_stored = []        #storage of the total voltage of the peak
n_elements_peak = []             #this will store the number of voltages I sum,
                #for each peak, for the error calc
delta_single_V_measurement = []         #storage of the error of the voltage 
                    #measurements, for the error calc
                    
voltage_peak_stored.append(sum_voltage)
n_elements_peak.append(len_peak)
delta_single_V_measurement.append(delta_V)


#Plot (debug)
#plt.plot([1e6 * x for x in time_stored[-1]],[1e3 * x for x in voltage_stored[-1]], 'bo-')    #-1 chooses last element, which is the

plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
plt.plot([1e6 * x for x in time_stored[3]],[1e3 * x for x in voltage_stored[3]], 'b.-')    
plt.plot([1e6 * x for x in time_peak],[1e3 * x for x in voltage_peak_neg], 'r.-')   
        #widht so that each bar touches each other!
plt.title("Waveform of Cs137 with LYSO", fontsize=22)           #title
plt.xlabel("time [us]", fontsize=14)                        #xlabel
plt.ylabel("voltage [mV]", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.savefig('Signal_LYSO_sum.png', format='png')

####CsI##########

#TO find the peak I could easily do:
peak = min(voltage_stored[4])                 #Peak value
index_peak = voltage_stored[4].index(peak)        #Index of the peak value

#using the index of the peak I see the interval by looking at the .csv file
#and the plot

index_min = 542
index_max = 2100#1921         #index that ends the peak (by looking the graph)

voltage_peak_neg = voltage_stored[4][index_min-1:index_max-1]  
time_peak = time_stored[4][index_min-1:index_max-1]
voltage_peak = [abs(x) for x in voltage_peak_neg]     #abs() because this
        #peak contains both > and ,0 values, so to add them in order to count
        #them, I have to put everything in the positive value
len_peak =  len(voltage_stored[4][index_min-1:index_max-1])           #len of the peak  
delta_V = 1/5 *20e-3                        #[V] 

sum_voltage = sum(voltage_peak)         #[V] total voltage of the peak

#Storing of the current

voltage_peak_stored.append(sum_voltage)
n_elements_peak.append(len_peak)
delta_single_V_measurement.append(delta_V)


#Plot (debug)
plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
plt.plot([1e6 * x for x in time_stored[4]],[1e3 * x for x in voltage_stored[4]], 'b.-')    
plt.plot([1e6 * x for x in time_peak],[1e3 * x for x in voltage_peak_neg], 'r.-')   
        #widht so that each bar touches each other!
plt.title("Waveform of Cs137 with CsI", fontsize=22)           #title
plt.xlabel("time [us]", fontsize=14)                        #xlabel
plt.ylabel("voltage [mV]", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.savefig('Signal_CsI_sum.png', format='png')

###BGO###########

#TO find the peak I could easily do:
peak = min(voltage_stored[5])                 #Peak value
index_peak = voltage_stored[5].index(peak)        #Index of the peak value

#using the index of the peak I see the interval by looking at the .csv file
#and the plot

index_min = 585
index_max = 1300#1060         #index that ends the peak (by looking the graph)

voltage_peak_neg = voltage_stored[5][index_min-1:index_max-1]  
time_peak = time_stored[5][index_min-1:index_max-1]
voltage_peak = [abs(x) for x in voltage_peak_neg]     #abs() because this
        #peak contains both > and ,0 values, so to add them in order to count
        #them, I have to put everything in the positive value
len_peak =  len(voltage_stored[5][index_min-1:index_max-1])           #len of the peak  
delta_V = 1/5 *20e-3                                     #[V] 

sum_voltage = sum(voltage_peak)         #[V] total voltage of the peak

#Storing of the current

voltage_peak_stored.append(sum_voltage)
n_elements_peak.append(len_peak)
delta_single_V_measurement.append(delta_V)


#Plot (debug)
plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
plt.plot([1e6 * x for x in time_stored[5]],[1e3 * x for x in voltage_stored[5]], 'b.-')    
plt.plot([1e6 * x for x in time_peak],[1e3 * x for x in voltage_peak_neg], 'r.-')   
        #widht so that each bar touches each other!
plt.title("Waveform of Cs137 with BGO", fontsize=22)           #title
plt.xlabel("time (us)", fontsize=14)                        #xlabel
plt.ylabel("voltage (mV)", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.savefig('Signal_BGO_sum.png', format='png')


##########################################################


#So, we have the charge (voltage), which is proportional to the GM cells that
#have been fired, which is also proportional to the number of photons that
#have reach the detector, i.e., the photons generated (yield). Since we do
#not know this proportionaly factor, we could just compute ratio of yields,
#since this will null out this proportionality factors.

delta_voltage_peak = np.multiply(delta_single_V_measurement, 
                                 [np.sqrt(x) for x in n_elements_peak])       
                    #[V] error of the voltage of the peak, the pseudo
                    #light yield


ratios = [voltage_peak_stored[0]/voltage_peak_stored[1],
          voltage_peak_stored[0]/voltage_peak_stored[2],
          voltage_peak_stored[2]/voltage_peak_stored[1],
          ]   #LYSO/CsI, LYSO/BGO, BGO/CsI

#error calculation of the ratio
auxiliar = [(delta_voltage_peak[0]/voltage_peak_stored[0])**2 + (delta_voltage_peak[1]/voltage_peak_stored[1])**2, 
            (delta_voltage_peak[0]/voltage_peak_stored[0])**2 + (delta_voltage_peak[2]/voltage_peak_stored[2])**2,
            (delta_voltage_peak[2]/voltage_peak_stored[2])**2 + (delta_voltage_peak[1]/voltage_peak_stored[1])**2] 
            #1st row: element 0 and 1
            #2nd row: element 0 and 2
            #3rd row: elememtn 1 and 2   (same order as the ratio!!!!)

                    #the thing to be square rooted
delta_ratios = np.multiply(ratios, [np.sqrt(x) for x in auxiliar])
                                             #this will be though

#Print:
print('Light yield ratio LYSO/CsI: ' + str(ratios[0]) + ' +/- ' + str(delta_ratios[0]))
print('Light yield ratio LYSO/BGO: ' + str(ratios[1]) + ' +/- ' + str(delta_ratios[1]))
print('Light yield ratio BGO/CsI: ' + str(ratios[2]) + ' +/- ' + str(delta_ratios[2]))
#LYSO/CsI, LYSO/BGO, CsI/BGO


plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
plt.title("Light yield ratio of the Cs137 peak for several scintillators", fontsize=22, wrap=True)           #title
plt.bar(['LYSO/CsI', 'LYSO/BGO', 'BGO/CsI'], ratios, yerr = delta_ratios, edgecolor="black")#, yerr = delta_light_yield)
#plt.xlabel("ADC channels", fontsize=10)                        #xlabel
plt.ylabel("Light yield ratio", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.savefig('Light_yield_ratio.png', format='png')


plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
plt.title("Pseudo Light yield of the Cs137 peak for several scintillators", fontsize=22, wrap=True)           #title
plt.bar(['LYSO', 'CsI', 'BGO'], voltage_peak_stored, edgecolor="black")#, yerr = delta_light_yield)
#plt.xlabel("ADC channels", fontsize=10)                        #xlabel
plt.ylabel("Pseudo Light yield", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.savefig('Pseudo_light_yield.png', format='png')



#CsI/BGO match,LYSO/CsI almost, LYSO/BGO do, should be greater than 1 an 
    #is lower than 1:(


#################################################
##############################################

#%% Rise and decay time

#Usando misma configuracion, mido el máximo, del pre, q sera el del fotopico. Pre mejor
#q raw pq el raw tenia mas pileups q el pre. Mediciones tanto para la raw signal como para
#la del pre!!!.

#Variables that will store(st) the results 
t_rise_st = []
t_decay_st = []
delta_t_rise_st = []
delta_t_decay_st = []

########LYSO RAW#######

#RAW

t_rise = 18.8e3                            #[ns]
t_decay = 77*1e3                      #[ns], 3.28us (micros )
delta_t_rise = 1/5*5e3                #[ns], error of the oscilloscope
                #it basically is 1/5 of the escale, since each square is divided
                #into 5 intervals, and the scale is the length of the sides of
                #the square
delta_t_decay = 1/5*10e3
  
#storing
t_rise_st.append(t_rise)
t_decay_st.append(t_decay)
delta_t_rise_st.append(delta_t_rise)
delta_t_decay_st.append(delta_t_decay)

########CsI RAW#######


t_rise = 5e3                            #[ns]
t_decay = 104*1e3                      #[ns], 3.28us (micros )
delta_t_rise = 1/5*25e3                   #[ns], error of the oscilloscope
delta_t_decay = 1/5*25e3
  
#storing
t_rise_st.append(t_rise)
t_decay_st.append(t_decay)
delta_t_rise_st.append(delta_t_rise)
delta_t_decay_st.append(delta_t_decay)

########BGO RAW#######

t_rise = 400                           #[ns]
t_decay = 52e3                     #[ns], 3.28us (micros )
delta_t_rise = 1/5*10e3                   #[ns], error of the oscilloscope
delta_t_decay = 1/5*10e3
  
#storing
t_rise_st.append(t_rise)
t_decay_st.append(t_decay)
delta_t_rise_st.append(delta_t_rise)
delta_t_decay_st.append(delta_t_decay)

##################################
#############PRE#################
################################
#This is the data from the new measurements, using the threshold to only measure
#the gamma peak. The old results are commented

########LYSO PRE#######

t_rise =  156  #28                           #[ns]
t_decay = 464  #640                      #[ns], 3.28us (micros )
delta_t_rise = 1/5 * 100 #1/5*100                   #[ns], error of the oscilloscope
delta_t_decay = 1/5 * 100#1/5*100
  
#storing
t_rise_st.append(t_rise)
t_decay_st.append(t_decay)
delta_t_rise_st.append(delta_t_rise)
delta_t_decay_st.append(delta_t_decay)



########CsI PRE#######

t_rise = 320 #360                           #[ns]
t_decay = 2.88e3    #4.52e3                      #[ns], 3.28us (micros )
delta_t_rise = 1/5 * 500#1/5*1e3                  #[ns], error of the oscilloscope
delta_t_decay = 1/5 * 500 #1/5*1e3
  
#storing
t_rise_st.append(t_rise)
t_decay_st.append(t_decay)
delta_t_rise_st.append(delta_t_rise)
delta_t_decay_st.append(delta_t_decay)



########BGO Pre#######

t_rise = 104#160                           #[ns]
t_decay = 1.15e3 #1.68e3                      #[ns], 3.28us (micros )
delta_t_rise = 1/5 * 100#1/5*1e3                  #[ns], error of the oscilloscope
delta_t_decay = 1/5 * 250#1/5*1e3
  
#storing
t_rise_st.append(t_rise)
t_decay_st.append(t_decay)
delta_t_rise_st.append(delta_t_rise)
delta_t_decay_st.append(delta_t_decay)


####Plot

#RAW

plt.figure(figsize=(13,6))  #width, heigh 6.4*4.8 inches by default
plt.subplot(1, 2, 1)
plt.suptitle("Rise and decay time of the raw signal of the Cs137 peak", fontsize=22, wrap=True)           #title
plt.bar(['LYSO', 'CsI', 'BGO'], t_rise_st[0:3], yerr = delta_t_rise_st[0:3], edgecolor="black")
#plt.xlabel("ADC channels", fontsize=10)                        #xlabel
plt.ylabel("Rise time (ns)", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 

plt.subplot(1, 2, 2)
plt.bar(['LYSO', 'CsI', 'BGO'], t_decay_st[0:3], yerr = delta_t_decay_st[0:3])
#plt.xlabel("ADC channels", fontsize=10)                        #xlabel
plt.ylabel("Decay time (ns)", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.savefig('Rise_decay_time_raw.png', format='png')


#Print:
print('Decay time of raw LYSO: (' + str(t_decay_st[0]) + ' +/- ' + str(delta_t_decay_st[0]) + ')ns')
print('Decay time of raw CsI: (' + str(t_decay_st[1]) + ' +/- ' + str(delta_t_decay_st[1]) + ')ns')
print('Decay time of raw BGO: (' + str(t_decay_st[2]) + ' +/- ' + str(delta_t_decay_st[2]) + ')ns'+"\n")
#################################################

#Pre

plt.figure(figsize=(13,6))  #width, heigh 6.4*4.8 inches by default
plt.subplot(1, 2, 1)
plt.suptitle("Rise and decay time of the signal (pre) of the Cs137 peak", fontsize=22, wrap=True)           #title
plt.bar(['LYSO', 'CsI', 'BGO'], t_rise_st[3:6], yerr = delta_t_rise_st[3:6], edgecolor="black")
#plt.xlabel("ADC channels", fontsize=10)                        #xlabel
plt.ylabel("Rise time (ns)", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 

plt.subplot(1, 2, 2)
plt.bar(['LYSO', 'CsI', 'BGO'], t_decay_st[3:6], yerr = delta_t_decay_st[3:6])
#plt.xlabel("ADC channels", fontsize=10)                        #xlabel
plt.ylabel("Decay time (ns)", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.savefig('Rise_decay_time_pre.png', format='png')


#Print:
print('Decay time of pre LYSO: (' + str(t_decay_st[3]) + ' +/- ' + str(delta_t_decay_st[3]) + ')ns')
print('Decay time of pre CsI: (' + str(t_decay_st[4]) + ' +/- ' + str(delta_t_decay_st[4]) + ')ns')
print('Decay time of pre BGO: (' + str(t_decay_st[5]) + ' +/- ' + str(delta_t_decay_st[5]) + ')ns'+"\n")
##############################################

#Lets now comput the ratio of the decay time, of both signals, raw and pre.
#0,1,2 ra ofLYSO, CsI, BGO respect
#3,4,5 pre of LYSO, CsI, BGO respect

ratio_decay = [t_decay_st[0]/t_decay_st[1], t_decay_st[0]/t_decay_st[2], t_decay_st[2]/t_decay_st[1],
              t_decay_st[3]/t_decay_st[4], t_decay_st[3]/t_decay_st[5], t_decay_st[5]/t_decay_st[4]
              ]  
                  #LYSO/CsI,  LYSO/BGO,  BGO/CsI, first raw and then pre

#lets now compute its error:
    
#error calculation of the ratio
auxiliar2 = [(delta_t_decay_st[0]/t_decay_st[0])**2 + (delta_t_decay_st[1]/t_decay_st[1])**2, 
            (delta_t_decay_st[0]/t_decay_st[0])**2 + (delta_t_decay_st[2]/t_decay_st[2])**2,
            (delta_t_decay_st[2]/t_decay_st[2])**2 + (delta_t_decay_st[1]/t_decay_st[1])**2,
            (delta_t_decay_st[3]/t_decay_st[3])**2 + (delta_t_decay_st[4]/t_decay_st[4])**2, 
            (delta_t_decay_st[3]/t_decay_st[3])**2 + (delta_t_decay_st[5]/t_decay_st[5])**2,
            (delta_t_decay_st[5]/t_decay_st[5])**2 + (delta_t_decay_st[4]/t_decay_st[4])**2,
            ] 
            #this are (delta_t1/t1)^2 +  (delta_t2/t2)^2, each row is this, in the order of the 
            #variable ratio_decay. This is basically the thing to be squared rooted to compute the 
            #error

delta_ratios = np.multiply(ratio_decay, [np.sqrt(x) for x in auxiliar2])
                                             #this will be though

print('Decay time ratio of raw LYSO/CsI: (' + str(ratio_decay[0]) + ' +/- ' + str(delta_ratios[0]) )
print('Decay time of raw LYSO/BGO: (' + str(ratio_decay[1]) + ' +/- ' + str(delta_ratios[1]) )
print('Decay time of raw BGO/CsI: (' + str(ratio_decay[2]) + ' +/- ' + str(delta_ratios[2]) +"\n")

print('Decay time ratio of pre LYSO/CsI: (' + str(ratio_decay[3]) + ' +/- ' + str(delta_ratios[3]))
print('Decay time of pre LYSO/BGO: (' + str(ratio_decay[4]) + ' +/- ' + str(delta_ratios[4]))
print('Decay time of pre BGO/CsI: (' + str(ratio_decay[5]) + ' +/- ' + str(delta_ratios[5]) +"\n")
#%% RESOLUTION
#forthis I use the common configuration that allows to plot a common histogram for every crystal
#similar as what is done by CAEN!!.

#########################4), Data loading #####################3
#The files to load are in txt. The best way to read is:

#All with same histogram, gain, bias, threshold, but different gate to optimize 
#each signal

#Firstly create the data by setting the led driver amplitude to 0, obtaining 
#extrange results, so will now unplug it, and closing the SiPM to see if now 
#the results are fine.

total_counts = []       #variable that will contain the total counts of all the
                        #spectras
counts_stored = []                          #storing variable of the counts
rate_stored = []                #storing variable of the count rate
time = [600, 200, 200]              #[s] duration time od the measurements
#CsI, BGO, LYSO

#All have the same channels, so only one variable will be saved

with open('Cs_137_1_elevacion_CsI_histo.txt') as file_object:
            
        lines = file_object.readlines()
        print('the number of lines of the Cs137_CsI is',len(lines))
        #This contains strings (have to be converted to numbers using int()
        #and \n, so the \n (salto de linea) have to be removed
        
        ADC_channel = []
        counts_help = []                    #mid variable, to be used to count
        for i in range(len(lines)):            
            ADC_channel.append(float(lines[i].split()[0])) #store 1st number of the
                            #column
            counts_help.append(float(lines[i].split()[1]))    #store 2nd number of the
                            #column
                            
        total_counts.append(sum(counts_help))  #total counts of the spectra   
        count_rate_help = [x/time[0] for x in counts_help] #count rate
        
        #Storing
        counts_stored.append(counts_help)
        rate_stored.append(count_rate_help)

with open('Cs_137_1_elevacion_BGO_histo.txt') as file_object:
            
        lines = file_object.readlines()
        print('the number of lines of the Cs137_BGO is',len(lines))
        #This contains strings (have to be converted to numbers using int()
        #and \n, so the \n (salto de linea) have to be removed
        
        counts_help = []
        for i in range(len(lines)):            
            counts_help.append(float(lines[i].split()[1]))    #store 2nd number of the
                            #column
       
        total_counts.append(sum(counts_help))  #total counts of the spectra         
        count_rate_help = [x/time[0] for x in counts_help] #count rate
        
        #Storing
        counts_stored.append(counts_help)
        rate_stored.append(count_rate_help)

with open('Cs_137_1_elevacion_LYSO_histo.txt') as file_object:
            
        lines = file_object.readlines()
        print('the number of lines of the Cs137_LYSO is',len(lines))
        #This contains strings (have to be converted to numbers using int()
        #and \n, so the \n (salto de linea) have to be removed
        
        counts_help = []
        for i in range(len(lines)):            
            counts_help.append(float(lines[i].split()[1]))    #store 2nd number of the
                            #column
       
        total_counts.append(sum(counts_help))  #total counts of the spectra       
        count_rate_help = [x/time[0] for x in counts_help] #count rate
        
        #Storing
        counts_stored.append(counts_help)
        rate_stored.append(count_rate_help)
        


#First column of counts stored is CsI, second GBo, 3rd LYSO

#    0.1. Representacion
            
#CsI
plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
plt.bar(ADC_channel,counts_stored[0], width = ADC_channel[1]-ADC_channel[0])     
        #widht so that each bar touches each other!
plt.title("Spectra of Cs^{137} with CsI", fontsize=24)           #title
plt.xlabel("ADC channels", fontsize=14)                        #xlabel
plt.ylabel("Counts", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.xlim(0,max(ADC_channel))                       #limits of x axis
#plt.ylim(0,11000)                            #limits of y axis


#BGO
plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
plt.bar(ADC_channel,counts_stored[1], width = ADC_channel[1]-ADC_channel[0])     
        #widht so that each bar touches each other!
plt.title("Spectra of Cs^{137} with BGO", fontsize=24)           #title
plt.xlabel("ADC channels", fontsize=14)                        #xlabel
plt.ylabel("Counts", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.xlim(0,max(ADC_channel))                       #limits of x axis
#plt.ylim(0,11000)                            #limits of y axis
#plt.yscale('log')  


#LYSO
plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
plt.bar(ADC_channel,counts_stored[2], width = ADC_channel[1]-ADC_channel[0])     
        #widht so that each bar touches each other!
plt.title("Spectra of Cs^{137} with LYSO", fontsize=24)           #title
plt.xlabel("ADC channels", fontsize=14)                        #xlabel
plt.ylabel("Counts", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.xlim(0,max(ADC_channel))                       #limits of x axis
#plt.ylim(0,)                        c    #limits of y axis



###Plot combined, as CAEN's
#plotting the count rate, since the measure time is differen!!
plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
plt.plot(ADC_channel, rate_stored[0], color='black', label = 'CsI')    
plt.plot(ADC_channel, rate_stored[1], color = 'red', label = 'BGO')
plt.plot(ADC_channel, rate_stored[2], color = 'blue', label = 'LYSO')      
plt.legend(['CsI', 'BGO', 'LYSO'], fontsize=10) 
plt.title("Cs137 spectra for several scintillators", fontsize=20)           #title
plt.xlabel("ADC channels", fontsize=10)                        #xlabel
plt.ylabel("Count rate [Hz]", fontsize=10)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=10)              #size of axis
plt.grid(True) 
plt.savefig('count_rate_Cs137_vs_scintillator_type.png', format='png')

#plt.xlim(0,800)                       #limits of x axis
#plt.ylim(0,11000)                            #limits of y axis


#%% ##########5) FIT #################################

#Lets do the gaussian fit to the gamma () of Cs137, to see the FWHM as a function
#of the scintillation crystal

def gaussian(x, a, b, c):       #Definition of the function to use to fit the data
    return a * np.exp(- (x-b)**2 / (2 * c**2)) 

        #this is a gaussian function (more general than normal distribution)
        
        #if using math.exp the fit gives error: 
        #(only size-1 arrays can be converted to Python scalars)
        #with numpy.exp everything fine.       
        #
        #FRIENDSHIP ENDED WITH math.exp(), NOW numpy.exp() IS MY 
        #BEST FRIEND
        #


#$$$$$$$$$$$$$$$$$$$ CsI $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#To find the intervals, I could easily find the index of the maximum, and from
#there move by looking at the .txt, or maybe say, peak me 100 hundred values above
#and below the max, etc.

#I could use here the counts, do not need to use the count rate, although it
#should give the same result.


peak = max(counts_stored[0])        #max value of the count rate, i.e., peak
peak_index = counts_stored[0].index(peak)

#But, since here the peak do not have the max value, I will have to serach it 
#by hand :()

import Gaussian_fit

fit = Gaussian_fit.Gaussian_fit(ADC_channel[4923-1:5741-1], 
                                   counts_stored[0][4923-1:5741-1])

#Storing of the relevant data, sigma and its error
sigma_stored = []
mean_stored = []
delta_mean_stored = []
delta_sigma_stored = []
FWHM_stored = []
delta_FWHM_stored = []

sigma_stored.append(fit['sigma'])
mean_stored.append(fit['mean'])
delta_mean_stored.append(fit['\Delta(mean)'])
delta_sigma_stored.append(fit['\Delta(sigma)'])
FWHM_stored.append(fit['FWHM'])
delta_FWHM_stored.append(fit['\Delta(FWHM)'])


#$$$$$$$$$$$$$$$$$$$ BGO $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

peak = max(counts_stored[1])        #max value of the count rate, i.e., peak
peak_index = counts_stored[1].index(peak) #1143

fit = Gaussian_fit.Gaussian_fit(ADC_channel[1036-1:1273-1], 
                                   counts_stored[1][1036-1:1273-1])

#Storing of the relevant data, sigma and its error

sigma_stored.append(fit['sigma'])
mean_stored.append(fit['mean'])
delta_mean_stored.append(fit['\Delta(mean)'])
delta_sigma_stored.append(fit['\Delta(sigma)'])
FWHM_stored.append(fit['FWHM'])
delta_FWHM_stored.append(fit['\Delta(FWHM)'])

#$$$$$$$$$$$$$$$$$$$ LYSO $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#here, again, the maximum peak is not the one of the gamma peak, so again have
#to find the peak by hand :)

fit = Gaussian_fit.Gaussian_fit(ADC_channel[2548-1:3191-1], 
                                   counts_stored[2][2548-1:3191-1])

#Storing of the relevant data, sigma and its error

sigma_stored.append(fit['sigma'])
mean_stored.append(fit['mean'])
delta_mean_stored.append(fit['\Delta(mean)'])
delta_sigma_stored.append(fit['\Delta(sigma)'])
FWHM_stored.append(fit['FWHM'])
delta_FWHM_stored.append(fit['\Delta(FWHM)'])


#%% ############### 6) Plot of R vs Scintillator ######################

#resolution = FWHM/<E> , <E> the centroid of the peak.

R_stored = np.multiply(FWHM_stored, [1/x for x in mean_stored])           #R
       
R_stored_100 = [100*x for x in R_stored]                                 #R[%]

#calc of delta R:
delta_FWHM_FWHM = np.multiply(delta_FWHM_stored,[1/x for x in FWHM_stored]) 
                                               #\Delta(FWHM)/FWHM
delta_mean_mean = np.multiply(delta_mean_stored,[1/x for x in mean_stored]) 
                                               #\Delta(FWHM)/FWHM
sum__ = np.multiply(delta_FWHM_FWHM, delta_FWHM_FWHM) + np.multiply(delta_mean_mean, delta_mean_mean)
        
sqrt_sum = [np.sqrt(x) for x in sum__]          #sqrt of the sum of relative errors

delta_R_stored = np.multiply(R_stored, sqrt_sum)              #delta(R)
delta_R_stored_100 = [100*x for x in delta_R_stored]              #delta_R[%]

#Plot

plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
plt.bar(['CsI', 'BGO', 'LYSO'], R_stored_100, yerr = delta_R_stored_100, edgecolor="black")
plt.title("Resolution of the Cs137 peak for several scintillators", fontsize=22, wrap=True)           #title
#plt.xlabel("ADC channels", fontsize=10)                        #xlabel
plt.ylabel("R (%)", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.savefig('Resolution_vs_scintillator.png', format='png')




###print

#Print:
print('FWHM CsI: ' + str(FWHM_stored[0]) + ' +/- ' + str(delta_FWHM_stored[0]))
print('FWHM BGO: ' + str(FWHM_stored[1]) + ' +/- ' + str(delta_FWHM_stored[1]))
print('FWHM LYSO: ' + str(FWHM_stored[2]) + ' +/- ' + str(delta_FWHM_stored[2])+"/n")

print('<channels> CsI: ' + str(mean_stored[0]) + ' +/- ' + str(delta_mean_stored[0]))
print('<channels> BGO: ' + str(mean_stored[1]) + ' +/- ' + str(delta_mean_stored[1]))
print('<channels> LYSO: ' + str(mean_stored[2]) + ' +/- ' + str(delta_mean_stored[2])+"/n")

print('R CsI: ' + str(R_stored_100[0]) + ' +/- ' + str(delta_R_stored_100[0]))
print('R BGO: ' + str(R_stored_100[1]) + ' +/- ' + str(delta_R_stored_100[1]))
print('R LYSO: ' + str(R_stored_100[2]) + ' +/- ' + str(delta_R_stored_100[2]))

#Lyso 2 BGO 1




#%% TRY

#Juanpa suggest to simply compute the ratio of the total number of counts to see
#whether this could be similar to the light yield ratio or no. So, come on!
#Since the time are different, I must compute the count/rate in order
#to compare the values

#0 csI, 1 BGO, 2 LYSO. 

one_slash_time = [1/x for x in time]                #1/time [s-1] to compute the ratio

total_rate = np.multiply(total_counts, one_slash_time) #[s-1] total counts/total time

ratio_total_counts = [total_rate[2]/total_rate[0], total_rate[2]/total_rate[1], 
                      total_rate[1]/total_rate[2]] #LYSO/CsI, LYSO/BGO, BGO/CsI

#Error

delta_total_counts = [np.sqrt(x) for x in total_counts]     #error of the total counts

#asuming the time has no error, the error of the ratio of total counts is:

aux = [ (np.sqrt(total_counts[2]) / total_counts[2])**2 + (np.sqrt(total_counts[0]) / total_counts[0])**2,
       (np.sqrt(total_counts[2]) / total_counts[2])**2 + (np.sqrt(total_counts[1]) / total_counts[1])**2,
       (np.sqrt(total_counts[1]) / total_counts[1])**2 + (np.sqrt(total_counts[2]) / total_counts[2])**2,
       ]      #(delta(c_1)/c_1)^2 + (delta(c_2)/c_2)^2), c = total counts. SAME ORDER THAN THE RATIO!!!!

aux = [np.sqrt(x) for x in aux]         #sqrt( (delta(c_1)/c_1)^2 + (delta(c_2)/c_2)^2)
                #c = total counts

delta_ratio_total_counts = np.multiply(ratio_total_counts, aux)
    

#ratio_total_counts = [total_counts[2]/total_counts[0], total_counts[2]/total_counts[1], 
#                      total_counts[1]/total_counts[2]] #LYSO/CsI, LYSO/BGO, BGO/CsI



plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
plt.title("Light yield ratio of the total counts of the spectra of the Cs137 (total count division)", 
          fontsize=22, wrap=True)           #title
plt.bar(['LYSO/CsI', 'LYSO/BGO', 'BGO/CsI'], ratio_total_counts, yerr = delta_ratio_total_counts, 
        edgecolor="black")#, yerr = delta_light_yield)
#plt.xlabel("ADC channels", fontsize=10)                        #xlabel
plt.ylabel("Light yield ratio", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.savefig('Ratio_total_counts_vs_scintillator.png', format='png')

#Dude, what the actual fuck, this is sick, totally match the light yield ratio.









#%% REDO OF THE RESOLUTION CALC, BUT THIS TIME ONLY MEASURING THE PEAK.

# #since this data correspond to the photopeak only, all the variables will contain
# #photo in the name


# total_counts_photo = []       #variable that will contain the total counts of all the
#                         #spectras
# ADC_channels_stored_photo = []                  #storing of the channels
#         #necause the .txt is fucked up, and several lines are repetaed, and so on,
#         #and spotting them by hand is though, so fuck it
        
# counts_stored_photo = []                          #storing variable of the counts
# rate_stored_photo = []                #storing variable of the count rate
# time_photo = [800, 200, 200]              #[s] duration time od the measurements
# #CsI, BGO, LYSO

# #All have the same channels, so only one variable will be saved

# with open('Cs137_CsI_photopeak_histo.txt') as file_object:
            
#         lines = file_object.readlines()
#         print('the number of lines of the Cs137_CsI (photo) is',len(lines))
#         #This contains strings (have to be converted to numbers using int()
#         #and \n, so the \n (salto de linea) have to be removed
        
#         ADC_help = []
#         counts_help = []                    #mid variable, to be used to count
#         for i in range(len(lines)):            
#             ADC_help.append(float(lines[i].split()[0])) #store 1st number of the
#                             #column
#             counts_help.append(float(lines[i].split()[1]))    #store 2nd number of the
#                             #column
                            
#         total_counts_photo.append(sum(counts_help))  #total counts of the spectra   
#         count_rate_help = [x/time_photo[0] for x in counts_help] #count rate
        
#         #Storing
#         counts_stored_photo.append(counts_help)
#         rate_stored_photo.append(count_rate_help)
#         ADC_channels_stored_photo.append(ADC_help)

# with open('Cs137_BGO_photopeak_histo.txt') as file_object:
            
#         lines = file_object.readlines()
#         print('the number of lines of the Cs137_BGO (photo) is',len(lines))
#         #This contains strings (have to be converted to numbers using int()
#         #and \n, so the \n (salto de linea) have to be removed
        
#         ADC_help = []
#         counts_help = []                    #mid variable, to be used to count
#         for i in range(len(lines)):            
#             ADC_help.append(float(lines[i].split()[0])) #store 1st number of the
#                             #column
#             counts_help.append(float(lines[i].split()[1]))    #store 2nd number of the
#                             #column
                            
#         total_counts_photo.append(sum(counts_help))  #total counts of the spectra   
#         count_rate_help = [x/time_photo[0] for x in counts_help] #count rate
        
#         #Storing
#         counts_stored_photo.append(counts_help)
#         rate_stored_photo.append(count_rate_help)
#         ADC_channels_stored_photo.append(ADC_help)


# with open('Cs137_LYSO_photopeak_histo.txt') as file_object:
            
#         lines = file_object.readlines()
#         print('the number of lines of the Cs137_LYSO (photo) is',len(lines))
#         #This contains strings (have to be converted to numbers using int()
#         #and \n, so the \n (salto de linea) have to be removed
        
#         ADC_help = []
#         counts_help = []                    #mid variable, to be used to count
#         for i in range(len(lines)):            
#             ADC_help.append(float(lines[i].split()[0])) #store 1st number of the
#                             #column
#             counts_help.append(float(lines[i].split()[1]))    #store 2nd number of the
#                             #column
                            
#         total_counts_photo.append(sum(counts_help))  #total counts of the spectra   
#         count_rate_help = [x/time_photo[0] for x in counts_help] #count rate
        
#         #Storing
#         counts_stored_photo.append(counts_help)
#         rate_stored_photo.append(count_rate_help)
#         ADC_channels_stored_photo.append(ADC_help)


# #First column of counts stored is CsI, second GBo, 3rd LYSO


# #    0.1. Representacion
            
# #CsI
# plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
# plt.bar(ADC_channels_stored_photo[0],counts_stored_photo[0], 
#         width = ADC_channels_stored_photo[0][1]-ADC_channels_stored_photo[0][0])     
#         #widht so that each bar touches each other!
# plt.title("Spectra (photo) of Cs^{137} with CsI", fontsize=24)           #title
# plt.xlabel("ADC channels", fontsize=14)                        #xlabel
# plt.ylabel("Counts", fontsize=14)              #ylabel
# # Set size of tick labels.
# plt.tick_params(axis='both', labelsize=14)              #size of axis
# plt.grid(True) 
# plt.xlim(0,max(ADC_channels_stored_photo[0]))                       #limits of x axis
# #plt.ylim(0,11000)                            #limits of y axis


# #BGO
# plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
# plt.bar(ADC_channels_stored_photo[1],counts_stored_photo[1], 
#         width = ADC_channels_stored_photo[0][1]-ADC_channels_stored_photo[0][0])    
#                 #widht so that each bar touches each other!
# plt.title("Spectra (photo) of Cs^{137} with BGO", fontsize=24)           #title
# plt.xlabel("ADC channels", fontsize=14)                        #xlabel
# plt.ylabel("Counts", fontsize=14)              #ylabel
# # Set size of tick labels.
# plt.tick_params(axis='both', labelsize=14)              #size of axis
# plt.grid(True) 
# plt.xlim(0,max(ADC_channels_stored_photo[0]))                       #limits of x axis


# #LYSO
# plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
# plt.bar(ADC_channels_stored_photo[2],counts_stored_photo[2], 
#         width = ADC_channels_stored_photo[0][1]-ADC_channels_stored_photo[0][0])  
#         #widht so that each bar touches each other!
# plt.title("Spectra (photo) of Cs^{137} with LYSO", fontsize=24)           #title
# plt.xlabel("ADC channels", fontsize=14)                        #xlabel
# plt.ylabel("Counts", fontsize=14)              #ylabel
# # Set size of tick labels.
# plt.tick_params(axis='both', labelsize=14)              #size of axis
# plt.grid(True) 
# plt.xlim(0,max(ADC_channels_stored_photo[0]))                       #limits of x axis


# #$$$$$$$$$$$$$$$$$$$ CsI $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# #watch out, the data is being overwritten, so data from the whole spectrum is being lost!!!!!



# peak = max(counts_stored_photo[0])        #max value of the count rate, i.e., peak
# peak_index = counts_stored_photo[0].index(peak)

# fit = Gaussian_fit.Gaussian_fit(ADC_channels_stored_photo[0][4932-1:5732-1], 
#                                    counts_stored_photo[0][4932-1:5732-1]) 

# #Storing of the relevant data, sigma and its error
# sigma_stored = []
# mean_stored = []
# delta_mean_stored = []
# delta_sigma_stored = []
# FWHM_stored = []
# delta_FWHM_stored = []

# sigma_stored.append(fit['sigma'])
# mean_stored.append(fit['mean'])
# delta_mean_stored.append(fit['\Delta(mean)'])
# delta_sigma_stored.append(fit['\Delta(sigma)'])
# FWHM_stored.append(fit['FWHM'])
# delta_FWHM_stored.append(fit['\Delta(FWHM)'])


# #$$$$$$$$$$$$$$$$$$$ BGO $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# peak = max(counts_stored[1])        #max value of the count rate, i.e., peak
# peak_index = counts_stored[1].index(peak) #1143

# # ####debug
# # plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
# # plt.bar(ADC_channels_stored_photo[1][626-1:865-1],counts_stored_photo[1][626-1:865-1], 
# #         width = ADC_channels_stored_photo[0][1]-ADC_channels_stored_photo[0][0]) 
# # plt.plot(ADC_channels_stored_photo[1][626-1:865-1], 
# #          gaussian(ADC_channels_stored_photo[1][626-1:865-1], a, mean_stored[-1], cc), 'ro', label = 'fit')

# # ####
# fit = Gaussian_fit.Gaussian_fit(ADC_channels_stored_photo[1][636-1:865-1], #865
#                                    counts_stored_photo[1][636-1:865-1])

# #Storing of the relevant data, sigma and its error

# sigma_stored.append(fit['sigma'])
# mean_stored.append(fit['mean'])
# delta_mean_stored.append(fit['\Delta(mean)'])
# delta_sigma_stored.append(fit['\Delta(sigma)'])
# FWHM_stored.append(fit['FWHM'])
# delta_FWHM_stored.append(fit['\Delta(FWHM)'])

# #$$$$$$$$$$$$$$$$$$$ LYSO $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# #here, again, the maximum peak is not the one of the gamma peak, so again have
# #to find the peak by hand :)

# fit = Gaussian_fit.Gaussian_fit(ADC_channels_stored_photo[2][1787-1:2346-1], #1787 2346
#                                    counts_stored_photo[2][1787-1:2346-1])

# #Storing of the relevant data, sigma and its error

# sigma_stored.append(fit['sigma'])
# mean_stored.append(fit['mean'])
# delta_mean_stored.append(fit['\Delta(mean)'])
# delta_sigma_stored.append(fit['\Delta(sigma)'])
# FWHM_stored.append(fit['FWHM'])
# delta_FWHM_stored.append(fit['\Delta(FWHM)'])


# ############### 6) Plot of R vs Scintillator ######################

# #resolution = FWHM/<E> , <E> the centroid of the peak.

# R_stored = np.multiply(FWHM_stored, [1/x for x in mean_stored])           #R
       
# R_stored_100 = [100*x for x in R_stored]                                 #R[%]

# #calc of delta R:
# delta_FWHM_FWHM = np.multiply(delta_FWHM_stored,[1/x for x in FWHM_stored]) 
#                                                #\Delta(FWHM)/FWHM
# delta_mean_mean = np.multiply(delta_mean_stored,[1/x for x in mean_stored]) 
#                                                #\Delta(FWHM)/FWHM
# sum__ = np.multiply(delta_FWHM_FWHM, delta_FWHM_FWHM) + np.multiply(delta_mean_mean, delta_mean_mean)
        
# sqrt_sum = [np.sqrt(x) for x in sum__]          #sqrt of the sum of relative errors

# delta_R_stored = np.multiply(R_stored, sqrt_sum)              #delta(R)
# delta_R_stored_100 = [100*x for x in delta_R_stored]              #delta_R[%]

# #Plot

# plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
# plt.bar(['CsI', 'BGO', 'LYSO'], R_stored_100, yerr = delta_R_stored_100, edgecolor="black")
# plt.title("Resolution of the Cs137 peak (photo) for several scintillators", fontsize=22, wrap=True)           #title
# #plt.xlabel("ADC channels", fontsize=10)                        #xlabel
# plt.ylabel("R (%)", fontsize=14)              #ylabel
# # Set size of tick labels.
# plt.tick_params(axis='both', labelsize=14)              #size of axis
# plt.grid(True) 
# plt.savefig('Resolution_vs_scintillator_gamma.png', format='png')




# ###print

# #Print:
# print('FWHM CsI photo: ' + str(FWHM_stored[0]) + ' +/- ' + str(delta_FWHM_stored[0]))
# print('FWHM BGO photo: ' + str(FWHM_stored[1]) + ' +/- ' + str(delta_FWHM_stored[1]))
# print('FWHM LYSO photo: ' + str(FWHM_stored[2]) + ' +/- ' + str(delta_FWHM_stored[2])+"/n")

# print('<channels> CsI photo: ' + str(mean_stored[0]) + ' +/- ' + str(delta_mean_stored[0]))
# print('<channels> BGO photo: ' + str(mean_stored[1]) + ' +/- ' + str(delta_mean_stored[1]))
# print('<channels> LYSO photo: ' + str(mean_stored[2]) + ' +/- ' + str(delta_mean_stored[2])+"/n")

# print('R CsI photo: ' + str(R_stored_100[0]) + ' +/- ' + str(delta_R_stored_100[0]))
# print('R BGO photo: ' + str(R_stored_100[1]) + ' +/- ' + str(delta_R_stored_100[1]))
# print('R LYSO photo: ' + str(R_stored_100[2]) + ' +/- ' + str(delta_R_stored_100[2]))

# #Lyso 2 BGO 1





























#%% RESIDUOS


###plotly

# data = [Bar(x=time, y=voltage)]
# x_axis_config = {'title': 'Time'}
# y_axis_config = {'title': 'Voltage'}
# my_layout = Layout(title='LYSO',
# xaxis=x_axis_config, yaxis=y_axis_config)
# offline.plot({'data': data, 'layout': my_layout}, filename='LYSO.html')
   
# #Dude, the plots are fucking horrible, so moving on to only plot the bars of the
# #voltage

# plt.figure(figsize=(10,5))  #width, heigh 6.4*4.8 inches by default
# plt.plot(list(range(len(voltage))),voltage)   
###############

# with open('TEK0000_a_mano.txt') as file_object:
            
#         lines = file_object.readlines()
#         print('the number of lines of the Cs137_CsI is',len(lines))
#         #This contains strings (have to be converted to numbers using int()
#         #and \n, so the \n (salto de linea) have to be removed
        
#         t = []
#         V = []                    #mid variable, to be used to count
#         for i in range(len(lines)):            
#             t.append(float(lines[i].split()[0])) #store 1st number of the
#                             #column
#             V.append(float(lines[i].split()[1]))    #store 2nd number of the
#                             #column
   
# plt.figure(figsize=(10,5))  #width, heigh 6.4*4.8 inches by default
# plt.plot(t,V, 'bo-')    
                     