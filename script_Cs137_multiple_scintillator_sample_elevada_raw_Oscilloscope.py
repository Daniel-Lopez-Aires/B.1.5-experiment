#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 10:00:44 2021

@author: dla

This script contains the calcs of the B.1.5. experiment using the signals from the oscilloscope. The signals are from the pre (divided with the disivor, to see the histogram) and directly from the detector.
The wrong way to compute the light yield ratio is also included, commented.

The sample is elevated, 1 cylinder only (plus the little ring placed just above the scintillator (see pictures in logbook))
"""

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


############TO DO

#1) MAKE EVERY ARRAY NUMPY ARRAY, SINCE NUMPY ARRAYS ARE ABLE TO DO ELEMENT WISE OPERATIONS
    #EASILY, WHICH IS SOMETHING I NEED A FUCKING LOT!
    



#%%
#########################1), Data loading #####################3
#The files to load are in txt. The best way to read is:


##Data of the natural signal, from the oscilloscope, in.csv. 

voltage_stored = []
time_stored = []


#################################
###########RAW##################3
################################


#with open('TEK0001.CSV') as file_object:            #LYSO, raw
with open('TEK0002_LYSO_raw.CSV') as file_object:            #LYSO, raw   NEW(13/5)         
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




#with open('TEK0005.CSV') as file_object:            #CsI, raw
with open('TEK0001_CsI_raw.CSV') as file_object:            #CsI, raw   NEW(13/5)                
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





#with open('TEK0007.CSV') as file_object:            #BGo, raw
with open('TEK0000_BGO_raw.CSV') as file_object:            #BGO, raw   NEW(13/5)              
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
plt.title("Waveform of Cs137 with LYSO", fontsize=22)           #title
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
plt.title("Waveform of Cs137 with CsI", fontsize=22)           #title
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
plt.title("Waveform of Cs137 with BGO", fontsize=22)           #title
plt.xlabel("time [us]", fontsize=14)                        #xlabel
plt.ylabel("voltage [mV]", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
#plt.xlim(0,max(ADC_channel))                       #limits of x axis     
plt.savefig('Signal_BGO.png', format='png')







#%% Rise and decay time######################################################

#Usando misma configuracion, mido el máximo, del pre, q sera el del fotopico. Pre mejor
#q raw pq el raw tenia mas pileups q el pre. Mediciones tanto para la raw signal como para
#la del pre!!!.

#Variables that will store(st) the results 
t_rise_st = []
t_decay_st = []
delta_t_rise_st = []
delta_t_decay_st = []

###################################LYSO RAW########################

#RAW

t_rise = 280#18.8e3                            #[ns]
t_decay = 55*1e3                      #[ns], 3.28us (micros )
delta_t_rise = 1/5*500                #[ns], error of the oscilloscope
                #it basically is 1/5 of the escale, since each square is divided
                #into 5 intervals, and the scale is the length of the sides of
                #the square
delta_t_decay = 1/5*10e3
  
#storing
t_rise_st.append(t_rise)
t_decay_st.append(t_decay)
delta_t_rise_st.append(delta_t_rise)
delta_t_decay_st.append(delta_t_decay)

############################CsI RAW#####################3#####


t_rise = 3.2e3                            #[ns]
t_decay = 145*1e3                      #[ns], 3.28us (micros )
delta_t_rise = 1/5*5e3                   #[ns], error of the oscilloscope
delta_t_decay = 1/5*25e3
  
#storing
t_rise_st.append(t_rise)
t_decay_st.append(t_decay)
delta_t_rise_st.append(delta_t_rise)
delta_t_decay_st.append(delta_t_decay)

####################################BGO RAW###################

t_rise = 1.1e3                           #[ns]
t_decay = 60.8e3                     #[ns], 3.28us (micros )
delta_t_rise = 1/5*2.5e3                #[ns], error of the oscilloscope
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

#################################LYSO PRE##################3####

t_rise =  156  #28                           #[ns]
t_decay = 464  #640                      #[ns], 3.28us (micros )
delta_t_rise = 1/5 * 100 #1/5*100                   #[ns], error of the oscilloscope
delta_t_decay = 1/5 * 100#1/5*100
  
#storing
t_rise_st.append(t_rise)
t_decay_st.append(t_decay)
delta_t_rise_st.append(delta_t_rise)
delta_t_decay_st.append(delta_t_decay)



###############################CsI PRE#########################

t_rise = 320 #360                           #[ns]
t_decay = 2.88e3    #4.52e3                      #[ns], 3.28us (micros )
delta_t_rise = 1/5 * 500#1/5*1e3                  #[ns], error of the oscilloscope
delta_t_decay = 1/5 * 500 #1/5*1e3
  
#storing
t_rise_st.append(t_rise)
t_decay_st.append(t_decay)
delta_t_rise_st.append(delta_t_rise)
delta_t_decay_st.append(delta_t_decay)



###############################3###BGO Pre############################

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
plt.suptitle("Rise and decay time of the signal of the Cs137 peak", fontsize=22, wrap=True)           #title
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




 #%% 3) Calcs of the light yield############################################

#Since the signals from the Pre are way better, will use them. Remember they are
#the last stored (6 total stored), so from 3 to 5 (0 the first). LYSO, CsI, BGO the order
#Note we are not sure about the light yield concept, but computing the amplitude will be good.


####################################LYSO##########################3

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

#Baseline.
    #to choose the baseline, the 1st approach was to do the mean between channel 0 and the
    #channel where the peak start to appear. But this is abd idea becaue the error is very high,
    #so the 2nd approach will simply be to choose a value, the one that appears most, by eyeseeking

#len_baseline_signal_points = len(voltage_stored[3][0:index_min-1])  #length of the voltages used
        #to compute the baseline. Will be neccesary for error calcs
#baseline_signal = sum(voltage_stored[3][0:index_min-1]) / len_baseline_signal_points    
        #[V] baseline, to compute the peak amplitude. I average between the initial value and
        #when the peak starts
baseline = 0.16         #[V], eyeseeking on the interval voltage_stored[3][0:index_min-1]


delta_V = 1/5 *100e-3           #[V] error of the voltage emasurements, from the scale of
                            #the oscilloscope, which come from the photos
        
sum_voltage = sum(voltage_peak)         #[V] total voltage of the peak

#Storing of the current

voltage_peak_stored = []                #storage of the total voltage of the peak

baseline_stored = []               #storage of the baseline
peak_stored = []                  #storage of the peak value, for the max amplitude
n_elements_peak = []                    #this will store the number of voltages I sum,
                #for each peak, for the error calc
delta_single_V_measurement = []        #storage of the error of the voltage 
                    #measurements, for the error calc
len_baseline_signal_points_stored = []                    
                    
voltage_peak_stored.append(sum_voltage)
baseline_stored.append(baseline)
peak_stored.append(peak)
n_elements_peak.append(len_peak)
delta_single_V_measurement.append(delta_V)
#len_baseline_signal_points_stored.append(len_baseline_signal_points)

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


###############################################CsI#####################

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

baseline = 0.1664         #[V], eyeseeking on the interval voltage_stored[4][0:index_min-1]

delta_V = 1/5 *20e-3                        #[V] 

sum_voltage = sum(voltage_peak)         #[V] total voltage of the peak

#Storing of the current

voltage_peak_stored.append(sum_voltage)
baseline_stored.append(baseline)
peak_stored.append(peak)
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

#######################################BGO##############################

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

baseline = 0.1656        #[V], eyeseeking on the interval voltage_stored[5][0:index_min-1]

delta_V = 1/5 *10e-3                                     #[V] 

sum_voltage = sum(voltage_peak)         #[V] total voltage of the peak

#Storing of the current

voltage_peak_stored.append(sum_voltage)
baseline_stored.append(baseline)
peak_stored.append(peak)
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

#The amplitude of the signal wil be the difference between the baseline and the
#peak. The peak have both signs, but the amplitude can be compute similarly, baseline-peak
#to do so, will create numpy arrays that can do element wise operations!!

peak_stored = np.array(peak_stored)
baseline_stored = np.array(baseline_stored)
amplitude_max_signals = baseline_stored - peak_stored

#Now lets compute the error:
delta_single_V_measurement = np.array(delta_single_V_measurement) #conversion to np array
#len_baseline_signal_points_stored = np.array(len_baseline_signal_points_stored) #conversion to np array

#delta_amplitude_max = 1e3 * ( delta_single_V_measurement * np.sqrt(len_baseline_signal_points_stored+1) ) #[mV]
delta_amplitude_max = delta_single_V_measurement * np.sqrt(2)           #[V] error of the amplitude

print('Peak LYSO (V) ' + str(peak_stored[0]) + ' +- ' + str(delta_single_V_measurement[0]))
print('Peak CsI (V) ' + str(peak_stored[1]) + ' +- ' + str(delta_single_V_measurement[1]))
print('Peak BGO (V)' + str(peak_stored[2]) + ' +- ' + str(delta_single_V_measurement[2]) + '\n')

print('Baseline LYSO (V) ' + str(baseline_stored[0]) + ' +- ' + str(delta_single_V_measurement[0]))
print('Baseline CsI (V) ' + str(baseline_stored[0]) + ' +- ' + str(delta_single_V_measurement[1]))
print('Baseline BGO (V) ' + str(baseline_stored[0]) + ' +- ' + str(delta_single_V_measurement[2]) + '\n')

print('Amplitude max LYSO (mV) = ' + str(amplitude_max_signals[0]*1e3) + ' +- ' + str(delta_amplitude_max[0]*1e3))
print('Amplitude max CsI (mV) =' + str(amplitude_max_signals[1]*1e3) + ' +- ' + str(delta_amplitude_max[1]*1e3))
print('Amplitude maxBGO (mV) = ' + str(amplitude_max_signals[2]*1e3) + ' +- ' + str(delta_amplitude_max[2]*1e3) +'\n')



#And the amplitude ratio:
    #BGO/ CsI & LYSO/CsI & LYSO/BGO

ratio_ampl = np.array([amplitude_max_signals[2] / amplitude_max_signals[1], 
                      amplitude_max_signals[0] / amplitude_max_signals[1],
                      amplitude_max_signals[0] / amplitude_max_signals[2]
                      ])                #ratio of amplitudes, order BGO/CsI, LYSO/CsI, LYSO/BGO


delta_ratio_ampl = ratio_ampl * np.sqrt(2) * np.sqrt( (delta_single_V_measurement/ratio_ampl)**2 + 
                                                     (delta_single_V_measurement/ratio_ampl)**2)
    
print('Amplitude ratio BGO/CsI = ' + str(ratio_ampl[0]) + ' +- ' + str(delta_ratio_ampl[0]))
print('Amplitude ratio LYSO/CsI =' + str(ratio_ampl[1]) + ' +- ' + str(delta_ratio_ampl[1]))
print('Amplitude ratio LYSO/BGO = ' + str(ratio_ampl[2]) + ' +- ' + str(delta_ratio_ampl[2]) +'\n')








##########Light yield calc#############


# #So, we have the charge (voltage), which is proportional to the GM cells that
# #have been fired, which is also proportional to the number of photons that
# #have reach the detector, i.e., the photons generated (yield). Since we do
# #not know this proportionaly factor, we could just compute ratio of yields,
# #since this will null out this proportionality factors.

# delta_voltage_peak = np.multiply(delta_single_V_measurement, 
#                                   [np.sqrt(x) for x in n_elements_peak])       
#                     #[V] error of the voltage of the peak, the pseudo
#                     #light yield


# ratios = [voltage_peak_stored[0]/voltage_peak_stored[1],
#           voltage_peak_stored[0]/voltage_peak_stored[2],
#           voltage_peak_stored[2]/voltage_peak_stored[1],
#           ]   #LYSO/CsI, LYSO/BGO, BGO/CsI

# #error calculation of the ratio
# auxiliar = [(delta_voltage_peak[0]/voltage_peak_stored[0])**2 + (delta_voltage_peak[1]/voltage_peak_stored[1])**2, 
#             (delta_voltage_peak[0]/voltage_peak_stored[0])**2 + (delta_voltage_peak[2]/voltage_peak_stored[2])**2,
#             (delta_voltage_peak[2]/voltage_peak_stored[2])**2 + (delta_voltage_peak[1]/voltage_peak_stored[1])**2] 
#             #1st row: element 0 and 1
#             #2nd row: element 0 and 2
#             #3rd row: elememtn 1 and 2   (same order as the ratio!!!!)

#                     #the thing to be square rooted
# delta_ratios = np.multiply(ratios, [np.sqrt(x) for x in auxiliar])
#                                               #this will be though

# #Print:
# print('Light yield ratio LYSO/CsI: ' + str(ratios[0]) + ' +/- ' + str(delta_ratios[0]))
# print('Light yield ratio LYSO/BGO: ' + str(ratios[1]) + ' +/- ' + str(delta_ratios[1]))
# print('Light yield ratio BGO/CsI: ' + str(ratios[2]) + ' +/- ' + str(delta_ratios[2]))
# #LYSO/CsI, LYSO/BGO, CsI/BGO


# plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
# plt.title("Light yield ratio of the Cs137 peak for several scintillators", fontsize=22, wrap=True)           #title
# plt.bar(['LYSO/CsI', 'LYSO/BGO', 'BGO/CsI'], ratios, yerr = delta_ratios, edgecolor="black")#, yerr = delta_light_yield)
# #plt.xlabel("ADC channels", fontsize=10)                        #xlabel
# plt.ylabel("Light yield ratio", fontsize=14)              #ylabel
# # Set size of tick labels.
# plt.tick_params(axis='both', labelsize=14)              #size of axis
# plt.grid(True) 
# plt.savefig('Light_yield_ratio.png', format='png')


# plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
# plt.title("Pseudo Light yield of the Cs137 peak for several scintillators", fontsize=22, wrap=True)           #title
# plt.bar(['LYSO', 'CsI', 'BGO'], voltage_peak_stored, edgecolor="black")#, yerr = delta_light_yield)
# #plt.xlabel("ADC channels", fontsize=10)                        #xlabel
# plt.ylabel("Pseudo Light yield", fontsize=14)              #ylabel
# # Set size of tick labels.
# plt.tick_params(axis='both', labelsize=14)              #size of axis
# plt.grid(True) 
# plt.savefig('Pseudo_light_yield.png', format='png')



#CsI/BGO match,LYSO/CsI almost, LYSO/BGO do, should be greater than 1 an 
    #is lower than 1:(


################################################
#############################################
