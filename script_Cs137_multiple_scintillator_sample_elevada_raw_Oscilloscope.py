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
from scipy import stats as stats     #to find the most commom element in a list

import Read_csv_oscilloscope
######3

#plt.close("all")


############TO DO

#1) MAKE EVERY ARRAY NUMPY ARRAY, SINCE NUMPY ARRAYS ARE ABLE TO DO ELEMENT WISE OPERATIONS
    #EASILY, WHICH IS SOMETHING I NEED A FUCKING LOT!                   DONE :))
     



#%% #########################################################
#########################1), Data loading #####################
#############################################################

#The files to load are in .csv. A reader function has been created, following Manu's advice
#to make the code as short as possible to improve the readability for error debugging, the real
#programmer work. Bro, from 43 to 283 to from 43 to 120. Lol! (I have also removes unnecesary plots!)


#Variables that will store the results
voltage_stored = np.array(np.array([]))
time_stored = np.array([])


#################################
###########RAW##################3
################################

###LYSO
load = Read_csv_oscilloscope.Read_csv_oscilloscope('TEK0002_LYSO_raw_24_5_both.CSV')    

#Storing of the values
voltage_stored = np.append(voltage_stored,load[1])
time_stored = np.append(time_stored,load[0])
###       
      
###CsI
load = Read_csv_oscilloscope.Read_csv_oscilloscope('TEK0004_CsI_raw_24_5_both.CSV')    

#Storing of the values
voltage_stored = np.column_stack((voltage_stored,load[1]))
time_stored = np.column_stack((time_stored,load[0]))
        #have to write column stack so it creates columns!
        # voltage_stored = np.column_stack((voltage_stored, voltage_help))
###  

###BGO
load = Read_csv_oscilloscope.Read_csv_oscilloscope('TEK0000_BGO_raw_24_5_both.CSV')    

#Storing of the values
voltage_stored = np.column_stack((voltage_stored,load[1]))
time_stored = np.column_stack((time_stored,load[0]))
###          


#################################
###########PRE##################3
################################

###LYSO
load = Read_csv_oscilloscope.Read_csv_oscilloscope('TEK0003_LYSO_pre_24_5_both.CSV')    

#Storing of the values
voltage_stored = np.column_stack((voltage_stored,load[1]))
time_stored = np.column_stack((time_stored,load[0]))
###       
      
###CsI
load = Read_csv_oscilloscope.Read_csv_oscilloscope('TEK0005_CsI_pre_24_5_both.CSV')    

#Storing of the values
voltage_stored = np.column_stack((voltage_stored,load[1]))
time_stored = np.column_stack((time_stored,load[0]))
###  

###BGO
load = Read_csv_oscilloscope.Read_csv_oscilloscope('TEK0001_BGO_pre_24_5_both.CSV')    

#Storing of the values
voltage_stored = np.column_stack((voltage_stored,load[1]))
time_stored = np.column_stack((time_stored,load[0]))
###        


        
# # #Plot

# plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
# plt.plot(1e6 *time_stored,1e3 * voltage_stored, 'bo-')    #-1 chooses last element, which is the
#         #one that have been added to the lsit the lastest ;)    
#         #widht so that each bar touches each other!
# plt.title("Raw Waveform of Cs137 with LYSO", fontsize=22)           #title
# plt.xlabel("time (us)", fontsize=14)                        #xlabel
# plt.ylabel("voltage (mV)", fontsize=14)              #ylabel
# # Set size of tick labels.
# plt.tick_params(axis='both', labelsize=14)              #size of axis
# plt.grid(True) 
# #plt.xlim(0,max(ADC_channel))                       #limits of x axis     
# plt.savefig('Raw_signal_LYSO.png', format='png')



####Plot of all the waves#######

#####Raw

#since not all the waveforms have the same baseline, for a good plot, they should all
#have the same baseline, so the idea will be to move (substract) ti akk the waveforms to 
#decrease them to the waveform with the minimum baseline

baseline_raw = [
            min(voltage_stored[:,0]), min(voltage_stored[:,1]), min(voltage_stored[:,2])
                ]
                #LYSO, CSI, BGO

baseline_raw_min = min(baseline_raw)                           #min value of the baseline (raw)
baseline_raw_min_index = np.where(baseline_raw == baseline_raw_min)        #index of the min value

#now that we have the index (note the right order of the basline indexes, 0, 1, 2) we can substract.

#plot
plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
plt.title("Raw Waveforms of Cs137", fontsize=22)           #title

plt.plot(1e6 *time_stored[:,0], voltage_stored[:,0] - (baseline_raw[0]- baseline_raw_min), 'k-')
plt.plot(1e6 *time_stored[:,1], voltage_stored[:,1] - (baseline_raw[1]- baseline_raw_min), 'b-')
plt.plot(1e6 *time_stored[:,2], voltage_stored[:,2] - (baseline_raw[2]- baseline_raw_min), 'r-') 

plt.xlabel("time (us)", fontsize=14)                        #xlabel
plt.ylabel("voltage (V)", fontsize=14)              #ylabel
plt.legend(['LYSO', 'CsI', 'BGO'], fontsize=10) 
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
#plt.xlim(0,max(ADC_channel))                       #limits of x axis     
plt.savefig('Raw_signals.png', format='png')



######Pre
plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
plt.title("Waveforms of Cs137", fontsize=22)           #title

plt.plot(1e6 *time_stored[:,3], voltage_stored[:,3], 'k-')
plt.plot(1e6 *time_stored[:,4], voltage_stored[:,4], 'b-')
plt.plot(1e6 *time_stored[:,5], voltage_stored[:,5], 'r-')

plt.xlabel("time (us)", fontsize=14)                        #xlabel
plt.ylabel("voltage (V)", fontsize=14)              #ylabel
plt.legend(['LYSO', 'CsI', 'BGO'], fontsize=10) 
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
#plt.xlim(0,max(ADC_channel))                       #limits of x axis     
plt.savefig('Signals.png', format='png')





#%% ###############################################################################
########################### Rise and decay time#####################################
###################################################################################

#Usando misma configuracion, mido el máximo, del pre, q sera el del fotopico. Pre mejor
#q raw pq el raw tenia mas pileups q el pre. Mediciones tanto para la raw signal como para
#la del pre!!!.

#Variables that will store(st) the results 
t_rise_st = np.array([])
t_decay_st = np.array([])
delta_t_rise_st = np.array([])
delta_t_decay_st = np.array([])

###################################LYSO RAW########################

#RAW

t_rise = 70  #280#18.8e3                            #[ns]
t_decay = 890   #55*1e3                      #[ns], 3.28us (micros )
delta_t_rise = 1/5*50  #500                #[ns], error of the oscilloscope
                #it basically is 1/5 of the escale, since each square is divided
                #into 5 intervals, and the scale is the length of the sides of
                #the square
delta_t_decay = 1/5*250  #10e3          #[ns]
  
#storing
t_rise_st = np.append(t_rise_st,t_rise)
t_decay_st = np.append(t_decay_st,t_decay)
delta_t_rise_st = np.append(delta_t_rise_st,delta_t_rise)
delta_t_decay_st = np.append(delta_t_decay_st,delta_t_decay)

############################CsI RAW#####################3#####


t_rise = 380  #3.2e3                            #[ns]
t_decay = 1e3* 4.32 #145                      #[ns], 3.28us (micros )
delta_t_rise = 1/5*250 #5e3                   #[ns], error of the oscilloscope
delta_t_decay = 1/5*1e3 #25e3
  
#storing
t_rise_st = np.append(t_rise_st,t_rise)
t_decay_st = np.append(t_decay_st,t_decay)
delta_t_rise_st = np.append(delta_t_rise_st,delta_t_rise)
delta_t_decay_st = np.append(delta_t_decay_st,delta_t_decay)

####################################BGO RAW###################

t_rise = 176 #1.1e3                           #[ns]
t_decay = 1.54e3 #60.8e3                     #[ns], 3.28us (micros )
delta_t_rise = 1/5*100 #2.5e3                #[ns], error of the oscilloscope
delta_t_decay = 1/5*500 #10e3
  
#storing
t_rise_st = np.append(t_rise_st,t_rise)
t_decay_st = np.append(t_decay_st,t_decay)
delta_t_rise_st = np.append(delta_t_rise_st,delta_t_rise)
delta_t_decay_st = np.append(delta_t_decay_st,delta_t_decay)

##################################
#############PRE#################
################################
#This is the data from the new measurements, using the threshold to only measure
#the gamma peak. The old results are commented

#################################LYSO PRE##################3####

t_rise =  64#94#156  #28                           #[ns]
t_decay = 690 #680#464  #640                      #[ns], 3.28us (micros )
delta_t_rise = 1/5 * 50#100 #1/5*100                   #[ns], error of the oscilloscope
delta_t_decay = 1/5 * 250#100#1/5*100
  

#new waveform, not weird (20_5_21):
    #t rise = 70ns, escala 100ns==> error 50/5ns = 10ns
    #t decay = 530ns, error de 1/5 * 100ns = 20ns
    
    
#wave, no plitter, 24_5 (the last one is not this, but 2 waves together!!!)
    #t dec: 680ns, error 100/5 ns
    #t rise: 94 ns, error 1/5 50ns

#storing
t_rise_st = np.append(t_rise_st,t_rise)
t_decay_st = np.append(t_decay_st,t_decay)
delta_t_rise_st = np.append(delta_t_rise_st,delta_t_rise)
delta_t_decay_st = np.append(delta_t_decay_st,delta_t_decay)


###############################CsI PRE#########################

t_rise = 360#310 #320 #360                                      #[ns]
t_decay = 3.34e3#6.52e3 #2.88e3    #4.52e3                      #[ns], 3.28us (micros )
delta_t_rise = 1/5 * 250#250 #500#1/5*1e3                       #[ns], error of the oscilloscope
delta_t_decay = 1/5 * 500#1e3#500 #1/5*1e3
  
#storing
t_rise_st = np.append(t_rise_st,t_rise)
t_decay_st = np.append(t_decay_st,t_decay)
delta_t_rise_st = np.append(delta_t_rise_st,delta_t_rise)
delta_t_decay_st = np.append(delta_t_decay_st,delta_t_decay)

#new waveform, not weird (21_5_21):
    #t rise = 310ns,  error 250/5ns
    #t decay = 4.08us, error de 1us/5  

#wave, no plitter, 24_5
    #t dec: 6.52us, error 1/5 us
    #t rise: 310nsm error 250*1/5ns

###############################3###BGO Pre############################

t_rise = 184 #165 #104#160                           #[ns]
t_decay = 1.48e3#1.52e3 #1.15e3 #1.68e3                      #[ns], 3.28us (micros )
delta_t_rise = 1/5 * 100#1/5*1e3                  #[ns], error of the oscilloscope
delta_t_decay = 1/5 * 100#500 #250#1/5*1e3
  
#storing
t_rise_st = np.append(t_rise_st,t_rise)
t_decay_st = np.append(t_decay_st,t_decay)
delta_t_rise_st = np.append(delta_t_rise_st,delta_t_rise)
delta_t_decay_st = np.append(delta_t_decay_st,delta_t_decay)

#new waveform, not weird (21_5_21):
    #t rise = 165ns,  error 100/5ns 
    #t decay = 1.52us, error de 500ns/5  
    
#wave, no plitter, 24_5
    #t dec: 1.62 us, error 500/5 ns
    #t rise: 132 ns, error 1/5 100ns 

    
####Plot

#have to reorder the data so that they are arranged from lowest to highers: LYSO-->BGO-->CsI
#the current order (storing order) is LYSO-->CsI-->BGO. I must not change the storing order since
#it would affect everyhing, so I will alter the plot order!!

#RAW

plt.figure(figsize=(13,6))  #width, heigh 6.4*4.8 inches by default
plt.subplot(1, 2, 1)
plt.suptitle("Rise and decay time of the raw signal of the Cs137 waveform", fontsize=22, wrap=True)           #title
plt.bar(['LYSO', 'BGO', 'CsI'], np.array([t_rise_st[0], t_rise_st[2], t_rise_st[1] ])*1e-3, 
        yerr = np.array([delta_t_rise_st[0], delta_t_rise_st[2], delta_t_rise_st[1] ])*1e-3, 
        edgecolor="black")
#plt.xlabel("ADC channels", fontsize=10)                        #xlabel
plt.ylabel("Rise time (us)", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 

plt.subplot(1, 2, 2)
plt.bar(['LYSO', 'BGO', 'CsI'], np.array([t_decay_st[0], t_decay_st[2], t_decay_st[1] ])*1e-3, 
        yerr = np.array([delta_t_decay_st[0], delta_t_decay_st[2], delta_t_decay_st[1] ])*1e-3, 
        edgecolor="black")
#plt.xlabel("ADC channels", fontsize=10)                        #xlabel
plt.ylabel("Decay time (us)", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.savefig('Rise_decay_time_raw.png', format='png')


#Print:
print('Decay time of raw LYSO: (' + str(t_decay_st[0]*1e-3) + ' +/- ' + str(delta_t_decay_st[0]*1e-3) + ')us')
print('Decay time of raw BGO: (' + str(t_decay_st[2]*1e-3) + ' +/- ' + str(delta_t_decay_st[2]*1e-3) + ')us')
print('Decay time of raw CsI: (' + str(t_decay_st[1]*1e-3) + ' +/- ' + str(delta_t_decay_st[1]*1e-3) + ')us' +"\n")

#################################################

#Pre

plt.figure(figsize=(13,6))  #width, heigh 6.4*4.8 inches by default
plt.subplot(1, 2, 1)
plt.suptitle("Rise and decay time of the signal of the Cs137 waveform", fontsize=22, wrap=True)    #title
plt.bar(['LYSO', 'BGO', 'CsI'], np.array([t_rise_st[3], t_rise_st[5], t_rise_st[4] ])*1e-3, 
        yerr = np.array([delta_t_rise_st[3], delta_t_rise_st[5], delta_t_rise_st[4] ])*1e-3, 
        edgecolor="black")
#plt.xlabel("ADC channels", fontsize=10)                        #xlabel
plt.ylabel("Rise time (us)", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 

plt.subplot(1, 2, 2)
plt.bar(['LYSO', 'BGO', 'CsI'], np.array([t_decay_st[3], t_decay_st[5], t_decay_st[4] ])*1e-3, 
        yerr = np.array([delta_t_decay_st[3], delta_t_decay_st[5], delta_t_decay_st[4] ])*1e-3, 
        edgecolor="black")
#plt.xlabel("ADC channels", fontsize=10)                        #xlabel
plt.ylabel("Decay time (us)", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.savefig('Rise_decay_time_pre.png', format='png')


#Print:
print('Decay time of pre LYSO: (' + str(t_decay_st[3]*1e-3) + ' +/- ' + str(delta_t_decay_st[3]*1e-3) + ')us')
print('Decay time of pre BGO: (' + str(t_decay_st[5]*1e-3) + ' +/- ' + str(delta_t_decay_st[5]*1e-3) + ')us')
print('Decay time of pre CsI: (' + str(t_decay_st[4]*1e-3) + ' +/- ' + str(delta_t_decay_st[4]*1e-3) + ')us' +"\n")
##############################################

#Lets now comput the ratio of the decay time, of both signals, raw and pre.

#0,1,2 raw of LYSO, CsI, BGO respect
#3,4,5 pre of LYSO, CsI, BGO respect

ratio_decay = np.array( [t_decay_st[0]/t_decay_st[1], t_decay_st[0]/t_decay_st[2], 
                         t_decay_st[2]/t_decay_st[1], t_decay_st[3]/t_decay_st[4],
                         t_decay_st[3]/t_decay_st[5], t_decay_st[5]/t_decay_st[4]
                         ]  )
                  #The order is: LYSO/CsI,  LYSO/BGO,  BGO/CsI; first raw and then pre

    
#Error calculation of the ratio
auxiliar2 = np.array( [(delta_t_decay_st[0]/t_decay_st[0])**2 + (delta_t_decay_st[1]/t_decay_st[1])**2, 
            (delta_t_decay_st[0]/t_decay_st[0])**2 + (delta_t_decay_st[2]/t_decay_st[2])**2,
            (delta_t_decay_st[2]/t_decay_st[2])**2 + (delta_t_decay_st[1]/t_decay_st[1])**2,
            (delta_t_decay_st[3]/t_decay_st[3])**2 + (delta_t_decay_st[4]/t_decay_st[4])**2, 
            (delta_t_decay_st[3]/t_decay_st[3])**2 + (delta_t_decay_st[5]/t_decay_st[5])**2,
            (delta_t_decay_st[5]/t_decay_st[5])**2 + (delta_t_decay_st[4]/t_decay_st[4])**2,
            ] )
            #this are (delta_t1/t1)^2 +  (delta_t2/t2)^2, each row is this, in the order of the 
            #variable ratio_decay. This is basically the thing to be squared rooted to compute the 
            #error

delta_ratios = ratio_decay * np.sqrt(auxiliar2)
                                             #this will be though

print('Decay time ratio of raw BGO/CsI: (' + str(ratio_decay[2]) + ' +/- ' + str(delta_ratios[2]))
print('Decay time ratio of raw LYSO/CsI: (' + str(ratio_decay[0]) + ' +/- ' + str(delta_ratios[0]) )
print('Decay time ratio of raw LYSO/BGO: (' + str(ratio_decay[1]) + ' +/- ' + str(delta_ratios[1]) +"\n" )

print('Decay time ratio of pre BGO/CsI: (' + str(ratio_decay[5]) + ' +/- ' + str(delta_ratios[5]))
print('Decay time ratio of pre LYSO/CsI: (' + str(ratio_decay[3]) + ' +/- ' + str(delta_ratios[3]))
print('Decay time ratio of pre LYSO/BGO: (' + str(ratio_decay[4]) + ' +/- ' + str(delta_ratios[4]) +"\n")




 #%% ########################################################################
 ###################3) Calcs of the light yield############################################
###############################################################################

#Since the signals from the Pre are way better, will use them. Remember they are
#the last stored (6 total stored), so from 3 to 5 (0 the first). LYSO, CsI, BGO the order
#Note we are not sure about the light yield concept, but computing the amplitude will be good.

#I have also implemented the integral of the curve. I will use trapz,which is the easiest way

     # numpy.trapz(y, x, dx=1.0, axis=-1)[source]

     #    Integrate along the given axis using the composite trapezoidal rule.
     #    Integrate y (x) along given axis.

####################################LYSO##########################3

#TO find the peak I could easily do:
peak = min(voltage_stored[:,3])                           #[V] Peak value
index_peak = np.where(voltage_stored[:,3] == peak)        #index

#using the index of the peak I see the interval by looking at the .csv file
#and the plot
index_min = 660#499         #index that starts the peak (by looking the graph)
index_max = 1140#1800#820#774         #index that ends the peak (by looking the graph)


voltage_peak_neg = voltage_stored[:,3][index_min-1:index_max-1]  
time_peak = time_stored[:,3][index_min-1:index_max-1]
voltage_peak = np.array(np.absolute(voltage_peak_neg))   #abs() because this
        #peak contains both > and ,0 values, so to add them in order to count
        #them, I have to put everything in the positive value
        
len_peak =  len(voltage_stored[:,3][index_min-1:index_max-1])           #len of the peak   

#Baseline.
    #to choose the baseline, the 1st approach was to do the mean between channel 0 and the
    #channel where the peak start to appear. But this is abd idea becaue the error is very high,
    #so the 2nd approach will simply be choose the most frquent value

#len_baseline_signal_points = len(voltage_stored[3][0:index_min-1])  #length of the voltages used
        #to compute the baseline. Will be neccesary for error calcs
#baseline_signal = sum(voltage_stored[3][0:index_min-1]) / len_baseline_signal_points    
        #[V] baseline, to compute the peak amplitude. I average between the initial value and
        #when the peak starts
        

baseline = stats.mode(voltage_stored[:,3])[0]       #[V] baseline voltage, the most common value
                    #[0]contains the value, [1] the frequency
                    
delta_V = 1/5 *1#100e-3           #[V] error of the voltage measurements, from the scale of
                            #the oscilloscope, which come from the photos (pre) (.csv)
 
delta_t = 1/5 * 250e-9                              #[s] error of the times in the .csv, 
                                            #from the oscilloscope

        
sum_voltage = sum(-voltage_peak + baseline)         #[V] total voltage of the peak, corrected with
                        #the baseline, so that this is the real voltage created!! Baseline-voltage because
                        #voltage of the amplitude is the greatest

integral = np.trapz(voltage_peak, time_peak)        #[V*t] area under the peak
delta_integral = 2 * np.sqrt( delta_V / (max(voltage_peak) - min(voltage_peak) ) 
                             + delta_t/ (max(time_peak) - min(time_peak) ) )    #overstimation of the
                #error of the integral. This is the error of the 
                #area of the rectangle (Vmax-Vmin)*(tmax-tmin)



#Storing

voltage_peak_st = np.array(np.array([]))                #storage of the total voltage of the peak,
                    #removing the baseline!!!!!!

baseline_st = np.array([])               #storage of the baseline
peak_st = np.array([])                  #storage of the peak value, for the max amplitude
n_elements_peak_st = np.array([])                        #this will store the number of voltages I sum,
                                            #for each peak, for the error calc
delta_single_V_measurement = np.array([])         #storage of the error of the voltage 
                                                #measurements, for the error calc
                                                
integral_st = np.array([])                      #peak integral
delta_integral_st = np.array([])                #error of the peak integral (overstimation)

#len_baseline_signal_points_stored = np.array([])                    
                    

voltage_peak_st = np.append(voltage_peak_st,sum_voltage)
baseline_st = np.append(baseline_st,baseline)
peak_st = np.append(peak_st,peak)
n_elements_peak_st = np.append(n_elements_peak_st,len_peak)
delta_single_V_measurement = np.append(delta_single_V_measurement,delta_V)
integral_st = np.append(integral_st, integral)
delta_integral_st = np.append(delta_integral_st, delta_integral)

#len_baseline_signal_points_stored.append(len_baseline_signal_points)

#Plot (debug)
#plt.plot([1e6 * x for x in time_stored[-1]],[1e3 * x for x in voltage_stored[-1]], 'bo-')    #-1 chooses last element, which is the

plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
plt.plot( 1e6 * time_stored[:,3], 1e3 * voltage_stored[:,3], 'b.-')    
plt.plot( 1e6 * time_peak , 1e3 * voltage_peak_neg, 'r.-')   
        #widht so that each bar touches each other!
plt.title("Waveform of Cs137 with LYSO", fontsize=22)           #title
plt.xlabel("time (us)", fontsize=14)                        #xlabel
plt.ylabel("voltage (mV)", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.savefig('Signal_LYSO_sum.png', format='png')


###############################################CsI#####################

#TO find the peak I could easily do:
peak = min(voltage_stored[:,4])                           #[V] Peak value
index_peak = np.where(voltage_stored[:,4] == peak)        #index

#using the index of the peak I see the interval by looking at the .csv file
#and the plot

index_min = 510#500#542
index_max = 2200#1900#2100#1921         #index that ends the peak (by looking the graph)

voltage_peak_neg = voltage_stored[:,4][index_min-1:index_max-1]  
time_peak = time_stored[:,4][index_min-1:index_max-1]
voltage_peak = np.array(np.absolute(voltage_peak_neg))   

len_peak =  len(voltage_stored[:,4][index_min-1:index_max-1])           #len of the peak  


baseline = stats.mode(voltage_stored[0:index_min-1,4])[0]       #[V] baseline voltage, the most common value
    #WHATCH OUT! Since here there is a baseline shift, the interval in which I should choose the most
    #common value should be until the peak appears

delta_V = 1/5 *500e-3#20e-3                        #[V] 
delta_t = 1/5 * 500e-9                             #[s] error of the times in the .csv, 
                                        #from the oscilloscope

sum_voltage = sum(-voltage_peak + baseline)         #[V] total voltage of the peak, corrected with
                        #the baseline, so that this is the real voltage created!!
integral = np.trapz(voltage_peak, time_peak)        #[V*t] area under the peak
delta_integral = 2 * np.sqrt( delta_V / (max(voltage_peak) - min(voltage_peak) ) 
                             + delta_t/ (max(time_peak) - min(time_peak) ) )    

#Storing

voltage_peak_st = np.append(voltage_peak_st,sum_voltage)
baseline_st = np.append(baseline_st,baseline)
peak_st = np.append(peak_st,peak)
n_elements_peak_st = np.append(n_elements_peak_st,len_peak)
delta_single_V_measurement = np.append(delta_single_V_measurement,delta_V)
integral_st = np.append(integral_st, integral)
delta_integral_st = np.append(delta_integral_st, delta_integral)


#Plot (debug)
plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
plt.plot( 1e6 * time_stored[:,4], 1e3 * voltage_stored[:,4], 'b.-')    
plt.plot( 1e6 * time_peak , 1e3 * voltage_peak_neg, 'r.-')   
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
peak = min(voltage_stored[:,5])                           #[V] Peak value
index_peak = np.where(voltage_stored[:,5] == peak)        #index

#using the index of the peak I see the interval by looking at the .csv file
#and the plot

index_min = 530#500#585
index_max = 1000#1450#1300#1060         #index that ends the peak (by looking the graph)

voltage_peak_neg = voltage_stored[:,5][index_min-1:index_max-1]  
time_peak = time_stored[:,5][index_min-1:index_max-1]
voltage_peak = np.array(np.absolute(voltage_peak_neg))

len_peak =  len(voltage_stored[:,5][index_min-1:index_max-1])           #len of the peak  

baseline = stats.mode(voltage_stored[0:index_min-1,5])[0]       #[V] baseline voltage, the most common value

delta_V = 1/5 *200e-3#10e-3                                     #[V] 
delta_t = 1/5 *500e-9                                                #[s]

sum_voltage = sum(-voltage_peak + baseline)          #[V] total voltage of the peak, corrected with
                        #the baseline, so that this is the real voltage created!!
integral = np.trapz(voltage_peak, time_peak)        #[V*t] area under the peak      
delta_integral = 2 * np.sqrt( delta_V / (max(voltage_peak) - min(voltage_peak) ) 
                             + delta_t/ (max(time_peak) - min(time_peak) ) )   

#Storing of the current

voltage_peak_st = np.append(voltage_peak_st,sum_voltage)
baseline_st = np.append(baseline_st,baseline)
peak_st = np.append(peak_st,peak)
n_elements_peak_st = np.append(n_elements_peak_st,len_peak)
delta_single_V_measurement = np.append(delta_single_V_measurement,delta_V)
integral_st = np.append(integral_st, integral)
delta_integral_st = np.append(delta_integral_st, delta_integral)


#Plot (debug)
plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
plt.plot( 1e6 * time_stored[:,5], 1e3 * voltage_stored[:,5], 'b.-')    
plt.plot( 1e6 * time_peak , 1e3 * voltage_peak_neg, 'r.-')     
        #widht so that each bar touches each other!
plt.title("Waveform of Cs137 with BGO", fontsize=22)           #title
plt.xlabel("time (us)", fontsize=14)                        #xlabel
plt.ylabel("voltage (mV)", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.savefig('Signal_BGO_sum.png', format='png')


##########################################################
#########################################################


#The amplitude of the signal wil be the difference between the baseline and the
#peak. The peak are negatives, while the baseline are positive, so will added them

max_voltage_signals = baseline_st + peak_st           #max value of the voltage of the peak

delta_max_voltage = delta_single_V_measurement * np.sqrt(2)           
            #[V] error of the amplitude = peak-baseline

delta_voltage_peak_st = np.sqrt(n_elements_peak_st) * delta_single_V_measurement 
                            #[V] error of the total voltage of the peak (integration), all voltages
                            #have the same error


#Amplitudes (voltage integration)
print('Voltage integration LYSO (baseline removed) (V) = ' + str(voltage_peak_st[0]) + ' +- ' + str(delta_voltage_peak_st[0]))
print('Voltage integration BGO (baseline removed)  (V) = ' + str(voltage_peak_st[2]) + ' +- ' + str(delta_voltage_peak_st[2]))
print('Voltage integration CsI (baseline removed)  (V) =' + str(voltage_peak_st[1]) + ' +- ' + str(delta_voltage_peak_st[1]) +'\n')

print('Amplitude (Vmax) LYSO (V) = ' + str(max_voltage_signals[0]) + ' +- ' + str(delta_max_voltage[0]))
print('Amplitude (Vmax) BGO (V) = ' + str(max_voltage_signals[2]) + ' +- ' + str(delta_max_voltage[2]))
print('Amplitude (Vmax) CsI (V) =' + str(max_voltage_signals[1]) + ' +- ' + str(delta_max_voltage[1]) +'\n')

print('Integral LYSO (baseline removed) (V*s) = ' + str(integral_st[0]) + ' +- ' + str(delta_integral_st[0]) )
print('Integral BGO (baseline removed) (V*s) =' + str(integral_st[2]) + ' +- ' + str(delta_integral_st[2]) )
print('Integral CsI (baseline removed) (V*s) =' + str(integral_st[1]) + ' +- ' + str(delta_integral_st[1]) + '\n')


#And the ratio of the maximum voltage:
    #BGO/ CsI & LYSO/CsI & LYSO/BGO

ratio_voltage_max = np.array([max_voltage_signals[2] / max_voltage_signals[1], 
                      max_voltage_signals[0] / max_voltage_signals[1],
                      max_voltage_signals[0] / max_voltage_signals[2]
                      ])                #ratio of amplitudes, order BGO/CsI, LYSO/CsI, LYSO/BGO


delta_ratio_voltage_max= ratio_voltage_max * np.sqrt( np.array(
    [(delta_max_voltage[2]/max_voltage_signals[2])**2 + (delta_max_voltage[1]/max_voltage_signals[1])**2,
      (delta_max_voltage[0]/max_voltage_signals[0])**2 + (delta_max_voltage[1]/max_voltage_signals[1])**2,
      (delta_max_voltage[0]/max_voltage_signals[0])**2 + (delta_max_voltage[2]/max_voltage_signals[2])**2,
              ]
    )
    )   #computing the error. The sqrt of the quotient have to be computed by hand :))


print('Amplitude ratio BGO/CsI = ' + str(ratio_voltage_max[0]) + ' +- ' + str(delta_ratio_voltage_max[0]))
print('Amplitude ratio LYSO/CsI =' + str(ratio_voltage_max[1]) + ' +- ' + str(delta_ratio_voltage_max[1]))
print('Amplitude ratio LYSO/BGO = ' + str(ratio_voltage_max[2]) + ' +- ' + str(delta_ratio_voltage_max[2]) +'\n')

#And the amplitude ratio, which should be related with the light yield, is:

ratio_total_ampl = np.array([voltage_peak_st[2] / voltage_peak_st[1], 
                      voltage_peak_st[0] / voltage_peak_st[1],
                     voltage_peak_st[0] / voltage_peak_st[2]
                      ])

delta_ratio_total_ampl = ratio_total_ampl * np.sqrt( np.array(
    [(delta_voltage_peak_st[2]/voltage_peak_st[2])**2 + (delta_voltage_peak_st[1]/voltage_peak_st[1])**2,
     (delta_voltage_peak_st[0]/voltage_peak_st[0])**2 + (delta_voltage_peak_st[1]/voltage_peak_st[1])**2,
     (delta_voltage_peak_st[0]/voltage_peak_st[0])**2 + (delta_voltage_peak_st[2]/voltage_peak_st[2])**2,
             ] ) )   #computing the error. The sqrt of the quotient have to be computed by hand :))


print('Total amplitude ratio BGO/CsI = ' + str(ratio_total_ampl[0]) + ' +- ' + str(delta_ratio_total_ampl[0]))
print('Total amplitude ratio LYSO/CsI =' + str(ratio_total_ampl[1]) + ' +- ' + str(delta_ratio_total_ampl[1]))
print('Total amplitude ratio LYSO/BGO = ' + str(ratio_total_ampl[2]) + ' +- ' + str(delta_ratio_total_ampl[2]) +'\n')

#Integral ratio

ratio_int = np.array([integral_st[2] / integral_st[1], 
                      integral_st[0] / integral_st[1],
                     integral_st[0] / integral_st[2]
                      ])

delta_ratio_int = ratio_total_ampl * np.sqrt( np.array(
    [(delta_integral_st[2]/integral_st[2])**2 + (delta_integral_st[1]/integral_st[1])**2,
     (delta_integral_st[0]/integral_st[0])**2 + (delta_integral_st[1]/integral_st[1])**2,
     (delta_integral_st[0]/integral_st[0])**2 + (delta_integral_st[2]/integral_st[2])**2,
             ] ) )   #computing the error. The sqrt of the quotient have to be computed by hand :))

print('Integral (no baseline) ratio BGO/CsI = ' + str(ratio_int[0]) + ' +- ' + str(delta_ratio_int[0]) )
print('Integral (no baseline) ratio  LYSO/CsI =' + str(ratio_int[1]) + ' +- ' + str(delta_ratio_int[1]) )
print('Integral (no baseline) ratio  LYSO/BGO = ' + str(ratio_int[2]) + ' +- ' + str(delta_ratio_int[2]) +'\n')
















#######################RESIDUO: Light yield calc######################################


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
