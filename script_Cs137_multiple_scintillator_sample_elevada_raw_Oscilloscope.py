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

import sys                   #to import functions from other folders!!
sys.path.insert(0, '/home/dla/Python/Functions_homemade')   #path where I have the functions

import Read_csv_oscilloscope
import Peak_analysis_oscilloscope
######3

#plt.close("all")




#%% #########################################################
#########################1), Data loading #####################
#############################################################

#The files to load are in .csv. A reader function has been created, following Manu's advice
#to make the code as short as possible to improve the readability for error debugging, the real
#programmer work. Bro, from 43 to 283 to from 43 to 120. Lol! (I have also removes unnecesary plots!)


#Variables that will store the results
voltage_stored = np.array(np.array([]))
time_stored = np.array([])



########LYSO RAW########
load = Read_csv_oscilloscope.Read_csv_oscilloscope('TEK0002.CSV')    

#Storing of the values
voltage_stored = np.append(voltage_stored,load[1])
time_stored = np.append(time_stored,load[0])     

  
#########CsI RAW##########
load = Read_csv_oscilloscope.Read_csv_oscilloscope('TEK0004.CSV')    

#Storing of the values
voltage_stored = np.column_stack((voltage_stored,load[1]))
time_stored = np.column_stack((time_stored,load[0]))
        #have to write column stack so it creates columns!
        # voltage_stored = np.column_stack((voltage_stored, voltage_help)) 


############BGO RAW#################
load = Read_csv_oscilloscope.Read_csv_oscilloscope('TEK0000.CSV')    

#Storing of the values
voltage_stored = np.column_stack((voltage_stored,load[1]))
time_stored = np.column_stack((time_stored,load[0]))       


#############LYSO PRE####################
load = Read_csv_oscilloscope.Read_csv_oscilloscope('TEK0003.CSV')    

#Storing of the values
voltage_stored = np.column_stack((voltage_stored,load[1]))
time_stored = np.column_stack((time_stored,load[0]))      
   
   
#########CsI PRE####################3
load = Read_csv_oscilloscope.Read_csv_oscilloscope('TEK0005.CSV')    

#Storing of the values
voltage_stored = np.column_stack((voltage_stored,load[1]))
time_stored = np.column_stack((time_stored,load[0]))

##########BGO PRE##################
load = Read_csv_oscilloscope.Read_csv_oscilloscope('TEK0001.CSV')    

#Storing of the values
voltage_stored = np.column_stack((voltage_stored,load[1]))
time_stored = np.column_stack((time_stored,load[0]))       


        
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



###############Plot of all the waves##################################

#####Raw####

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
plt.title("Raw Waveforms of the Cs137 photopeak", fontsize=22)           #title

plt.plot(1e6 *time_stored[:,0], voltage_stored[:,0] - (baseline_raw[0]- baseline_raw_min), 'k-')
plt.plot(1e6 *time_stored[:,1], voltage_stored[:,1] - (baseline_raw[1]- baseline_raw_min), 'b-')
plt.plot(1e6 *time_stored[:,2], voltage_stored[:,2] - (baseline_raw[2]- baseline_raw_min), 'r-') 

plt.xlabel("time (us)", fontsize=14)                        #xlabel
plt.ylabel("voltage (V)", fontsize=14)              #ylabel
plt.legend(['LYSO', 'CsI', 'BGO'], fontsize=16) 
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
#plt.xlim(0,max(ADC_channel))                       #limits of x axis     
plt.savefig('Raw_signals.png', format='png')



######Pre####

plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
plt.title("Waveforms of the Cs137 photopeak", fontsize=22)           #title

plt.plot(1e6 *time_stored[:,3], voltage_stored[:,3], 'k-')
plt.plot(1e6 *time_stored[:,4], voltage_stored[:,4], 'b-')
plt.plot(1e6 *time_stored[:,5], voltage_stored[:,5], 'r-')

plt.xlabel("time (us)", fontsize=14)                        #xlabel
plt.ylabel("voltage (V)", fontsize=14)              #ylabel
plt.legend(['LYSO', 'CsI', 'BGO'], fontsize=16) 
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
#plt.xlim(0,max(ADC_channel))                       #limits of x axis     
plt.savefig('Signals.png', format='png')



###############Plot of both waves for each crystal######################

#LYSO
plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
plt.title("Waveforms of the Cs137 photopeak using LYSO", fontsize=22)           #title

plt.plot(1e6 *time_stored[:,0], voltage_stored[:,0] )
plt.plot(1e6 *time_stored[:,3], voltage_stored[:,3] -.15)
plt.xlabel("time (us)", fontsize=14)                        #xlabel
plt.ylabel("voltage (V)", fontsize=14)              #ylabel
plt.legend(['from detector', 'from SP5600'], fontsize=14, loc = 'lower right') 
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=16)              #size of axis
plt.grid(True) 
#plt.xlim(0,max(ADC_channel))                       #limits of x axis    
plt.savefig('Waves_LYSO.png', format='png')


#CsI
plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
plt.title("Waveforms of the Cs137 photopeak using CsI", fontsize=22)           #title

plt.plot(1e6 *time_stored[:,1], voltage_stored[:,1] -.8 )
plt.plot(1e6 *time_stored[:,4], voltage_stored[:,4] )
plt.xlabel("time (us)", fontsize=14)                        #xlabel
plt.ylabel("voltage (V)", fontsize=14)              #ylabel
plt.legend(['from detector', 'from SP5600'], fontsize=16, loc = 'lower right') 
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
#plt.xlim(0,max(ADC_channel))                       #limits of x axis   
plt.savefig('Waves_CsI.png', format='png')


#BGO
plt.figure(figsize=(10,6))  #width, heigh 6.4*4.8 inches by default
plt.title("Waveforms of the Cs137 photopeak using BGO", fontsize=22)           #title

plt.plot(1e6 *time_stored[:,2], voltage_stored[:,2] -1.2 )
plt.plot(1e6 *time_stored[:,5], voltage_stored[:,5] )
plt.xlabel("time (us)", fontsize=14)                        #xlabel
plt.ylabel("voltage (V)", fontsize=14)              #ylabel
plt.legend(['from detector', 'from SP5600'], fontsize=16, loc = 'lower right')  
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
#plt.xlim(0,max(ADC_channel))                       #limits of x axis   
plt.savefig('Waves_BGO.png', format='png')



 #%% ########################################################################
 ###################3) Pulse analysis ############################################
###############################################################################

#Since the signals from the Pre are way better, will use them. Remember they are
#the last stored (6 total stored), so from 3 to 5 (0 the first). LYSO, CsI, BGO the order
#Note we are not sure about the light yield concept, but computing the amplitude will be good.

#I have also implemented the integral of the curve. I will use trapz,which is the easiest way

     # numpy.trapz(y, x, dx=1.0, axis=-1)[source]

     #    Integrate along the given axis using the composite trapezoidal rule.
     #    Integrate y (x) along given axis.

voltage_peak_st = np.array(np.array([]))                #storage of the total voltage of the peak,
                    #removing the baseline!!!!!!

baseline_st = np.array([])               #storage of the baseline
peak_st = np.array([])                  #storage of the peak value, for the max amplitude
n_elements_peak_st = np.array([])                        #this will store the number of voltages I sum,
                                            #for each peak, for the error calc
delta_single_V_measurement = np.array([])         #storage of the error of the voltage 
                                                #measurements, for the error calc
amplitude_st = np.array( [] )               #storage of the amplitudes                                           
delta_amplitude_st = np.array([])           #storage of \Delta(amplitude)                                                
integral_st = np.array([])                      #peak integral
delta_integral_st = np.array([])                #error of the peak integral (overstimation)

t_rise_st = np.array([])
t_decay_st = np.array([])
delta_t_rise_st = np.array([])
delta_t_decay_st = np.array([])



#note that the function that analyses the peak are about 50 lines. Multiplying buy 6, we have
#reduced the code is 300 lines!!!!!!!!!!


##################
#################
#################
#%% 
####################################LYSO PRE ##########################3

LYSO_pre = Peak_analysis_oscilloscope.Peak_analysis_oscillo(voltage_stored,time_stored, 'pre', 3, 495, 1350, 
                                                 1/5 *1, 1/5 * 250e-9 )

#Debug
#t_decayy = time_stored[1140,3]*1e6 - time_stored[703,3] *1e6 
#     			#time when the peak ends - time when the max value is reached
# t_rise = -time_stored[660,3]*1e6 + time_stored[703,3]*1e6  #peak - start
# 		# -time when the peak starts + time when the max value is reached


#Storing

baseline_st = np.append(baseline_st, LYSO_pre['baseline[V]'] )  
peak_st = np.append(peak_st, LYSO_pre['|peak|[V]'] )
n_elements_peak_st = np.append(n_elements_peak_st, LYSO_pre['N_peak'] )
delta_single_V_measurement = np.append(delta_single_V_measurement, LYSO_pre['\Delta(V[V])'] )
integral_st = np.append(integral_st, LYSO_pre['integral[V*s]'] )
delta_integral_st = np.append(delta_integral_st, LYSO_pre['\Delta(integral[V*s])'] )
amplitude_st = np.append(amplitude_st, LYSO_pre['amplitude[V]'])
delta_amplitude_st = np.append(delta_amplitude_st, LYSO_pre['\Delta(amplitude[V])'] )
t_rise_st = np.append(t_rise_st,LYSO_pre['t_rise[s]'])
t_decay_st = np.append(t_decay_st,LYSO_pre['t_decay[s]'])
delta_t_rise_st = np.append(delta_t_rise_st,LYSO_pre['\Delta(t_rise[s])'])
delta_t_decay_st = np.append(delta_t_decay_st,LYSO_pre['\Delta(t_decay[s])'])


###############################################CsI PRE   #####################

CsI_pre = Peak_analysis_oscilloscope.Peak_analysis_oscillo(voltage_stored,time_stored,'pre',4, 490, 1900, 
                                                 1/5 *500e-3, 1/5 * 500e-9 )

#Storing

baseline_st = np.append(baseline_st, CsI_pre['baseline[V]'] )  
peak_st = np.append(peak_st, CsI_pre['|peak|[V]'] )
n_elements_peak_st = np.append(n_elements_peak_st, CsI_pre['N_peak'] )
delta_single_V_measurement = np.append(delta_single_V_measurement, CsI_pre['\Delta(V[V])'] )
integral_st = np.append(integral_st, CsI_pre['integral[V*s]'] )
delta_integral_st = np.append(delta_integral_st, CsI_pre['\Delta(integral[V*s])'] )
amplitude_st = np.append(amplitude_st, CsI_pre['amplitude[V]'])
delta_amplitude_st = np.append(delta_amplitude_st, CsI_pre['\Delta(amplitude[V])'] )
t_rise_st = np.append(t_rise_st,CsI_pre['t_rise[s]'])
t_decay_st = np.append(t_decay_st,CsI_pre['t_decay[s]'])
delta_t_rise_st = np.append(delta_t_rise_st,CsI_pre['\Delta(t_rise[s])'])
delta_t_decay_st = np.append(delta_t_decay_st,CsI_pre['\Delta(t_decay[s])'])




#######################################BGO PRE ##############################

BGO_pre = Peak_analysis_oscilloscope.Peak_analysis_oscillo(voltage_stored,time_stored,'pre',5, 550, 1500, 
                                                 1/5 *200e-3, 1/5 * 500e-9 )

#Storing

baseline_st = np.append(baseline_st, BGO_pre['baseline[V]'] )  
peak_st = np.append(peak_st, BGO_pre['|peak|[V]'] )
n_elements_peak_st = np.append(n_elements_peak_st, BGO_pre['N_peak'] )
delta_single_V_measurement = np.append(delta_single_V_measurement, BGO_pre['\Delta(V[V])'] )
integral_st = np.append(integral_st, BGO_pre['integral[V*s]'] )
delta_integral_st = np.append(delta_integral_st, BGO_pre['\Delta(integral[V*s])'] )
amplitude_st = np.append(amplitude_st, BGO_pre['amplitude[V]'])
delta_amplitude_st = np.append(delta_amplitude_st, BGO_pre['\Delta(amplitude[V])'] )
t_rise_st = np.append(t_rise_st,BGO_pre['t_rise[s]'])
t_decay_st = np.append(t_decay_st,BGO_pre['t_decay[s]'])
delta_t_rise_st = np.append(delta_t_rise_st,BGO_pre['\Delta(t_rise[s])'])
delta_t_decay_st = np.append(delta_t_decay_st,BGO_pre['\Delta(t_decay[s])'])



######################
################RAW####
#######################

####################################LYSO RAW  ##########################3

LYSO_raw = Peak_analysis_oscilloscope.Peak_analysis_oscillo(voltage_stored,time_stored,'raw',0, 480, 1350, 
                                                 1/5 *50e-3, 1/5 * 250e-9 )

#Storing

baseline_st = np.append(baseline_st, LYSO_raw['baseline[V]'] )  
peak_st = np.append(peak_st, LYSO_raw['|peak|[V]'] )
n_elements_peak_st = np.append(n_elements_peak_st, LYSO_raw['N_peak'] )
delta_single_V_measurement = np.append(delta_single_V_measurement, LYSO_raw['\Delta(V[V])'] )
integral_st = np.append(integral_st, LYSO_raw['integral[V*s]'] )
delta_integral_st = np.append(delta_integral_st, LYSO_raw['\Delta(integral[V*s])'] )
amplitude_st = np.append(amplitude_st, LYSO_raw['amplitude[V]'])
delta_amplitude_st = np.append(delta_amplitude_st, LYSO_raw['\Delta(amplitude[V])'] )
t_rise_st = np.append(t_rise_st,LYSO_raw['t_rise[s]'])
t_decay_st = np.append(t_decay_st,LYSO_raw['t_decay[s]'])
delta_t_rise_st = np.append(delta_t_rise_st,LYSO_raw['\Delta(t_rise[s])'])
delta_t_decay_st = np.append(delta_t_decay_st,LYSO_raw['\Delta(t_decay[s])'])



####################################CsI RAW ##########################3

CsI_raw = Peak_analysis_oscilloscope.Peak_analysis_oscillo(voltage_stored,time_stored,'raw',1, 480, 2100, 
                                                 1/5 *20e-3, 1/5 * 1e-6 )

#Storing

baseline_st = np.append(baseline_st, CsI_raw['baseline[V]'] )  
peak_st = np.append(peak_st, CsI_raw['|peak|[V]'] )
n_elements_peak_st = np.append(n_elements_peak_st, CsI_raw['N_peak'] )
delta_single_V_measurement = np.append(delta_single_V_measurement, CsI_raw['\Delta(V[V])'] )
integral_st = np.append(integral_st, CsI_raw['integral[V*s]'] )
delta_integral_st = np.append(delta_integral_st, CsI_raw['\Delta(integral[V*s])'] )
amplitude_st = np.append(amplitude_st, CsI_raw['amplitude[V]'])
delta_amplitude_st = np.append(delta_amplitude_st, CsI_raw['\Delta(amplitude[V])'] )
t_rise_st = np.append(t_rise_st,CsI_raw['t_rise[s]'])
t_decay_st = np.append(t_decay_st,CsI_raw['t_decay[s]'])
delta_t_rise_st = np.append(delta_t_rise_st,CsI_raw['\Delta(t_rise[s])'])
delta_t_decay_st = np.append(delta_t_decay_st,CsI_raw['\Delta(t_decay[s])'])



####################################BGO RAW##########################3

BGO_raw = Peak_analysis_oscilloscope.Peak_analysis_oscillo(voltage_stored,time_stored,'raw',2, 540, 1500, 
                                                 1/5 *10e-3, 1/5 * 500e-9 )

#Storing

baseline_st = np.append(baseline_st, BGO_raw['baseline[V]'] )  
peak_st = np.append(peak_st, BGO_raw['|peak|[V]'] )
n_elements_peak_st = np.append(n_elements_peak_st, BGO_raw['N_peak'] )
delta_single_V_measurement = np.append(delta_single_V_measurement, BGO_raw['\Delta(V[V])'] )
integral_st = np.append(integral_st, BGO_raw['integral[V*s]'] )
delta_integral_st = np.append(delta_integral_st, BGO_raw['\Delta(integral[V*s])'] )
amplitude_st = np.append(amplitude_st, CsI_raw['amplitude[V]'])
delta_amplitude_st = np.append(delta_amplitude_st, CsI_raw['\Delta(amplitude[V])'] )
t_rise_st = np.append(t_rise_st,BGO_raw['t_rise[s]'])
t_decay_st = np.append(t_decay_st,BGO_raw['t_decay[s]'])
delta_t_rise_st = np.append(delta_t_rise_st,BGO_raw['\Delta(t_rise[s])'])
delta_t_decay_st = np.append(delta_t_decay_st,BGO_raw['\Delta(t_decay[s])'])


##############################################################################
##################PLOT OF TIMES####################################################
###################################################################################


#have to reorder the data so that they are arranged from lowest to highers: LYSO-->BGO-->CsI
#the current order (storing order) is LYSO-->CsI-->BGO. I must not change the storing order since
#it would affect everyhing, so I will alter the plot order!!

#RAW

plt.figure(figsize=(13,6))  #width, heigh 6.4*4.8 inches by default
plt.subplot(1, 2, 1)
plt.suptitle("Rise and decay time of the signals coming from the crystals", fontsize=22, wrap=True)           #title
plt.bar(['LYSO', 'BGO', 'CsI'], np.array([t_rise_st[0], t_rise_st[2], t_rise_st[1] ])*1e6, 
        yerr = np.array([delta_t_rise_st[0], delta_t_rise_st[2], delta_t_rise_st[1] ])*1e6, 
        edgecolor="black")
#plt.xlabel("ADC channels", fontsize=10)                        #xlabel
plt.ylabel("Rise time (us)", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 

plt.subplot(1, 2, 2)
plt.bar(['LYSO', 'BGO', 'CsI'], np.array([t_decay_st[0], t_decay_st[2], t_decay_st[1] ])*1e6, 
        yerr = np.array([delta_t_decay_st[0], delta_t_decay_st[2], delta_t_decay_st[1] ])*1e6, 
        edgecolor="black")
#plt.xlabel("ADC channels", fontsize=10)                        #xlabel
plt.ylabel("Decay time (us)", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.savefig('Rise_decay_time_pre.png', format='png')


#Print:
print('Decay time of pre LYSO: (' + str(t_decay_st[0]*1e6) + ' +/- ' + str(delta_t_decay_st[0]*1e6) + ')us')
print('Decay time of pre BGO: (' + str(t_decay_st[2]*1e6) + ' +/- ' + str(delta_t_decay_st[2]*1e6) + ')us')
print('Decay time of pre CsI: (' + str(t_decay_st[1]*1e6) + ' +/- ' + str(delta_t_decay_st[1]*1e6) + ')us' +"\n")

#################################################

#Pre

plt.figure(figsize=(13,6))  #width, heigh 6.4*4.8 inches by default
plt.subplot(1, 2, 1)
plt.suptitle("Rise and decay time of the raw signals coming from the crystals", fontsize=22, wrap=True)    #title
plt.bar(['LYSO', 'BGO', 'CsI'], np.array([t_rise_st[3], t_rise_st[5], t_rise_st[4] ])*1e6, 
        yerr = np.array([delta_t_rise_st[3], delta_t_rise_st[5], delta_t_rise_st[4] ])*1e6, 
        edgecolor="black")
#plt.xlabel("ADC channels", fontsize=10)                        #xlabel
plt.ylabel("Rise time (us)", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 

plt.subplot(1, 2, 2)
plt.bar(['LYSO', 'BGO', 'CsI'], np.array([t_decay_st[3], t_decay_st[5], t_decay_st[4] ])*1e6, 
        yerr = np.array([delta_t_decay_st[3], delta_t_decay_st[5], delta_t_decay_st[4] ])*1e6, 
        edgecolor="black")
#plt.xlabel("ADC channels", fontsize=10)                        #xlabel
plt.ylabel("Decay time (us)", fontsize=14)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=14)              #size of axis
plt.grid(True) 
plt.savefig('Rise_decay_time_raw.png', format='png')


#Print:
print('Decay time of raw LYSO: (' + str(t_decay_st[3]*1e6) + ' +/- ' + str(delta_t_decay_st[3]*1e6) + ')us')
print('Decay time of raw BGO: (' + str(t_decay_st[5]*1e6) + ' +/- ' + str(delta_t_decay_st[5]*1e6) + ')us')
print('Decay time of raw CsI: (' + str(t_decay_st[4]*1e6) + ' +/- ' + str(delta_t_decay_st[4]*1e6) + ')us' +"\n")
##############################################



####################################################################
##################DEcay time ratios#############################
##################################################################
#ratio between the decay time from the pre and from the detector

ratio_t_decay_pre_raw = t_decay_st[3:6]/t_decay_st[0:3] #correct indexes, yes
print('Ratio t_decay pre/raw: (LYSO, CsI, BGO) ' + str(ratio_t_decay_pre_raw) + '\n')

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



############################################################
##################Amplitudes ratios#############################
##################################################################
#ratio between the amplitudes from the pre and from the detector


#Amplitudes (voltage integration)
print('Amplitude (Vmax) LYSO (V) = ' + str(amplitude_st[0]) + ' +- ' + str(delta_amplitude_st[0]))
print('Amplitude (Vmax) BGO (V) = ' + str(amplitude_st[2]) + ' +- ' + str(delta_amplitude_st[2]))
print('Amplitude (Vmax) CsI (V) =' + str(amplitude_st[1]) + ' +- ' + str(delta_amplitude_st[1]) +'\n')

print('Amplitude (Vmax) raw LYSO (V) = ' + str(amplitude_st[3]) + ' +- ' + str(delta_amplitude_st[3]))
print('Amplitude (Vmax) raw BGO (V) = ' + str(amplitude_st[5]) + ' +- ' + str(delta_amplitude_st[5]))
print('Amplitude (Vmax) raw CsI (V) =' + str(amplitude_st[4]) + ' +- ' + str(delta_amplitude_st[4]) +'\n')


#Amplitudes ratios, pre/raw
    #order: LYSO,CsI, BGO, 0,1,2 pre, 3,4,5 raw

ratio_ampl_pre_raw = amplitude_st[0:3]/amplitude_st[3:6]  #yup, correct indexes
delta_ratio_ampl_pre_raw = ratio_ampl_pre_raw * np.sqrt( 
    np.array( [ (delta_amplitude_st[0]/amplitude_st[0])**2 + (delta_amplitude_st[3]/amplitude_st[3])**2, 
    (delta_amplitude_st[1]/amplitude_st[1])**2 + (delta_amplitude_st[4]/amplitude_st[4])**2, 
    (delta_amplitude_st[2]/amplitude_st[2])**2 + (delta_amplitude_st[5]/amplitude_st[5])**2   
    ] ) )                               #delta of the amplitudes, following the order defined in the variable

#print('Ratio amplitude pre/raw: (LYSO, CsI, BGO) ' + str(ratio_ampl_pre_raw) + '\n')

print('Amplitude ratio LYSO pre/raw = ' + str(ratio_ampl_pre_raw[0]) + ' +- ' + str(delta_ratio_ampl_pre_raw[0]) )
print('Amplitude ratio BGO pre/raw = ' + str(ratio_ampl_pre_raw[2]) + ' +- ' + str(delta_ratio_ampl_pre_raw[2]) )
print('Amplitude ratio CsI pre/raw = ' + str(ratio_ampl_pre_raw[1]) + ' +- ' + str(delta_ratio_ampl_pre_raw[1]) + '\n')



print('Integral LYSO (baseline removed) (V*s) = ' + str(integral_st[0]) + ' +- ' + str(delta_integral_st[0]) )
print('Integral BGO (baseline removed) (V*s) =' + str(integral_st[2]) + ' +- ' + str(delta_integral_st[2]) )
print('Integral CsI (baseline removed) (V*s) =' + str(integral_st[1]) + ' +- ' + str(delta_integral_st[1]) + '\n')


#And the ratio of the maximum voltage:
    #BGO/ CsI & LYSO/CsI & LYSO/BGO

ratio_voltage_max = np.array([amplitude_st[2] / amplitude_st[1], 
                      amplitude_st[0] / amplitude_st[1],
                      amplitude_st[0] / amplitude_st[2]
                      ])                #ratio of amplitudes, order BGO/CsI, LYSO/CsI, LYSO/BGO


delta_ratio_voltage_max= ratio_voltage_max * np.sqrt( np.array(
    [(delta_amplitude_st[2]/amplitude_st[2])**2 + (delta_amplitude_st[1]/amplitude_st[1])**2,
      (delta_amplitude_st[0]/amplitude_st[0])**2 + (delta_amplitude_st[1]/amplitude_st[1])**2,
      (delta_amplitude_st[0]/amplitude_st[0])**2 + (delta_amplitude_st[2]/amplitude_st[2])**2,
              ]
    )
    )   #computing the error. The sqrt of the quotient have to be computed by hand :))


print('Amplitude ratio BGO/CsI = ' + str(ratio_voltage_max[0]) + ' +- ' + str(delta_ratio_voltage_max[0]))
print('Amplitude ratio LYSO/CsI =' + str(ratio_voltage_max[1]) + ' +- ' + str(delta_ratio_voltage_max[1]))
print('Amplitude ratio LYSO/BGO = ' + str(ratio_voltage_max[2]) + ' +- ' + str(delta_ratio_voltage_max[2]) +'\n')


#Integral ratio

ratio_int = np.array([integral_st[2] / integral_st[1], 
                      integral_st[0] / integral_st[1],
                     integral_st[0] / integral_st[2]
                      ])

delta_ratio_int = ratio_int * np.sqrt( np.array(
    [(delta_integral_st[2]/integral_st[2])**2 + (delta_integral_st[1]/integral_st[1])**2,
     (delta_integral_st[0]/integral_st[0])**2 + (delta_integral_st[1]/integral_st[1])**2,
     (delta_integral_st[0]/integral_st[0])**2 + (delta_integral_st[2]/integral_st[2])**2,
             ] ) )   #computing the error. The sqrt of the quotient have to be computed by hand :))
        #FOUNDED ERROR HERE!!!, WAS NOT MULTIPLYING BU THE INTEGRAL RATOI BUT WITH THE TOTAL VOLTAGE RATIO
        #(NOW DELETED SINCE IT IS USELESS)
        #
        #DONT MATTER ACTUALLY, SINCE THE ERROR IS SO HIGH I AM NOT USGIN IT Xdddddd
        
        
print('Integral (no baseline) ratio BGO/CsI = ' + str(ratio_int[0]) + ' +- ' + str(delta_ratio_int[0]) )
print('Integral (no baseline) ratio  LYSO/CsI =' + str(ratio_int[1]) + ' +- ' + str(delta_ratio_int[1]) )
print('Integral (no baseline) ratio  LYSO/BGO = ' + str(ratio_int[2]) + ' +- ' + str(delta_ratio_int[2]) +'\n')








######REsiduo: Rise and decay time del osciloscopio#####################################
###################################################################################

# #This is the data from the new measurements, using the threshold to only measure
# #the gamma peak. The old results are commented

# #Variables that will store(st) the results 


# ###################################LYSO RAW########################
# #the error  is basically is 1/5 of the escale, since each square is divided
#                 #into 5 intervals, and the scale is the length of the sides of
#                 #the square
# #RAW

# t_rise, delta_t_rise = 70, 1/5*50        #[ns] rise time and its error
# t_decay, delta_t_decay = 890, 1/5*250   #[ns] decay time and its error
  
# #storing


# ############################CsI RAW#####################3#####

# t_rise, delta_t_rise = 380, 1/5*250        #[ns] rise time and its error
# t_decay, delta_t_decay = 4.32e3, 1/5*1e3   #[ns] decay time and its error
  
# #storing
# t_rise_st = np.append(t_rise_st,t_rise)
# t_decay_st = np.append(t_decay_st,t_decay)
# delta_t_rise_st = np.append(delta_t_rise_st,delta_t_rise)
# delta_t_decay_st = np.append(delta_t_decay_st,delta_t_decay)

# ####################################BGO RAW###################

# t_rise, delta_t_rise = 176, 1/5*100        #[ns] rise time and its error
# t_decay, delta_t_decay = 1.54e3, 1/5*500   #[ns] decay time and its error
  
# #storing
# t_rise_st = np.append(t_rise_st,t_rise)
# t_decay_st = np.append(t_decay_st,t_decay)
# delta_t_rise_st = np.append(delta_t_rise_st,delta_t_rise)
# delta_t_decay_st = np.append(delta_t_decay_st,delta_t_decay)


# #################################LYSO PRE##################3####

# t_rise, delta_t_rise = 64, 1/5*50        #[ns] rise time and its error
# t_decay, delta_t_decay = 690, 1/5*250   #[ns] decay time and its error
  

# #new waveform, not weird (20_5_21):
#     #t rise = 70ns, escala 100ns==> error 50/5ns = 10ns
#     #t decay = 530ns, error de 1/5 * 100ns = 20ns
    
    
# #wave, no plitter, 24_5 (the last one is not this, but 2 waves together!!!)
#     #t dec: 680ns, error 100/5 ns
#     #t rise: 94 ns, error 1/5 50ns

# #storing
# t_rise_st = np.append(t_rise_st,t_rise)
# t_decay_st = np.append(t_decay_st,t_decay)
# delta_t_rise_st = np.append(delta_t_rise_st,delta_t_rise)
# delta_t_decay_st = np.append(delta_t_decay_st,delta_t_decay)


# ###############################CsI PRE#########################

# t_rise, delta_t_rise = 360, 1/5*250        #[ns] rise time and its error
# t_decay, delta_t_decay = 3.34e3, 1/5*500   #[ns] decay time and its error
  
# #storing
# t_rise_st = np.append(t_rise_st,t_rise)
# t_decay_st = np.append(t_decay_st,t_decay)
# delta_t_rise_st = np.append(delta_t_rise_st,delta_t_rise)
# delta_t_decay_st = np.append(delta_t_decay_st,delta_t_decay)

# #new waveform, not weird (21_5_21):
#     #t rise = 310ns,  error 250/5ns
#     #t decay = 4.08us, error de 1us/5  

# #wave, no plitter, 24_5
#     #t dec: 6.52us, error 1/5 us
#     #t rise: 310nsm error 250*1/5ns

# ###############################3###BGO Pre############################

# t_rise, delta_t_rise = 184, 1/5*100        #[ns] rise time and its error
# t_decay, delta_t_decay = 1.48e3, 1/5*100   #[ns] decay time and its error
  
# #storing
# t_rise_st = np.append(t_rise_st,t_rise)
# t_decay_st = np.append(t_decay_st,t_decay)
# delta_t_rise_st = np.append(delta_t_rise_st,delta_t_rise)
# delta_t_decay_st = np.append(delta_t_decay_st,delta_t_decay)

# #new waveform, not weird (21_5_21):
#     #t rise = 165ns,  error 100/5ns 
#     #t decay = 1.52us, error de 500ns/5  
    
# #wave, no plitter, 24_5
#     #t dec: 1.62 us, error 500/5 ns
#     #t rise: 132 ns, error 1/5 100ns 

    










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
