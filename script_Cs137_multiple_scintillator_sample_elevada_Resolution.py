#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 10:28:13 2021

@author: dla


This contains the calcs of the B.1.5. experiment with the spectras. Spectras obtained with a similar configuration within each other, changing only the cristal.

The sample is elevated, 1 cylinder only (plus the little ring placed just above the scintillator (see pictures in logbook))
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

#%% ###############################################################
#########################0), Data loading #####################3
################################################################

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


#%% #################################################################
############################1) FIT #################################
###########################################################

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
sigma_stored = np.array([])
mean_stored = np.array([])
delta_mean_stored = np.array([])
delta_sigma_stored = np.array([])
FWHM_stored = np.array([])
delta_FWHM_stored = np.array([])


sigma_stored = np.append(sigma_stored, fit['sigma'])
mean_stored = np.append(mean_stored, fit['mean'])
delta_mean_stored = np.append(delta_mean_stored, fit['\Delta(mean)'])
delta_sigma_stored = np.append(delta_sigma_stored, fit['\Delta(sigma)'])
FWHM_stored = np.append(FWHM_stored, fit['FWHM'])
delta_FWHM_stored = np.append(delta_FWHM_stored, fit['\Delta(FWHM)'])


#$$$$$$$$$$$$$$$$$$$ BGO $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

peak = max(counts_stored[1])        #max value of the count rate, i.e., peak
peak_index = counts_stored[1].index(peak) #1143

fit = Gaussian_fit.Gaussian_fit(ADC_channel[1036-1:1273-1], 
                                   counts_stored[1][1036-1:1273-1])


sigma_stored = np.append(sigma_stored, fit['sigma'])
mean_stored = np.append(mean_stored, fit['mean'])
delta_mean_stored = np.append(delta_mean_stored, fit['\Delta(mean)'])
delta_sigma_stored = np.append(delta_sigma_stored, fit['\Delta(sigma)'])
FWHM_stored = np.append(FWHM_stored, fit['FWHM'])
delta_FWHM_stored = np.append(delta_FWHM_stored, fit['\Delta(FWHM)'])

#$$$$$$$$$$$$$$$$$$$ LYSO $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#here, again, the maximum peak is not the one of the gamma peak, so again have
#to find the peak by hand :)

fit = Gaussian_fit.Gaussian_fit(ADC_channel[2548-1:3191-1], 
                                   counts_stored[2][2548-1:3191-1])


sigma_stored = np.append(sigma_stored, fit['sigma'])
mean_stored = np.append(mean_stored, fit['mean'])
delta_mean_stored = np.append(delta_mean_stored, fit['\Delta(mean)'])
delta_sigma_stored = np.append(delta_sigma_stored, fit['\Delta(sigma)'])
FWHM_stored = np.append(FWHM_stored, fit['FWHM'])
delta_FWHM_stored = np.append(delta_FWHM_stored, fit['\Delta(FWHM)'])


#%% #########################################3############################
#################### 2) Plot of R vs Scintillator #######################
##########################################################################

#resolution = FWHM/<E> , <E> the centroid of the peak.

R_stored = FWHM_stored / mean_stored           #channel Resolution 
R_stored_100 = 100 * R_stored                                 #R[%]

#calc of delta R:
    
sqrt_sum = np.sqrt( (delta_FWHM_stored / FWHM_stored)**2 + (delta_mean_stored / mean_stored)**2 ) 
                                         #sqrt of the sum of relative errors

delta_R_stored = R_stored * sqrt_sum                                #delta(R)
delta_R_stored_100 = 100 * delta_R_stored                           #delta_R[%]



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
print('FWHM LYSO: ' + str(FWHM_stored[2]) + ' +/- ' + str(delta_FWHM_stored[2])+"\n")

print('<channels> CsI: ' + str(mean_stored[0]) + ' +/- ' + str(delta_mean_stored[0]))
print('<channels> BGO: ' + str(mean_stored[1]) + ' +/- ' + str(delta_mean_stored[1]))
print('<channels> LYSO: ' + str(mean_stored[2]) + ' +/- ' + str(delta_mean_stored[2])+"\n")

print('R CsI: ' + str(R_stored_100[0]) + ' +/- ' + str(delta_R_stored_100[0]))
print('R BGO: ' + str(R_stored_100[1]) + ' +/- ' + str(delta_R_stored_100[1]))
print('R LYSO: ' + str(R_stored_100[2]) + ' +/- ' + str(delta_R_stored_100[2]) + '\n')

#Lyso 2 BGO 1




#%%########################################################################
####################### 3)Ratio of total count rate of the spectra#############
######################################################################

#Juanpa suggest to simply compute the ratio of the total number of counts to see
#whether this could be similar to the light yield ratio or no. So, come on!
#Since the time are different, I must compute the count/rate in order
#to compare the values

#0 csI, 1 BGO, 2 LYSO. 

one_slash_time = [1/x for x in time]                #1/time [s-1] to compute the ratio

total_rate = np.multiply(total_counts, one_slash_time) #[s-1] total counts/total time

ratio_total_counts = [total_rate[2]/total_rate[0], 
                      total_rate[2]/total_rate[1], 
                      total_rate[1]/total_rate[2] ] #LYSO/CsI, LYSO/BGO, BGO/CsI

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
plt.title("Light yield ratio of the spectra of  Cs137", 
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

print('Ratio LYSO/CsI: (' + str(ratio_total_counts[0]) + ' +/- ' + str(delta_ratio_total_counts[0]) )
print('Ratio LYSO/BGO: (' + str(ratio_total_counts[1]) + ' +/- ' + str(delta_ratio_total_counts[1]) )
print('Ratio BGO/CsI: (' + str(ratio_total_counts[2]) + ' +/- ' + str(delta_ratio_total_counts[2]) + '\n' )









#%%#########################################################################
##############################3) Peka position ratio#######################
############################################################################

#This is basically the ratio of the position of the peak for each scintillator. We assume this is the p
#position in the hist!!!.

#we have the mean_stored of the gaussian fit, so it is asically using it  #0 csI, 1 BGO, 2 LYSO. 

peak_position_ratio = np.array( [ mean_stored[1] / mean_stored[0] , 
         mean_stored[2] / mean_stored[0], 
         mean_stored[2] / mean_stored[1] ])
                            #BGO/ CsI, LYSO/CsI, LYSO/BGO

aux = [np.sqrt( ( delta_mean_stored[1]/ mean_stored[1] )**2 + ( delta_mean_stored[0]/ mean_stored[0] )**2 ),
       np.sqrt( ( delta_mean_stored[2]/ mean_stored[2] )**2 + ( delta_mean_stored[0]/ mean_stored[0] )**2 ),
       np.sqrt( ( delta_mean_stored[2]/ mean_stored[2] )**2 + ( delta_mean_stored[1]/ mean_stored[1] )**2 )
       ]            #order of the peak position ratio
aux = np.array(aux)

delta_peak_position_ratio = peak_position_ratio * aux
                
                                        #since the order of the peak position ratio is not the one of
                                        #the stored mean values, have to create the arrayof sqrt by hand

print('Peak position ratio BGO/CsI: (' + str(peak_position_ratio[0]) + ' +/- ' + str(delta_peak_position_ratio[0]) )
print('Peak position ratio LYSO/CsI: (' + str(peak_position_ratio[1]) + ' +/- ' + str(delta_peak_position_ratio[1]) )
print('Peak position ratio LYSO/BGO: (' + str(peak_position_ratio[2]) + ' +/- ' + str(delta_peak_position_ratio[2]) )












#%% RESIDUO: REDO OF THE RESOLUTION CALC, BUT THIS TIME ONLY MEASURING THE PEAK.

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
                     
