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
import scipy
import numpy as np
    #np contain linspaces as np.linspace(a,b,N)
from scipy.stats import norm  #to do gaussian fits

        
######3

plt.close("all")

#%%
#########################1), Data loading #####################3
#The files to load are in txt. The best way to read is:

#All with same histogram, gain, bias, threshold, but different gate to optimize 
#each signal

#Firstly create the data by setting the led driver amplitude to 0, obtaining 
#extrange results, so will now unplug it, and closing the SiPM to see if now 
#the results are fine.

total_counts = []       #variable that will contain the total counts of all the
                        #spectras
counts_stored = []           #storing variable of the counts
time = 180                              #[s] duration time od the measurements

#All have the same channels, so only one variable will be saved

with open('Cs137_CsI_new_histo.txt') as file_object:
            
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
        #count_rate_CsI = [x/time_CsI for x in counts_CsI] #count rate
        
        counts_stored.append(counts_help)

with open('Cs137_BGO_new_histo.txt') as file_object:
            
        lines = file_object.readlines()
        print('the number of lines of the Cs137_BGO is',len(lines))
        #This contains strings (have to be converted to numbers using int()
        #and \n, so the \n (salto de linea) have to be removed
        
        counts_help = []
        for i in range(len(lines)):            
            counts_help.append(float(lines[i].split()[1]))    #store 2nd number of the
                            #column
       
        total_counts.append(sum(counts_help))  #total counts of the spectra         
        counts_stored.append(counts_help)           #stored of counts

with open('Cs137_LYSO_new_histo.txt') as file_object:
            
        lines = file_object.readlines()
        print('the number of lines of the Cs137_LYSO is',len(lines))
        #This contains strings (have to be converted to numbers using int()
        #and \n, so the \n (salto de linea) have to be removed
        
        counts_help = []
        for i in range(len(lines)):            
            counts_help.append(float(lines[i].split()[1]))    #store 2nd number of the
                            #column
       
        total_counts.append(sum(counts_help))  #total counts of the spectra       
        counts_stored.append(counts_help)
        
#First column of counts stored is CsI, second GBo, 3rd LYSO

#%%    0.1. Representacion
            
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
#plt.ylim(0,)                            #limits of y axis



###Plot combined, as CAEN's

plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
plt.plot(ADC_channel, counts_stored[0], color='black', label = 'CsI')    
plt.plot(ADC_channel, counts_stored[1], color = 'red', label = 'BGO')
plt.plot(ADC_channel, counts_stored[2], color = 'blue', label = 'LYSO')      
plt.legend(['CsI', 'BGO', 'LYSO'], fontsize=10) 
plt.title("Cs137 for several scintillators (similar config)", fontsize=20)           #title
plt.xlabel("ADC channels", fontsize=10)                        #xlabel
plt.ylabel("Counts", fontsize=10)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=10)              #size of axis
plt.grid(True) 
plt.savefig('count_rate_Cs137_vs_scintillator_type.png', format='png')

#plt.xlim(0,800)                       #limits of x axis
#plt.ylim(0,11000)                            #limits of y axis


#%% ##########2) FIT #################################

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

import Gaussian_fit

ajuste = Gaussian_fit.Gaussian_fit(ADC_channel[5551-1:6741-1], 
                                   counts_stored[0][5551-1:6741-1])


#Creation of the array needed to do the fit
#indexes found by hand. Look af peakUtils, may be useful to do this automatically,
#saving time!

x_data = np.array(ADC_channel[5551-1:6741-1])  #-1 because python starts at 0
#and the indexes were found at the .txt reader, that starts in line 1      
y_data = np.array(counts_stored[0][5551-1:6741-1])


initial = [max(y_data), x_data[0], (x_data[1] - x_data[0]) * 5]
                #initial guesses for the fit. If None, this does not work, so this
                #is very important when having an offset! Thank you 
                #Lucas Hermann Negri (PeakUtils)
gaussian_fit = scipy.optimize.curve_fit(gaussian, x_data, y_data, initial)

opt_values = gaussian_fit[0]   #optimal values of the function to fit the data
cov_of_opt_val = gaussian_fit[1]            #covariances of the optimal values
    #the diagonal are the variance of the parameter to estimate.
    
a = opt_values[0]  
b = opt_values[1]
cc = opt_values[2]
        #similar values as the ones given by the fit function in Matlab :)

perr = np.sqrt(np.diag(cov_of_opt_val))        #standard deviation error (el 
                                                #error de toa la via vamos)
    #source: 
    #https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html


sigma = cc/np.sqrt(2)                   #standard deviation of the gaussian fit
Delta_sigma = perr[2]/np.sqrt(2)        #error of the standar deviation
print('sigma: ' + str(sigma) + ' +/- ' + str(Delta_sigma) + ' MeV')

FWHM = 2 * np.sqrt(2 * np.exp(2)) * sigma                   #FWHM of the peak
Delta_FWHM = 2 * np.sqrt(2 * np.exp(2)) * Delta_sigma     #error of the FWHM
print('FWHM: ' + str(FWHM) + ' +/- ' + str(Delta_FWHM) + ' MeV')

##Plot of the fit
plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
plt.plot(x_data, y_data, label = 'data')        #original data
plt.plot(x_data, gaussian(x_data, a, b, cc), 'ro', label = 'fit')
plt.title('Gaussian fit to the peak, CsI', fontsize=20)          #title
plt.xlabel("E (MeV)", fontsize=10)                                    #xlabel
plt.ylabel("Cuentas", fontsize=10)                                    #ylabel
plt.tick_params(axis='both', labelsize=10)            #size of tick labels  
plt.grid(True)                                              #show grid
#plt.xlim(5.35,5.55)                                         #limits of x axis
##good enough, so moving on xD

#Storing of the relevant data, sigma and its error
sigma_stored = []
delta_sigma_stored = []
FWHM_stored = []
delta_FWHM_stored = []

sigma_stored.append(sigma)
delta_sigma_stored.append(Delta_sigma)
FWHM_stored.append(FWHM)
delta_FWHM_stored.append(Delta_FWHM)


#$$$$$$$$$$$$$$$$$$$ BGO $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#Creation of the array needed to do the fit
x_data = np.array(ADC_channel[972-1:1336-1])  #-1 because python starts at 0
#and the indexes were found at the .txt reader, that starts in line 1      
y_data = np.array(counts_stored[1][972-1:1336-1])


initial = [max(y_data), x_data[0], (x_data[1] - x_data[0]) * 5]
                #initial guesses for the fit. If None, this does not work, so this
                #is very important when having an offset! Thank you 
                #Lucas Hermann Negri (PeakUtils)
gaussian_fit = scipy.optimize.curve_fit(gaussian, x_data, y_data, initial)

opt_values = gaussian_fit[0]   #optimal values of the function to fit the data
cov_of_opt_val = gaussian_fit[1]            #covariances of the optimal values
    #the diagonal are the variance of the parameter to estimate.
    
a = opt_values[0]  
b = opt_values[1]
cc = opt_values[2]
        #similar values as the ones given by the fit function in Matlab :)

perr = np.sqrt(np.diag(cov_of_opt_val))        #standard deviation error (el 
                                                #error de toa la via vamos)
sigma = cc/np.sqrt(2)                   #standard deviation of the gaussian fit
Delta_sigma = perr[2]/np.sqrt(2)        #error of the standar deviation
print('sigma: ' + str(sigma) + ' +/- ' + str(Delta_sigma) + ' MeV')

FWHM = 2 * np.sqrt(2 * np.exp(2)) * sigma                   #FWHM of the peak
Delta_FWHM = 2 * np.sqrt(2 * np.exp(2)) * Delta_sigma     #error of the FWHM
print('FWHM: ' + str(FWHM) + ' +/- ' + str(Delta_FWHM) + ' MeV')

##Plot of the fit
plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
plt.plot(x_data, y_data, label = 'data')        #original data
plt.plot(x_data, gaussian(x_data, a, b, cc), 'ro', label = 'fit')
plt.title('Gaussian fit to the peak, BGO', fontsize=20)          #title
plt.xlabel("E (MeV)", fontsize=10)                                    #xlabel
plt.ylabel("Cuentas", fontsize=10)                                    #ylabel
plt.tick_params(axis='both', labelsize=10)            #size of tick labels  
plt.grid(True)                                              #show grid
#plt.xlim(5.35,5.55)                                         #limits of x axis
##good enough, so moving on xD

#Storing of the relevant data, sigma and its error

sigma_stored.append(sigma)
delta_sigma_stored.append(Delta_sigma)
FWHM_stored.append(FWHM)
delta_FWHM_stored.append(Delta_FWHM)

#$$$$$$$$$$$$$$$$$$$ LYSO $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#Creation of the array needed to do the fit
x_data = np.array(ADC_channel[2572-1:3362-1])  #-1 because python starts at 0
#and the indexes were found at the .txt reader, that starts in line 1      
y_data = np.array(counts_stored[2][2572-1:3362-1])


initial = [max(y_data), x_data[0], (x_data[1] - x_data[0]) * 5]
                #initial guesses for the fit. If None, this does not work, so this
                #is very important when having an offset! Thank you 
                #Lucas Hermann Negri (PeakUtils)
gaussian_fit = scipy.optimize.curve_fit(gaussian, x_data, y_data, initial)

opt_values = gaussian_fit[0]   #optimal values of the function to fit the data
cov_of_opt_val = gaussian_fit[1]            #covariances of the optimal values
    #the diagonal are the variance of the parameter to estimate.
    
a = opt_values[0]  
b = opt_values[1]
cc = opt_values[2]
        #similar values as the ones given by the fit function in Matlab :)

perr = np.sqrt(np.diag(cov_of_opt_val))        #standard deviation error (el 
                                                #error de toa la via vamos)
sigma = cc/np.sqrt(2)                   #standard deviation of the gaussian fit
Delta_sigma = perr[2]/np.sqrt(2)        #error of the standar deviation
print('sigma: ' + str(sigma) + ' +/- ' + str(Delta_sigma) + ' MeV')

FWHM = 2 * np.sqrt(2 * np.exp(2)) * sigma                   #FWHM of the peak
Delta_FWHM = 2 * np.sqrt(2 * np.exp(2)) * Delta_sigma     #error of the FWHM
print('FWHM: ' + str(FWHM) + ' +/- ' + str(Delta_FWHM) + ' MeV')

##Plot of the fit
plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
plt.plot(x_data, y_data, label = 'data')        #original data
plt.plot(x_data, gaussian(x_data, a, b, cc), 'ro', label = 'fit')
plt.title('Gaussian fit to the peak, LYSO', fontsize=20)          #title
plt.xlabel("E (MeV)", fontsize=10)                                    #xlabel
plt.ylabel("Cuentas", fontsize=10)                                    #ylabel
plt.tick_params(axis='both', labelsize=10)            #size of tick labels  
plt.grid(True)                                              #show grid
#plt.xlim(5.35,5.55)                                         #limits of x axis
##good enough, so moving on xD

#Storing of the relevant data, sigma and its error

sigma_stored.append(sigma)
delta_sigma_stored.append(Delta_sigma)
FWHM_stored.append(FWHM)
delta_FWHM_stored.append(Delta_FWHM)






#%% ############### 4) Plot of FWHM vs Scintillator ######################

plt.figure(figsize=(8,5))  #width, heigh 6.4*4.8 inches by default
plt.bar(['CsI', 'BGO', 'LYSO'], FWHM_stored, yerr = delta_FWHM_stored)
plt.title("FWHM of the Cs137 peak for several scintillators (similar config)", fontsize=20, wrap=True)           #title
#plt.xlabel("ADC channels", fontsize=10)                        #xlabel
plt.ylabel("FWHM", fontsize=10)              #ylabel
# Set size of tick labels.
plt.tick_params(axis='both', labelsize=10)              #size of axis
plt.grid(True) 
plt.savefig('FWHM_vs_scintillator.png', format='png')
