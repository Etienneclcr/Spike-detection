# -*- coding: utf-8 -*-
"""
Created on Fri May 21 13:39:34 2021

@author: Etienne.CLAUSS
"""

import pandas as pd, winsound
from ToolKit_AP import TraceManager_AP, Analysor_AP, PlotFactory_AP

folder_exp = r"C:\Etienne.CLAUSS\Lab\Exp\ephy\I_C\A5\LA"

'''                   Define parameters of the analysis                   '''        
#Define parameters of the analysis & create empty lists
drug, event_type, phenotype = (folder_exp.split('\\')[-1], 
                               folder_exp.split('\\')[-2], 
                               folder_exp.split('\\')[-3])
drug_arrival, exp_duration, pol, bin_10s, bin_1s = 120, 600, -1, 10, 1
chosen_wind, wind_bl, wind_dr, wind_wsh, quantif_wind = True, 60, 120, 480, 60
Paired = True

colorz = ('blue', 'limegreen', 'red', 'darkorange', 'lawngreen')

#Create folders to stock the analysis
FM = TraceManager_AP.FileManager(folder_exp)
(folder_exp_analysis, individual_traces, stock , writer) = FM.FileMaker(
                                                            quantif_wind, drug)

#Find the files containing the data
folder_exp_list, n_cell = FM.FileSeeker()

#Create the Output time course dictionnary
Stock = {x:[] for x in ('Cell ID', 'Binary', 'Spike Total', 'Bin 10s', 
                        'Bin 1s', 'bl_spikes', 'Ratio 10s', 
                        'Ratio 1s', 'Indexes', 'Ecdf')}

#Create a loop to analyses the cells one by one
for file_name in folder_exp_list:
    file = f'{folder_exp}\{file_name}'
    print('='*30, '\n'*2, file_name, '\n'*2, '='*30)
    
    '''                          Pre-Analysis part                     '''
    Trace = TraceManager_AP.Trace(file)
    
    #Apply a bandpass filter and downsample the data
    T = Trace.Filter(file, pol, drug_arrival, 
                     exp_duration, event_type)
    
    (sampling_Hz, X_label, Y_label, time, rec_duration, data) = T
     
    '''                             Analysis                               '''
    A = Analysor_AP.Analysor(rec_duration, sampling_Hz, drug_arrival)
    
    #Find Events
    threshold, indexes = A.PeakFinder(data)
    
    #Create a binary list of event and its cumulative function
    peak_binary, spike_tot, ecdf, bl_spikes = A.Binary()

    #Calculate the frequency in each bin of time
    (bin_cell_10s, ratio_max_10s, 
     bin_cell_1s, ratio_max_1s)  = A.BinCalculator(bin_10s, bin_1s)
    
    '''                       Plot section                          '''
    Plot = PlotFactory_AP.PlotFactory(data, threshold, indexes, file_name, 
                                      exp_duration, rec_duration, drug_arrival)
    
    #Create an overlay of the curent detected to check their kinetics
    Plot.TracePloter(individual_traces, sampling_Hz)
   
    '''                           Stock Section                           '''            
    for k, v in zip(Stock, (file_name, peak_binary, spike_tot, bin_cell_10s, 
                            bin_cell_1s, bl_spikes, ratio_max_10s, ratio_max_1s,
                            indexes/sampling_Hz, ecdf)):
        Stock[k].append(v)
    
    All_Mighty = pd.DataFrame.from_dict(Stock)
    All_Mighty = All_Mighty.set_index('Cell ID')
      
'''                Create the final Plots                '''
FP = PlotFactory_AP.FinalPlots(exp_duration, rec_duration, drug, colorz, 
                               folder_exp_analysis, drug_arrival, phenotype)
    
# Show the time course of the frequency variations
FP.TimeCourses(Stock['Bin 10s'], Stock['bl_spikes'], Stock['Ratio 10s'], Stock['Spike Total'], 
               wind_bl/10, wind_dr/10, wind_wsh/10, quantif_wind/10)

# Create histograms
indiv_val, val_stat = FP.Histo(Stock['Bin 1s'], wind_bl, wind_dr, wind_wsh, quantif_wind, Paired)

df = pd.DataFrame(val_stat)
writer_exp = pd.ExcelWriter(f'{folder_exp_analysis}/excel.xlsx')
df.to_excel(writer_exp)
writer_exp.save()


winsound.Beep(200, 1000)