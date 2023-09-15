# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 08:51:36 2021

@author: Etienne.CLAUSS
"""
import matplotlib.pyplot as plt
import numpy as np, seaborn as sns, scipy.stats as Stat
from ToolKit_PSC.IntraGrpStat import IntraGrpStat


class PlotFactory():
    def __init__(self, data, threshold, indexes, file_name, 
                 exp_duration, rec_duration, drug_arrival):
        self.data, self.threshold = data, threshold
        self.indexes = indexes
        self.file_name = file_name.split('.')[0]
        self.exp_duration = exp_duration
        self.rec_duration = rec_duration
        self.time2 = np.linspace(0, self.exp_duration, self.rec_duration)
        self.drug_arrival = drug_arrival

    def TracePloter(self, individual_traces, sampling_Hz, show_fig = False):
        
        plt.figure().suptitle(self.file_name)
        X = np.linspace(0, self.exp_duration, len(self.data))
        plt.plot(X, self.data, c='b', lw=.5, label='data')
        plt.axhline(self.threshold)
        plt.axvline(self.drug_arrival, c='gold', lw=2)
        for index in self.indexes:
            plt.scatter(X[index], self.data[index], c='r',marker='o')
        plt.legend(loc='upper right')
        plt.savefig(rf'{individual_traces}/{self.file_name}.pdf')
        if not show_fig: plt.close()
     
        
class FinalPlots():
    def __init__(self, exp_duration, rec_duration, drug, colorz, 
                 folder, drug_arrival, phenotype):
        self.exp_duration = exp_duration
        self.rec_duration = rec_duration
        self.sampling_Hz = rec_duration/exp_duration
        self.x_ax = np.linspace(0, self.exp_duration, self.rec_duration)
        self.drug = drug
        self.colorz = colorz
        self.folder = folder
        self.drug_arrival = drug_arrival
        self.phenotype = phenotype
 
    def TimeCourses(self, bin_10s, bl_spikes, ratio_max, spike_tot, wind_bl, 
                    wind_dr, wind_wsh, quantif_wind, show_fig = False):
       
        #Time course normalised to the total spikes in the recording
        plt.figure(), plt.title(f'{self.drug} {self.phenotype} norm to tot spikes')
        plt.xlabel('Bin 10s'), plt.ylabel('% activity')
        bin_norm = [[(x / spike_tot[y]) for x in bin_10s[y]] 
                    for y, _ in enumerate(bin_10s)]
        bin_norm, x_ax = np.asarray(bin_norm), np.arange(len(bin_norm[0]))
        mean_bin = np.nanmean(bin_norm, axis = 0)
        
        sem_bin = [Stat.sem(bin_norm[:,i]) for i in x_ax]
        plt.fill_between(x_ax, mean_bin-sem_bin, mean_bin+sem_bin)
        plt.axvline( self.drug_arrival/10, c='gold', lw = 2)
        
        plt.axvspan(wind_bl/10, wind_bl/10+quantif_wind/10,
            color='yellow', alpha=.1)
        plt.axvspan(wind_dr/10, wind_dr/10+quantif_wind/10,
            color='b', alpha=.1)
        plt.axvspan(wind_wsh/10, wind_wsh/10+quantif_wind/10,
            color='green', alpha=.1)
        
        plt.plot(x_ax, mean_bin, c='r', zorder=2)      
        plt.savefig(rf'{self.folder}\Time course tot.pdf') 
        if not show_fig: plt.close()

        # Individual time course normalised to the total spikes
        for i,x in enumerate(bin_norm):
            plt.figure(), plt.title(f'{self.drug} {self.phenotype} cell {i+1}')
            plt.xlabel('Bin 10s'), plt.ylabel('% activity')
            plt.axvline( self.drug_arrival/10, c='gold', lw = 2)
            
            plt.axvspan(wind_bl/10, wind_bl/10+quantif_wind/10,
                color='yellow', alpha=.1)
            plt.axvspan(wind_dr/10, wind_dr/10+quantif_wind/10,
                color='b', alpha=.1)
            plt.axvspan(wind_wsh/10, wind_wsh/10+quantif_wind/10,
                color='green', alpha=.1)
            plt.plot(x, c = 'black')
            plt.savefig(rf'{self.folder}\Time course spike tot Freq cell {i+1}.pdf') 
            if not show_fig: plt.close()
           
    
    #Histos
    def Histo(self, bin_1s, wind_bl, wind_dr, wind_wsh, 
              quantif_wind, Paired, show_fig = False):        
        x_ax = (0, 0.1, 0.2)     
        bl = [bin_1s[i][wind_bl:wind_bl+quantif_wind] for i,_ in enumerate(bin_1s)]        
        dr = [bin_1s[i][wind_dr:wind_dr+ quantif_wind] 
                  for i,_ in enumerate(bin_1s)]        
        dr2 = [bin_1s[i][wind_dr + quantif_wind:wind_dr + 2*quantif_wind] 
                  for i,_ in enumerate(bin_1s)]        
        wsh = [bin_1s[i][wind_wsh:wind_wsh + quantif_wind] 
                  for i,_ in enumerate(bin_1s)]        

        bl = [np.mean(bl[i]) for i,_ in enumerate(bl)]
        dr = [np.mean(dr[i]) for i,_ in enumerate(dr)]
        dr2 = [np.mean(dr2[i]) for i,_ in enumerate(dr2)]
        wsh = [np.mean(wsh[i]) for i,_ in enumerate(wsh)]
        # val_stat = np.asarray((bl, dr, dr2, wsh)) 
        val_stat = np.asarray((bl, dr, wsh)) 
        indiv_val = val_stat.transpose()
        # meanz = np.asarray((np.mean(bl), np.mean(dr), np.mean(dr2), np.mean(wsh)))
        meanz = np.asarray((np.mean(bl), np.mean(dr), np.mean(wsh)))
        
        # bl_sem, dr_sem, dr2_sem, wsh_sem = (Stat.sem(bl), Stat.sem(dr), 
        #                                     Stat.sem(dr2),Stat.sem(wsh)) 
        bl_sem, dr_sem, wsh_sem = (Stat.sem(bl), Stat.sem(dr),Stat.sem(wsh)) 
        # semz = np.asarray((np.mean(bl_sem), np.mean(dr_sem), 
        #                     np.mean(dr2_sem), np.mean(wsh_sem)))
        semz = np.asarray((np.mean(bl_sem), np.mean(dr_sem), np.mean(wsh_sem)))

        
        plt.figure(), plt.xticks(x_ax, ['Baseline', self.drug, 'Wash']), 
        plt.title(f'{self.drug} {self.phenotype}')
        plt.ylabel('Frequency (Hz)') 
        plt.bar(x_ax, meanz, width = 0.1,
                yerr=semz, color=self.colorz[1], capsize=10, zorder=0)
        plt.plot(1,0)
        for x in indiv_val:
            [plt.scatter(x_ax[i], x, s=20, c='w', marker='o',
                            edgecolor='k', zorder=2) for i,x in enumerate(x)]  
            plt.plot(x_ax, x, c='k', lw=0.5, zorder=1)     
    
        stat, pval, ph_out = IntraGrpStat(val_stat, Paired)
        
        if pval > 0.05:
            plt.text(x_ax[-1]+0.1, meanz[2], f'{pval} \n {ph_out["Test"]}')
        else:
            plt.text(x_ax[-1]+0.1, meanz[2], f'{pval} \n {ph_out["Test"]} \n {ph_out["Time 0 vs Time 1"]} \n {ph_out["Time 0 vs Time 2"]} \n {ph_out["Time 1 vs Time 2"]}')
                          
        plt.savefig(rf'{self.folder}\Histo_wind_fixed_3.pdf')
        if not show_fig: plt.close()
        
        return indiv_val, val_stat
        