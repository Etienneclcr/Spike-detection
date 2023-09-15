# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 11:54:34 2021

@author: Angel.BAUDON
"""
from scipy.signal import find_peaks
import numpy as np


class Analysor():
    def __init__(self, rec_duration, sampling_Hz, drug_arrival):
        self.rec_duration = rec_duration
        self.sampling_Hz = sampling_Hz
        self.exp_duration = int(rec_duration/sampling_Hz)
        self.drug_arrival = drug_arrival    
    
    def PeakFinder(self, data, show_fig=False):
        threshold = 30
        
        ind, _ = find_peaks(data, height = 10, prominence = 6*np.std(data))
        self.indexes = ind
        return threshold, self.indexes

    def Binary(self):
        #Binary list of events
        self.peak_binary = np.zeros(self.rec_duration)
        for index in self.indexes: self.peak_binary[index] += 1
        spike_tot = np.sum(self.peak_binary)
        if spike_tot == 0: spike_tot = 1 
        ecdf = np.cumsum(self.peak_binary)
        ecdf = ecdf/ecdf[-1]
        ecdf = [x for i,x in enumerate(ecdf) if not i%100]
        
        splitted_frame = np.split(self.peak_binary, self.exp_duration)
        mini_bin = [np.nansum(x) for x in splitted_frame]
        bl_spikes = np.sum(mini_bin[:self.drug_arrival]) 

        return self.peak_binary, spike_tot, ecdf, bl_spikes

        
    def BinCalculator(self, bin_10s, bin_1s):
        #Split the data into bins & attribute a binary value to each bin
        #10s bins
        n_bin_10s = int(self.rec_duration/(bin_10s*self.sampling_Hz))
        splitted_cell_binary_10s = np.split(self.peak_binary, n_bin_10s)
        self.bin_cell_10s = [np.nansum(scb) for scb in splitted_cell_binary_10s]  
        ratio_10s = max(self.bin_cell_10s)   
        if ratio_10s == 0: ratio_10s = 1
    
        #1s bins
        n_bin_1s = int(self.rec_duration/(bin_1s*self.sampling_Hz))
        splitted_cell_binary_1s = np.split(self.peak_binary, n_bin_1s)
        self.bin_cell_1s = [np.nansum(scb) for scb in splitted_cell_binary_1s]  
        ratio_1s = max(self.bin_cell_1s)
        if ratio_1s == 0: ratio_1s = 1

        return self.bin_cell_10s, ratio_10s, self.bin_cell_1s, ratio_1s
