# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 11:55:24 2021

@author: Etienne.CLAUSS
"""
import numpy as np, pyabf, pyabf.filter
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
import os, pandas as pd


class FileManager():
    def __init__(self, folder_exp):
        self.folder_exp = folder_exp


    def FileMaker(self, wind_length, drug):
        """ Create 3 files to stock the analysis """
        analysis = f'{self.folder_exp}\Analysis_{wind_length}s'
        traces = f'{analysis}\Indiv traces'
        stock = f'{analysis}\Stock'
        [os.makedirs(path) for path in [analysis, traces, stock] 
         if not os.path.exists(path)]
        writer = pd.ExcelWriter('{}/{}.xlsx'.format(analysis, drug))
        return analysis, traces, stock, writer
    
    def FileSeeker(self):
        folder_exp_list = os.listdir(self.folder_exp)
        folder_exp_list = [x for x in folder_exp_list 
                           if x.split('.')[-1] == 'abf']
        return folder_exp_list, len(folder_exp_list)
    
    
class Trace():

    def __init__(self, file):
        self.file = file
    
    def Filter(self, *args):
        file, pol, drug_arrival, exp_duration, event_type = args

        abf = pyabf.ABF(self.file)

        #Extract info
        reduction_factor = 10 if abf.dataRate == 20000 else 20
        sampling_Hz = int(abf.dataRate/reduction_factor)
        X_label, Y_label = abf.sweepLabelX, abf.sweepLabelY
        self.rec_duration = exp_duration*sampling_Hz    
        
        #Downsample the trace to reduce handling cost
        time, data  = (abf.sweepX, pol*abf.sweepY)
        
        bl_est = np.asarray(savgol_filter(data, 1001, 1))
        
        data = data - bl_est
        
        self.dat = [x for i, x in enumerate(data) if not i%reduction_factor]
        self.data = [x for x in self.dat if str(x)!='nan'][:self.rec_duration]
        self.tim = [x for i, x in enumerate(time) if not i%reduction_factor]
        self.time = [x for x in self.tim if str(x)!='nan'][:self.rec_duration]

        
        return (sampling_Hz, X_label, Y_label, self.time, 
                self.rec_duration, self.data)

    
