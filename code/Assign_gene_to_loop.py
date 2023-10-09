from __future__ import  division
from itertools import  islice
from scipy import optimize
import numpy as np
import pandas as pd 
# import bisect, cPickle
import pysam, sys, os

def Loading_Loops(loop_fil):
    """
    """
    Loop_type = np.dtype({'names':['chr','start','end','M_IF','P_IF','P_value'],
                          'formats':['S4',np.int,np.int,np.float,np.float,np.float]})
    
    Loops = []
    with open(loop_fil,'r') as f:
        for line in islice(f,1,None):
            line = line.strip().split()
            if line[-1] == 'NA':
                continue
            else:
                tmp = (line[0],line[1],line[2],line[3],line[4],line[5])
                Loops.append(tmp)
    
    Loops = np.array(Loops,dtype = Loop_type)
    
    return Loops

def Loading_max_Loops(loop_fil):
    """
    """
    Loop_type = np.dtype({'names':['chr','start','end','M_IF','P_IF','P_value','ss','ee'],
                          'formats':['S4',np.int,np.int,np.float,np.float,np.float,np.int,np.int]})
    
    Loops = []
    with open(loop_fil,'r') as f:
        for line in islice(f,1,None):
            line = line.strip().split()
            if line[-1] == 'NA':
                continue
            else:
                tmp = (line[0],line[1],line[2],line[3],line[4],line[5],line[11],line[12])
                Loops.append(tmp)
    
    Loops = np.array(Loops,dtype = Loop_type)
    
    return Loops


def Loading_Allel_gene(fil):
    """
    """
    dtype = np.dtype({'names':['seqnames','start','end','type','gene_name'],
                      'formats':['S4',np.int,np.int,'S4','S8']})
    
    Peaks = np.loadtxt(fil,dtype = dtype,usecols = (0,1,2,3,4),skiprows = 1)
    
    return Peaks

def Integrate_Allelic_Peak_Loops(Peak_fil,Loop_fil):
    """
    """
    Peaks = Loading_Allel_gene(Peak_fil)
    #Loops = Loading_Loops(Loop_fil)
    Loops = Loading_max_Loops(Loop_fil)
    Loops = Loops[Loops['P_value'] < 0.1]
    
    Peak_To_Loop = {}
    for lp in Loops:
        chro = lp['chr']
        start = lp['start'] - 500000
        end = lp['end'] + 500000
        # start = lp['ss']
        # end = lp['ee']
        tmp_peaks = Peaks[Peaks['seqnames'] == chro]
        mask = ((tmp_peaks['start'] > start) & (tmp_peaks['start'] < end)) | ((tmp_peaks['end'] > start) & (tmp_peaks['end'] < end))
        peak = tmp_peaks[mask]
        if peak.size == 0:
            continue
        else:
            Peak_To_Loop[tuple(lp)] = peak
    
    return Peak_To_Loop


Integrate_Allelic_Peak_Loops('/public/home/shidetong/projects/lxx/RNA-seq/haplotype/dif/padj0.1/gene_sit.csv','/public/home/shidetong/wedata/test/Hic/Allilc_loop/Allelic_input/loop/pvalue0.1/cut/pvalue0.1_sort_cut_noX.csv')