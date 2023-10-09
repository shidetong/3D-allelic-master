from __future__ import  division
from itertools import  islice
from scipy import optimize
import numpy as np
import pandas as pd 
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


def Loading_Allel_Peak(fil):
    """
    """
    dtype = np.dtype({'names':['chr','start','end','B6M_C','B6P_C','CastM_C','CastP_C'],
                      'formats':['S4',np.int,np.int,np.int,np.int,np.int,np.int]})
    
    Peaks = np.loadtxt(fil,dtype = dtype,usecols = (0,1,2,3,4,5,6))
    
    return Peaks


def Integrate_Allelic_Peak_Loops(Peak_fil,Loop_fil):
    """
    """
    Peaks = Loading_Allel_Peak(Peak_fil)
    Loops = Loading_Loops(Loop_fil)
    Loops = Loops[Loops['P_value'] < 0.1]
    
    Peak_To_Loop = {}
    for lp in Loops:
        chro = lp['chr']
        start = lp['start']
        end = lp['end']
        tmp_peaks = Peaks[Peaks['chr'] == chro]
        mask1 = ((tmp_peaks['start'] > start) & (tmp_peaks['end'] < start + 40000)) 
        #mask2 = ((tmp_peaks['start'] > end) & (tmp_peaks['end'] < end + 40000))
        
        peak = tmp_peaks[mask1]
        #peak = tmp_peaks[mask2]
        if peak.size == 0:
            continue
        else:
            Peak_To_Loop[tuple(lp)] = peak
    
    return Peak_To_Loop


b = Integrate_Allelic_Peak_Loops('/public/home/shidetong/wedata/test/Allelic/chip_allelic/rmdup_peak/peak_signal_350bp/allpeak/classify/B6.csv' ,'/public/home/shidetong/wedata/test/Hic/Allilc_loop/Allelic_input/loop/pvalue0.1/cut/pvalue0.1_sort_cut_B6.csv')
b = Integrate_Allelic_Peak_Loops('/public/home/shidetong/wedata/test/Allelic/chip_allelic/rmdup_peak/peak_150bp/depth/B6bias.csv' ,'/public/home/shidetong/wedata/test/Hic/Allilc_loop/Allelic_input/loop/pvalue0.1/cut/pvalue0.1_sort_cut_B6.csv')

bias = [[],[]]
peak = []
for key,values in b.items():
    k = key
    p = values
    w = pd.DataFrame(p)

    #m = w[w['B6Score']==w['B6Score'].max()]
    q = np.array(w)
    bias[0].append(k)
    bias[1].append(q)
dic = zip(bias[0],bias[1])
dicc = dict(dic)

i = []
for key,values in dicc.items():
    for v in values:
        q = v
        q = tuple(q)
        i.append(q)
    
print(v)    

ss = np.array(i)

r = pd.DataFrame(ss)

r[0] = 'chr'+ r[0]
r.columns = ['chr','start','end','B6M','B6P','CastM','CastP','B6score','Castscore']
mm = r.sort_values(by=['chr'],ascending=True)

mm.to_csv('/public/home/shidetong/wedata/meme/bins_all_350bp/Cast/peak/loc_2.csv',sep = '\t',index = None)

