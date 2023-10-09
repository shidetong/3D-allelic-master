import pandas as pd
import numpy as np
from __future__ import  division
from itertools import  islice
from scipy import optimize
import numpy as np
import bisect, cPickle
import pysam, sys, os  


peak = pd.read_csv('/public/home/shidetong/wedata/test/Allelic/chip_allelic/rmdup_peak/peak_signal_350bp/allpeak/alldepth_norchr.csv',sep = '\t')

B6 = peak[(peak[3]/(peak[3]+peak[4])>0.6)&(peak[6]/(peak[5]+peak[6])>0.6)]
Cast = peak[(peak[4]/(peak[3]+peak[4])>0.6)&(peak[5]/(peak[5]+peak[6])>0.6)]
B6.columns = ['chr','satrt','end','B6M_depth','B6P_depth','CastM_depth','CastP_depth']
Cast.columns = ['chr','satrt','end','B6M_depth','B6P_depth','CastM_depth','CastP_depth']
B6['B6_score'] = B6['B6M_depth']/(B6['B6M_depth']+B6['B6P_depth'])
B6['Cast_score'] = B6['CastP_depth']/(B6['CastM_depth']+B6['CastP_depth'])
Cast['Cast_score'] = Cast['CastM_depth']/(Cast['CastM_depth']+Cast['CastP_depth'])
Cast['B6_score'] = Cast['B6P_depth']/(Cast['B6M_depth']+Cast['B6P_depth'])

"""
B6M = peak[peak[3]/(peak[3]+peak[4])>0.6]
B6P = peak[peak[4]/(peak[3]+peak[4])>0.6]
CastM = peak[peak[5]/(peak[5]+peak[6])>0.6]
CastP = peak[peak[6]/(peak[5]+peak[6])>0.6]

B6 = pd.concat([B6M,B6P],axis = 0)
Cast = pd.concat([CastM,CastP],axis = 0)
B6.to_csv('/public/home/shidetong/wedata/test/Allelic/chip_allelic/rmdup_peak/peak_signal_350bp/allpeak/classify/B6.csv',sep = '\t',index = None,header = None)
Cast.to_csv('/public/home/shidetong/wedata/test/Allelic/chip_allelic/rmdup_peak/peak_signal_350bp/allpeak/classify/Cast.csv',sep = '\t',index = None,header = None)
"""
#增加一列score

B6 = pd.read_csv('/public/home/shidetong/wedata/test/Allelic/chip_allelic/rmdup_peak/peak_signal_350bp/allpeak/classify/B6.csv',sep = '\t',header = None)
Cast = pd.read_csv('/public/home/shidetong/wedata/test/Allelic/chip_allelic/rmdup_peak/peak_signal_350bp/allpeak/classify/Cast.csv',sep = '\t',header = None)
B6[7] = B6[3]/(B6[3]+B6[4])
Cast[7] = Cast[5]/(Cast[5]+Cast[6])
B6.to_csv('/public/home/shidetong/wedata/test/Allelic/chip_allelic/rmdup_peak/peak_signal_350bp/allpeak/classify/B6_score.csv',index = None,sep = '\t',header = None)
Cast.to_csv('/public/home/shidetong/wedata/test/Allelic/chip_allelic/rmdup_peak/peak_signal_350bp/allpeak/classify/Cast_score.csv',index = None,sep = '\t',header = None)
def add_score(file,outfile):
    f = pd.read_csv(file,sep = '\t',header = None)
    f[7] = f[3]/(f[3]+f[4])
    f[8] = f[5]/(f[5]+f[6])
    f.to_csv(outfile,index = None,sep = '\t',header = None)

#peak to loop
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
    dtype = np.dtype({'names':['chr','start','end','M_C','P_C','Score'],
                      'formats':['S4',np.int,np.int,np.int,np.int,np.float]})
    
    Peaks = np.loadtxt(fil,dtype = dtype,usecols = (0,1,2,5,6,-1))
    
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
        mask = ((tmp_peaks['start'] > start) & (tmp_peaks['end'] < start + 40000)) | \
               ((tmp_peaks['start'] > end) & (tmp_peaks['end'] < end + 40000))
        peak = tmp_peaks[mask]
        if peak.size == 0:
            continue
        else:
            Peak_To_Loop[tuple(lp)] = peak
    
    return Peak_To_Loop

## 可以将loop的位置拆分为两部分，分别看其位置的peak
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
        #mask1 = ((tmp_peaks['start'] > start) & (tmp_peaks['end'] < start + 40000)) 
        mask2 = ((tmp_peaks['start'] > end) & (tmp_peaks['end'] < end + 40000))
        
        #peak = tmp_peaks[mask1]
        peak = tmp_peaks[mask2]
        if peak.size == 0:
            continue
        else:
            Peak_To_Loop[tuple(lp)] = peak
    
    return Peak_To_Loop    

ptl = Integrate_Allelic_Peak_Loops(Peak_fil,Loop_fil)

#筛选每个loop位置中差异最大的峰
"""
输入的ptl  为 peak_to_loop 
"""
bias = [[],[]]
for key,values in ptl.items():
    k = key
    p = values
    w = pd.DataFrame(p)
    m = w[w['Score']==w['Score'].max()]
    q = np.array(m)
    bias[0].append(k)
    bias[1].append(q)
dic = zip(bias[0],bias[1])
dicc = dict(dic)

def peak_max(Peak_To_Loop):
    for key,values in Peak_To_Loop.items():
        k = key
        p = values
        w = pd.DataFrame(p)
        m = w[w['Score']==w['Score'].max()]
        q = np.array(m)
        bias[0].append(k)
        bias[1].append(q)
    dic = zip(bias[0],bias[1])
    dicc = dict(dic)



def get_loc1_loc2(outfile):
    """
    列出染色质环的loc的两个位置
    """
    b = Integrate_Allelic_Peak_Loops(Peak_fil,Loop_fil)
    kk = []
    for key,values in b.items():
        k =key
        kk.append(k)
    kkk = np.array(kk)
    q = pd.DataFrame(kkk)
    q[1] = q[1].astype('int')
    q[2] = q[2].astype('int')

    q[6] = q[1] + 40000
    q[7] = q[2] + 40000
    q[0] = 'chr' + q[0]

    q.columns = ['chr','ss','es','B6M_C','B6P_C','pvalue','se','ee']

    q.to_csv(outfile,index = None,sep = '\t')
    return q 


##列出染色质环的loc的两个位置。
kk = []
for key,values in ptl.items():
    k = key
    kk.append(k)

kkk = np.array(kk)
q = pd.DataFrame(kkk)

q[1] = q[1].astype('int')
q[6] = q[1] + 40000


w = np.array(q)

nn = []
for i in w[:,[1,2,3]] :
    nn.append(i)

pp = []
for i in nn:
    r = i.tolist()
    pp.append(r)

result = list([tuple(e) for e in pp])


###peak_max
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
        ####取loc_1位置或者loc_2位置
        mask1 = ((tmp_peaks['start'] > start) & (tmp_peaks['end'] < start + 40000)) 
        #mask2 = ((tmp_peaks['start'] > end) & (tmp_peaks['end'] < end + 40000))
        
        peak = tmp_peaks[mask1]
        #peak = tmp_peaks[mask2]
        if peak.size == 0:
            continue
        else:
            Peak_To_Loop[tuple(lp)] = peak 
    return Peak_To_Loop

b = Integrate_Allelic_Peak_Loops('/public/home/shidetong/wedata/test/Allelic/chip_allelic/rmdup_peak/peak_signal_350bp/allpeak/classify/allpeak.csv' ,'/public/home/shidetong/wedata/test/Hic/Allilc_loop/Allelic_input/loop/pvalue0.1/cut/pvalue0.1_sort_cut_noX.csv')

def peak_max(peak,loop,out):
    """
    获得peak_to_loop中peak位置的等位特异性差异最大的peak
    """
    dicc = Integrate_Allelic_Peak_Loops(peak,loop)  ##可以取loc_1和loc_2
    bias = [[],[]]
    peak = []
    for key,values in dicc.items():
        k = key
        p = values
        w = pd.DataFrame(p)
        m = w[w['B6Score']==w['B6Score'].max()]   #取在某个loc中peak的等位特异性最大的peak
        q = np.array(m)
        bias[0].append(k)
        bias[1].append(q)
    dic = zip(bias[0],bias[1])   #压缩为字典
    dicc = dict(dic)

    i = []
    for key,values in dicc.items():
        v = values[0]
        o = v.tolist()
        e = tuple(v)   #将list转为tuple
        i.append(e)

    ss = np.array(i)
    r = pd.DataFrame(ss)
    r[0] = 'chr'+ r[0]   #第一列添加chr
    r.columns = ['chr','start','end','B6M','B6P','CastM','CastP','B6score','Castscore']
    mm = r.sort_values(by=['chr'],ascending=True) #排序
    mm.read_csv(out,sep = '\t',index = None)



def all_peak(peak,loop,out):
    """
    将loc位置中的所有peak均取出来，但是要注意判断B6还是Cast
    """
    b = Integrate_Allelic_Peak_Loops(peak,loop)  ##可以取loc_1和loc_2
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
    ss = np.array(i)
    r = pd.DataFrame(ss)
    r[0] = 'chr'+ r[0]
    r.columns = ['chr','start','end','B6M','B6P','CastM','CastP','B6score','Castscore']
    mm = r.sort_values(by=['chr'],ascending=True)
    mm.read_csv(out,sep = '\t',index = None)



        









