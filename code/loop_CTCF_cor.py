#%%
from __future__ import division
import pandas as pd 
import numpy as np
#from tadlib.calfea.analyze import getmatrix
import matplotlib
# Use a non-interactive backend
# matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import os
from scipy.special import ndtr
from scipy import optimize

#--------------------------------------------------------------------------
## Matplotlib Settings
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])
my_cmap.set_bad('#2672a1')

#%%
def f_1(x, A, B):
    return A * x + B


def curve_fitting(x , y):
    A1, B1 = optimize.curve_fit(f_1, x, y)[0]
    x1 = np.arange(min(x), max(x), 0.001)
    y1 = A1 * x1 + B1
    return x1, y1



def TF_Loop_Allelic_corr(Loop , signal):
    '''
    Parameters
    ----------
    Loop : Allelic_Loop_array
        DESCRIPTION.
    signal : Allelic_TF_array
        DESCRIPTION.

    Returns
    -------
    bias : list 
        list of Loop allelic_bias and TF Signal allelic_bias.

    '''
    bias = [[],[]]
    n=0
    for i in Loop:
        g = i['chr']
        loop_ss = i['start'] 
        loop_se = i['start'] +  res 
        loop_es = i['end'] 
        loop_ee = i['end'] +  res
        loop_M = i['M_C']
        loop_P = i['P_C']
        tmp = signal[signal['chr'] == g]
        mask = (tmp['start'] <= loop_se) & (tmp['end'] >= loop_ss) | (tmp['start'] <= loop_ee) & (tmp['end'] >= loop_es) 
        overlap = tmp[mask]

        
        if overlap.size != 0:
            
            n += 1
            loop_bias = loop_M / (loop_M + loop_P)
            CTCF_b =[]
            for j in overlap:
                b = j['M_C'] /(j['M_C'] + j['P_C'])
                CTCF_b.append(b)

            if loop_M > loop_P:
                CTCF_bias = max(CTCF_b)
            else:
                CTCF_bias = min(CTCF_b)
            # CTCF_bias = (overlap1['M_C'].sum() + overlap2['M_C'].sum()) / (overlap1['M_C'].sum() + overlap1['P_C'].sum() + overlap2['M_C'].sum() + overlap2['P_C'].sum())
             
            bias[0].append(loop_bias)
            bias[1].append(CTCF_bias)
            if (CTCF_bias > 0.6) and (CTCF_bias < 0.4):
                print (i , overlap1 , overlap2)
    print(n)
                
    return bias 


def Bias_scatter_Plot(bias , signal_name):
    cor = np.round(np.corrcoef(bias[0] , bias[1])[0][1] , 5)
    x1 , y1 = curve_fitting(bias[0], bias[1])       
    fig = plt.figure(figsize = (11, 10))
    ax = fig.add_axes([0.15 , 0.15 , 0.7, 0.7])
    ax.scatter(bias[0] , bias[1] , c = 'red')
    ax.text(0.7 , 1 , 'PCC = ' + str(cor) , size = 20 )
    ax.set_xlim((-0.1 , 1.1))
    ax.set_ylim((-0.1 , 1.1))
    ax.set_xlabel('Loop Maternal/(Maternal + Paternal)_num:' + str(len(bias[0]))+'Loop',size = 20)
    ax.set_ylabel(signal_name + ' Maternal/(Maternal + Paternal)',size = 20)
    ax.set_title('Parent RNA correlation of loop',size = 30)
    
    ax.plot(x1, y1, "blue")
    
    return fig

#%%
##Loop Data----------------------------------------------------------------

res = 40000


loop_type = np.dtype({'names':['chr' , 'start' , 'end' , 'M_C' , 'P_C' , 'p_value'],
                      'formats':['U4' , np.int , np.int , np.float , np.float , 'U64']})
loop_type_1 = np.dtype({'names':['chr' , 'start' , 'end' , 'M_C' , 'P_C' , 'p_value'],
                        'formats':['U4' , np.int , np.int , np.float , np.float , np.float]})                      
                      
GM_Loop = np.loadtxt('/public/home/shidetong/wedata/test/Hic/Allilc_loop/Allelic_input/loop/pvalue0.1/cut/pvalue0.1_sort_cut_noX.csv',
                  dtype = loop_type , usecols = (0,1,2,5,6,-1) , skiprows = 1)
                         
GM_Loop = np.array(list(GM_Loop[GM_Loop['p_value'] != 'NA']) , dtype = loop_type_1)

GM_Loop = GM_Loop[GM_Loop['p_value'] <= 0.1]



GM_Loop = GM_Loop[GM_Loop['end'] - GM_Loop['start'] >= 10 * res]

Loop_M = GM_Loop[GM_Loop['M_C'] >1.2* GM_Loop['P_C']]

Loop_P = GM_Loop[GM_Loop['M_C'] *1.2  <  GM_Loop['P_C']]


GM_Loop = np.hstack((Loop_M , Loop_P))

'''
loop = pd.read_csv('/public/home/shidetong/projects/sdt/hic-seq/haplotype/loop_number/loop_B6_Cast_hap/allloop/adjust/adjustloop.csv',sep = '\t')
strainbias = loop[(loop['B6_adjustP']<=0.1)&(loop['Cast_adjustP']<=0.1)]
strainbias = np.array(list(strainbias.iloc[:,0:6]))
'''


#%%
sig_type = np.dtype({'names':['chr' , 'start' , 'end' , 'M_C' , 'P_C'],
                     'formats':['U4' , np.int , np.int , np.int , np.int]})
sig_type1 = np.dtype({'names':['chr' , 'start' , 'end' ],
                     'formats':['U4' , np.int , np.int ]})

CTCF = np.loadtxt('/public/home/shidetong/projects/sdt/chip-seq/haplotype/diffbind/number4/depth/depth/al_peak_nochr_var.csv' , dtype = sig_type)
Rad21 = np.loadtxt('/public/home/shidetong/wedata/test/Allelic/RNA_allelic/Cast_depth.csv' , dtype = sig_type)

signal = {'CTCF':CTCF , 'Rad21':Rad21}

#%%
bias = {'CTCF':[[],[]] , 'Rad21':[[],[]]}

bias['CTCF'] = TF_Loop_Allelic_corr(GM_Loop, signal['CTCF'])
bias['Rad21'] = TF_Loop_Allelic_corr(GM_Loop, signal['Rad21'])

fig1 = Bias_scatter_Plot(bias['CTCF'], 'RNA')
fig2 = Bias_scatter_Plot(bias['Rad21'], 'CBF1_Cast')

# pp = PdfPages('/public/home/shidetong/wedata/test/Allelic/Allelic_Loop&BCF1_CTCF_scatter.pdf')
# pp.savefig(fig1)
# pp.close()

# pp = PdfPages('/public/home/shidetong/wedata/test/chip/chip_ALLelic/Allelic_Loop&Rad21_scatter.pdf')
# pp.savefig(fig2)
# pp.close()
 # %%
