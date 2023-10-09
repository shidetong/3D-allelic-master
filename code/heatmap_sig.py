%matplotlib inline
from __future__ import division
import numpy as np
#from tadlib.calfea.analyze import getmatrix
import matplotlib
# Use a non-interactive backend
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import os
import pickle
import pyBigWig
from itertools import islice
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
# my_cmap = LinearSegmentedColormap.from_list('interaction',
#                                             ['#CD0000','#191970'])
# my_cmap.set_bad('#2672a1')


def binbias(up, down):
    bias = 0
    zeromask = (up != 0) & (down != 0)
    if  zeromask.sum() <= 1:
        return bias
        
    upmean = up.mean(); downmean = down.mean()
    SD_1 = np.sum((up - upmean) ** 2) / (up.size * (up.size - 1))
    SD_2 = np.sum((down - downmean) ** 2) / (down.size * (down.size - 1))
    SD_pool = np.sqrt(SD_1 + SD_2)
    if SD_pool != 0:
        bias = (upmean - downmean) / SD_pool

    return bias
    
    
def CalDI(matrix):
    n = 0
    DIs = []
    shape = matrix.shape
    for j in matrix[:shape[0]]:
        if n <20:
            DIs.append(0)
        elif n >= 20 and n <= int(shape[0])-21:
            if len(j[j!=0])/int(shape[0]) < 0.05:
                bias = 0
            else:
                up = j[n-20:n] 
                down = j[n+1:n+21]
                bias = binbias(up,down) 
            DIs.append(bias)
        else:
            DIs.append(0)
        n += 1
    return DIs


def caxis_H(ax):
    """
    Axis Control for HeatMaps.
    
    """
    ax.yaxis.set_ticks_position('left')
    #ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis = 'both', bottom = False, top = False, left = True,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = True, labelright = False , length = 5 , labelsize = 10)
                   

    
def caxis_colorbar(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.tick_params(axis = 'both', bottom = True, top = False, left = False,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = False, labelright = False , labelsize = 25)
    
    

def caxis_S_horizontal(ax, color):
    """
    Axis Control for PCA plots.
    """
    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis = 'both', bottom = False, top = False, left = True,
                   right = False, labelbottom = False, labeltop = False,
                   labelleft = True, labelright = False , labelsize = 15)
    ax.spines['left'].set_lw(1.5)
    ax.spines['left'].set_color(color)
    ax.spines['left'].set_alpha(0)
    ax.spines['left'].set_linestyle('dotted')

    
def caxis_DI(ax, color):
    """
    Axis Control for PCA plots.
    """
    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis = 'both', bottom = True, top = False, left = False,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = False, labelright = False)
    ax.spines['bottom'].set_lw(1.5)
    ax.spines['bottom'].set_color(color)
    ax.spines['bottom'].set_alpha(0.9)
    ax.spines['bottom'].set_linestyle('dotted')
    

def properU(pos):
    """
    Express a genomic position in a proper unit (KB, MB, or both).
    
    """
    i_part = int(pos) // 1000000 # Integer Part
    d_part = (int(pos) % 1000000) // 1000 # Decimal Part
    
    if (i_part > 0) and (d_part > 0):
        return ''.join([str(i_part), 'M', str(d_part), 'K'])
    elif (i_part == 0):
        return ''.join([str(d_part), 'K'])
    else:
        return ''.join([str(i_part), 'M'])
    
def DI_Calling(HiC_fil):
    """
    """
    HiC_Data = np.load(HiC_fil,allow_pickle=True,encoding='latin1')
    Lib = HiC_Data['arr_0'][()]
    DI_Data = {}
    for chro in Lib.keys():
        Matrix = Lib[chro]
        DI_score = CalDI(Matrix)
        DI_Data[chro] = DI_score
    
    return DI_Data



        
def pca_To_20K(fil):
    """
    """
    pca_type = np.dtype({'names':['chr' , 'PCA'],
                     'formats':['S4' , np.float]})
    PCA_Data = np.loadtxt(fil , dtype=pca_type)
    
    chroms = set(PCA_Data['chr'])
    New_Data = {}
    for c in chroms:
        New_Data[c] = {}
        tmp_data = PCA_Data[PCA_Data['chr'] == c]
        New_Data[c] = []
        for i in tmp_data:
            New_Data[c].extend([i['PCA']] * 10)
            
    
    return New_Data



def standard_axes_lim(N,*figs):
    """
    """    
    lim = []
    lim_1 = []
    for s_f in figs:
        lim.append(s_f.axes[N].get_ylim()[1])
        lim_1.append(s_f.axes[N].get_ylim()[0])
    
    for s_f in figs:
        s_f.axes[N].set_ylim(np.round(min(lim_1), 2 ) - 0.01,np.round(max(lim) , 2) + 0.01)
        
        
        
def UpdateDI(DI):
    """
    """
    New_DI = []
    New_index = []

    for index in range(len(DI) - 1):
        if DI[index] * DI[index + 1] < 0:
            New_DI.append(DI[index])
            New_DI.append(0)
            New_index.append(index)
            New_index.append(index + 0.5)
        else:
            New_DI.append(DI[index])
            New_index.append(index)
    
    return np.array(New_index), np.array(New_DI)

        
def Get_loops(LoopSource):
    """
    cell : String ,The name of cells  e.g.'fESC' , 'ESC' , 'NT5' , 'NT6' , 'CCS'
    get the loops' location as a dictionary
    """
    loops = []
    Loop = np.loadtxt(LoopSource, usecols = (0,1,2) , dtype = loop_type, skiprows = 1)
    for i in Loop:
        if i['end'] - i['start'] >= 40000:
            loops.append(i) 
        else:
            continue
    loops = np.array(loops , dtype = loop_type)
    return loops


def get_loc(loopspp):
    loops = []
    Loop = np.loadtxt(loopspp, usecols = (0,11,12) , dtype = loop_type, skiprows = 1)
    for i in Loop:
        if i['end'] - i['start'] >= 40000:
            loops.append(i) 
        else:
            continue
    loops = np.array(loops , dtype = loop_type)
    loops = loops.tolist()
    return loops 
    
    
# def Get_loops(LoopSource):
#     """
#     cell : String ,The name of cells  e.g.'fESC' , 'ESC' , 'NT5' , 'NT6' , 'CCS'
#     get the loops' location as a dictionary
#     """
 
#     Loop = np.loadtxt(LoopSource, usecols = (0,1,2) , dtype = loop_type, skiprows = 1)
#     loops = np.array(loop , dtype = loop_type)
#     return loops  


def Load_gtf(gtfil):
    gtf_type = np.dtype({'names':['gene_id' , 'gene_name' , 'chr' , 'strand' , 'start' , 'end'],
                     'formats':['S64' , 'S64' , 'S8' , 'S4' , np.int , np.int]})
    gtf = open(gtfil , 'r')
    gtf_1 = []
    for i in islice(gtf , 5 , None):
        a = i.strip().split()
        if a[2] == 'gene':
            gene_id = i.strip().split('\"')[1]
            gene_name = i.strip().split('\"')[5]
            chro = a[0]
            strand = a[6]
            start = a[3]
            end = a[4]
            gtf_1.append((gene_id , gene_name , chro , strand , start , end))
    gtf = np.array(gtf_1 , dtype = gtf_type)
    return gtf
    
    
# def bigwig_10bp_Plot(fil,chro,start,end,fig,location,color,label):
#     """
#     """
#     sig_type = np.dtype({'names':['start' , 'end' , 'score'],
#                       'formats':[np.int , np.int , np.float]})
#     bw = pyBigWig.open(fil)
#     bw = bw.intervals(chro, start, end)
#     tmp_data = np.array(list(bw) , dtype = sig_type)
#     bin_size = (end - start) // 100 + 1
#     sig_data = np.zeros((bin_size,))
#     for line in tmp_data:
#         s = line['start'] // 100 - start // 100
#         e = line['end'] // 100 - start // 100
#         for i in range(s,e):
#             if i >= 0 and i < bin_size: 
#                 sig_data[i] += line['score']
#             else:
#                 pass
#     ax = fig.add_axes(location)
#     ax.fill_between(np.arange(len(sig_data)),sig_data, facecolor = color, edgecolor = 'none')
# #     ax.set_xlim((0 , len(sig_data)))
#     ticks = list(np.linspace(0 , len(sig_data) , 5).astype(float))
#     pos = [((startHiC + t) * R) for t in ticks]
#     labels = [properU(p) for p in pos[:5]]

#     cxlim = ax.get_xlim()
#     ax.set_xticks(ticks)
#     ax.set_xticklabels(labels)
#     ax.set_ylim((0 , 200))
#     ax.set_ylabel(label,fontsize = 15,rotation = 'horizontal',labelpad = 50)
#     caxis_S_horizontal(ax,color)

def bigwig_10bp_Plot(fil,chro,start,end,fig,location,color,label):
    """
    """
    sig_type = np.dtype({'names':['start' , 'end' , 'score'],
                      'formats':[np.int , np.int , np.float]})
    bw = pyBigWig.open(fil)
    bw = bw.intervals(chro, start, end)
    tmp_data = np.array(list(bw) , dtype = sig_type)
    bin_size = (end - start)//100 +1
    sig_data = np.zeros((bin_size,))
    for line in tmp_data:
        s = line['start'] //100 - start //100
        e = line['end'] // 100 - start // 100
        for i in range(s,e):
            if i >= 0 and i < bin_size: 
                sig_data[i] += line['score']
            else:
                pass
    ax = fig.add_axes(location)
    ax.fill_between(np.arange(len(sig_data)),sig_data, facecolor = color, edgecolor = 'none')
#     ax.set_xlim((0 , len(sig_data)))
    ticks = list(np.linspace(0 , len(sig_data) , 5).astype(float))
    pos = [((startHiC + t) * R) for t in ticks]
    labels = [properU(p) for p in pos[:5]]

    cxlim = ax.get_xlim()
#     ax.set_xticks(ticks)
#     ax.set_xticklabels(labels)
    ax.set_xticklabels([])
    ax.set_ylim((0 , 200))
    ax.set_ylabel(label,fontsize = 15,rotation = 'horizontal',labelpad = 50)
    caxis_S_horizontal(ax,color)
    


def bigwig_1bp_Plot(fil,chro,start,end,fig,location,color,label):
    """
    """
    sig_type = np.dtype({'names':['start' , 'end' , 'score'],
                      'formats':[np.int , np.int , np.float]})
    bw = pyBigWig.open(fil)
    bw = bw.intervals(chro, start, end)
    tmp_data = np.array(list(bw) , dtype = sig_type)
    bin_size = (end - start)//100 +1
    sig_data = np.zeros((bin_size,))
    for line in tmp_data:
        s = line['start'] //100 - start //100
        e = line['end'] // 100 - start // 100
        for i in range(s,e):
            if i >= 0 and i < bin_size: 
                sig_data[i] += line['score']
            else:
                pass
    ax = fig.add_axes(location)
    ax.fill_between(np.arange(len(sig_data)),sig_data, facecolor = color, edgecolor = 'none')
#     ax.set_xlim((0 , len(sig_data)))
    ticks = list(np.linspace(0 , len(sig_data) , 5).astype(float))
    pos = [((startHiC + t) * R) for t in ticks]
    labels = [properU(p) for p in pos[:5]]

    cxlim = ax.get_xlim()
#     ax.set_xticks(ticks)
#     ax.set_xticklabels(labels)
    ax.set_xticklabels([])
    ax.set_ylim((0 , 300))
    ax.set_ylabel(label,fontsize = 15,rotation = 'horizontal',labelpad = 50)
    caxis_S_horizontal(ax,color)



    
def Sig_Plot(data,start,end,chro,fig,location,color,label):
    """
    data: signal data 
    
    """
    tmp = data[chro]
    # sig_start = start // R - 50
    # sig_end = end // R + 50
    sig_data = tmp[start:end]
    max_ = sig_data.max()
    

    ax = fig.add_axes(location)
    ax.fill_between(np.arange(len(sig_data)),sig_data, facecolor = color, edgecolor = 'none')

    ax.set_ylabel(label,fontsize = 15,rotation = 'horizontal',labelpad = 50)
    ax.set_xlim((0,len(sig_data)))
    # ax.set_xlabel(chro,fontsize = 5,rotation = 'horizontal',labelpad = 20)
    
    caxis_S(ax)
    return max_



cell = ['B6','Cast','B6_M' , 'B6_P' , 'Cast_M' , 'Cast_P' ]
cell1 = ['B6','Cast']
cell2 = ['B6_M' , 'B6_P' , 'Cast_M' , 'Cast_P' ]
cells = {'B6':0,'Cast':1,'B6_M':2 , 'B6_P':3 , 'Cast_M':4 , 'Cast_P':5 }



size = (12, 12)
# Left = 0.25 ; HB = 0.1 ; width = 0.6 ; HH = 0.6
Left = 0.1 ; HB = 0.1 ; width = 0.8 ; HH = 0.8

tad_type = np.dtype({'names':['chr' , 'start' , 'end'],
                     'formats':['S4' , np.int , np.int]})
di_type = np.dtype({'names':['chr' , 'DI'],
                     'formats':['S4' , np.float]})

sig_type = np.dtype({'names':['chr','start' , 'end' , 'score'],
                      'formats':['S4',np.int , np.int , np.float]})
    
loop_type = np.dtype({'names':['chr' , 'start' , 'end'],
                      'formats':['S4' , np.int , np.int]})

    
    
f = open('/public/home/shidetong/projects/sdt/snp/B6_M/genome/genomeSize' , 'r')
#f = open('/public/home/shidetong/projects/sdt/snp/Cast_M/genome/genomeSize' , 'r')
mm ={}
for i in f:
    i = i.split()
    mm[i[0]] = int(i[1])
tmp = {}




#单倍型与未拆分
B6 = np.load('/public/home/shidetong/projects/sdt/hic-seq/haplotype/B6/B6_Matrix_40K/Cooler_40K/cooltonpz/Merged_Traditional_Multi.npz',allow_pickle=True)
B6_Lib = B6['arr_0'][()]
Cast = np.load('/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_Matrix_40K/Cooler/cooltonpz/Merged_Traditional_Multi.npz',allow_pickle=True)
Cast_Lib = Cast['arr_0'][()]
B6_M_Lib = np.load('/public/home/shidetong/projects/sdt/hic-seq/haplotype/B6/B6_Matrix_40K/Cooler_40K/cooltonpz/Merged_Imputated_Haplotype_Multi_M.npz',allow_pickle=True)
B6_M_Lib = B6_M_Lib['arr_0'][()]
B6_P_Lib = np.load('/public/home/shidetong/projects/sdt/hic-seq/haplotype/B6/B6_Matrix_40K/Cooler_40K/cooltonpz/Merged_Imputated_Haplotype_Multi_P.npz',allow_pickle=True)
B6_P_Lib = B6_P_Lib['arr_0'][()]
Cast_M_Lib = np.load('/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_Matrix_40K/Cooler/cooltonpz/Merged_Imputated_Haplotype_Multi_M.npz',allow_pickle=True)
Cast_M_Lib = Cast_M_Lib['arr_0'][()]
Cast_P_Lib = np.load('/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_Matrix_40K/Cooler/cooltonpz/Merged_Imputated_Haplotype_Multi_P.npz',allow_pickle=True)
Cast_P_Lib = Cast_P_Lib['arr_0'][()]

HiC_Data = {'B6':B6_Lib,'Cast':Cast_Lib,'B6_M':B6_M_Lib,'B6_P':B6_P_Lib,'Cast_M':Cast_M_Lib,'Cast_P':Cast_P_Lib}





#传统hic
B6_M_Lib = np.load('/public/home/shidetong/projects/sdt/hic-seq/HiC_workspace/cooler2np/B61.npz',allow_pickle=True)
B6_M_Lib = B6_M_Lib['arr_0'][()]
B6_P_Lib = np.load('/public/home/shidetong/projects/sdt/hic-seq/HiC_workspace/cooler2np/B62.npz',allow_pickle=True)
B6_P_Lib = B6_P_Lib['arr_0'][()]
Cast_M_Lib = np.load('/public/home/shidetong/projects/sdt/hic-seq/HiC_workspace/cooler2np/Cast1.npz',allow_pickle=True)
Cast_M_Lib = Cast_M_Lib['arr_0'][()]
Cast_P_Lib = np.load('/public/home/shidetong/projects/sdt/hic-seq/HiC_workspace/cooler2np/Cast2.npz',allow_pickle=True)
Cast_P_Lib = Cast_P_Lib['arr_0'][()]
sample = ['B6_1','B6_2','Cast_1','Cast_2']
HiC_Data = {'B6_1':B6_M_Lib,'B6_2':B6_P_Lib,'Cast_1':Cast_M_Lib,'Cast_2':Cast_P_Lib}



selected_interval = [('1',7000000,11000000)]
R =40000
res = '40K'

pvalue = 0.05
pp1 = PdfPages('/public/home/shidetong/projects/sdt/hic-seq/haplotype/plot_matrix/Pbias/B6_M.pdf')
pp2 = PdfPages('/public/home/shidetong/projects/sdt/hic-seq/haplotype/plot_matrix/Pbias/B6_P.pdf')
pp3 = PdfPages('/public/home/shidetong/projects/sdt/hic-seq/haplotype/plot_matrix/Pbias/Cast_M.pdf')
pp4 = PdfPages('/public/home/shidetong/projects/sdt/hic-seq/haplotype/plot_matrix/Pbias/Cast_P.pdf')
pp5 = PdfPages('/public/home/shidetong/projects/sdt/hic-seq/haplotype/plot_matrix/Pbias/B6.pdf')


for i in selected_interval:
    print (i)
    g = i[0]
    startHiC = i[1] // R 
    endHiC = i[2] // R 
    for c in cell:
#     for c in sample:
        lib = HiC_Data[c][g]
#        DIData = DI_Data[c][g]
        matrix = lib[startHiC:endHiC , startHiC:endHiC]
        nonzero = matrix[np.nonzero(matrix)]
        vmax = np.percentile(nonzero, 95)
        # DI = np.array(DIData[startHiC : endHiC])
        # PCAData = PCA_Data[c][g]
        # PCA = np.array(PCAData[startHiC:endHiC])
#         loops = Loop_Data[c][Loop_Data[c]['chr'] == g]
        
        
        #=============HeatMap + colorbar + DI=================================
        
        
        fig = plt.figure(figsize = size)
        ax1 = fig.add_axes([Left  , HB , width , HH])
        sc = ax1.imshow(matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                       extent = (0, len(matrix), 0, len(matrix)), vmax =vmax, origin = 'lower')
        cxlim = ax1.get_xlim()
        cylim = ax1.get_ylim()
        ## Ticks and Labels
        ticks = list(np.linspace(0 , len(matrix) , 5).astype(float))
        pos = [((startHiC + t) * R) for t in ticks]
        labels = [properU(p) for p in pos[:5]]
        ax1.set_xticks(ticks)
        ax1.set_xticklabels(labels)
        ax1.set_yticks(ticks)
        ax1.set_yticklabels(labels)
#         ax1.set_yticklabels(labels, rotation = 'horizontal')
#         ax1.set_xlabel('chr'+g , fontsize = 30 )
#         plt.savefig('/public/home/shidetong/'+c+'.png')
#         ax1.set_xlabel(c , fontsize = 30 )
#         plt.xticks([])  #去掉横坐标值
#         plt.yticks([])
#         plt.savefig('/public/home/shidetong/projects/sdt/hic-seq/haplotype/plot_matrix/'+c+'_chr2_13M_14M_500max.png')
        
#         # loop data
#         mask = (loops['start'] >= startHiC * R) & (loops['end'] < endHiC * R)
#         extract = loops[mask]
#         for p in extract:
#             y = p['start'] // R - startHiC
#             x = p['end'] // R - startHiC
#             ax1.scatter(y + 0.5 ,x + 0.5 , color = '', edgecolors = 'black', s =2500)
            
#         ax1.set_xlim(cxlim)
#         ax1.set_ylim(cylim)                    
#         caxis_H(ax1)
#         Colorbar
        ax2 = fig.add_axes([Left + width + 0.015 , HB , 0.035 , 0.1])
#         fig.colorbar(sc,cax = ax2, orientation='vertical' , ticks = [matrix.min() , matrix.max()])
        fig.colorbar(sc,cax = ax2, orientation='vertical' , ticks = [0 , vmax])
#         plt.savefig('/public/home/shidetong/projects/sdt/hic-seq/haplotype/plot_matrix/'+c+'bar.png')
# #         ##chip_2
#         location = [Left + width + 0.03,HB,0.035,0.1]
#         ax6 = fig.add_axes(location)
#         bigwig_10bp_Plot(Chip_Data[c] , chro , i[1] , i[2] , fig , location , 'mediumblue' , 'CTCF')   
#         caxis_H(ax6)
        

        ##gene Tracks
#         ax3 = fig.add_axes([Left , HB + HH +0.15, width , 0.04])
#         ax3.fill_between('rr','GM16432','ee',facecolor = '#E47833')
#         caxis_S_horizontal(ax3, 'black')
        
#         if c == 'B6_M' or c== 'B6_P' or c== 'Cast_M' or c== 'Cast_P':
#             chro = 'chr' + g
#         else:
#             chro = g 
        #RNA Tracks
#         location = [Left , HB + HH + 0.03+0.1, width , 0.075]
#         ax4 = fig.add_axes(location)
#         bigwig_1bp_Plot(RNA_Data[c] , chro , i[1] , i[2] , fig , location , 'fuchsia' , 'RNA')
#         caxis_H(ax4)
#         ##Chip-input Tracks
#         location = [Left , HB + HH + 0.03 , width , 0.075]
#         ax5 = fig.add_axes(location)
#         bigwig_10bp_Plot(Chip_Data_input[c] , chro , i[1] , i[2] , fig , location , 'green' , 'Chip-input')
        
        #Chip Tracks
#         location = [Left , HB + HH + 0.03 , width , 0.075]
#         ax5 = fig.add_axes(location)
#         bigwig_10bp_Plot(Chip_Data[c] , chro , i[1] , i[2] , fig , location , 'blue' , 'CTCF')   
#         caxis_H(ax5)
        
#         tmp[c] = fig
    
    
#     ccs_fig = tmp['B6_M']
#     fesc_fig = tmp['B6_P']
#     nt5_fig = tmp['Cast_M']
#     nt6_fig = tmp['Cast_P']
#     B6 = tmp['B6']
#     Cast = tmp['Cast']
 
    
#     figs = [ccs_fig,fesc_fig,nt5_fig,nt6_fig]
#     standard_axes_lim(0,*figs)
# #     standard_axes_lim(2,*figs)
# #     standard_axes_lim(3,*figs)
  

#     pp1.savefig(ccs_fig)
#     pp2.savefig(fesc_fig)
#     pp3.savefig(nt5_fig)
#     pp4.savefig(nt6_fig)
#     pp5.savefig(B6)
#     pp6.savefig(Cast)
#     plt.close(ccs_fig)
#     plt.close(fesc_fig)
#     plt.close(nt5_fig)
#     plt.close(nt6_fig)
#     plt.close(B6)
#     plt.close(Cast)

# pp1.close()
# pp2.close()
# pp3.close()
# pp4.close()
# pp5.close()
# pp6.close()