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
            New_Data[c].extend([i['PCA']])
#             New_Data[c].extend([i['PCA']] * 10)
            
    
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
 

# def bigwig_100bp_Plot(fil,chro,start,end,fig,location,color,label):
#     """
#     """
#     sig_type = np.dtype({'names':['start' , 'end' , 'score'],
#                       'formats':[np.int , np.int , np.float]})
#     bw = pyBigWig.open(fil)
#     bw = bw.intervals(chro, start, end)
#     tmp_data = np.array(list(bw) , dtype = sig_type)
#     bin_size = (end - start)*5000 +1
#     sig_data = np.zeros((bin_size,))
#     for line in tmp_data:
#         s = line['start'] *5000 - start *5000
#         e = line['end'] *5000 - start *5000
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
# #     ax.set_xticks(ticks)
# #     ax.set_xticklabels(labels)
#     ax.set_xticklabels([])
#     ax.set_ylim((0 , 200))
#     ax.set_ylabel(label,fontsize = 15,rotation = 'horizontal',labelpad = 50)
#     caxis_S_horizontal(ax,color)
    


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
    ax.set_ylim((0 , 500))
    ax.set_ylabel(label,fontsize = 15,rotation = 'horizontal',labelpad = 50)
#     caxis_S_horizontal(ax,color)




cell = ['B6','Cast','B6_M' , 'B6_P' , 'Cast_M' , 'Cast_P' ]
cell1 = ['B6','Cast']
cell = ['B6_M' , 'B6_P' , 'Cast_M' , 'Cast_P' ]
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



B6M_PCA = pca_To_20K('/public/home/shidetong/projects/sdt/hic-seq/haplotype/B6/B6_Matrix_200K/Cooler/Maternal_PC/Maternal_PC_Compartment_200K.txt')
B6P_PCA = pca_To_20K('/public/home/shidetong/projects/sdt/hic-seq/haplotype/B6/B6_Matrix_200K/Cooler/Paternal_PC/Paternal_PC_Compartment_200K.txt')
CastM_PCA = pca_To_20K('/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_Matrix_200k/Cooler/Maternal_PC/Maternal_PC_Compartment_200K.txt')
CastP_PCA = pca_To_20K('/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_Matrix_200k/Cooler/Paternal_PC/Paternal_PC_Compartment_200K.txt')
PCA_Data = {'B6_M':B6M_PCA,
            'B6_P':B6P_PCA,
            'Cast_M':CastM_PCA,
            'Cast_P':CastP_PCA}



B6_M_RNA = '/public/home/shidetong/projects/lxx/RNA-seq/haplotype/B6_M/bw/B6_RNA_M_merge_100bp.bw'
B6_P_RNA = '/public/home/shidetong/projects/lxx/RNA-seq/haplotype/B6_M/bw/B6_RNA_P_merge_100bp.bw'
Cast_M_RNA = '/public/home/shidetong/projects/lxx/RNA-seq/haplotype/Cast_M/bw/Cast_RNA_M_merge_100bp.bw'
Cast_P_RNA = '/public/home/shidetong/projects/lxx/RNA-seq/haplotype/Cast_M/bw/Cast_RNA_P_merge_100bp.bw'
RNA_Data = {'B6_M':B6_M_RNA,'B6_P':B6_P_RNA,'Cast_M':Cast_M_RNA,'Cast_P':Cast_P_RNA}




pvalue = 0.05
pp1 = PdfPages('/public/home/shidetong/projects/sdt/hic-seq/haplotype/PC_A2B/PCsit/PC_classify/PC_gene/plot/Castbias_B6_M.pdf')
pp2 = PdfPages('/public/home/shidetong/projects/sdt/hic-seq/haplotype/PC_A2B/PCsit/PC_classify/PC_gene/plot/Castbias_B6_P.pdf')
pp3 = PdfPages('/public/home/shidetong/projects/sdt/hic-seq/haplotype/PC_A2B/PCsit/PC_classify/PC_gene/plot/Castbias_Cast_M.pdf')
pp4 = PdfPages('/public/home/shidetong/projects/sdt/hic-seq/haplotype/PC_A2B/PCsit/PC_classify/PC_gene/plot/Castbias_Cast_P.pdf')


selected_interval = [('5',130400000, 131000000)]
R =200000



for i in selected_interval:
    print (i)
    g = i[0]
    startHiC = i[1] // R 
    endHiC = i[2] // R
    n = 0
    max_rna = []
    for c in cell:
#         lib = HiC_Data[c][g]
#         DIData = DI_Data[c][g]
#         matrix = lib[startHiC:endHiC , startHiC:endHiC]
#         nonzero = matrix[np.nonzero(matrix)]
#         vmax = np.percentile(nonzero, 95)
#         DI = np.array(DIData[startHiC : endHiC])
        PCAData = PCA_Data[c][g]
        PCA = np.array(PCAData[startHiC:endHiC])
        
                
        fig = plt.figure(figsize = size)
#         ax1 = fig.add_axes([Left  , HB , width , HH])
#         sc = ax1.imshow(matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
#                        extent = (0, len(matrix), 0, len(matrix)), vmax =30, origin = 'lower')
#         cxlim = ax1.get_xlim()
#         cylim = ax1.get_ylim()
#         ## Ticks and Labels
#         ticks = list(np.linspace(0 , len(matrix) , 5).astype(float))
#         pos = [((startHiC + t) * R) for t in ticks]
#         labels = [properU(p) for p in pos[:5]]
#         ax7.set_xticks(ticks)
#         ax7.set_xticklabels(labels)
#         ax1.set_yticks(ticks)
#         ax1.set_yticklabels(labels)
# #         ax1.set_yticklabels(labels, rotation = 'horizontal')
#         ax1.set_xlabel('chr'+g , fontsize = 30 )
    
        PCA_index , PCA = UpdateDI(PCA)
        ax7 = fig.add_axes([Left, HB+0.1 , width ,0.075])
        ax7.fill_between(PCA_index , PCA , where = PCA >= 0 , facecolor = '#E47833' , edgecolor = 'none' )
        ax7.fill_between(PCA_index , PCA , where = PCA <= 0 , facecolor = '#7093DB' , edgecolor = 'none' )
        ax7.set_xlim(0 , PCA_index.max())
        ax7.set_ylim((-0.01 , 0.01))
#         ax7.xtick(fontsize=5)
        ax7.set_ylabel('PC1',fontsize = 15,rotation = 'horizontal',labelpad = 50)
        caxis_S_horizontal(ax7, 'black')
        
#         if c == 'B6_M' or c== 'B6_P' or c== 'Cast_M' or c== 'Cast_P':
#             chro = 'chr' + g
#         else:
#             chro = g 
        
#         location = [Left , HB , width , 0.075]
#         ax4 = fig.add_axes(location)
#         bigwig_1bp_Plot(RNA_Data[c] , chro , i[1] , i[2] , fig , location , 'black' , 'RNA')
#         caxis_S_horizontal(ax4, 'black')
# # #         caxis_H(ax4)
        
# #         location = [Left , HB + HH + 0.05++0.05+0.075+0.075, width , 0.075]
# #         ax4 = fig.add_axes(location)
# #         bigwig_10bp_Plot(RNA_Data['B6_M'] , chro , i[1] , i[2] , fig , location , 'fuchsia' , 'RNA')
# # #         caxis_H(ax4)
        
        tmp[c] = fig
    
    
    ccs_fig = tmp['B6_M']
    fesc_fig = tmp['B6_P']
    nt5_fig = tmp['Cast_M']
    nt6_fig = tmp['Cast_P']
# #     B6 = tmp['B6']
# #     Cast = tmp['Cast']
 
    
    figs = [ccs_fig,fesc_fig,nt5_fig,nt6_fig]
    standard_axes_lim(0,*figs)
#     standard_axes_lim(1,*figs)
#     standard_axes_lim(2,*figs)
  

    pp1.savefig(ccs_fig)
    pp2.savefig(fesc_fig)
    pp3.savefig(nt5_fig)
    pp4.savefig(nt6_fig)
# #     pp5.savefig(B6)
# #     pp6.savefig(Cast)
    plt.close(ccs_fig)
    plt.close(fesc_fig)
    plt.close(nt5_fig)
    plt.close(nt6_fig)
# #     plt.close(B6)
# #     plt.close(Cast)

pp1.close()
pp2.close()
pp3.close()
pp4.close()
        
        