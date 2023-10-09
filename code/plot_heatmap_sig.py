

#-- coding:UTF-8 --
from __future__ import division
import numpy as np
#from tadlib.calfea.analyze import getmatrix
import matplotlib
# Use a non-interactive backend
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import os
import pyBigWig
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
                   right = False, labelbottom = False, labeltop = False,
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
                   labelleft = True, labelright = False , labelsize = 23)
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
        if i['end'] - i['start'] >= 300000:
            loops.append(i) 
        else:
            continue
    loops = np.array(loops , dtype = loop_type)
    return loops
    
def bigwig_10bp_Plot(fil,chro,start,end,fig,location,color,label):
    """
    """
    sig_type = np.dtype({'names':['start' , 'end' , 'score'],
                      'formats':[np.int , np.int , np.float]})
    bw = pyBigWig.open(fil)
    bw = bw.intervals(chro, start, end)
    tmp_data = np.array(list(bw) , dtype = sig_type)
    bin_size = (end - start) // 100 + 1
    sig_data = np.zeros((bin_size,))
    for line in tmp_data:
        s = line['start'] // 100 - start // 100
        e = line['end'] // 100 - start // 100
        for i in range(s,e):
            if i >= 0 and i < bin_size: 
                sig_data[i] += line['score']
            else:
                pass
    ax = fig.add_axes(location)
    ax.fill_between(np.arange(len(sig_data)),sig_data, facecolor = color, edgecolor = 'none')
    ax.set_xlim((0 , len(sig_data)))
    ax.set_ylabel(label,fontsize = 15,rotation = 'horizontal',labelpad = 50)
    caxis_S_horizontal(ax,color)


OutFolder = '/public/home/shidetong/projects/sdt/hic-seq/haplotype/plot_heatmap/strain'

cells = {'B6_M':0 , 'B6_P':1 , 'Cast_M':2 , 'Cast':3 }
cell = ['B6_M' , 'B6_P' , 'Cast_M' , 'Cast_P']
R = 40000
res = '40K'

#selected_interval = [('3' , 33760000 , 35680000 , 'Sox2') , ('4' , 55200000 , 57000000 , 'Klf4') , 
#                     ('17' , 35280000 , 36360000 , 'Oct4') , ('6' , 121400000 , 123320000 , 'Nanog') , 
#                     ('12' , 86080000 , 87120000 , 'Esrrb') , ('2' , 106900000 , 108000000 , 'Fig2')]
selected_interval = [('1',170000000,17200000),('1',176000000,178000000),('1',177000000,179000000),('1',30000000,32000000),
                    ('2',54000000,56000000),('2',164000000,165000000),('3',90000000,92000000),('4',143000000,145000000),
                    ('4',121000000,124000000),('5',29000000,32000000),('5',145000000,147000000),('5',3000000,5000000),
                    ('5',5000000,7000000),('6',59000000,61000000),('6',125000000,127000000),('6',69000000,71000000),
                    ('6',123000000,125000000),('7',5000000,7000000),('7',24000000,26000000),('7',26000000,28000000),
                    ('7',41000000,44000000),('7',59000000,61000000),('7',85000000,87000000),('7',17000000,19000000),
                    ('7',41000000,43000000),('8',19000000,21000000),('9',35000000,37000000),('10',3000000,5000000),
                    ('11',3000000,6000000),('11',73000000,75000000),('12',93000000,95000000),('11',101000000,103000000),
                    ('13',12000000,14000000),('15',89000000,92000000),('17',6000000,8000000),('19',7000000,9000000),
                    ('19',33000000,35000000),('X',91000000,93000000)]

#selected_interval = [('11' , 84000000 , 88000000 )]

# genomeSize = {'1' : 195471971,'2' : 182113224,'3' : 160039680,
#                 '4' : 156508116,'5' : 151834684,'6' : 149736546,
#                 '7' : 145441459,'8' : 129401213,'9' : 124595110,
#                 '10' : 130694993,'11' : 122082543,'12' : 120129022,
#                 '13' : 120421639,'14' : 124902244,'15' : 104043685,
#                 '16' : 98207768,'17' : 94987271,'18' : 90702639,
#                 '19' : 61431566,'X' : 171031299} 

# genome_size = [('1' , 0 , 195471971),('2', 0 , 182113224) , ('3' , 0 , 160039680),
#                  ('4',0 , 156508116),('5',0,151834684),('6',0,149736546),
#                  ('7',0,145441459),('8',0,129401213),('9',0,124595110),
#                  ('10',0,130694993),('11',0,122082543),('12',0,120129022),
#                  ('13',0,120421639),('14',0,124902244),('15',0,104043685),
#                  ('16',0,98207768),('17',0,94987271),('18',0,90702639),
#                  ('19',0,61431566),('X',0,171031299)] 
# sele =[]
# for i in genome_size :
#     chrom = i[0]
#     start = i[1]
#     end = i[1] + 4000000
#     sele.append((chrom,start,end))
#     while end <= i[2]-3999999:
#         start = end + 1 
#         end = start + 3999999
#         sele.append((chrom,start,end))
#     if end <= i[2]:
#         start = end +1 
#         end = i[2]
#         sele.append((chrom,start,end))



size = (12, 12)
Left = 0.25 ; HB = 0.05 ; width = 0.6 ; HH = 0.6

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


# B6_M_DI = DI_Calling('/public/home/shidetong/projects/sdt/hic-seq/haplotype/B6/new_B6_Matrix/Cooler/cooltonpz/Merged_Imputated_Haplotype_Multi_M.npz')
# B6_P_DI = DI_Calling('/public/home/shidetong/projects/sdt/hic-seq/haplotype/B6/new_B6_Matrix/Cooler/cooltonpz/Merged_Imputated_Haplotype_Multi_P.npz')
# Cast_M_DI = DI_Calling('/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_Matrix/Cooler/cooltpnpz/Merged_Imputated_Haplotype_Multi_M.npz')
# Cast_P_DI = DI_Calling('/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_Matrix/Cooler/cooltpnpz/Merged_Imputated_Haplotype_Multi_P.npz')
# DI_Data = {'B6_M':B6_M_DI,'B6_P':B6_P_DI,'Cast_M':Cast_M_DI,'Cast_P':Cast_P_DI}

B6_M_Lib = np.load('/public/home/shidetong/projects/sdt/hic-seq/haplotype/B6/new_B6_Matrix/Cooler/cooltonpz/Merged_Imputated_Haplotype_Multi_M.npz',allow_pickle=True)
B6_M_Lib = B6_M_Lib['arr_0'][()]
B6_P_Lib = np.load('/public/home/shidetong/projects/sdt/hic-seq/haplotype/B6/new_B6_Matrix/Cooler/cooltonpz/Merged_Imputated_Haplotype_Multi_P.npz',allow_pickle=True)
B6_P_Lib = B6_P_Lib['arr_0'][()]
Cast_M_Lib = np.load('/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_Matrix/Cooler/cooltpnpz/Merged_Imputated_Haplotype_Multi_M.npz',allow_pickle=True)
Cast_M_Lib = Cast_M_Lib['arr_0'][()]
Cast_P_Lib = np.load('/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_Matrix/Cooler/cooltpnpz/Merged_Imputated_Haplotype_Multi_P.npz',allow_pickle=True)
Cast_P_Lib = Cast_P_Lib['arr_0'][()]
HiC_Data = {'B6_M':B6_M_Lib,'B6_P':B6_P_Lib,'Cast_M':Cast_M_Lib,'Cast_P':Cast_P_Lib}

# B6_M_PCA = pca_To_20K('/public/home/shidetong/projects/sdt/hic-seq/haplotype/B6/new_B6_Matrix/Cooler/Maternal_PC/Maternal_PC_Compartment_500K.txt')
# B6_P_PCA = pca_To_20K('/public/home/shidetong/projects/sdt/hic-seq/haplotype/B6/new_B6_Matrix/Cooler/Paternal_PC/Paternal_PC_Compartment_500K.txt')
# Cast_M_PCA = pca_To_20K('/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_Matrix/Cooler/Maternal_PC/Maternal_PC_Compartment_500K.txt')
# Cast_P_PCA = pca_To_20K('/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_Matrix/Cooler/Paternal_PC/Paternal_PC_Compartment_500K.txt')
# PCA_Data = {'B6_M':B6_M_PCA,'B6_P':B6_P_PCA,'Cast_M':Cast_M_PCA,'Cast_P':Cast_P_PCA}

# B6_M_Loop = Get_loops('/public/home/shidetong/projects/sdt/hic-seq/haplotype/B6/new_B6_Matrix/Cooler/Maternal_Loops/Cluster_Maternal_Loops_Loops_40K.txt')
# B6_P_Loop = Get_loops('/public/home/shidetong/projects/sdt/hic-seq/haplotype/B6/new_B6_Matrix/Cooler/Paternal_Loops/Cluster_Paternal_Loops_Loops_40K.txt')
# Cast_M_Loop = Get_loops('/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_Matrix/Cooler/Maternal_Loops/Cluster_Maternal_Loops_Loops_40K.txt')
# Cast_P_Loop = Get_loops('/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_Matrix/Cooler/Paternal_Loops/Cluster_Paternal_Loops_Loops_40K.txt')
# Loop_Data = {'B6_M':B6_M_Loop,'B6_P':B6_P_Loop,'Cast_M':Cast_M_Loop,'Cast_P':Cast_P_Loop}

B6_M_RNA = '/public/home/shidetong/projects/lxx/RNA-seq/haplotype/B6_M/bw/B6_RNA_M_merge_100bp.bw'
B6_P_RNA = '/public/home/shidetong/projects/lxx/RNA-seq/haplotype/B6_M/bw/B6_RNA_P_merge_100bp.bw'
Cast_M_RNA = '/public/home/shidetong/projects/lxx/RNA-seq/haplotype/Cast_M/bw/Cast_RNA_M_merge_100bp.bw'
Cast_P_RNA = '/public/home/shidetong/projects/lxx/RNA-seq/haplotype/Cast_M/bw/Cast_RNA_P_merge_100bp.bw'
RNA_Data = {'B6_M':B6_M_RNA,'B6_P':B6_P_RNA,'Cast_M':Cast_M_RNA,'Cast_P':Cast_P_RNA}

B6_M_Chip = '/public/home/shidetong/projects/sdt/chip-seq/haplotype/B6_M/bw/B6_Chip_CTCF_M_merge_100bp.bw'
B6_P_Chip = '/public/home/shidetong/projects/sdt/chip-seq/haplotype/B6_M/bw/B6_Chip_CTCF_P_merge_100bp.bw'
Cast_M_Chip = '/public/home/shidetong/projects/sdt/chip-seq/haplotype/Cast_M/bw/Cast_Chip_CTCF_M_merge_100bp.bw'
Cast_P_Chip = '/public/home/shidetong/projects/sdt/chip-seq/haplotype/Cast_M/bw/Cast_Chip_CTCF_P_merge_100bp.bw'
Chip_Data = {'B6_M':B6_M_Chip,'B6_P':B6_P_Chip,'Cast_M':Cast_M_Chip,'Cast_P':Cast_P_Chip}

#chip-input
# B6_M_Chip_input = '/public/home/shidetong/projects/sdt/chip-seq/haplotype/B6_M/bw/Speci_Unique_B6_batch2_Chip_CTCF_1_input_M_100bp.bw'
# B6_P_Chip_input = '/public/home/shidetong/projects/sdt/chip-seq/haplotype/B6_M/bw/Speci_Unique_B6_batch2_Chip_CTCF_1_input_P_100bp.bw'
# Cast_M_Chip_input = '/public/home/shidetong/projects/sdt/chip-seq/haplotype/Cast_M/bw/Speci_Unique_Cast_batch2_Chip_CTCF_1_input_M_100bp.bw'
# Cast_P_Chip_input = '/public/home/shidetong/projects/sdt/chip-seq/haplotype/Cast_M/bw/Speci_Unique_Cast_batch2_Chip_CTCF_1_input_P_100bp.bw'
# Chip_Data_input = {'B6_M':B6_M_Chip_input,'B6_P':B6_P_Chip_input,'Cast_M':Cast_M_Chip_input,'Cast_P':Cast_P_Chip_input}

pp1 = PdfPages('/public/home/shidetong/projects/sdt/hic-seq/haplotype/plot_heatmap/strain/B6_M_Genetic_transmission.pdf')
pp2 = PdfPages('/public/home/shidetong/projects/sdt/hic-seq/haplotype/plot_heatmap/strain/B6_P_Genetic_transmission.pdf')
pp3 = PdfPages('/public/home/shidetong/projects/sdt/hic-seq/haplotype/plot_heatmap/strain/Cast_M_Genetic_transmission.pdf')
pp4 = PdfPages('/public/home/shidetong/projects/sdt/hic-seq/haplotype/plot_heatmap/strain/Cast_P_Genetic_transmission.pdf')

for i in selected_interval:
    print (i)
    g = i[0]
    startHiC = i[1] // R 
    endHiC = i[2] // R 
    for c in cell:
        lib = HiC_Data[c][g]
        #DIData = DI_Data[c][g]
        matrix = lib[startHiC:endHiC , startHiC:endHiC]
        nonzero = matrix[np.nonzero(matrix)]
        vmax = np.percentile(nonzero, 95)
        # DI = np.array(DIData[startHiC : endHiC])
        # PCAData = PCA_Data[c][g]
        # PCA = np.array(PCAData[startHiC:endHiC])
        #loops = Loop_Data[c][Loop_Data[c]['chr'] == g]
        
        
        #=============HeatMap + colorbar + DI=================================
        
        fig = plt.figure(figsize = size)
        ax1 = fig.add_axes([Left  , HB , width , HH])
        sc = ax1.imshow(matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                       extent = (0, len(matrix), 0, len(matrix)), vmax = vmax, origin = 'lower')
        cxlim = ax1.get_xlim()
        cylim = ax1.get_ylim()
        ## Ticks and Labels
        ticks = list(np.linspace(0 , len(matrix) , 5).astype(float))
        pos = [((startHiC + t) * R) for t in ticks]
        labels = [properU(p) for p in pos[:5]]
        ax1.set_xticks(ticks)
        ax1.set_xticklabels(labels)
        ax1.set_yticks(ticks)
        ax1.set_yticklabels(labels, rotation = 'horizontal')
        ax1.set_xlabel('chr'+g , fontsize = 30 )
        
        ## loop data
        # mask = (loops['start'] >= startHiC * R) & (loops['end'] < endHiC * R)
        # extract = loops[mask]
        # for p in extract:
        #     x = p['start'] // R - startHiC
        #     y = p['end'] // R - startHiC
        #     ax1.scatter(y + 0.5 ,x + 0.5 , color = '', edgecolors = 'black', s = 200)
            
        ax1.set_xlim(cxlim)
        ax1.set_ylim(cylim)                    
        #caxis_H(ax1)
        ## Colorbar
        ax2 = fig.add_axes([Left + width + 0.015 , HB , 0.035 , 0.1])
        fig.colorbar(sc,cax = ax2, orientation='vertical' , ticks = [matrix.min() , vmax])

        ##gene Tracks
        #ax3 = fig.add_axes([Left , HB + HH , width , 0.04])
        #caxis_S_horizontal(ax3, 'black')
        
        if c == 'B6_M' or c== 'B6_P' or c== 'Cast_M' or c== 'Cast_P':
            chro = 'chr' + g
        else:
            chro = g 
        ##RNA Tracks
        location = [Left , HB + HH + 0.03, width , 0.075]
        ax4 = fig.add_axes(location)
        bigwig_10bp_Plot(RNA_Data[c] , chro , i[1] , i[2] , fig , location , 'fuchsia' , 'RNA')
        caxis_H(ax4)
        # ##Chip-input Tracks
        # location = [Left , HB + HH + 0.03 + 0.075, width , 0.075]
        # ax5 = fig.add_axes(location)
        # bigwig_10bp_Plot(Chip_Data_input[c] , chro , i[1] , i[2] , fig , location , 'green' , 'Chip-input')
        # ##Chip Tracks
        location = [Left , HB + HH + 0.03 + 0.075, width , 0.075]
        ax5 = fig.add_axes(location)
        bigwig_10bp_Plot(Chip_Data[c] , chro , i[1] , i[2] , fig , location , 'mediumblue' , 'CTCF')   
        caxis_H(ax5)
        # ##DI Tracks
        # ax6 = fig.add_axes([Left , HB + HH + 0.05 + 0.075 + 0.075, width , 0.075])
        # DI_index , DI = UpdateDI(DI)
        # ax6.fill_between(np.arange(DI.size) , DI , where = DI >= 0 , facecolor = '#E47833' , edgecolor = 'none' )
        # ax6.fill_between(np.arange(DI.size) , DI , where = DI <= 0 , facecolor = '#7093DB' , edgecolor = 'none' )
        # ax6.set_ylabel('DI',fontsize = 15,rotation = 'horizontal',labelpad = 50)
        # ax6.set_xlim((cxlim[0] , cxlim[1] - 1))
        # caxis_S_horizontal(ax6, 'black')
        # ##PCA Tracks
        # PCA_index , PCA = UpdateDI(PCA)
        # ax7 = fig.add_axes([Left, HB + width + 0.05 + 0.075 + 0.075 + 0.075, width , 0.075])
        # ax7.fill_between(PCA_index , PCA , where = PCA >= 0 , facecolor = '#E47833' , edgecolor = 'none' )
        # ax7.fill_between(PCA_index , PCA , where = PCA <= 0 , facecolor = '#7093DB' , edgecolor = 'none' )
        # ax7.set_xlim(0 , PCA_index.max())
        # ax7.set_ylabel('PC1',fontsize = 15,rotation = 'horizontal',labelpad = 50)
        # caxis_S_horizontal(ax7, 'black')

        tmp[c] = fig
       


    B6_M_fig = tmp['B6_M']
    B6_P_fig = tmp['B6_P']
    Cast_M_fig = tmp['Cast_M']
    Cast_P_fig = tmp['Cast_P']
   
    
    figs = [B6_M_fig,B6_P_fig,Cast_M_fig,Cast_P_fig]
    standard_axes_lim(3,*figs)
    # standard_axes_lim(4,*figs)
    # standard_axes_lim(5,*figs)
    

    pp1.savefig(B6_M_fig)
    pp2.savefig(B6_P_fig)
    pp3.savefig(Cast_M_fig)
    pp4.savefig(Cast_P_fig)
   
    plt.close(B6_M_fig)
    plt.close(B6_P_fig)
    plt.close(Cast_M_fig)
    plt.close(Cast_P_fig)
   

pp1.close()
pp2.close()
pp3.close()
pp4.close()

# %%
