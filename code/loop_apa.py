%matplotlib inline
from __future__ import division
import numpy as np
#from tadlib.calfea.analyze import getmatrix
import matplotlib
# Use a non-interactive backend
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import os
from scipy.special import ndtr
import math
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

def apa_submatrix(M, pos, w=5):
    
    Len = M.shape[0]

    apa = []
    for i, j in pos:
        if (i-w>=0) and (i+w+1<=Len) and (j-w>=0) and (j+w+1<=Len):
            tmp = M[i-w:i+w+1, j-w:j+w+1]
            if tmp.mean()==0:
                continue
            mask = np.isnan(tmp)
            if mask.sum() > 0:
                continue
            tmp = tmp / tmp.mean()
            apa.append(tmp)
    
    return apa

def apa_analysis(apa, w=5, cw=3):
    
    avg = apa.mean(axis=0)
    lowerpart = avg[-cw:,:cw]
    upperpart = avg[:cw,-cw:]
    maxi = upperpart.mean() * 5
    ## APA score
    score = avg[w,w] / lowerpart.mean()
    ## z-score
    z = (avg[w,w] - lowerpart.mean()) / lowerpart.std()
    p = 1 - ndtr(z)
    
    return avg, score, z, p, maxi
    

def Get_nan_zero_Matrix(HiC_Lib):
    '''
    '''
    Lib_new = {}
    for g in HiC_Lib:
        tmp = HiC_Lib[g]
        tmp[    (tmp)] = 0
        Lib_new[g] = tmp
    return Lib_new
                
    
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
   
    
def Correct_VC(X, alpha):
    x = np.array(X,float)
    s1 = np.sum(x, axis = 1)
    s1 = s1 ** alpha
#    s /= np.mean(s[s!=0])
    s1[s1 == 0] = 1
    s2 = np.sum(x, axis = 0)
    s2 = s2 ** alpha
#    s2 /= np.mean(s2[s2 != 0])
    s2[s2 == 0] = 1
    return x / (s2[None, :] * s1[:, None])


def Normal_VC_Correct(NPZ):
    """
    """
    Raw_Lib = np.load(NPZ)
    Nor_Lib = {}
    for c in Raw_Lib.keys():
        Raw_M = Raw_Lib[c]
        Nor_M = Correct_VC(Raw_M, 2/3)
        
        #Recall
        Re_factor = Raw_M.mean() / Nor_M.mean()
        Nor_M = Nor_M * Re_factor
        Nor_Lib[c] = Nor_M
    
    return Nor_Lib


chros = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X']
res = 40000
initial_distance = 2 * math.sqrt(2)

data_type = np.dtype({'names':['chr' , 'start' , 'end'] ,
                      'formats':['S8' , np.int , np.int]})
cells = ['B6_M','B6_P','Cast_M','Cast_P']


#HiC Data Process
B6_M_Lib = np.load('/public/home/shidetong/projects/sdt/hic-seq/haplotype/B6/B6_Matrix_40K/Cooler_40K/cooltonpz/Merged_Imputated_Haplotype_Multi_M.npz' , allow_pickle=True)
B6_M_Lib = B6_M_Lib['arr_0'][()]
B6_P_Lib = np.load('/public/home/shidetong/projects/sdt/hic-seq/haplotype/B6/B6_Matrix_40K/Cooler_40K/cooltonpz/Merged_Imputated_Haplotype_Multi_P.npz' , allow_pickle=True)
B6_P_Lib = B6_P_Lib['arr_0'][()]
Cast_M_Lib = np.load('/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_Matrix_40K/Cooler/cooltonpz/Merged_Imputated_Haplotype_Multi_M.npz' , allow_pickle=True)
Cast_M_Lib = Cast_M_Lib['arr_0'][()]
Cast_P_Lib = np.load('/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_Matrix_40K/Cooler/cooltonpz/Merged_Imputated_Haplotype_Multi_P.npz' , allow_pickle=True)
Cast_P_Lib = Cast_P_Lib['arr_0'][()]
HiC_Data = {'B6_M':B6_M_Lib,'B6_P':B6_P_Lib,'Cast_M':Cast_M_Lib,'Cast_P':Cast_P_Lib}




B6_M_speci_Loops = np.loadtxt('/public/home/shidetong/projects/sdt/hic-seq/haplotype/B6/new_B6_Matrix/Cooler/Maternal_Loops/Cluster_Maternal_Loops_Loops_40K.txt' , 
                             dtype = data_type ,skiprows=1 , usecols = (0 , 1 , 2))
B6_P_speci_Loops = np.loadtxt('/public/home/shidetong/projects/sdt/hic-seq/haplotype/B6/new_B6_Matrix/Cooler/Paternal_Loops/Cluster_Paternal_Loops_Loops_40K.txt' , 
                             dtype = data_type ,skiprows=1 , usecols = (0 , 1 , 2))
Cast_M_speci_Loops = np.loadtxt('/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_Matrix/Cooler/Maternal_Loops/Cluster_Maternal_Loops_Loops_40K.txt' , 
                             dtype = data_type ,skiprows=1 , usecols = (0 , 1 , 2))
Cast_P_speci_Loops = np.loadtxt('/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_Matrix/Cooler/Paternal_Loops/Cluster_Paternal_Loops_Loops_40K.txt' , 
                             dtype = data_type ,skiprows=1 , usecols = (0 , 1 , 2))
Loops = {'B6_M':B6_M_speci_Loops ,
        'B6_P':B6_P_speci_Loops,
        'Cast_M':Cast_M_speci_Loops,
        'Cast_P':Cast_P_speci_Loops}



# Left = 0.2 ; HB = 0.2 ; width = 0.6 ; HH = 0.6                     

for cl in ['Cast_P']:
    for c in ['Cast_P']:
        P_Loops = Loops[cl]
        apa = []                    
        for g in chros:
            peaks = P_Loops[P_Loops['chr'] == g]
            pos = []
            M = HiC_Data[c][g]
            for p in peaks:
                x, y = p[1], p[2]
                if abs(y-x) < 15 * res:
                    continue
                s_l = range(p[1]//res, int(np.ceil((p[1]+40000)/float(res))))
                e_l = range(p[2]//res, int(np.ceil((p[2]+40000)/float(res))))
                si, ei = None, None
                for st in s_l:
                    for et in e_l:
                        if (st < M.shape[0]) and (et < M.shape[0]):
                            if si is None:
                                si, ei = st, et
                            else:
                                if M[st,et]>M[si,ei]:
                                    si,ei = st,et
                if not si is None:
                    if si < ei:
                        pos.append((si,ei))
                    else:
                        pos.append((ei,si))
            tmp = apa_submatrix(M, pos)
            apa.extend(tmp)
        apa = np.r_[apa]
        avg,score,z,p_value,maxi = apa_analysis(apa)
        
        fig = plt.figure(figsize = (8,8))
        plt.tick_params(axis='both', bottom=False, top=False, left=False, right=False,
                        labelbottom=False, labeltop=False, labelleft=False, labelright=False)
        ax = fig.add_axes([Left  , HB , width , HH])
        sc = ax.imshow(avg, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                       extent = (0, len(avg), 0, len(avg)), vmin = 0 , vmax = 3, origin = 'lower')
        ax.set_xticks([0 , 5 , 10])
        ax.set_xticklabels(['-100K' , 'Loop'  , '100K'])
        ax.set_yticks([0 , 5 , 10])
        ax.set_yticklabels(['-100K' , 'Loop'  , '100K'])
        ax.set_xlabel('Cast_P', fontsize = 15)
        ax.set_title('APA score = {0:.3g}, p-value = {1:.3g} , loop_num: {2:.0f}'.format(score, p_value, len(P_Loops)))
        ax = fig.add_axes([Left + width + 0.03 , HB , 0.035 , 0.1])
        cbar = fig.colorbar(sc,cax = ax, orientation='vertical')
        cbar.set_ticks([0 , 3])

