from numpy import histogram, histogram_bin_edges
import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn as sns
import numpy as np




H6C7 = pd.read_csv('/public/home/shidetong/projects/sdt/hic-seq/HiC_workspace/B6_batch3_HiC_3_2/HiCHap_workspace/B6_batch3_HiC_3_2/Matrix/Cooler/Traditional_PC/Traditional_PC_Compartment_2M.txt',sep = '\t',header = None) 
BXPC = pd.read_csv('/public/home/shidetong/projects/sdt/hic-seq/HiC_workspace/B6_batch3_HiC_5_1/Matrix/Cooler/Traditional_PC/Traditional_PC_Compartment_2M.txt',sep = '\t',header = None)
PANC = pd.read_csv('/public/home/shidetong/projects/sdt/hic-seq/HiC_workspace/Cast_batch2_HiC_2_1/Matrix/Cooler/Traditional_PC/Traditional_PC_Compartment_2M.txt',sep = '\t',header = None)
P = pd.read_csv('/public/home/shidetong/projects/sdt/hic-seq/HiC_workspace/Cast_batch2_HiC_4_3/Matrix/Cooler/Traditional_PC/Traditional_PC_Compartment_2M.txt',sep = '\t',header = None)


H6C7 = pd.read_csv('/public/home/shidetong/projects/sdt/hic-seq/haplotype/B6/B6_Matrix_40K/Cooler_40K/Maternal_PC/Maternal_PC_Compartment_500K.txt',sep = '\t',header = None) 
BXPC = pd.read_csv('/public/home/shidetong/projects/sdt/hic-seq/haplotype/B6/B6_Matrix_40K/Cooler_40K/Paternal_PC/Paternal_PC_Compartment_500K.txt',sep = '\t',header = None)
PANC = pd.read_csv('/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_Matrix_40K/Cooler/Maternal_PC/Maternal_PC_Compartment_500K.txt',sep = '\t',header = None)
P = pd.read_csv('/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_Matrix_40K/Cooler/Paternal_PC/Paternal_PC_Compartment_500K.txt',sep = '\t',header = None)


pc = pd.concat([H6C7,BXPC,PANC,P],axis=1)

pc.columns = ['a','b','c','d','e','f','g','h']

IF = pc.loc[:,['b','d','f','h']]
IF.columns = ['H6C7','BXPC','PANC','P']

cor = IF.corr()
cor = pd.read_csv('/public/home/shidetong/projects/yf/ATAC/bw/merge_plot/pearsonCorr_readCounts.tab')
# ax = sns.heatmap(cor,vmax=0.4,annot=True)
ax = sns.heatmap(cor,vmax=0.9,vmin = 0.85,annot=False,cmap='GnBu',square=True)
# sns.clustermap(cor,vmax=0.7,annot=True)
# sns.set_style('whitegrid', {'font.sans-serif': ['simhei','FangSong']})
ax.set_title("Hi-C pearson correlation",size = 20)

plt.savefig('/public/home/shidetong/projects/sdt/hic-seq/HiC_workspace/plot/hic_cor.pdf')



H6C7 = np.load('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/H6C7/npz/20200518H6C7_40K.npz')
BXPC = np.load('/public/home/shidetong/projects/yf/hic/HiCHap_workspace/Matrix/Bxpc/npz/20200318Bxpc_40K.npz')

cor=np.vstack(H6C7,BXPC)