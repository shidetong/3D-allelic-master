
from __future__ import division
import math
import numpy as np
import csv , copy
# import xlrd
import pandas as pd
from itertools import islice
from sklearn.cluster import KMeans
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster import hierarchy
from scipy import cluster 
import seaborn as sns
import copy
import scipy
import scipy.cluster.hierarchy as sch
from itertools import islice  
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap





    
#代码实现：已知基因组序列的两个位置，看这两个位置之间有哪些gene

#重新定义gtf文件，
def Load_gtf(gtfil):
    gtf_type = np.dtype({'names':['gene_id' , 'gene_name' , 'chr' , 'strand' , 'start' , 'end'],
                     'formats':['U64' , 'U64' , 'U8' , 'U4' , np.int , np.int]})
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

#找到位置
def Get_interval_genes(interval , gtf ):
    '''
    '''
    genes = []
    for i in interval:
        g = i[0]
        start = i[1] 
        end = i[2]
        tmp = gtf[gtf['chr'] == "chr"+ g]
        mask = (tmp['start'] <= end) & (tmp['end'] >= start)
        overlap = tmp[mask]
        if overlap.size != 0:
            for j in overlap:
                genes.append(j)
            
            
    genes = np.array(genes , dtype = gtf.dtype)
    return genes



# gtf = Load_gtf('E:\\gtf\\gencode.v37lift37.annotation.gtf')
gtf = Load_gtf('/public/home/shidetong/wedata/test/genenum/gencode.v37lift37.annotation.gtf')



#read file
f = open ('/public/home/shidetong/projects/yf/hic/test/hic_breakfinder/PANC/panc.assemblies.txt') 

file = f.readlines()
#len(file)



for i in range(len(file)):
    line = file[i]

    chr1 = line.split("\t")[1].split(",")[1]
    chr2 = line.split("\t")[1].split(",")[4]
    chr1_start = line.split("\t")[1].split(",")[2]
    chr1_end = line.split("\t")[2].split(",")[1]
    chr2_start = line.split("\t")[1].split(",")[5]
    chr2_end = line.split("\t")[3].split(",")[1]
    num = line.split("\t")[0]
    
#判断start与end大小  
#list=[1,2]

#for i in list:


    if chr1_start > chr1_end:
        chr1_start = chr1_start
        chr1_end = chr1_end
    else:
        chr1_start , chr1_end = chr1_end , chr1_start
    
    
    if chr2_start > chr2_end:
        chr2_start = chr2_start
        chr2_end = chr2_end
    else:
        chr2_start , chr2_end = chr2_end , chr2_start
    


 
    interval_1 =[(chr1 , int(chr1_end) ,int(chr1_start) )]
    genes_1 = Get_interval_genes(interval_1 , gtf)
    
    
    interval_2 =[(chr2 , int(chr2_end) ,int(chr2_start) )]
    genes_2 = Get_interval_genes(interval_2 , gtf)
    
    
    
    
    out_file =  open (f"/public/home/shidetong/projects/yf/hic/test/hic_breakfinder/PANC/gene/{num}_1_chr{chr1}_gene.txt","w")
    for i in genes_1:
        
        out_file.writelines('\t'.join([str(x) for x in i]) + '\n')
    out_file.close()
   
    
    out_file =  open (f"/public/home/shidetong/projects/yf/hic/test/hic_breakfinder/PANC/gene/{num}_2_chr{chr2}_gene.txt","w")
    for i in genes_2:
        
        out_file.writelines('\t'.join([str(x) for x in i]) + '\n')
    out_file.close()
    
    
    
    
#----------------------------------------------------------#    
    out_file =  open (f"E:\\neo\\gene_find\\bxpc\\{num}_1_chr{chr1}_gene.txt","w")
    out_file.write('"gene_id","gene_name","chr","strand","start","end"\n' )  #添加列名
    out_file.write(str(genes_1))
    out_file.close()
    
    out_file =  open (f"E:\\neo\\gene_find\\bxpc\\{num}_2_chr{chr2}_gene.txt","w")
    out_file.write('"gene_id","gene_name","chr","strand","start","end"\n' )  #添加列名
    out_file.write(str(genes_2))
    out_file.close()