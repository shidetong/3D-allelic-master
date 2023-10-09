from itertools import  islice
from scipy import optimize
import numpy as np
import pandas as pd 
import pysam, sys, os


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
            strand = a[2]
            start = a[3]
            end = a[4]
            gtf_1.append((gene_id , gene_name , chro , strand , start , end))
    gtf = np.array(gtf_1 , dtype = gtf_type)
    return gtf


a = Load_gtf('/public/home/shidetong/wedata/test/genenum/gencode.v37lift37.annotation.gtf')

r = pd.DataFrame(a)
b = r.drop_duplicates(subset=['gene_name'],keep='first')
b.columns = ['gene_id' , 'symbol' , 'chr' , 'strand' , 'start' , 'end']
b.to_csv('/public/home/shidetong/projects/lxx/RNA-seq/haplotype/dif/all_gene_sit.csv',sep = '\t',index = None)