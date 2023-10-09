#单倍体
import numpy as np
import cooler
from scipy import sparse
import sys
import os

def cooler2npz(DataFile , OutFile , res , maternal):
    Matrix = cooler.Cooler(DataFile + '::' + str(res))
    
    tmp = {}
    for g in chroms:
        if g[0] == maternal: 
            matrix = Matrix.matrix(balance=False).fetch(g)
            tmp[g.lstrip(maternal)] = matrix

    np.savez_compressed(OutFile,tmp)
    
chromsM = ['M' + i for i in map(str,range(1,20))] + ['MX']
chromsP = ['P' + i for i in map(str,range(1,20))] + ['PX']
chroms = chromsM + chromsP

cooler2npz('/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_Matrix_40K/Cooler/cooltonpz/Merged_Imputated_Haplotype_Multi.cool' , '/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_Matrix_40K/Cooler/cooltonpz/Merged_Imputated_Haplotype_Multi_P.npz' ,  40000 , maternal='P')
cooler2npz('/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_Matrix_40K/Cooler/cooltonpz/Merged_Imputated_Haplotype_Multi.cool' , '/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_Matrix_40K/Cooler/cooltonpz/Merged_Imputated_Haplotype_Multi_M.npz' ,  40000 , maternal='M')


#二倍体
import numpy as np
import cooler
from scipy import sparse
import sys
import os
def cooler2npz_diploid(DataFile , OutFile , res , cor):
    Matrix = cooler.Cooler(DataFile + '::' + str(res))
    chroms = [i for i in map(str,range(1,20))] + ['X']
    # chroms = ['7']
    
    tmp = {}
    if cor == True:
        for g in chroms:
            matrix = Matrix.matrix(balance=True).fetch(g)
            tmp[g] = matrix
    else:
        for g in chroms:
            matrix = Matrix.matrix(balance=False).fetch(g)
            tmp[g] = matrix

    np.savez_compressed(OutFile, tmp)

cooler2npz_diploid('/public/home/shidetong/projects/sdt/hic-seq/haplotype/B6/B6_Matrix_40K/Cooler_40K/cooltonpz/Merged_Traditional_Multi.cool' , '/public/home/shidetong/projects/sdt/hic-seq/haplotype/B6/B6_Matrix_40K/Cooler_40K/cooltonpz/Merged_Traditional_Multi.npz' , 40000 , True)
cooler2npz_diploid('/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_Matrix_40K/Cooler/cooltonpz/Merged_Traditional_Multi.cool' , '/public/home/shidetong/projects/sdt/hic-seq/haplotype/Cast/Cast_Matrix_40K/Cooler/cooltonpz/Merged_Traditional_Multi.npz' , 40000 , True)




Matrix()