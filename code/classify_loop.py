import pandas as pd 


#cut之后的
loop = pd.read_csv('/public/home/shidetong/projects/sdt/hic-seq/haplotype/loop_number/loop_B6_Cast_hap/allloop/loop/adjust/adjustloop.csv',sep = '\t')
#adjust取小于0.1
"""
因为adjust取小于0.1的值，loop均为strain类型的，parent另外在找条件
"""
adjust = loop[(loop['B6_adjustP']<=0.1)&(loop['Cast_adjustP']<=0.1)]


B6 = adjust[(adjust['B6M_IF']>adjust['B6P_IF'])&(adjust['CastM_IF']<adjust['CastP_IF'])]
Cast = adjust[(adjust['B6M_IF']<adjust['B6P_IF'])&(adjust['CastM_IF']>adjust['CastP_IF'])]

#parent
maternal = loop[(loop['B6M_IF']>loop['B6P_IF'])&(loop['CastM_IF']>loop['CastP_IF'])]
paternal = loop[(loop['B6M_IF']<loop['B6P_IF'])&(loop['CastM_IF']<loop['CastP_IF'])]


#pvalue<0.05
pvalue = loop[(loop['B6_pvalue']<=0.05)&(loop['Cast_pvalue']<=0.05)]

#pvalue<0.1
pvalue = loop[(loop['B6_pvalue']<=0.1)&(loop['Cast_pvalue']<=0.1)]


##***----------------------------***##

##可视化之后筛选出的染色质环
#找出空白行
 



##------------------------------##
def Classify_loop(file):
    loop = pd.read_csv(file,sep ='\t')
    maternal = loop[(loop['B6M_IF']>loop['B6P_IF'])&(loop['CastM_IF']>loop['CastP_IF'])]
    paternal = loop[(loop['B6M_IF']<loop['B6P_IF'])&(loop['CastM_IF']<loop['CastP_IF'])]
    B6 = loop[(loop['B6M_IF']>loop['B6P_IF'])&(loop['CastM_IF']<loop['CastP_IF'])]
    Cast = loop[(loop['B6M_IF']<loop['B6P_IF'])&(loop['CastM_IF']>loop['CastP_IF'])]
    return B6, Cast, maternal, paternal


    




