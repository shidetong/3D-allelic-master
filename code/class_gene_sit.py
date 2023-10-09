import pandas as pd 


classify = pd.read_csv('/public/home/shidetong/projects/lxx/RNA-seq/haplotype/dif/classify/Castbias.csv',sep = '\t')

genesit = pd.read_csv('/public/home/shidetong/projects/lxx/RNA-seq/haplotype/dif/all_gene_sit.csv',sep = '\t')

classify['symbol'] = classify['symbol'].str.upper()
endfil  =pd.merge(classify,genesit,on=['symbol'],how = 'inner')

a = endfil.loc[:,['chr','start','end','symbol','log2FC_BCF1','pvalue_BCF1','log2FC_CBF1','pvalue_CBF1']]
a.to_csv('/public/home/shidetong/projects/lxx/RNA-seq/haplotype/dif/classify/gene_sit/Castbias.csv',sep = '\t',index = None)