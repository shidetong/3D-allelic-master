import pandas as pd 
B6_1 = pd.read_csv('/public/home/shidetong/projects/lxx/RNA-seq/20210324/ballgown/B6-RNA-1A/B6-RNA-1A.tab',sep = '\t')
B6_2 = pd.read_csv('/public/home/shidetong/projects/lxx/RNA-seq/20210324/ballgown/B6-RNA-3A/B6-RNA-3A.tab',sep = '\t')
Cast_1 = pd.read_csv('/public/home/shidetong/projects/lxx/RNA-seq/20210324/ballgown/Cast-RNA-1A/Cast-RNA-1A.tab',sep = '\t')
Cast_2 = pd.read_csv('/public/home/shidetong/projects/lxx/RNA-seq/20210324/ballgown/Cast-RNA-3A/Cast-RNA-3A.tab',sep = '\t')
B6_1 = B6_1.iloc[:,[1,7]]
B6_1 = B6_1.drop_duplicates(subset=['Gene Name'],keep = False)
B6_2 = B6_2.iloc[:,[1,7]]
B6_2 = B6_2.drop_duplicates(subset=['Gene Name'],keep = False)
Cast_1 = Cast_1.iloc[:,[1,7]]
Cast_1 = Cast_1.drop_duplicates(subset=['Gene Name'],keep = False)
Cast_2 = Cast_2.iloc[:,[1,7]]
Cast_2 = Cast_2.drop_duplicates(subset=['Gene Name'],keep = False)
B6 = pd.merge(B6_1,B6_2,on = ['Gene Name'])
Cast = pd.merge(Cast_1,Cast_2,on=['Gene Name'])
all_fpkm =pd.merge(B6,Cast,on = ['Gene Name']) 
all_fpkm.columns = ['genename','B6_1','B6_2','Cast_1','Cast_2']
all_fpkm.to_csv('/public/home/shidetong/projects/lxx/RNA-seq/20210324/ballgown/fpkm.csv',sep = '\t',index = None)

