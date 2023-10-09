import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn as sns

fl = pd.read_csv('/public/home/shidetong/projects/lxx/RNA-seq/20210324/ballgown/count/test/gene_count_matrix.csv',sep = ',')
fl = fl.loc[:,['B6-RNA-1A','B6-RNA-3A','Cast-RNA-1A','Cast-RNA-3A']]

# fl = pd.read_csv('/public/home/shidetong/projects/yf/RNA-seq/20200924/count/count_new.csv',sep = '\t')
# fl = fl.iloc[:,1:10]
 
cor = fl.corr()
 
# ax = sns.heatmap(cor,vmax=0.79,annot=True)
# ax.set_title("RNA pearson correlation",size = 20)
# plt.savefig('/public/home/shidetong/projects/yf/RNA-seq/20200924/count/RNA_cor.png')


ax = sns.heatmap(cor,vmax=1,vmin = 0.95,annot=False,cmap='GnBu',square=True)
ax.set_title("RNA pearson correlation",size = 20)
plt.savefig('/public/home/shidetong/projects/lxx/RNA-seq/20210324/plot/RNA_cor_new.pdf')




