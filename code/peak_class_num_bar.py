import matplotlib.pyplot as plt 
import seaborn as sns  



peak = [879,799,4,4,2577,1315,1214,1565,1652]
loop = [83,27]
loop = [174,206,1,6,514,119,61,45,116]
gene = [40,94]
# compartment = [11,20,2916,961,985,115,137,74,59]
compartment = [11,20,400,261,285,115,137,74,59]

###peak
plt.figure(figsize=(7,5))
plt.bar(range(len(peak)),peak,color=['#FB8402','#219EBC','#90C9E6'])
plt.tick_params(labelsize=18)
x = [0,1,2,3,4,5,6,7,8]
lab = ['1','2','3','4','5','6','7','8','9']
plt.xticks(x,lab)
plt.savefig('/public/home/shidetong/projects/sdt/chip-seq/haplotype/plot/peakclassify.png',dpi=300)

###loop
plt.figure(figsize=(7,5))
# plt.bar(range(len(loop)),loop,color=['#FB8402','#219EBC','#90C9E6'],width = 0.5)
plt.bar(range(len(loop)),loop,color=['#FB8402','#219EBC','#90C9E6'])
plt.tick_params(labelsize=18)
# x = [0,1]
# lab = ['1','2']
x = [0,1,2,3,4,5,6,7,8]
lab = ['1','2','3','4','5','6','7','8','9']
plt.xticks(x,lab)
plt.savefig('/public/home/shidetong/projects/sdt/hic-seq/haplotype/plot/loopnum_new.png',dpi=300)


#RNA
plt.figure(figsize=(3,5))
plt.bar(range(len(gene)),gene,color=['#FF0000','#0000FF','#90C9E6'],width = 0.5)
plt.tick_params(labelsize=18)
x = [0,1]
lab = ['1','2']
plt.xticks(x,lab)
plt.savefig('/public/home/shidetong/projects/lxx/RNA-seq/haplotype/plot/difgenenum_red.png',dpi=300)

#compartment
plt.figure(figsize=(7,5))
plt.bar(range(len(compartment)),compartment,color=['#FB8402','#219EBC','#90C9E6'])
plt.tick_params(labelsize=18)
x = [0,1,2,3,4,5,6,7,8]
lab = ['1','2','3','4','5','6','7','8','9']
plt.xticks(x,lab)
plt.savefig('/public/home/shidetong/projects/sdt/hic-seq/haplotype/PC_A2B/PC_classify/compartment_classify.png',dpi=300)

