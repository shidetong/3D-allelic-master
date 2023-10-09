import cooler
B6_T = cooler.Cooler('Merged_Traditional_Multi.cool::40000')
B6_T.matrix(balance = False).fetch('1') #Get the raw chromosome 1 Matrix
B6_T.matrix(balance = True).fetch('1') #Get the balanced chromosome 1 Matrix
#Haplotype-resolved Matrix
B6_Haplotype = cooler.Cooler('Merged_Imputated_Haplotype_Multi.cool::40000')
B6_Haplotype.matrix(balance = False).fetch('M1') #Get chromosome 1 Maternal Matrix
B6_Haplotype.matrix(balance = False).fetch('P1') #Get chromosome 1 Paternal Matrix


from HiCHap.StructureFind import StructureFind

#========== === Compartment==============
## For traditional Hi-C
GM_T_PC = StructureFind(cooler_fil = 'Merged_Traditional_Multi.cool', Res = 500000, Allelic = False, genomeSize_fil = 'genomeSize')
GM_T_PC.run_Compartment(OutPath = 'Traditional_PC', plot = True, MS = 'IF', SA = False)

## For haplotype-resolved Hi-C
GM_M_PC = StructureFind(cooler_fil = 'Merged_Imputated_Haplotype_Multi.cool', Res = 500000, Allelic = 'Maternal', genomeSize_fil = 'genomeSize')
GM_M_PC.run_Compartment(OutPath = 'Maternal_PC', plot = True, MS = 'IF', SA = False, Tranditional_PC_file = 'Traditional_PC/Traditional_PC_Compartment_500K.txt')

GM_P_PC = StructureFind(cooler_fil = 'Merged_Imputated_Haplotype_Multi.cool', Res = 500000, Allelic = 'Paternal', genomeSize_fil = 'genomeSize')
GM_P_PC.run_Compartment(OutPath = 'Paternal_PC', plot = True, MS = 'IF', SA = False, Tranditional_PC_file = 'Traditional_PC/Traditional_PC_Compartment_500K.txt')


#============= TADs calling=============
## For traditional Hi-C
GM_tads_T = StructureFind(cooler_fil = 'Merged_Traditional_Multi.cool', Res = 40000, Allelic = False , genomeSize_fil = 'genomeSize')
GM_tads_T.run_TADs(OutPath = 'Traditional_TADs', plot = True)

## For haplotype-resolved Hi-C
GM_tads_M = StructureFind(cooler_fil = 'Merged_Imputated_Haplotype_Multi.cool', Res = 40000, Allelic ='Maternal', genomeSize_fil = 'genomeSize')
GM_tads_M.run_TADs(OutPath = 'Maternal_TADs', plot = True)

GM_tads_P = StructureFind(cooler_fil = 'Merged_Imputated_Haplotype_Multi.cool', Res = 40000, Allelic ='Paternal', genomeSize_fil = 'genomeSize')
GM_tads_P.run_TADs(OutPath = 'Paternal_TADs', plot = True)


#============= Loops calling=============
## For traditonal Hi-C
GM_Loop_T = StructureFind(cooler_fil = 'Merged_Traditional_Multi.cool', Res = 40000, Allelic = False , genomeSize_fil = 'genomeSize')
GM_Loop_T.run_Loops(OutPath = 'Traditional_Loops', plot = True)

## For haplotype-resolved Hi-C
GM_Loop_M = StructureFind(cooler_fil = 'Merged_Imputated_Haplotype_Multi.cool', Res = 40000, Allelic = 'Maternal', GapFile = 'Merged_Imputated_Gap.npz' , genomeSize_fil = 'genomeSize')
GM_Loop_M.run_Loops(OutPath = 'Maternal_Loops', plot = True)

GM_Loop_P = StructureFind(cooler_fil = 'Merged_Imputated_Haplotype_Multi.cool', Res = 40000, Allelic = 'Paternal', GapFile = 'Merged_Imputated_Gap.npz' , genomeSize_fil = 'genomeSize')
GM_Loop_P.run_Loops(OutPath = 'Paternal_Loops', plot = True)
