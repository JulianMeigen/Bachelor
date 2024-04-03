import pybedtools
import numpy as np

# path_01 = "/sybig/home/jme/Bachelorarbeit/Collect_Data/Test_BED_files/ARNT.bed"
# path_02 = "Collect_Data/Test_BED_files/ASCL1.bed"

# a = pybedtools.BedTool(path_01)
# b = pybedtools.BedTool(path_02)

# a_and_b = a.intersect(b)

#a.head()
#print()
#b.head()
#print()
#a_and_b.head()


####

raw_UniBind_path = "/sybig/projects/GeneRegulation/data/jme/Bachelorarbeit/raw_data/hg38_UniBind_allTFBSs.bed"

raw_UniBind = pybedtools.BedTool(raw_UniBind_path)

# Filter the raw data and only save chr1-22 and chrX and chrY in BED file

chr_ordered = np.array(["chr1", "chr2", "chr3", "chr4","chr5","chr6","chr7","chr8", "chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19", "chr20","chr21","chr22","chrX","chrY"])
output = "/sybig/projects/GeneRegulation/data/jme/Bachelorarbeit/raw_data/UniBind_allTFBS_filtered.bed"

data_unibind = raw_UniBind.filter(lambda x: np.isin(x.chrom, chr_ordered)).saveas(output)

print("Output erstellt.\n")
print("Der Output hat anstatt 97492844 Einträgen nur:\n")
print(data_unibind.count())


