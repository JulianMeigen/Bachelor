import pybedtools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


tfbs_all = pybedtools.BedTool("/sybig/projects/GeneRegulation/data/jme/Bachelorarbeit/data/tfbs/GTEx_prom_200bp_with_all_TFBS/output_data.bed")
gtex_prom = pybedtools.BedTool("/sybig/projects/GeneRegulation/data/jme/Bachelorarbeit/data/GTEx_Promotors.bed")


intersect = gtex_prom.intersect(tfbs_all, wa=True, wb=True).sort()

def unique_BedTool_field(BedTool, field_number):
    """
    function that get unique intervall names of given BedTool object. Return a list with unique names.
    """
    map_name = map(lambda x: x.fields[field_number], BedTool)
    names = np.array(list(map_name))
    names_unique = np.unique(names)
    return list(names_unique)



def merge_bed_considering_field(BedTool, field_number):
    """
    function that will use BedTool.merge on BedTool entrys with same BedTool.field[field_number]
    """

    BedTool_name_unique = unique_BedTool_field(BedTool, field_number)
    
    merged_tfbs_tmp = pybedtools.BedTool(())
    for tfbs_name in BedTool_name_unique:
        single_tfbs_merged = BedTool.filter(lambda x: x.fields[7] == tfbs_name).merge(s=True, c=[4,5,6,7,8,9,10,  14,12,13,16, 18, 14, 14], o=["distinct","distinct","distinct","distinct","distinct","distinct","distinct",  "collapse","collapse","collapse","collapse", "distinct", "count", "count_distinct"])
        merged_tfbs_tmp = merged_tfbs_tmp.cat(single_tfbs_merged, postmerge=False)
    
    merged_tfbs = merged_tfbs_tmp.sort()

    return merged_tfbs


merged_outout = merge_bed_considering_field(intersect,7)
merged_outout.saveas("/sybig/projects/GeneRegulation/data/jme/Bachelorarbeit/data/Promotor_with_TFBS/prom_with_tfbs_02.bed")



