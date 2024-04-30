"""
The task of the script is to create a csv file for each pair, which contains the most important information for this pair and which can be used later for further processing. 

The input file is a BED file with the following format:

Promotor Genomic Location   	geneID 	    TFBS    	TFBS distance to TSS 	    TFBS strand 	*
chr1 17436 17636 - 	            ENSGXXX 	MYCN 	        74 86 	                   - 	        ***

And each output file should have the following format:
chr,geneID,TFBS_close,TFBS_dist,closest_TSS_distance,furthest_TSS_distance,orientation_of_pair,distance_between_pair,tfbs_count,tfbs_count_unique
chr7,ENSG00000003147.17,CTCF,REST,-10.0,103.0,convergent,73.0,47.0,35.0

Each pair can thus be successfully described.



"""

import pybedtools
import pandas as pd
import numpy as np
import argparse

from timebudget import timebudget
from multiprocessing import Pool

@timebudget
def func_for_tfbs_subset(BedTool_Interval, tfbs_lst):
    """
    Function that returns TRUE, if all of the tfbs in tfbs_lst are in Interval.fields[6]
    It also checks if tfbs in tfbs_lst occur more often then tfbs in Interval.fields[6]. 
    Therefore if tfbs_lst contains the same tfbs twice, it will only return True if Interval.fields[6] also contains the tfbs at least twice. 
    """

    # Processing tfbs_lst --> unique with counts
    input_tfbs_unique, input_tfbs_counts = np.unique(np.array(tfbs_lst), return_counts=True)
    input_tfbs_dict = dict(zip(input_tfbs_unique, input_tfbs_counts))
    
    # Processing tfbs in interval --> split tfbs, distance, strand into lists
    tfbs_arr = np.array(BedTool_Interval.fields[6].split(","))

    # get unique tfbs and counts and save in dict
    tfbs_unique, tfbs_counts = np.unique(tfbs_arr, return_counts=True)
    tfbs_dict = dict(zip(tfbs_unique, tfbs_counts))

    # Compare unique values and counts (equal or more (?))

    # Checks if all unique TFBS in tfbs_lst occur in tfbs_unique.
    if np.all(np.isin(input_tfbs_unique, tfbs_unique)):
        
        # Checks how often tfbs occur in input list, compared to interval tfbs. If tfbs count in input list is higher then in interval it will return False.
        bool_lst = []
        for tfbs, input_count in input_tfbs_dict.items():
            tfbs_count = tfbs_dict[tfbs]
            if input_count > tfbs_count:
                bool_lst.append(False)
            else:
                bool_lst.append(True)
        if np.all(bool_lst):
            return True
        else:
            return False
        
        
    else:
        return False
@timebudget   
def get_tfbs_subset(BedTool, tfbs_lst):
    """
    Filtering BedTool, so that only Promotors/GeneIDs remain, that contain every TFBS in tfbs_lst at least once.
    """
    return BedTool.filter(func_for_tfbs_subset, tfbs_lst)


def main():
    data = pybedtools.BedTool("/sybig/projects/GeneRegulation/data/jme/Bachelorarbeit/data/Promotor_with_TFBS/All_GTEx_Prom_with_TFBS.bed")
    func_for_tfbs_subset(data[0], ("CTCF", "REST"))
main()