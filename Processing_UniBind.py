"""
This Script processes the BED File from the UniBind robust Collection, which contains all TFBS (in the human Genome).
The Data looks like this:
chr1	10143	10157	EXP040021_K562--myelogenous-leukemia-_NR2F6_MA0677.1	0	-	10143	10157	198,244,48
chr1	10143	10158	EXP040021_K562--myelogenous-leukemia-_NR2F6_MA0728.1	0	-	10143	10158	198,244,48
chr1	10175	10192	EXP000883_T47D--invasive-ductal-carcinoma-_ESR1_MA0112.3	0	-	10175	10192	209,252,66

Important is that the "name" column of the BED file contain the information about ”ChipSeq-ID_Zelllinie_TF-name_JASPAR-ID”, 
because this information will be necasarry to process the file based on the TF-Name.

The Processing of the file has these steps:

Step 1: Filtering the genomic location by a second BED File, which contains the important Regions (e.g. Promotor regions etc.)

Step 2: (Optional) Filtering the genomic location based on the chromosomes. (default= chr1-22 and X,Y)

Step 3: Refining the Data and changing columns: Information from "name" column will be extracted and they'll replace the unnecesarry columns.
        (If possible with .each(func) and not with Dataframe. In case the Bedfile will be too big for beeing saved in a DataFrame.)

Step 4: Merge repetitive Entrys


"""


import pybedtools
import numpy as np
import pandas as pd
import argparse


def reduce_bed_by_intersect_wa(bed_file_path, intersect_file_path):
    """
    This function loads two BED files as BedTool Objects.
    The genomic regions of the bed_file, who intersect with a genomic region in the intersect_file, 
    will be filtered and returned as new BedTool object.
    """
    # Load input Bed file into BedTools object
    raw_data_bed = pybedtools.BedTool(bed_file_path)

    intersect_data = pybedtools.BedTool(intersect_file_path)

    reduced_data = raw_data_bed.intersect(intersect_data, wa=True)
    
    return reduced_data

def filtering_chromosomes(BedTool, chr_lst=["chr1", "chr2", "chr3", "chr4","chr5","chr6","chr7","chr8", "chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19", "chr20","chr21","chr22","chrX","chrY"]):
    """
    This function is filtering a BedTool Object by its chromosomes. 
    As default the filtered object will contain chr1-22 and chrX, chrY, but any list of chromosome names can be handed over.
    """
    filtered_data = BedTool.filter(lambda x: np.isin(x.chrom, np.array(chr_lst)))

    return filtered_data

def refine_name_column(BedTool, )

def merge_bed_considering_name(bedtool_data, )

