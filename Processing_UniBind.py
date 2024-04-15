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



def refine_BedTool_intervall(BedTool_intervall):
    """
    function that processes single entry and rename columns (fields)  based on ”ChipSeq-ID_Zelllinie_TF-name_JASPAR-ID” information from UniBind column.
    """
    # split name
    split_name = BedTool_intervall.fields[3].split("_")
    chipseq_id = split_name[0]
    tissue = split_name[1]
    tf_name = split_name[2]
    jaspar_id = split_name[3]

    # rename columns (fields)
    BedTool_intervall[3] = tf_name
    BedTool_intervall[6] = chipseq_id
    BedTool_intervall[7] = tissue
    BedTool_intervall[8] = jaspar_id

    return BedTool_intervall


def refine_BedTool(BedTool):
    """
    function that processes each entry of BedTool object and rename the columns (fields) based on refine_BedTool_intervall.
    """
    refined_BedTool = BedTool.each(refine_BedTool_intervall)
    return refined_BedTool
    

def unique_BedTool_names(BedTool):
    """
    function that get unique intervall names of given BedTool object. Return a list with unique names.
    """
    map_name = map(lambda x: x.name, BedTool)
    names = np.array(list(map_name))
    names_unique = np.unique(names)
    return list(names_unique)



def merge_bed_considering_name(BedTool):
    """
    function that will use BedTool.merge on BedTool entrys with same BedTool.name
    """

    BedTool_name_unique = unique_BedTool_names(BedTool)
    
    merged_tfbs_tmp = pybedtools.BedTool(())
    for tfbs_name in BedTool_name_unique:
        single_tfbs_merged = BedTool.filter(lambda x: x.name == tfbs_name).merge(s=True, c=[4,5,6,7,8,9,1], o=["distinct","sum","distinct","distinct", "distinct", "distinct","count"])
        merged_tfbs_tmp = merged_tfbs_tmp.cat(single_tfbs_merged, postmerge=False)
    
    merged_tfbs = merged_tfbs_tmp.sort()

    return merged_tfbs



def main():
    """
    Step 1: Filtering the genomic location by a second BED File, which contains the important Regions (e.g. Promotor regions etc.)
    """
    print("\n Running Script ...\n")
    reduced_data = reduce_bed_by_intersect_wa(bed_file_path=args.filename, intersect_file_path=args.intersect).saveas()
    print("Step 1 successfully completed. Only UniBind TFBS which intersect with Genomic Regions of -fb BED file remain.\n")

    # Optional save:
    if args.all_output:
        reduced_data.saveas(f"{args.output}/reduced_data.bed")
        print("--- reduced_data.bed generated")

    """
    Step 2: (Optional) Filtering the genomic location based on the chromosomes. (default= chr1-22 and X,Y)
    """
    # if no -c argument is used --> default = ""
    if args.chromosome_list == "":
        print("No filtering by chromosomes ...\n")
        filtered_data = reduced_data
        print("Step 2 skipped\n")
    # if -c argument is used. If -c "default" is given, chr1-22 and chrX and chrY will be filtered. Otherwise the -c argument need to be a list with suitable chromosome names.  
    else:
        print("... Filtering by chromosomes ...\n")
        if args.chromosome_list == "deault":
            filtered_data = filtering_chromosomes(reduced_data).saveas()
            print("Step 2 successfully completed. Only chr1-22 and chrX and chrY remain ...\n")
        else:
            filtered_data = filtering_chromosomes(reduced_data, args.chromosome_list).saveas()
            print(f"Step 2 successfully completed. Only {args.chromosome_list} remain ...\n")

        # Optional save:
        if args.all_output:
            filtered_data.saveas(f"{args.output}/filtered_data.bed")
            print("--- filtered_data generated.\n")


    """
    Step 3: Refining the Data and changing columns: Information from "name" column will be extracted and they'll replace the unnecesarry columns.
            (If possible with .each(func) and not with Dataframe. In case the Bedfile will be too big for beeing saved in a DataFrame.)
    """
    refined_data = refine_BedTool(filtered_data).saveas()

     # Optional save:
    if args.all_output:
        refined_data.saveas(f"{args.output}/refined_data.bed")
        print("--- refined_data generated.\n")

    """
    Step 4: Merge repetitive Entrys
    """
    merged_data = merge_bed_considering_name(refined_data).saveas(f"{args.output}/output_data.bed")

    print("--- output_data.bed generated.\n")
    print("... Script complete.\n")





if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog = 'Processes the BED File from the UniBind robust Collection by filterinng with a second BED File.',
                    description = """
                    This Script processes the BED File from the UniBind robust Collection, which contains all TFBS (in the human Genome). 
                    STEP 1: Filtering the genomic location by a second BED File, which contains the important Regions (e.g. Promotor regions etc.). 
                    STEP 2: (Optional) Filtering the genomic location based on the chromosomes. (default= chr1-22 and X,Y)
                    STEP 3: Refining the Data and changing columns: Information from "name" column will be extracted and they'll replace the unnecesarry columns.
                    STEP 4: Merge repetitive Entrys and save into a new BED file)
                    """   )
    # Required
    # Input and Output of the files:
    parser.add_argument('-f', '--filename', required=True,
            help='give "UniBind robust Collection" file path')
    parser.add_argument('-fb', '--intersect', required=True,
            help='give BED file for intersection')
    parser.add_argument('-out', '--output', required=True,
            help='give output path for folder where new BED file(s) were saved.')
    # Optional
    parser.add_argument('-c', '--chromosome_list', required=False, default="",
            help='give list with chromosomes. They have to match with BedTool.chrom entrys. If you Type  "default", chr1-22 and chrX and Y are used.')   
    parser.add_argument('-all', '--all_output', required=False, action="store_true",
            help='give True, default = False. Also saves files during the process as BED file.')

    args = parser.parse_args()
    
    main()
