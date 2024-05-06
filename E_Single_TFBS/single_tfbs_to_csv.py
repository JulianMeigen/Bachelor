





import pybedtools
import pandas as pd
import numpy as np
import argparse

from collections import Counter

from timebudget import timebudget
from multiprocessing import Pool
from itertools import repeat

####################################   Functions to Filter the BedTool  ####################################

def Check_if_input_in_BedTool_column(BedTool, column_number, input):
    # Creating an array with unique entrys of BedTool column
    BedTool_column = np.unique(np.array(list(map(lambda x: x.fields[column_number], BedTool))))
    # Checking if input is in column
    if np.isin(np.array(input), BedTool_column):
        return True
    else:
        return False


def filter_bed(BedTool, column_number, filter_content):
    """
    Function that filters a bed file by a specific content in a specific column. 
    It Checks first, if filter_content is in BedTool column.
    Bsp. for project-Data
    Filter for GeneType: filter_bed(BedTool, 10, "protein_coding")
    Filter for Chromosome: filter_bed(BedTool, 0, "chr22")
    """

    if Check_if_input_in_BedTool_column(BedTool, column_number, filter_content):
        filtered_BedTool = BedTool.filter(lambda x: x.fields[column_number] == filter_content)
    else:
        print(f"\nWARNING: '{filter_content}' could not be found in column {column_number} of given BedTool. Therefore the BedTool object could not be filtered!\n")
        print("The Original BedTool was returned.")
        filtered_BedTool = BedTool
    
    return filtered_BedTool


####################################  Function process BedTool Interval  ####################################


def get_tfbs_arrays_from_BedTool_Interval(BedTool_Interval):
    """
    Function that is splitting the infomration about the TFBS, such as TSS-distance and strand oriantation and saving it in diffrent numpy arrays.
    These numpy arrays can be used in other functions to create a DataFrame for the whole BedTool Object. 
    """
    # Splitting infomration into numpy arrays
    tfbs_arr = np.array(BedTool_Interval.fields[6].split(","))
    tfbs_close_tss_arr = np.array(BedTool_Interval.fields[7].split(","), dtype="int")
    tfbs_dist_tss_arr = np.array(BedTool_Interval.fields[8].split(","), dtype="int")

    # Get information if tfbs is on same strand as gene or not (True, False)
    tfbs_strand = np.array(BedTool_Interval.fields[9].split(","))
    tfbs_strand_orientation = np.array(["nT" if i==BedTool_Interval.fields[5] else "T" for i in tfbs_strand])

    # Save GeneID, chr necessary for additional functions. 
    geneID = BedTool_Interval.fields[3]
    chr = BedTool_Interval.fields[0]

    # Generate a numpy array, that contains the count of all tfbs. Its a array with the same length as the other, but all containting the same tfbs count.
    tfbs_count = len(tfbs_arr)
    # Generate a numpy array, that contains the count of all unqiue tfbs.(s.o)
    tfbs_unique_count = len(np.unique(tfbs_arr))

    return chr, geneID, tfbs_arr, tfbs_close_tss_arr, tfbs_dist_tss_arr, tfbs_strand_orientation, tfbs_count, tfbs_unique_count




####################################  Changed for single TFBS to process BedTool Interval ####################################
def process_prom_single_tfbs(BedTool_Interval):
    chr, geneID, tfbs_arr, tfbs_close_tss_arr, tfbs_dist_tss_arr, tfbs_strand_orientation, tfbs_count, tfbs_unique_count = get_tfbs_arrays_from_BedTool_Interval(BedTool_Interval)
    
    # generate an array containing the homotypic counts of the tfbs 
    homotypic_count = np.array([Counter(tfbs_arr)[tfbs] for tfbs in tfbs_arr])
    # Summarize TFBS-Informations as concatenated string
    summarized_arrays = np.array([f"{chr},{geneID},{tfbs_arr[i]},{tfbs_close_tss_arr[i]},{tfbs_dist_tss_arr[i]},{tfbs_strand_orientation[i]},{homotypic_count[i]},{tfbs_count},{tfbs_unique_count}\n" for i in range(len(tfbs_arr))])
        
    return tfbs_arr, summarized_arrays



def get_GeneExpr_for_geneId(gtex_df, geneId):
        
    red_df = gtex_df[gtex_df.Name == geneId]
    if len(red_df) > 0 :
        expr = red_df.to_numpy()[0][1:]
    else:
        print(f"For {geneId} was no Expression Data found.")
        expr = np.empty((1,len(gtex_df.columns)-1))
        expr[:] = np.nan  
    return expr


def process_prom_single_tfbs_with_Expression(BedTool_Interval, gtex_df):
    """
    Same as process_prom_single_tfbs but extracting Genexpression Data from a Dataframe, which contains the GeneID and the correspponding GeneExpression for the Region.
    """
    chr, geneID, tfbs_arr, tfbs_close_tss_arr, tfbs_dist_tss_arr, tfbs_strand_orientation, tfbs_count, tfbs_unique_count = get_tfbs_arrays_from_BedTool_Interval(BedTool_Interval)
    
    # Extract Geneexpression from gtex_df and append it to the csv file.
    GeneExpr = get_GeneExpr_for_geneId(gtex_df, geneID)
    GeneExpr_str = ",".join([str(i) for i in GeneExpr])

    # generate an array containing the homotypic counts of the tfbs 
    homotypic_count = np.array([Counter(tfbs_arr)[tfbs] for tfbs in tfbs_arr])
    # Summarize TFBS-Informations as concatenated string
    summarized_arrays = np.array([f"{chr},{geneID},{tfbs_arr[i]},{tfbs_close_tss_arr[i]},{tfbs_dist_tss_arr[i]},{tfbs_strand_orientation[i]},{homotypic_count[i]},{tfbs_count},{tfbs_unique_count},{GeneExpr_str}\n" for i in range(len(tfbs_arr))])
        
    return tfbs_arr, summarized_arrays



#################################### Function to write single Interval into csv File ####################################

def write_and_append_csv_files_single_tfbs(tfbs_arr, summarized_arrays, output_folder):
    for idx, tf in enumerate(tfbs_arr):
        file_name = f"{output_folder}/{tf}.csv"
        with open(file_name, "a") as f:
            f.write(summarized_arrays[idx])


def write_csv_file_for_interval(BedTool_Interval, output_folder):
        if process_prom_single_tfbs(BedTool_Interval) is not None:
            tfbs_arr, summarized_array = process_prom_single_tfbs(BedTool_Interval)
            write_and_append_csv_files_single_tfbs(tfbs_arr, summarized_array, output_folder)

def write_csv_file_for_interval_with_GeneExpr(BedTool_Interval, output_folder, gtex_df):
        if process_prom_single_tfbs_with_Expression(BedTool_Interval, gtex_df) is not None:
            tfbs_arr, summarized_array = process_prom_single_tfbs_with_Expression(BedTool_Interval, gtex_df)
            write_and_append_csv_files_single_tfbs(tfbs_arr, summarized_array, output_folder)

@timebudget
def main():
    output_folder = args.output
    data = pybedtools.BedTool(args.filename)

    if len(args.filter_GeneType)>0:
        data = filter_bed(data, 10, args.filter_GeneType)
    if len(args.filter_chromosome)>0:
        data = filter_bed(data, 0, args.filter_chromosome)

    # Specify Geneeypression Data
    if len(args.add_GeneExpression_Data) > 0:
        gtex_df = pd.read_csv(args.add_GeneExpression_Data, sep=",")  
        with Pool(4) as pool:
            pool.starmap(write_csv_file_for_interval_with_GeneExpr, zip(data, repeat(output_folder), repeat(gtex_df)))

    else:
        with Pool(4) as pool:
            pool.starmap(write_csv_file_for_interval, zip(data, repeat(output_folder)))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                        prog = 'This Script returns a csv file with all important informations about all single TFBS.',
                        description = """Each Promotor in the input BedFile will be processed and information about each tf will be saved in seperate csv files.
                        These files containing all important informations about that tf and can be used for further processing. """   )
    # Required
    # Input and Output of the files:
    parser.add_argument('-f', '--filename', required=True,
            help='give BED file path containing genomic regions of Promotors, the GeneID, all TFBS etc. It has to be in the right format (see Documentation)')
    parser.add_argument('-out', '--output', required=True,
            help='give output folder path for all csv files')
    
    # Optional
    # Specify GeneType or chromosome for filtering
    parser.add_argument('-geneType', '--filter_GeneType', required=False, default="",
            help='Specify GeneType for filtering given BedTool. ')  
    parser.add_argument('-chr', '--filter_chromosome', required=False, default="",
            help='Specify chromosome for filtering given BedTool. ') 
    
    parser.add_argument('-GeneExpr', '--add_GeneExpression_Data', required=False, default="",
            help='Add file path for an csv file, containing the Geneexpression Data for the GeneIDs (for format specifiactions see Documentation)') 

    args = parser.parse_args()
    
    main()
