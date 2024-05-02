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



def calculate_pairs(tf_s, close_s, dist_s, ori_s):
    """
    Function that calculates tfbs-pairs and associated information in numpy arays of tuples. 
    """
    tf_pairs = np.array([(a, b) for idx, a in enumerate(tf_s) for b in tf_s[idx + 1:]])
    close_pairs = np.array([(a, b) for idx, a in enumerate(close_s) for b in close_s[idx + 1:]])
    dist_pairs = np.array([(a, b) for idx, a in enumerate(dist_s) for b in dist_s[idx + 1:]])
    ori_pairs = np.array([(a, b) for idx, a in enumerate(ori_s) for b in ori_s[idx + 1:]])
    return tf_pairs, close_pairs, dist_pairs, ori_pairs



def sort_tfbs_by_order(tf_pairs, close_pairs, dist_pairs, ori_pairs):
    sorted_idx = np.argsort(close_pairs, axis=1)
    tf_pairs_order = np.take_along_axis(tf_pairs, sorted_idx, axis=1)
    close_pairs_order = np.take_along_axis(close_pairs, sorted_idx, axis=1)
    dist_pairs_order = np.take_along_axis(dist_pairs, sorted_idx, axis=1)
    ori_pairs_order = np.take_along_axis(ori_pairs, sorted_idx, axis=1)
    return tf_pairs_order, close_pairs_order, dist_pairs_order, ori_pairs_order


def get_pair_orientaion(close_ori, dist_ori):
    """
    Function that calculates the orientation of the pair based on the strand orientation of the tfbs that is closest to the tss and the tfbs that is further away from the tss.
    """
    if (close_ori == "T") and (dist_ori == "T"):
        Orientation_of_pair = "both-T" 
    elif (close_ori == "nT") and (dist_ori == "nT"):
        Orientation_of_pair = "both-nT"
    elif (close_ori == "T") and (dist_ori == "nT"):
        Orientation_of_pair = "convergent"
    elif (close_ori == "nT") and (dist_ori == "T"):
        Orientation_of_pair = "divergent"
    
    return Orientation_of_pair


def get_infos_to_pairs(tf_pairs_order, close_pairs_order, dist_pairs_order, ori_pairs_order):
    # Save TFBS-order
    close_tfbs = tf_pairs_order[:,0]
    dist_tfbs = tf_pairs_order[:,1]
    # Calculate closest and furthest distance to tss, which is the smallest number of the closest tfbs and the biggest number of the furthest tfbs.
    closest_dist_to_TSS = close_pairs_order[:,0]
    furthest_dist_to_TSS = dist_pairs_order[:,1]
    # Calculate the distance between pairs
    pair_dist = close_pairs_order[:,1] - dist_pairs_order[:,0]
    # Calculate the Orientation of the pair. Possible Oreintations are both-T, both-nT, konvergent (T,nT) and divergent (nT,T)
    pair_ori = np.array([get_pair_orientaion(close_ori= tfbs_ori[0], dist_ori=tfbs_ori[1]) for tfbs_ori in ori_pairs_order])

    return close_tfbs, dist_tfbs, closest_dist_to_TSS, furthest_dist_to_TSS, pair_dist, pair_ori


def process_prom(BedTool_Interval):
    chr, geneID, tfbs_arr, tfbs_close_tss_arr, tfbs_dist_tss_arr, tfbs_strand_orientation, tfbs_count, tfbs_unique_count = get_tfbs_arrays_from_BedTool_Interval(BedTool_Interval)
    
    #Check if the region has more then one tfbs to form a pair.
    if len(tfbs_arr)>1:

        # Sorting tfbs arrays according to tfbs_name and then to tss distance.
        sorted_idx = np.lexsort((tfbs_close_tss_arr, tfbs_arr))
        tf_s, close_s, dist_s, ori_s = tfbs_arr[sorted_idx], tfbs_close_tss_arr[sorted_idx], tfbs_dist_tss_arr[sorted_idx], tfbs_strand_orientation[sorted_idx]
        # Calculating Pairs of these sorted arrays
        tf_pairs, close_pairs, dist_pairs, ori_pairs = calculate_pairs(tf_s, close_s, dist_s, ori_s)
        # Sorting Pairs by tfbs_order
        tf_pairs_order, close_pairs_order, dist_pairs_order, ori_pairs_order = sort_tfbs_by_order(tf_pairs, close_pairs, dist_pairs, ori_pairs)

        # Get Information about TFBS-pairs
        close_tfbs, dist_tfbs, closest_dist_to_TSS, furthest_dist_to_TSS, pair_dist, pair_ori = get_infos_to_pairs(tf_pairs_order, close_pairs_order, dist_pairs_order, ori_pairs_order)
        # Summarize TFBS-Informations as concatenated string
        summarized_array = np.array([f"{chr},{geneID},{close_tfbs[i]},{dist_tfbs[i]},{closest_dist_to_TSS[i]},{furthest_dist_to_TSS[i]},{pair_dist[i]},{pair_ori[i]},{tfbs_count},{tfbs_unique_count}\n" for i in range(len(tf_pairs_order))])
        
        return tf_pairs, summarized_array
    else:
        print(f"The Region {geneID} does only contain {tfbs_arr[0]}.")
        return None


#################################### Function to write single Interval into csv File ####################################

def write_and_append_csv_files(tf_pairs, summarized_array, output_folder):
    for idx, pair in enumerate(tf_pairs):
        file_name = f"{output_folder}/{pair[0]}_{pair[1]}.csv"
        with open(file_name, "a") as f:
            f.write(summarized_array[idx])


def write_csv_file_for_interval(BedTool_Interval, output_folder):
        if process_prom(BedTool_Interval) is not None:
            tf_pairs, summarized_array = process_prom(BedTool_Interval)
            write_and_append_csv_files(tf_pairs, summarized_array, output_folder)



@timebudget
def main():
    output_folder = args.output
    data = pybedtools.BedTool(args.filename)

    if len(args.filter_GeneType)>0:
        data = filter_bed(data, 10, args.filter_GeneType)
    if len(args.filter_chromosome)>0:
        data = filter_bed(data, 0, args.filter_chromosome)

    
    with Pool(4) as pool:
        pool.starmap(write_csv_file_for_interval, zip(data, repeat(output_folder)))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                        prog = 'This Script returns a csv file with all important informations about all TFBS pairs',
                        description = """Each Promotor in the input BedFile will be processed and all possible pairs will be saved in seperate csv files.
                        These files containing all important informations about that tfbs pair and can be used for further processing. """   )
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

    args = parser.parse_args()
    
    main()
