'''
Calculate the spearman correlation between the tfbs count and GeneExpression
'''



import pybedtools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
import timebudget


###################################    Previous Functions    ###################################

##### Functions to filter for GeneType or Chromosome

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


##### Functions to filter for specific tfbs


#Function to filter BED file for specific 
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
 
def get_tfbs_subset(BedTool, tfbs_lst):
    """
    Filtering BedTool, so that only Promotors/GeneIDs remain, that contain every TFBS in tfbs_lst at least once.
    """
    return BedTool.filter(func_for_tfbs_subset, tfbs_lst)

### Extract BedTool informations

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



###################################    New Functions    ###################################


def homotypic_count_interval(BedTool_Interval, tf_name):
    chr, geneID, tfbs_arr, tfbs_close_tss_arr, tfbs_dist_tss_arr, tfbs_strand_orientation, tfbs_count, tfbs_unique_count= get_tfbs_arrays_from_BedTool_Interval(BedTool_Interval)
    # tfbs mask
    tfbs_bool = np.array([True if tfbs==tf_name else False for tfbs in tfbs_arr ])
    # get mean tss distance
    mean_tfbs_dist = tfbs_dist_tss_arr - tfbs_close_tss_arr
    tf_dist_all = mean_tfbs_dist[tfbs_bool]
    # get strand
    tf_ori = tfbs_strand_orientation[tfbs_bool]
    # get homotypic count
    homotypic_count = len(tf_dist_all)

    # generate 2d array for multiple tfbs in one prom
    homotypic_information = np.array([[chr,geneID,tf_dist,tf_ori[i], homotypic_count] for i, tf_dist in enumerate(tf_dist_all)])

    return homotypic_information


def get_homotypic_count_for_tf(BedTool, tf_name):
    arr = np.empty((0,5))
    for i in BedTool:
        homotypic_information = homotypic_count_interval(i, tf_name)
        arr = np.vstack((arr, homotypic_information))
    return arr



def get_GeneExpr_for_geneIds(gtex_df, geneIds):
    gene_expr = np.empty((0,len(gtex_df.columns)-1))
    for id in geneIds:
        red_df = gtex_df[gtex_df.Name == id]
        if len(red_df) > 0 :
         expr = red_df.to_numpy()[0][1:]
        else:
            print(f"For {id} was no Expression Data found.")
            expr = np.zeros((1,len(gtex_df.columns)-1))

        gene_expr = np.vstack((gene_expr , expr))    
    return gene_expr




def main(BedTool, geneexpr, tf_name):

    # Filter BedTool for GeneType for "protein_coding" and for specific TFBS
    data = filter_bed(BedTool, 10, "protein_coding")
    tfbs_subset = get_tfbs_subset(data, [tf_name])

    # Extract important information out of every interval
    tfbs_info_2d = get_homotypic_count_for_tf(tfbs_subset, tf_name)
    geneIds = tfbs_info_2d[:,1]
    counts = tfbs_info_2d[:,-1].astype(int)

    # Extract Geneexpresson for every GeneId out of geneexpr Dataframe
    tfbs_gene_expr = get_GeneExpr_for_geneIds(gtex_df, geneIds)

    # Expand every tfbs count for 54 times and flatten the geneexpression data in order to calculate the correlations between all tissues.
    counts_copy = np.repeat(counts, np.shape(tfbs_gene_expr)[1])
    tfbs_gene_expr_flatten = tfbs_gene_expr.flatten()

    # Calculate Spearman correlation
    r, p_value = spearmanr(counts_copy, tfbs_gene_expr_flatten)

    return r, p_value


if __name__ == "__main__":

    data = pybedtools.BedTool("/sybig/projects/GeneRegulation/data/jme/Bachelorarbeit/data/Promotor_with_TFBS/New_TFBS_BED/Prom_with_TFBSs.bed")
    gtex_df = pd.read_csv("/sybig/projects/GeneRegulation/data/jme/Bachelorarbeit/data/GTEx_GenExpr_ucsc.csv", sep=",")


    main(BedTool=data, geneexpr=gtex_df, tf_name="ESR1")