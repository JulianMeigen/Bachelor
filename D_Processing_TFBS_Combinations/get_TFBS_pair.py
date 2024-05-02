"""
This script gets a BedTool and a list with multiple TFBS as input and returns a processed DataFrame with informations about the Order, Orientation and TSS distances of the TFBS-pair.

"""

import pybedtools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

from timebudget import timebudget



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



def get_information_for_singel_tfbs_set(BedTool_Interval):
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

    # Generate a numpy array, that contains the GeneID, necessary for additional functions. Its a array with the same length as the other, but all containting the same GeneID.
    geneID = np.array([BedTool_Interval.fields[3] for i in range(len(tfbs_arr))])
    chr = np.array([BedTool_Interval.fields[0] for i in range(len(tfbs_arr))])

    # Generate a numpy array, that contains the count of all tfbs. Its a array with the same length as the other, but all containting the same tfbs count.
    tfbs_count = np.array([len(tfbs_arr)+1 for i in range(len(tfbs_arr))])
    # Generate a numpy array, that contains the count of all unqiue tfbs.(s.o)
    tfbs_unique_count = np.array([len(np.unique(tfbs_arr))+1 for i in range(len(tfbs_arr))])

    return chr, geneID, tfbs_arr, tfbs_close_tss_arr, tfbs_dist_tss_arr, tfbs_strand_orientation, tfbs_count, tfbs_unique_count


@timebudget
def expand_tfbs_information_to_df(BedTool):
    """
    Function that iterates over all BedTool-Intervals and extracting inforations about TFBS by using the get_information_for_singel_tfbs_set function.
    The processed arrays for every BedTool entry will be appended to one big Dataframe column.
    The function returns a DataFrame with an entry for every TFBS and their coressponding information such as tss distance and strand orientation in all Promotors.
    """

    # Information for Promotor/Gene. Will be the same for all TFBS in one Region
    chr_all = np.array([])
    geneID_all = np.array([])

    # Information about specific TFBS, such as name, tss-distance and orientation
    tfbs_all = np.array([])
    tfbs_close_tss_all = np.array([])
    tfbs_dist_tss_all = np.array([])
    tfbs_orientation_all = np.array([])

    # Meta-Information about all TFBSs in one promotor region. Will be the same for each geneID.
    tfbs_count_all = np.array([])
    tfbs_unique_count_all =np.array([])

    for BedTool_Interval in BedTool:
        # Extracting information about one Promotor. Each Promotor returning array with the length of the tfbs-count.
        chr, geneID, tfbs_arr, tfbs_close_tss_arr, tfbs_dist_tss_arr, tfbs_orientation, tfbs_count, tfbs_unique_count = get_information_for_singel_tfbs_set(BedTool_Interval)

        # Appending  the _all arrays. Every array will be a DataFrame column later on.
        chr_all = np.append(chr_all, chr)
        geneID_all = np.append(geneID_all, geneID)
        tfbs_all = np.append(tfbs_all, tfbs_arr)
        tfbs_close_tss_all = np.append(tfbs_close_tss_all, tfbs_close_tss_arr)
        tfbs_dist_tss_all = np.append(tfbs_dist_tss_all, tfbs_dist_tss_arr)
        tfbs_orientation_all = np.append(tfbs_orientation_all, tfbs_orientation)
        tfbs_count_all = np.append(tfbs_count_all, tfbs_count)
        tfbs_unique_count_all = np.append(tfbs_unique_count_all, tfbs_unique_count)

    # Generating DataFrame with specific column names.
    df = pd.DataFrame({"chr":chr_all, "geneID":geneID_all, "TFBS": tfbs_all, "Distance_to_TSS_close":tfbs_close_tss_all, "Distance_to_TSS_dist":tfbs_dist_tss_all, "Orientation": tfbs_orientation_all, "TFBS_count":tfbs_count_all, "TFBS_count_unique": tfbs_unique_count_all})
    return df


def filter_homotypic(df, tfbs_lst, get_homotypic=False):
    """"
    Function that filters out every Promotor/GeneID, which contains more then one TFBS of the TFs, specified in tfbs_lst.
    It returns a DataFrame which contains every TFBS in TFBS_lst exaclty once per Promotor/GeneID.
    If get_homotypic=True, it also returns filtered DataFrame AND the DataFrame containing the Promotors with the homotypic TFBS.
    """
    error_GeneIDs = []
    for tf in tfbs_lst:
        tf_homo = df[df["TFBS"]==tf]
        geneID, tf_homo_counts = np.unique(tf_homo["geneID"], return_counts=True)
        error_geneID = geneID[tf_homo_counts != 1]
        error_GeneIDs= np.append(error_GeneIDs, error_geneID)
    
    GeneIDs_not_in_error_GeneIDS = ~df["geneID"].isin(np.unique(error_GeneIDs))
    
    filtered_df = df[GeneIDs_not_in_error_GeneIDS]

    # If get_homotypic, it returns also the Dataframe containing homotypic TFBS.
    if get_homotypic:
        homotypic_df = df[df["geneID"].isin(np.unique(error_GeneIDs))]
        return filtered_df, homotypic_df

    return filtered_df
    

@timebudget
def get_orientation_for_pair(pair_df_sorted):
    """
    Function which is used in get_orientation_and_order_for_pair.
    It gets an DataFrame as input, which only consist of two rows, containing the two TFBS in one Promotor.
    It returns the Orientation of the TFBS-pair
    """

    tfbs_close_orientation = pair_df_sorted.iloc[0]["Orientation"]
    tfbs_dist_orientation = pair_df_sorted.iloc[1]["Orientation"]
    
    if (tfbs_close_orientation == "T") and (tfbs_dist_orientation == "T"):
        Orientation_of_pair = "both-T" 
    elif (tfbs_close_orientation == "nT") and (tfbs_dist_orientation == "nT"):
        Orientation_of_pair = "both-nT"
    elif (tfbs_close_orientation == "T") and (tfbs_dist_orientation == "nT"):
        Orientation_of_pair = "convergent"
    elif (tfbs_close_orientation == "nT") and (tfbs_dist_orientation == "T"):
        Orientation_of_pair = "divergent"

    return Orientation_of_pair



@timebudget
def get_orientation_and_order_for_pair(df):
    """
     Function that gets an DataFrame as input, which should only contain two entrys per Promotor/GeneID, because it will only consider the first two rows.
     If thats not the case you should use filter_homotypic() first.

     It returns a Dataframe, where each entry represents one Promotor/GeneID and the information about the TFBS-pair within these Promotors.
     Explanation for the diffrent columns:

     chr:                   Chromosome on which the TFBS pair is located
     geneID:                The GenID of the associated gene (for the determination of gene expression)
     TFBS_close             The name of the TF closest to the TSS
     TFBS_dist              The name of the TF that is furthest away from the TSS
     closest_TSS_distance   The next distance to the TSS. In other words, from the TF that is closest to the TSS, the side of the binding side that is closest.
     furthest_TSS_distance  The furthest distance to the TSS. In other words, from the TF that is furthest away from the TSS, the side of the binding point that is furthest away.
     orientation_of_pair    The strand orientation of the pair, depending on the gen. There are four different possibilities. Both TFBS are on template or non-template strand. Or they are on different strands, whereby the TFBS can be convergent or divergent to each other.
     distance_between_pair  The distance between the TFBS. If this is negative, the TFBS overlap
     tfbs_count             Information on the number of TFBSs present in the promoter region
     tfbs_count_unique      Information on the number of unique TFBSs present in the promoter region.

    """

    geneIDs = np.unique(df["geneID"])

    lst_2D = np.array(["chr", "geneID", "TFBS_close", "TFBS_dist", "closest_TSS_distance", "furthest_TSS_distance", "orientation_of_pair", "distance_between_pair", "tfbs_count", "tfbs_count_unique"])
    for geneID in geneIDs:
        pair_df_sorted = df[df["geneID"]==geneID].sort_values(by=["Distance_to_TSS_close"])

        # Columns for new Dataframe
        # get additional information about geneID such as chr, tfbs_count and tfbs_count_unique
        chr = np.unique(pair_df_sorted["chr"])[0]
        tfbs_count = np.unique(pair_df_sorted["TFBS_count"])[0]
        tfbs_count_unique = np.unique(pair_df_sorted["TFBS_count_unique"])[0]
        
        # Because the df is sorted by the tf closest to TSS, the first entry gets tf close to TSS, and the second entry the tf further away. 
        TFBS_close = pair_df_sorted.iloc[0]["TFBS"]
        TFBS_dist = pair_df_sorted.iloc[1]["TFBS"]

        # Get TFBS-Range (closest position to TSS and furthest position)
        closest_TSS_distance = pair_df_sorted.iloc[0]["Distance_to_TSS_close"]
        furthest_TSS_distance = pair_df_sorted.iloc[1]["Distance_to_TSS_dist"]

        # Get Orientation of the pair. It can be both-Template, both-not-Template, convergent (TSS-T-nT), Divergent(TSS-nT-T)
        orientation_of_pair = get_orientation_for_pair(pair_df_sorted)

        # Get Distance between pair. 
        distance_between_pair = pair_df_sorted.iloc[1]["Distance_to_TSS_close"] - pair_df_sorted.iloc[0]["Distance_to_TSS_dist"]
 

        values_per_col = [chr, geneID, TFBS_close, TFBS_dist, closest_TSS_distance, furthest_TSS_distance, orientation_of_pair, distance_between_pair, tfbs_count, tfbs_count_unique]
        lst_2D = np.vstack((lst_2D, np.array(values_per_col)))

    # Generating DataFrame
    df_pair = pd.DataFrame(data = lst_2D[1:,:], 
                      columns= lst_2D[0,:])
    df_pair["closest_TSS_distance"] = pd.to_numeric(df_pair["closest_TSS_distance"])
    df_pair["furthest_TSS_distance"] = pd.to_numeric(df_pair["furthest_TSS_distance"])
    df_pair["distance_between_pair"] = pd.to_numeric(df_pair["distance_between_pair"])
    df_pair["tfbs_count"] = pd.to_numeric(df_pair["tfbs_count"])
    df_pair["tfbs_count_unique"] = pd.to_numeric(df_pair["tfbs_count_unique"])

    return df_pair




@timebudget
def get_pair(BedTool, tfbs_lst, filter_by_geneType="", filter_by_chr=""):
    """
    Function that combines the previous functions to get the DataFrame for a Pair.
    """
    # Step 1: Filtering BedTool, so that only Promotors/GeneIDs remain, that contain every TFBS in tfbs_lst at least once.
    tbfs_lst_BedTool = BedTool.filter(func_for_tfbs_subset, tfbs_lst)

    # (Optional) Step 2: Filtering by geneType or chromosome
    if len(filter_by_geneType)>0:
        tfbs_lst_BedTool = filter_bed(tbfs_lst_BedTool, 10, filter_by_geneType)
    if len(filter_by_chr)>0:
        tfbs_lst_BedTool = filter_bed(tbfs_lst_BedTool, 0, filter_by_chr)
    
    # Step 3: 
    # Generate a BedTool object with information for every TFBS in each Promotor region
    all_tfbs_df = expand_tfbs_information_to_df(tbfs_lst_BedTool)
    # Reduce this DataFrame to the TFBS specified in tfbs_lst
    tfbs_pair_df = all_tfbs_df[all_tfbs_df["TFBS"].isin(tfbs_lst)]

    # Step 4: Filtering the Promotor/GeneIDs, which contain multiple TFBS of the same TF.
    tfbs_pair_hetero, tfbs_pair_homo = filter_homotypic(tfbs_pair_df, tfbs_lst, get_homotypic=True)
    print(f"Note: {len(tfbs_pair_homo)}/{len(tfbs_pair_hetero)+len(tfbs_pair_homo)} TFBS were removed because one of the TFBS from the TFBS pair appeared several times in the promoter region. {len(tfbs_pair_hetero)} TFBS remain.")
    
    # Step 5: Generating DataFrame, which contains important information about Order, Orientaion and TSS-Distance of TFBS-pair.
    tfbs_pair_ord_ori_df = get_orientation_and_order_for_pair(tfbs_pair_hetero)

    return  tfbs_pair_ord_ori_df


def get_possible_pairs(all_tfbs):
    all_tfbs.sort()
    pairs_as_tuple = [(a, b) for idx, a in enumerate(all_tfbs) for b in all_tfbs[idx + 1:]]
    return pairs_as_tuple


def get_csv_single_pair(data, tfbs_pair):
    df = get_pair(data, tfbs_pair, filter_by_geneType=args.filter_GeneType, filter_by_chr=args.filter_chromosome)
    df.to_csv(f"{args.output}/{tfbs_pair[0]}_{tfbs_pair[1]}.csv", index=False)
    return df

def main():
    data = pybedtools.BedTool(args.filename)

    # Set tfbs_lst to UniBind TFBS if -return_uniBind is specified.
    if args.return_all_uniBind_tfbs_pairs:
        tfbs_lst = ['AR', 'ARNT', 'ARNTL', 'ASCL1', 'ATF2', 'ATF3', 'ATF4', 'ATF7', 'BACH1', 'BACH2', 'BATF', 'BATF3', 'BCL6', 'BCL6B', 'BHLHE22', 'BHLHE40', 'CDX2', 'CEBPA', 'CEBPB', 'CEBPD', 'CEBPG', 'CLOCK', 'CREB1', 'CREM', 'CTCF', 'CTCFL', 'CUX1', 'DUX4', 'E2F1', 'E2F4', 'E2F6', 'E2F7', 'E2F8', 'EBF1', 'EBF3', 'EGR1', 'EGR2', 'EGR3', 'EHF', 'ELF1', 'ELF3', 'ELF4', 'ELF5', 'ELK1', 'ELK4', 'EOMES', 'ERF', 'ERG', 'ESR1', 'ESR2', 'ESRRA', 'ETS1', 'ETV1', 'ETV4', 'ETV5', 'ETV6', 'FLI1', 'FOS', 'FOSL1', 'FOSL2', 'FOXA1', 'FOXA2', 'FOXA3', 'FOXJ2', 'FOXK1', 'FOXK2', 'FOXO1', 'FOXP1', 'GABPA', 'GATA1', 'GATA2', 'GATA3', 'GATA4', 'GATA6', 'GFI1', 'GFI1B', 'GLIS1', 'GLIS2', 'GLIS3', 'GRHL2', 'HAND2', 'HIC1', 'HIF1A', 'HINFP', 'HLF', 'HMBOX1', 'HNF1A', 'HNF1B', 'HNF4A', 'HNF4G', 'HOXA9', 'HOXB13', 'HSF1', 'HSF2', 'IRF1', 'IRF2', 'IRF3', 'IRF4', 'ISL1', 'JUN', 'JUNB', 'JUND', 'KLF1', 'KLF11', 'KLF12', 'KLF15', 'KLF16', 'KLF4', 'KLF5', 'KLF9', 'LEF1', 'LHX2', 'LHX9', 'MAF', 'MAFB', 'MAFF', 'MAFK', 'MAX', 'MECOM', 'MEF2A', 'MEF2B', 'MEF2C', 'MEF2D', 'MEIS1', 'MEIS2', 'MGA', 'MITF', 'MLX', 'MNT', 'MXI1', 'MYB', 'MYBL2', 'MYC', 'MYCN', 'MYF5', 'MYOD1', 'MYOG', 'NEUROD1', 'NEUROG2', 'NFE2', 'NFE2L2', 'NFIA', 'NFIB', 'NFIC', 'NFIL3', 'NFKB1', 'NFKB2', 'NFYA', 'NFYB', 'NFYC', 'NKX2-5', 'NKX3-1', 'NR1H2', 'NR1H3', 'NR1H4', 'NR2C1', 'NR2C2', 'NR2F1', 'NR2F2', 'NR2F6', 'NR3C1', 'NR4A1', 'NR5A1', 'NR5A2', 'NRF1', 'OCT4', 'ONECUT1', 'ONECUT2', 'OSR2', 'OTX2', 'PAX5', 'PAX6', 'PBX1', 'PBX2', 'PBX3', 'PDX1', 'PHOX2B', 'PKNOX1', 'PLAG1', 'POU2F1', 'POU2F2', 'POU5F1', 'PPARG', 'PRDM1', 'PRDM4', 'RARA', 'RARG', 'RBPJ', 'REL', 'RELA', 'RELB', 'REST', 'RFX1', 'RFX2', 'RFX3', 'RFX5', 'RUNX1', 'RUNX2', 'RUNX3', 'RXRA', 'RXRB', 'SCRT1', 'SCRT2', 'SIX1', 'SIX2', 'SMAD2', 'SMAD3', 'SMAD4', 'SNAI2', 'SOX10', 'SOX13', 'SOX17', 'SOX2', 'SOX4', 'SOX6', 'SOX9', 'SP1', 'SP2', 'SP3', 'SP4', 'SPDEF', 'SPI1', 'SPIB', 'SREBF1', 'SRF', 'STAT1', 'STAT2', 'STAT3', 'STAT4', 'STAT5A', 'STAT5B', 'STAT6', 'T', 'TAL1', 'TBX21', 'TBX5', 'TCF12', 'TCF3', 'TCF4', 'TCF7', 'TCF7L1', 'TCF7L2', 'TEAD1', 'TEAD3', 'TEAD4', 'TFAP2A', 'TFAP2C', 'TFAP4', 'TFDP1', 'TFE3', 'TFEB', 'THAP11', 'THRB', 'TP53', 'TP63', 'TP73', 'TWIST1', 'USF1', 'USF2', 'VDR', 'WT1', 'XBP1', 'YY1', 'YY2', 'ZFX', 'ZNF143', 'ZNF263', 'ZNF740']
    else:
        tfbs_lst = args.tfbs_list

    # get csv for all pairs
    print("\nCalculate Script for each pair...\n")
    pairs_lst = get_possible_pairs(tfbs_lst)
    print(f"A file will be created for the following pair(s):\n")
    [print(pair) for pair in pairs_lst]

    number_pairs = len(pairs_lst)
    current_pair_number = 0
    print("----------------------   0 %   ----------------------")
    for idx, tfbs_pair in enumerate(pairs_lst):
        print(f"{tfbs_pair}:\n ... Generating csv file ...\n")
        get_csv_single_pair(data, tfbs_pair)
        print(f"... File has been created ...\n")
        
        # print status.
        current_pair_number =+ idx
        print(f"----------------------   {round(((current_pair_number+1)/number_pairs)*100, 2)} %   ----------------------")
        






if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog = 'This Script returns a csv file with all important informations about a given TFBS pair.',
                    description = """The input BedFile will be filtered for the tfbs pair subset. Each Promotor/GenID containing multiple copies of one od the TFBS will be filtered.
                     A csv file will be generated, containing all important informations about that tfbs pair and can be used for further processing. """   )
    # Required
    # Input and Output of the files:
    parser.add_argument('-f', '--filename', required=True,
            help='give BED file path containing genomic regions of Promotors, the GeneID, all TFBS etc. It has to be in the right format (see Documentation)')
    parser.add_argument('-tfbs_lst', '--tfbs_list', nargs="*", type=str, required=True, 
            help='give a list, containing TFBS names. You can specify more than two TFBSs, and the script will be run for each possible pair of TFBSs.')
    parser.add_argument('-out', '--output', required=True,
            help='give output path for csv file')
    
    # Optional
    # Specify GeneType or chromosome for filtering
    parser.add_argument('-geneType', '--filter_GeneType', required=False, default="",
            help='Specify GeneType for filtering given BedTool. ')  
    parser.add_argument('-chr', '--filter_chromosome', required=False, default="",
            help='Specify chromosome for filtering given BedTool. ') 
    
    # Apply script on every possible tfbs_pair.
    parser.add_argument('-return_uniBind', '--return_all_uniBind_tfbs_pairs', action="store_true",
        help='If -return_all is used, all possible TFBS pairs are calculated from the TFs specified in tfbs_pair. The original script will then be applied to all TFBS pairs. ') 

    args = parser.parse_args()
    
    
    data = pybedtools.BedTool("/sybig/projects/GeneRegulation/data/jme/Bachelorarbeit/data/Promotor_with_TFBS/All_GTEx_Prom_with_TFBS.bed")

    df = get_pair(data, ("CTCF", "REST"))
    print(df)






