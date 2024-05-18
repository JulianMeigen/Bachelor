import numpy as np
import pandas as pd

#import matplotlib as mpl
import matplotlib.pyplot as plt

import seaborn as sns
#import seaborn.objects as so

from scipy.stats import pearsonr, spearmanr
#import glob



#p
# path = "/sybig/projects/GeneRegulation/data/jme/Bachelorarbeit/data/Single_TFBS/All_TFBS.csv"
path = r"C:\Users\julia\MyPython\BA\All_TFBS.csv"
df_combined = pd.read_csv(path)
df_combined = df_combined[df_combined["chr"] != "chrM"]

# df_combined[df_combined["tf"]=="HMBOX1"]

# tfbs_multiple_df = df_combined[df_combined["homotypic_count"] > 1]
# tfbs_for_tf = tfbs_multiple_df[tfbs_multiple_df["tf"] == "HMBOX1"]
# tfbs_for_tf
# tfbs_multiple_df
# # Calculate Strand-Percent


# # change orientaion to 1 for nT and -1 for T strand in DataFrame
# df_change = tfbs_multiple_df.copy()
# df_change["strand_orientation"] = np.array([1 if i == "nT" else -1 for i in tfbs_multiple_df["strand_orientation"]])
# df_change

# # Calculate pct for orientation! in same len like df_change etc!
# test = df_change.groupby(by=["tf", "geneID", "homotypic_count"])["strand_orientation"].transform("mean")
# pct, count = np.unique(test, return_counts=True)
# for i in range(len(pct)):
#     print((pct[i], count[i]))


# test

def add_strand_orientation_pct(df):
    """
    Calculates the fraction for TFBS occurences for each TFBS and adds a new column at the end of the DataFrame.
    1 = all non-Template Strand 
    -1 =  all Template strand
    """
    df_change = df.copy()
    df_change["strand_orientation"] = np.array([1 if i == "nT" else -1 for i in df["strand_orientation"]])
    mean_orientation_arr = df_change.groupby(by=["tf", "geneID", "homotypic_count"])["strand_orientation"].transform("mean")
    df_new = df.insert(loc=6, column="homotypic_strand_orientation", value=mean_orientation_arr)
    return df_new
    
df_new = add_strand_orientation_pct(df_combined)
df_new 