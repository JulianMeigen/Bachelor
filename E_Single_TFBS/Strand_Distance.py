import numpy as np
import pandas as pd

#import matplotlib as mpl
import matplotlib.pyplot as plt

import seaborn as sns
#import seaborn.objects as so

from scipy.stats import pearsonr, spearmanr
#import glob



#path = "/sybig/projects/GeneRegulation/data/jme/Bachelorarbeit/data/Single_TFBS/All_TFBS.csv"
path = r"H:\Uni - iCloud - Alte Semester\BA\All_TFBS.csv"
df_combined = pd.read_csv(path)
df_combined = df_combined[df_combined["chr"] != "chrM"]
df_combined[df_combined["tf"]=="HMBOX1"]

tfbs_multiple_df = df_combined[df_combined["homotypic_count"] > 1]
tfbs_for_tf = tfbs_multiple_df[tfbs_multiple_df["tf"] == "HMBOX1"]
tfbs_for_tf
tfbs_multiple_df
# Calculate Strand-Percent


# change orientaion to 1 for nT and -1 for T strand in DataFrame
df_change = tfbs_multiple_df.copy()
df_change["strand_orientation"] = np.array([1 if i == "nT" else -1 for i in tfbs_multiple_df["strand_orientation"]])
df_change

# Calculate pct for orientation! in same len like df_change etc!
test = df_change.groupby(by=["tf", "geneID", "homotypic_count"])["strand_orientation"].transform("mean")
pct, count = np.unique(test, return_counts=True)
for i in range(len(pct)):
    print((pct[i], count[i]))