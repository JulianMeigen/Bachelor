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
df_combined

tfbs_multiple_df = df_combined[df_combined["homotypic_count"] > 1]
tfbs_for_tf = tfbs_multiple_df[tfbs_multiple_df["tf"] == "HMBOX1"]
tfbs_for_tf

