import numpy as np
import pandas as pd

#import matplotlib as mpl
#import matplotlib.pyplot as plt

#import seaborn as sns
#import seaborn.objects as so

#from scipy.stats import pearsonr, spearmanr
#import glob



def add_homotypic_TFBS_based_on_minimal_TFBS_len(df_combined):
    grouped_tf_ratio_TFBS = np.floor(df_combined.groupby(["tf"])[["dist_tss","close_tss"]].apply(lambda x: (x.dist_tss-x.close_tss)/(x.dist_tss-x.close_tss).min()))
    tf_ratio_TFBS = grouped_tf_ratio_TFBS.reset_index(level=[0]).sort_index()[0].astype(int)
    
    df_new = pd.DataFrame(columns=df_combined.columns)
    for idx in df_combined.index:
        for repeat in range(tf_ratio_TFBS[idx]):
            df_new = pd.concat([df_new, pd.DataFrame(df_combined.iloc[idx]).T], ignore_index=True)
    return df_new


if __name__ == "__main__":
    
     # LOAD DATA
    path = r"C:\Users\julia\MyPython\BA\All_TFBS_new.csv"
    df_combined = pd.read_csv(path)

    # PROCESS DATAFRAME
    df_add_homotypic =  add_homotypic_TFBS_based_on_minimal_TFBS_len(df_combined)

    # SAVE new DATAFRAME
    path_homotypic = r"C:\Users\julia\MyPython\BA\Homotypic_All_TFBS_new.csv"
    df_add_homotypic.to_csv(path_homotypic, index=False)