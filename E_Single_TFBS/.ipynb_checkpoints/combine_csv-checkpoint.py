import pybedtools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import pearsonr, spearmanr
import glob

#import seaborn as sns


def get_df_for_single_tfbs(csv_file):
    # Specify column names
    columns_01 = ["chr","geneID","tf","close_tss","dist_tss","strand_orientation","homotypic_count","all_tfbs_count","all_tfbs_unique_count"]
    columns_02 = ['Adipose - Subcutaneous', 'Adipose - Visceral (Omentum)',
       'Adrenal Gland', 'Artery - Aorta', 'Artery - Coronary',
       'Artery - Tibial', 'Bladder', 'Brain - Amygdala',
       'Brain - Anterior cingulate cortex (BA24)',
       'Brain - Caudate (basal ganglia)', 'Brain - Cerebellar Hemisphere',
       'Brain - Cerebellum', 'Brain - Cortex', 'Brain - Frontal Cortex (BA9)',
       'Brain - Hippocampus', 'Brain - Hypothalamus',
       'Brain - Nucleus accumbens (basal ganglia)',
       'Brain - Putamen (basal ganglia)', 'Brain - Spinal cord (cervical c-1)',
       'Brain - Substantia nigra', 'Breast - Mammary Tissue',
       'Cells - Cultured fibroblasts', 'Cells - EBV-transformed lymphocytes',
       'Cervix - Ectocervix', 'Cervix - Endocervix', 'Colon - Sigmoid',
       'Colon - Transverse', 'Esophagus - Gastroesophageal Junction',
       'Esophagus - Mucosa', 'Esophagus - Muscularis', 'Fallopian Tube',
       'Heart - Atrial Appendage', 'Heart - Left Ventricle', 'Kidney - Cortex',
       'Kidney - Medulla', 'Liver', 'Lung', 'Minor Salivary Gland',
       'Muscle - Skeletal', 'Nerve - Tibial', 'Ovary', 'Pancreas', 'Pituitary',
       'Prostate', 'Skin - Not Sun Exposed (Suprapubic)',
       'Skin - Sun Exposed (Lower leg)', 'Small Intestine - Terminal Ileum',
       'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina',
       'Whole Blood']
    columns = columns_01 + list(columns_02)
    # Read single csv file
    tfbs_df = pd.read_csv(csv_file, names=columns)
    return tfbs_df

def flatten_tissues_expand_df(tfbs_df):
    tfbs_part_repeat = tfbs_df.iloc[:,:9]
    tfbs_part_flatten = tfbs_df.iloc[:,9:]
    
    repeat_df =  pd.DataFrame(np.repeat(tfbs_part_repeat, len(tfbs_part_flatten.columns), axis=0))
    flatten_arr = tfbs_part_flatten.to_numpy().flatten()
    
    repeat_df.columns = tfbs_part_repeat.columns
    repeat_df["All_tissues"] = flatten_arr

    return repeat_df

def homotyic_pearson(tfbs_df):
    # First the Dataframe will be exoanded to include every tissue.
    expand_df = flatten_tissues_expand_df(tfbs_df)

    homotyic_count = expand_df.homotypic_count.to_numpy()
    geneexpr = expand_df.All_tissues.to_numpy()

    r,p_value = pearsonr(geneexpr, homotyic_count)

    return r,p_value


####################  Pearson etc ######################

def get_pearson_for_single_tf(csv_file):
    tfbs_df = get_df_for_single_tfbs(csv_file)
    tf_name = tfbs_df.tf.unique()[0]
    tfbs_df_flatten = flatten_tissues_expand_df(tfbs_df)
    r, p_value = homotyic_pearson(tfbs_df_flatten)
    return tf_name, r, p_value

def get_pearson_for_all_tf(csv_folder):
    tf_paths = glob.glob(f"{csv_folder}/*.csv")

    tf_dict = dict()
    #p_value_dict = dict()
    for tf_path in tf_paths:
        tf_name, r, p_value = get_pearson_for_single_tf(tf_path)
        tf_dict[tf_name] = r, p_value
        # p_value_dict[tf_name] = p_value
    return tf_dict

#####################  Gene Expr etc ###################

def get_flat_Expr_for_single_tfbs(csv_file):
    tfbs_df = get_df_for_single_tfbs(csv_file)
    tf_name = tfbs_df.tf.unique()[0]
    tfbs_df_flatten = flatten_tissues_expand_df(tfbs_df)

    expr = tfbs_df_flatten.All_tissues.to_numpy()
    
    return tf_name, expr

def get_flat_Expr_all_tfbs(csv_folder):
    tf_paths = glob.glob(f"{csv_folder}/*.csv")

    tf_dict = dict()
    for tf_path in tf_paths:
        tf_name, expr = get_flat_Expr_for_single_tfbs(tf_path)
        tf_dict[tf_name] = expr
    return tf_dict

#############################################################

def get_df_for_all_tfbs(csv_folder):
    # Specify column names
    columns_01 = ["chr","geneID","tf","close_tss","dist_tss","strand_orientation","homotypic_count","all_tfbs_count","all_tfbs_unique_count"]
    columns_02 = ['Adipose - Subcutaneous', 'Adipose - Visceral (Omentum)',
       'Adrenal Gland', 'Artery - Aorta', 'Artery - Coronary',
       'Artery - Tibial', 'Bladder', 'Brain - Amygdala',
       'Brain - Anterior cingulate cortex (BA24)',
       'Brain - Caudate (basal ganglia)', 'Brain - Cerebellar Hemisphere',
       'Brain - Cerebellum', 'Brain - Cortex', 'Brain - Frontal Cortex (BA9)',
       'Brain - Hippocampus', 'Brain - Hypothalamus',
       'Brain - Nucleus accumbens (basal ganglia)',
       'Brain - Putamen (basal ganglia)', 'Brain - Spinal cord (cervical c-1)',
       'Brain - Substantia nigra', 'Breast - Mammary Tissue',
       'Cells - Cultured fibroblasts', 'Cells - EBV-transformed lymphocytes',
       'Cervix - Ectocervix', 'Cervix - Endocervix', 'Colon - Sigmoid',
       'Colon - Transverse', 'Esophagus - Gastroesophageal Junction',
       'Esophagus - Mucosa', 'Esophagus - Muscularis', 'Fallopian Tube',
       'Heart - Atrial Appendage', 'Heart - Left Ventricle', 'Kidney - Cortex',
       'Kidney - Medulla', 'Liver', 'Lung', 'Minor Salivary Gland',
       'Muscle - Skeletal', 'Nerve - Tibial', 'Ovary', 'Pancreas', 'Pituitary',
       'Prostate', 'Skin - Not Sun Exposed (Suprapubic)',
       'Skin - Sun Exposed (Lower leg)', 'Small Intestine - Terminal Ileum',
       'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina',
       'Whole Blood']
    columns = columns_01 + list(columns_02)

    # Get a list of all CSV files in a directory
    csv_files = glob.glob(f"{csv_folder}/*.csv")

    # Create an empty dataframe to store the combined data
    combined_df = pd.DataFrame()

    # Loop through each CSV file and append its contents to the combined dataframe
    for csv_file in csv_files:
        df = pd.read_csv(csv_file)
        df.columns = columns
        combined_df = pd.concat([combined_df, df])
    combined_df.columns = columns
    return combined_df

        



csv_folder = "/sybig/projects/GeneRegulation/data/jme/Bachelorarbeit/data/Single_TFBS/Protein_Region_single_TFBS_with_GTEx"


tfbs = get_df_for_all_tfbs(csv_folder)

tfbs.to_csv(f"/sybig/projects/GeneRegulation/data/jme/Bachelorarbeit/data/Single_TFBS/All_TFBS.csv", index=False)

