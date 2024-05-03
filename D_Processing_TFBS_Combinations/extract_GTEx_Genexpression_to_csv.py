"""
Extracts the gene expression data for the GTEx file and the GeneID and Expression data are saved as a csv file.
By default it will expect the processed GTEx-TFBS-file where the GeneID = x.fields[3] and GeneExpr = x.fields[11]
But if -ucsc [TRUE]. It will expect the ucsc structure where GeneID = x.fields[6] and GeneExpr = x.fields[9]

"""

import pybedtools
import argparse

def extract_GTEx_data(BedTool_Interval, c_geneID, c_geneExpr):
    GeneID = BedTool_Interval.fields[c_geneID]
    GeneExpr = BedTool_Interval.fields[c_geneExpr][:-1]
    return GeneID, GeneExpr



def main():

    GTEx_BedTool = pybedtools.BedTool(args.gtex_file)

    with open(f"{args.output}/GTEx_GenExpr.csv", "w") as f:
        # Header
        f.write("Name,Adipose - Subcutaneous,Adipose - Visceral (Omentum),Adrenal Gland,Artery - Aorta,Artery - Coronary,Artery - Tibial,Bladder,Brain - Amygdala,Brain - Anterior cingulate cortex (BA24),Brain - Caudate (basal ganglia),Brain - Cerebellar Hemisphere,Brain - Cerebellum,Brain - Cortex,Brain - Frontal Cortex (BA9),Brain - Hippocampus,Brain - Hypothalamus,Brain - Nucleus accumbens (basal ganglia),Brain - Putamen (basal ganglia),Brain - Spinal cord (cervical c-1),Brain - Substantia nigra,Breast - Mammary Tissue,Cells - Cultured fibroblasts,Cells - EBV-transformed lymphocytes,Cervix - Ectocervix,Cervix - Endocervix,Colon - Sigmoid,Colon - Transverse,Esophagus - Gastroesophageal Junction,Esophagus - Mucosa,Esophagus - Muscularis,Fallopian Tube,Heart - Atrial Appendage,Heart - Left Ventricle,Kidney - Cortex,Kidney - Medulla,Liver,Lung,Minor Salivary Gland,Muscle - Skeletal,Nerve - Tibial,Ovary,Pancreas,Pituitary,Prostate,Skin - Not Sun Exposed (Suprapubic),Skin - Sun Exposed (Lower leg),Small Intestine - Terminal Ileum,Spleen,Stomach,Testis,Thyroid,Uterus,Vagina,Whole Blood\n")
        
        # Write csv for each interval
        for BedTool_Interval in GTEx_BedTool:
            # Check if file is in ucsc format
            if args.ucsc_format:
                GeneID, GeneExpr = extract_GTEx_data(BedTool_Interval, 6, 9)
            else:
                GeneID, GeneExpr = extract_GTEx_data(BedTool_Interval, 3, 11)
            
            f.write(f"{GeneID},{GeneExpr}\n")








if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog = '',
                    description = """"""   )
    # Required
    # Input and Output of the files:
    parser.add_argument('-f', '--gtex_file', required=True,
            help='give GTEx BED file path')
    parser.add_argument('-out', '--output', required=True,
            help='give output path')
    
    parser.add_argument('-ucsc', '--ucsc_format', action="store_true",
            help='specify file format. If -ucsc [TRUE], the file format of the GTEx ucsc file will be used. Otherwise the fileformat of the processed file will be used.')   

    
    args = parser.parse_args()

    main()