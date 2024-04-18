"""
This Script assign all TFBS in a BED file looking like this:
chr1	17510	17522	MYCN	0	-
chr1	629638	629650	OTX2	0	-
chr1	634195	634206	EGR2	0	+
chr1	758332	758348	JUN	    0	-

To a Genomic Region in another BED file (e.g. Promotor regions etc.) Important is that this file contains the genomic locations (chrom, start, end, strand) and the geneID.
z.b. 
chr1	11668	11868	ENSG00000223972.5	0	+

A special case is if a GTEx File is passed, which look like this:
chr1	11668	11868	DDX11L1	11	+	ENSG00000223972.5	transcribed_unprocessed_pseudogene	54	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.166,0,0,0,0,

In general the :::: Function have to be changed, in order to merge the regions together correctly, if other input files were used.

Step 1: 
- Loading tfbs BED file and Promotor BED file into a BedTool object
- Intersecting thesese files using intersect -wa -wb
        --> The genomic regions of the promotor remains. Every single TFBS intersecting with the promotor will be appended to the file as new columns.
    --> chr1	11668	11868	ENSG00000223972.5	0	+ ***prom_fields*** chr1	17510	17522	MYCN	0	- ***TFBS fields***

Step 2:
- Merging of the repetitive Promotor regions and remaining important informations.
    - Input of tfbs columns that should remain and associated operation
    - Input of prom columns that should remain and associated operation
    - Optional: Input of prom column that should be used to merge only region with same field content.
- Defining standarized columns and operations.


"""


import pybedtools
import numpy as np
import pandas as pd
import argparse

def unique_BedTool_field(BedTool, field_number):
    """
    function that get unique intervall names of given BedTool object. Return a list with unique names.
    """
    map_name = map(lambda x: x.fields[field_number], BedTool)
    names = np.array(list(map_name))
    names_unique = np.unique(names)
    return list(names_unique)



def merge_bed_considering_field(BedTool, field_number, c_for_merge, o_for_merge):
    """
    function that will use BedTool.merge on BedTool entrys with same BedTool.field[field_number]
    """

    BedTool_name_unique = unique_BedTool_field(BedTool, field_number)
    
    merged_tfbs_tmp = pybedtools.BedTool(())
    for tfbs_name in BedTool_name_unique:
        single_tfbs_merged = BedTool.filter(lambda x: x.fields[7] == tfbs_name).merge(s=True, c=c_for_merge, o=o_for_merge)
        merged_tfbs_tmp = merged_tfbs_tmp.cat(single_tfbs_merged, postmerge=False)
    
    merged_tfbs = merged_tfbs_tmp.sort()

    return merged_tfbs



def merge_bed(BedTool, c_for_merge, o_for_merge):
    return BedTool.merge(s=True, c=c_for_merge, o=o_for_merge).sort()




def generate_input_for_merge(prom, c_tfbs, c_prom, o_tfbs, o_prom):
    """
    Generating c_for_merge and o_for_merge by considering the input columns for the tfbs BED file and the prom BED file.
    First, the length of the promotor BED is considered to change the column index of the tfbs file to the new index position in the intersected file.
    Second the prom lists will be appended by the tfbs lists 
    """
    # Changing TFBS columns to new intersect columns
    prom_field_len = len(prom[0].fields)
    c_tfbs_new = [ci + prom_field_len for ci in c_tfbs]

    c_for_merge = c_prom + c_tfbs_new
    o_for_merge = o_prom + o_tfbs
    return c_for_merge, o_for_merge








def main():
    # Input files
    tfbs_all = pybedtools.BedTool(args.tfbs_file)
    prom = pybedtools.BedTool(args.intersect_file)

    # Intersect
    intersect= prom.intersect(tfbs_all, wa=True, wb=True).sort()

    # If GTEx file is used
    if args.gtex:
        c_for_merge=[4,5,6,7,8,9,10,  14,12,13,16, 18, 14, 14]
        o_for_merge=["distinct","distinct","distinct","distinct","distinct","distinct","distinct",  "collapse","collapse","collapse","collapse", "distinct", "count", "count_distinct"]
        merged_BedTool = merge_bed_considering_field(intersect,7)
    
    # Default
    else:
        c_for_merge, o_for_merge = generate_input_for_merge(prom, args.column_in_tfbs, args.column_in_prom, args.operation_for_tfbs, args.operation_for_prom)

        # Considering a column in the Promotor file to merge only Promotors with same content.
        if args.consider_field_in_prom != 0:
            merged_BedTool = merge_bed_considering_field(BedTool=intersect, field_number=args.consider_field_in_prom, c_for_merge=c_for_merge, o_for_merge=o_for_merge)

        else:
            merged_BedTool = merge_bed(BedTool=intersect, c_for_merge=c_for_merge, o_for_merge=o_for_merge)

    # Saving the intermediate BED Tool
    merged_BedTool.saveas(f"{args.output}/Prom_with_TFBSs_raw.bed")
    







if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog = 'This Script assigns all TFBS in a BED file to genomic locations of another BED file. Specific conditions for merging can be specified.',
                    description = """All TFBS in first file will be intersected with region in second file. The regions will be merged under user-defined conditions. 
                     By default specific columns of the tfbs and promotor remain intact. """   )
    # Required
    # Input and Output of the files:
    parser.add_argument('-f', '--tfbs_file', required=True,
            help='give TFBS BED file path')
    parser.add_argument('-fb', '--intersect_file', required=True,
            help='give BED file for intersection. Containing genomic regions e.g. promotors')
    parser.add_argument('-out', '--output', required=True,
            help='give output path')
    
    # Optional
    # Define if known GTEx file
    parser.add_argument('-gtex', '--gtex_as_input_file', required=False, action="store_true",
            help='If -gtex, then the -fb file will be considered as GTEx File and the input (see below) is not necessary to process it correctly.')  

    # Define if it should be merged by a specific column in Prom BED file
    parser.add_argument('-consider_field', '--consider_field_in_prom', required=False, default=0,
                help='Specify columns from the input prom file, where the content must be the same for the promoter to be merged. Caution: RunTime, The unique values in the column may not be too many. By default (0), no column will be considerd. ')

    # Input for specific columns and operation for input files, during merging.
    parser.add_argument('-c_tfbs', '--column_in_tfbs', required=False, default=[3, 1, 2, 5],
            help='Specify columns from the input tfbs file to operate upon (see -o_tfbs option, below). Multiple columns can be specified in a comma-delimited list. (see bedtool merge arguments)')   
    parser.add_argument('-o_tfbs', '--operation_for_tfbs', required=False, default=["collapse","collapse","collapse","collapse"],
            help='Specify the operation that should be applied to -c. Multiple operations can be specified in a comma-delimited list. (see bedtool merge arguments)')
    
    parser.add_argument('-c_prom', '--column_in_prom', required=False, default=[3, 4, 5],
            help='Specify columns from the input prom file to operate upon (see -o_prom option, below). Multiple columns can be specified in a comma-delimited list. (see bedtool merge arguments)')   
    parser.add_argument('-o_prom', '--operation_for_prom', required=False, default=["distinct", "distinct", "distinct"],
            help='Specify the operation that should be applied to -c_prom. Multiple operations can be specified in a comma-delimited list. (see bedtool merge arguments)')
    


    
    args = parser.parse_args()
    
    main()
