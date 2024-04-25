"""
This script filters a BED file by a specific content in a specific column. By default it will filter for the protein_coding GeneType in a GTEx BED file.
"""

import pybedtools
import argparse

def filter_bed(BedTool, column_number, filter_content):
    """
    Function that filters a bed file by a specific content in a specific column.
    """
    filtered_BedTool = BedTool.filter(lambda x: x.fields[column_number] == filter_content)
    return filtered_BedTool

    
    

def main():
    # input files
    BedTool = pybedtools.BedTool(args.filename)
    filtered_BedTool = filter_bed(BedTool, args.column_number, args.filter_content)

    filtered_BedTool.saveas(f"{args.output}/Filtered_promotors_by_{args.filter_content}.bed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser('This script filters a BED file by a specific content in a specific column. By default it will filter for the protein_coding GeneType in a GTEx BED file.',
                    description = """Uses pybedtools '.filter(FUNC) as input. FUNC is a function that is True if .fields[column_number] == filter_content. """   )
    
    parser.add_argument('-f', '--filename', required=True,
            help='give BED file path')
    parser.add_argument('-out', '--output', required=True,
            help='give output path for folder where new BED file will be saved.')
    
    parser.add_argument('-c', '--column_number', required=False, default=7,
            help='give column (field) number in which the "filter_content" input is. By default its 7, because in GTEx files it the column for the geneType.')

    parser.add_argument('-filtered_by', '--filter_content', required=False, default="protein_coding",
            help='Specify the content to be filtered. By default it is "protein_coding", because in GTEx files that indicates protein coding genes.')
     
    args = parser.parse_args()

    main()
