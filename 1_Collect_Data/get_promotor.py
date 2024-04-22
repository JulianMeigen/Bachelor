"""
This script gets Genes and returns the genomic locations of the promotor. The promotor is defined as a region 200bp upstream of the TSS.
The length of the promotor can be changed.
"""

import pybedtools
import argparse

def get_promotor(BedTool, prom_len):
    """
    Function that gets BedTool as input and returns the region "prom_len" upstream of these regions.
    """
    promotor_BedTool = BedTool.flank(genome= "hg38", l=prom_len, r=0, s=True).saveas()
    return promotor_BedTool

def main():
    # input files
    genes_BedTool = pybedtools.BedTool(args.filename)
    promotor_BedTool = get_promotor(genes_BedTool, args.promotor_length)

    promotor_BedTool.saveas(f"{args.output}/Promotors_{args.promotor_length}.bed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog = 'This script gets Genes and returns the genomic locations of the promotor.',
                    description = """Uses BedTools 'flank -g "hg38" -l 200 -s' as input. The length of the promotor region can be changed, by default its 200bp. """   )
    
    parser.add_argument('-f', '--filename', required=True,
            help='give BED file path')
    parser.add_argument('-out', '--output', required=True,
            help='give output path for folder where new BED file will be saved.')
    
    parser.add_argument('-l', '--promotor_length', required=False, default= 200,
            help='give length of Promotor region in bp. By default 200 bp.')
    
    args = parser.parse_args()

    main()
