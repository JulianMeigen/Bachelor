'''
This script needs genomic locations of genes (BED/GTF etc.) as input and returns the promotor region, defined as x bp upstream of the gene (depends on strand if gene).
'''

import pybedtools
import argparse


def get_promotor(BedTool, prom_len):
    promotors_bed = BedTool.flank(genome= "hg38", l=prom_len, r=0, s=True).saveas(f"{args.output}/Promotor_{prom_len}bp.bed")
    return promotors_bed

def main():
    genes = pybedtools.BedTool(args.filename)
    get_promotor(genes, args.promotor_length)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog = 'Generating BED file with Promotor Rewgion of given genomic locations in input file.',
                    description = """Flanking genomic locations in BED file with BedTool function 'flank -g "hg38" -l 200 -r 0 -s' """   )
    # Required
    # Input and Output of the files:
    parser.add_argument('-f', '--filename', required=True,
            help='give path for BED file with genomic locations of genes.')
    parser.add_argument('-out', '--output', required=True,
            help='give output path for folder.')
    # Optional
    parser.add_argument('-len', '--promotor_length', required=False, default=200
            help='give length of the Promotor region. By default 200bp.')
    
    main()