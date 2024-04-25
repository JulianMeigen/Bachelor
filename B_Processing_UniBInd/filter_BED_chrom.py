import pybedtools
import numpy as np
import argparse

def filter_BED_chrom(bed_file_path, output_path, lst_chrom):
    
    # Load input Bed file into BedTools object
    raw_data = pybedtools.BedTool(bed_file_path)

    # Filter the raw data and only save chr1-22 and chrX and chrY in BED file

    chr_ordered = np.array(lst_chrom)

    data = raw_data.filter(lambda x: np.isin(x.chrom, chr_ordered)).saveas(output_path)
    print("Output generated.\n")
    return data




if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog = 'Reduce BED file to common chromosomes',
                    description = 'Reduce given BED file to specific chromosomes. In default chr1-22 and chrX and chrY will be saved in a new BED file output.',)
    # Required
    # Input and Output of the files:
    parser.add_argument('-f', '--filename', required=True,
                    help='give BED file path')
    parser.add_argument('-out', '--output', required=True,
                help='give output path for new BED file')
    # Optional
    # List of chromosomes
    parser.add_argument('-chrom', '--chromosomes', required=False, default=["chr1", "chr2", "chr3", "chr4","chr5","chr6","chr7","chr8", "chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19", "chr20","chr21","chr22","chrX","chrY"],
                    help='give a list of chromosomes to which the BED file should be filtered. The chromosomes should be in BED file format such as chr1, chrX. On defalut chr1-22, chrX and chrY will be used.')

    args = parser.parse_args()
    
    filter_BED_chrom(bed_file_path=args.filename, output_path=args.output, lst_chrom=args.chromosomes)
