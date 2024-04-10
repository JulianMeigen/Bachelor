import pybedtools
import numpy as np
import argparse

def reduce_bed_by_intersect_wa(bed_file_path, intersect_file_path, output_path):
    
    # Load input Bed file into BedTools object
    raw_data_bed = pybedtools.BedTool(bed_file_path)

    intersect_data = pybedtools.BedTool(intersect_file_path)

    reduced_data = raw_data_bed.intersect(intersect_data, wa=True).saveas(output_path)

    print("Output generated.\n")
    return reduced_data




if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog = 'Reduce BED file to common chromosomes',
                    description = 'Reduce given BED file to specific chromosomes. In default chr1-22 and chrX and chrY will be saved in a new BED file output.',)
    # Required
    # Input and Output of the files:
    parser.add_argument('-fa', '--filename', required=True,
                    help='give BED file path')
    parser.add_argument('-fb', '--intersect', required=True,
            help='give BED file for intersection')
    
    parser.add_argument('-out', '--output', required=True,
                help='give output path for new BED file')
    # Optional
    args = parser.parse_args()
    
    reduce_bed_by_intersect_wa(bed_file_path=args.filename, intersect_file_path=args.intersect, output_path=args.output)