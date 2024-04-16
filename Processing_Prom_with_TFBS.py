
"""
This Script will assign to each region in a given BED file all TFBS from another BED file.

Step 1: Assign TFBS to Promotor region by using BedTools "intersect -wa -wb" method

Step 2: Merging Promotor regions together and changing important columns

"""



import pybedtools
import numpy as np
import pandas as pd
import argparse