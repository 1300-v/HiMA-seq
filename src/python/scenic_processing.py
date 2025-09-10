# conda activate pyscenic
import os
import sys
import loompy as lp
import numpy as np
import scanpy as sc
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True, help="input file")
ap.add_argument("-o", "--output", required=True, help="output file")
args = vars(ap.parse_args())

input_file = args["input"]  # counts
output_file = args["output"]  # sample.loom

counts = sc.read_csv(input_file)
row_attrs = {"Gene": np.array(counts.var_names),}
col_attrs = {"CellID": np.array(counts.obs_names)}
lp.create(output_file, counts.X.transpose(), row_attrs, col_attrs)
