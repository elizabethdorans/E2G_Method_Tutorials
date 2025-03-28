import pandas as pd
from glob import glob
import os
import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--scent_predictions_dir",
                    help = "Folder containing *.tsv files with SCENT predictions")
parser.add_argument("--output_file",
                    help = "Name of file to save concatenated SCENT predictions")

args = parser.parse_args()

scent_predictions_dir = args.scent_predictions_dir
output_file = args.output_file

out_dir = os.path.dirname(output_file)
if not os.path.exists(out_dir):
    print("Creating %s!" % out_dir)
    os.makedirs(out_dir)

pgl_files = glob("%s/*.tsv" % scent_predictions_dir)

merged = pd.DataFrame(columns = ["gene", "peak", "beta", "se", "z", "p", "boot_basic_p"])
for file in pgl_files:
    with open(file) as f:
        tmp = pd.read_csv(file, sep = "\t")
        merged = pd.concat([merged, tmp])
        
merged = merged.drop(columns = ["se", "z", "p"])
merged = merged.drop_duplicates()

merged.to_csv(output_file, sep = "\t", index = False)
print(merged)
print(output_file)