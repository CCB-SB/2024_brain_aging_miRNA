import pandas as pd
import sys

# Read command line arguments
isomir_expr_path = sys.argv[1]
mirna_expr_path = sys.argv[2]
mapping_info_path = sys.argv[3]
annot_path = sys.argv[4]
output_path = sys.argv[5]

# Load the tables
isomir_expr = pd.read_csv(isomir_expr_path, sep="\t", index_col=["miRNA", "iso_type"])
mirna_expr = pd.read_csv(mirna_expr_path, sep="\t")
annot = pd.read_csv(annot_path, sep='\t')
mapping = pd.read_csv(mapping_info_path, sep='\t')
mapping = mapping[mapping["Mismatches"] == 1]

# Merge annot sample names to mapping
mapping = pd.merge(mapping, annot[["SampleID", "fastq_name"]], on="fastq_name", suffixes=("", "_drop"))
mapping.set_index("SampleID", inplace=True)

# RPMM normalization
isomir_expr_rpmm = isomir_expr / mapping["reads_aligned"] * 1e6

# reset multi-index
isomir_expr_rpmm.reset_index(inplace=True)

# Get list of miRNAs from mirna_expr
mirnas_to_keep = mirna_expr["miRNA"].tolist()

# Filter rows in isomir_expr where miRNA is in the list
isomir_expr_rpmm_filtered = isomir_expr_rpmm[isomir_expr_rpmm["miRNA"].isin(mirnas_to_keep)]

# Keep only the 'miRNA' column and any other columns that also appear in mirna_expr
columns_to_keep = ["miRNA", "iso_type"] + [col for col in isomir_expr_rpmm_filtered.columns if (col in mirna_expr.columns) and (col not in["miRNA", "iso_type"])]
isomir_expr_rpmm_filtered = isomir_expr_rpmm_filtered[columns_to_keep]

# Save to output
isomir_expr_rpmm_filtered.to_csv(output_path, sep="\t", index=False)