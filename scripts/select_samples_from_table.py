import pandas as pd
import numpy as np
from pathlib import Path
import argparse

def select_samples_from_table(data_set, source_folder_parent, target_folder_parent, target_annot_path, rna_types=["miRNA"]):

    # add columns from small to big
    # the source folder (we will process all csv files of the fiven rna types)
    #source_folder_parent = "/local/23-04_brain_aging_sncrna_data_analysis/data_CA1/without_mid/sncrna-pipeline-input_CA1/"
    #source_folder_parent = "/local/23-04_brain_aging_sncrna_data_analysis_aengel/data_external/ROSMAP_COMPASS_MIRNAS_human/all_feature_filtering_disease_status/"
    ##source_folder_parent = "/local/23-04_brain_aging_sncrna_data_analysis/data_TMS/sncrna-pipeline-input_TMS_brain/"

    # path to save folder where new files will be saved
    #target_folder_parent = "/local/23-04_brain_aging_sncrna_data_analysis_aengel/data_external/ROSMAP_COMPASS_MIRNAS_human/all_healthy/"
    #target_folder_parent = "/local/23-04_brain_aging_sncrna_data_analysis/data_TMS/sncrna-pipeline-input_TMS_brain_age=3m_12m_15m/"

    # parent folder of expression folder
    #annotation_folder_parent = "/local/23-04_brain_aging_sncrna_data_analysis_aengel/data_external/ROSMAP_COMPASS_MIRNAS_human/all_healthy/"
    #annotation_folder_parent = "/local/23-04_brain_aging_sncrna_data_analysis/data_TMS/sncrna-pipeline-input_TMS_brain_age=3m_12m_15m/"

    # which RNA types should be considered, set to False if everything should be considered
    #rna_types = ["circRNA", "lncRNA", "piRNA", "rRNA", "scaRNA", "snoRNA", "snRNA", "tRNA", "ncRNA"]
    #rna_types = ["ncRNA", "tRNA"]
    #rna_types = ["miRNA"]

    annotation_folder_parent = target_annot_path

    source_folder_parent_path = Path(source_folder_parent)
    target_folder_parent_path = Path(target_folder_parent)
    annotation_path_parent = Path(annotation_folder_parent)

    # load annotation
    annotation_filename = f"annotation_{data_set}_filtered.csv"
    annotation_path = annotation_path_parent / annotation_filename
    print(annotation_path)
    annotation = pd.read_csv(annotation_path, sep="\t")

    # find the columns to copy from annotation file
    copy_columns = annotation["SampleID"].to_list()
    print("++++++++++++++++++++++++")
    print(f"number of columns to copy {len(copy_columns)}")
    # print(f"columns to copy {copy_columns}")

    # build source folder path
    source_folder_path = source_folder_parent_path
    target_folder_path = target_folder_parent_path

    # create target target_folder_path if necessary
    target_folder_path.mkdir(exist_ok=True)

    print(source_folder_path)

    # do this for all csv files
    for source_file_path in source_folder_path.glob("*.csv"): #_complete_quantification_raw.csv
        # ignore hidden files (starting with a ".")
        if str(source_file_path.name).startswith("."):
            continue
        # only process some RNA types
        if rna_types != False:
            valid_type = False
            for rna_type in rna_types:
                if str(source_file_path.name).startswith(rna_type):
                    valid_type = True
            if valid_type == False:
                continue
        # make path to target file
        target_file_path = target_folder_path / source_file_path.name
        # if we want to change the name:
        # target_file_path = target_folder_path / (source_file_path.stem + "_without_age=28m.csv")

        print("------------------------")
        print(f"source file {source_file_path}")
        # print(f"target file {target_file_path}")

        # load both with 
        df_source = pd.read_csv(source_file_path, sep="\t")

        # get the name of the first columns (miRNA or RNA)
        first_col_name = df_source.columns.to_list()[0]
        print(f"name of the first column {first_col_name}")

        # copy columns from complete dataframe
        df_removed = df_source[[first_col_name] + copy_columns]
        print(f"shape of source dataframe {df_source.shape}")
        print(f"shape of dataframe with removed columns {df_removed.shape}")

        # save
        print(f"saving df_removed to {target_file_path}")
        df_removed.to_csv(target_file_path, sep="\t", na_rep="NA", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process 4 input arguments")
    parser.add_argument("data_set", help="First input: data_set")
    parser.add_argument("source_folder_parent", help="Second input: source_folder_parent")
    parser.add_argument("target_folder_parent", help="Third input: target_folder_parent")
    parser.add_argument("target_annot_path", help="Fourth input: target_annot_path")
    parser.add_argument("rna_types", nargs="*", default=["miRNA"], help="RNA classes")

    args = parser.parse_args()

    select_samples_from_table(
        args.data_set,
        args.source_folder_parent,
        args.target_folder_parent,
        args.target_annot_path,
        args.rna_types
        )
