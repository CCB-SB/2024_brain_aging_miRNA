from os.path import join
from pathlib import Path

from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("5.4.0")


##### load config #####

configfile: "config_preprocessing.yaml"

result_files = [join(config["results_folder"], "steps/{}.done".format(m)) for m in config['modules'] if m != "common"]

##### target rules #####
rule all:
    input: result_files

#-------------------- preprocessing --------------------
rule preprocessing:
    input: 
        join(config["results_folder"], "steps/rpmmm_norm.done"),
        join(config["results_folder"], "steps/log_transform.done"),
        join(config["results_folder"], "steps/create_raw_detection_matrix.done"),
        join(config["results_folder"], "steps/filter_samples_number_alignment.done"),
        join(config["results_folder"], "steps/filter_expression_files_per_group.done"),
    output: join(config["results_folder"], "steps/preprocessing.done")
    shell:
        "touch {output};"

rule rpmm_norm_mirna:
    input: mapping_info=config["mapping_info"],
    output: result=join(config["results_folder"], "steps/rpmm_norm.done")
    params: 
        data_input=config["common"].get("input", None),
        id_col=config["preprocessing"].get("id_col", None),
    run:
        import pandas as pd
        data_input = params.data_input
        for rna_class in data_input["rna_class_list"]:    
            #rna_class = "miRNA"
            if rna_class != "miRNA":
                rna_class_key = "RNA"
            else:
                rna_class_key = "miRNA"
            #tbl = pd.read_csv(f"{data_input['preprocesing_data_folder']}_{data_input['data_sub_set']}/quantification/{rna_class}_complete_quantification_raw.csv", sep='\t', index_col=["precursor", "miRNA"])
            # get the second part of the data sub_set
            if data_input['data_sub_set'] in ("microglia", "ROSMAP"):
                suffix = data_input['data_sub_set']
            else:
                suffix = "_".join(data_input['data_sub_set'].split("_")[1:])
            tbl = pd.read_csv(f"{data_input['preprocesing_data_folder']}_{data_input['data_sub_set']}/{rna_class}_expression_raw_{suffix}.tsv", sep='\t', index_col=rna_class_key)
            # merge rows by miRNA name
            #tbl = tbl.groupby(rna_class).mean()
            mapping = pd.read_csv(input.mapping_info, sep='\t')
            mapping = mapping[mapping["Mismatches"] == 1]
            # merge annot sample names to mapping
            annot = pd.read_csv(f"{data_input['annotation_data_folder']}/annotation_{data_input['data_sub_set']}_filtered.csv", sep='\t')
            mapping = pd.merge(mapping, annot[[params.id_col, "fastq_name"]], on="fastq_name", suffixes=("", "_drop"))
            mapping.set_index(params.id_col, inplace=True)
            res = tbl / mapping["reads_aligned"] * 1e6
            # change colnames
            #annot = pd.read_csv(f"{data_input['preprocesing_data_folder']}_{data_input['data_sub_set']}/annotation_{data_input['data_sub_set']}.csv", sep='\t', index_col="Sample")
            #sample_2_name = annot.to_dict()["Name"]
            # delete precurser col and change colnames
            #cols_2_keep = [c for c in res.columns if c in sample_2_name or c in annot.index or c in {"miRNA", "RNA", "precursor"}]
            #res = res[cols_2_keep]
            #res.rename(columns=sample_2_name, inplace=True)
            res.to_csv(f"{data_input['preprocesing_data_folder']}_{data_input['data_sub_set']}/{rna_class}_complete_quantification_rpmm_norm.csv", sep='\t', index=True)
        with open(output[0], "w") as file:
            file.write("")

rule rpmmm_norm:
    input: 
    output: result=join(config["results_folder"], "steps/rpmmm_norm.done")
    params: data_input=config["common"].get("input", None),
    run:
        import pandas as pd
        from pathlib import Path
        data_input = params.data_input
        for rna_class in data_input["rna_class_list"]:    
            # load expression
            if rna_class != "miRNA":
                continue
            else:
                rna_class_key = "miRNA"
            # get the second part of the data sub_set
            if data_input['data_sub_set'] in ("microglia", "ROSMAP"):
                suffix = data_input['data_sub_set']
            else:
                suffix = "_".join(data_input['data_sub_set'].split("_")[1:])
            tbl = pd.read_csv(f"{data_input['preprocesing_data_folder']}_{data_input['data_sub_set']}/{rna_class}_expression_raw_{suffix}.tsv", sep='\t')
            tbl.iloc[:,1:] = tbl.iloc[:,1:] / tbl.groupby(rna_class_key).mean().sum() * 1e6
            tbl.to_csv(f"{data_input['preprocesing_data_folder']}_{data_input['data_sub_set']}/{rna_class}_complete_quantification_rpmmm_norm.csv", sep='\t', index=False)
        with open(output[0], "w") as file:
            file.write("")

rule log_transform:
    input: 
        result_rpmm=join(config["results_folder"], "steps/rpmm_norm.done"),
        result_rpmmm=join(config["results_folder"], "steps/rpmmm_norm.done")
    output: result=join(config["results_folder"], "steps/log_transform.done")
    params: 
        data_input=config["common"].get("input", None),
        log_list = config["preprocessing"].get("log_list", None),
    conda: "envs/analysis_packages.yml"
    script: "scripts/log_transform.R"

rule create_raw_detection_matrix:
    input: result=join(config["results_folder"], "steps/log_transform.done")
    output: result=join(config["results_folder"], "steps/create_raw_detection_matrix.done")
    params: 
        data_input=config["common"].get("input", None),
        min_detection_list = config["preprocessing"].get("min_detection_list", None),
    conda: "envs/analysis_packages.yml"
    script: "scripts/create_raw_detection_matrix.R"

rule filter_samples_number_alignment:
    input: 
        mapping_info=config["mapping_info"],
        result=join(config["results_folder"], "steps/create_raw_detection_matrix.done"),
    output: result=join(config["results_folder"], "steps/filter_samples_number_alignment.done")
    params: 
        data_input=config["common"].get("input", None),
        number_of_reads=config["preprocessing"].get("number_of_reads", None),
        id_col=config["preprocessing"].get("id_col", None),
        fastq_filename_col=config["preprocessing"].get("fastq_filename_col", None),
    conda: "envs/analysis_packages.yml"
    script: "scripts/filter_samples_number_alignment.R"

rule filter_expression_files_per_group:
    input: result=join(config["results_folder"], "steps/filter_samples_number_alignment.done")
    output: result=join(config["results_folder"], "steps/filter_expression_files_per_group.done")
    params:
        data_input=config["common"].get("input", None),
        min_detection_list = config["preprocessing"].get("min_detection_list", None),
        detection_rate_list=config["preprocessing"].get("detection_rate_list", None),
        filter_variable=config["preprocessing"].get("filter_variable", None),
    conda: "envs/analysis_packages.yml"
    script: "scripts/filter_expression_files_per_group.R"


