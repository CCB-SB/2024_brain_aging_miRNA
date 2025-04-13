from os.path import join

from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("5.4.0")


##### load config #####

configfile: "config.yaml"

def get_steps_folder(step_name):
    step_folder_path = f"{config['common'].get('input', None)['rna_class']}_{config['common'].get('input', None)['detection_rate']}/results_{config['common'].get('input', None)['data_sub_set']}/steps/{step_name}.done"
    return(step_folder_path)

#result_files = [join(config["results_folder"], "steps/{}.done".format(m)) for m in config['modules'] if m != "common"]
result_files = [join(config["results_folder"], get_steps_folder(m)) for m in config['modules'] if m != "common"]

# find which diff_exp_CA1 should be run
# we search for it in the modules
result_file_diff_exp = ""
for m in config['modules']:
    if m.startswith("diff_exp_CA1"):
        result_file_diff_exp = get_steps_folder(m)

##### target rules #####
rule all:
    input: result_files


#------------------ plots and tables ------------------
rule create_binned_human_age:
    input:
    output: result=join(config["results_folder"], get_steps_folder("create_binned_human_age"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        property_to_be_binned=config["create_binned_human_age"].get("property_to_be_binned", None),
        number_of_bins=config["create_binned_human_age"].get("number_of_bins", None),
        results_folder=config["results_folder"],
    conda: "envs/analysis_packages.yml"
    script: "scripts/create_binned_human_age.R"

rule metadata_overview_tables:
    input:
    output: result=join(config["results_folder"], get_steps_folder("metadata_overview_tables"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        combinations=config["metadata_overview_tables"].get("combinations", None),
        results_folder=config["results_folder"],
    conda: "envs/analysis_packages.yml"
    script: "scripts/metadata_overview_tables.R"
    
rule diff_exp_human:
    input: result=join(config["results_folder"], get_steps_folder("create_binned_human_age"))
    output: result=join(config["results_folder"], get_steps_folder("diff_exp_human"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        comparisons_supp_table=config["diff_exp_human"].get("comparisons", None),
        comparisons_volcano_age=config["volcano_scatter_plot"].get("comparisons_human", None),
        log_list=config["diff_exp_human"].get("log_list", None),
        adjustment=config["diff_exp_human"]["adjustment"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/diff_exp.R"


rule deregulation_mirna_human_heatmap_plot:
    input: result_local=join(config["results_folder"], get_steps_folder("diff_exp_human"))
    output: result_local=join(config["results_folder"], get_steps_folder("deregulation_mirna_human_heatmap_plot"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        bin_params=config["deregulation_mirna_human_heatmap_plot"].get("bin_params", None),
        comparisons=config["deregulation_mirna_human_heatmap_plot"].get("comparisons", None),
        params=config["deregulation_mirna_human_heatmap_plot"].get("params", None),
        plots_props=config["common"].get("plots", None),
        xticks_names=config["xticks_names"],
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/deregulation_mirna_human_heatmap_plot.R"

rule correlation_mirna_human_with_age_heatmap_plot:
    input: 
    output: result_local=join(config["results_folder"], get_steps_folder("correlation_mirna_human_with_age_heatmap_plot"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        feature=config["correlation_mirna_human_with_age_heatmap_plot"].get("feature", None),
        miRNA_human_path=config["correlation_mirna_human_with_age_heatmap_plot"].get("miRNA_human_path", None),
        split_prop=config["correlation_mirna_human_with_age_heatmap_plot"].get("split_prop", None),
        corr_prop=config["correlation_mirna_human_with_age_heatmap_plot"].get("corr_prop", None),
        corr_params=config["correlation_mirna_human_with_age_heatmap_plot"].get("corr_params", None),
        plots_props=config["common"].get("plots", None),
        xticks_names=config["xticks_names"],
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/correlation_mirna_human_with_age_heatmap_plot.R"

rule top_expressed_mirnas:
    input: 
    output: result=join(config["results_folder"],get_steps_folder("top_expressed_mirnas"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        results_folder=config["results_folder"],
        top_list=config["umap"].get("top_list", None),
    conda: "envs/analysis_packages.yml"
    script: "scripts/top_expressed_mirnas.R"

rule umap_analysis:
    input: result=join(config["results_folder"], get_steps_folder("top_expressed_mirnas"))
    output: result=join(config["results_folder"], get_steps_folder("umap_analysis"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        results_folder=config["results_folder"],
        top_list=config["umap"].get("top_list", None),
        n_neighbors=config["umap"].get("neighbors", None),
        min_dist=config["umap"].get("min_dist", None),
        init=config["umap"].get("init", None),
        metric=config["umap"].get("metric", None),
        norm=config["umap"].get("norm", None),
        seed=config["umap"].get("seed", None),
    conda: "envs/umap_analysis.yml"
    script: "scripts/umap_analysis.py"

rule umap:
    input: result=join(config["results_folder"], get_steps_folder("umap_analysis"))
    output: result=join(config["results_folder"], get_steps_folder("umap"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        results_folder=config["results_folder"],
        props=config["umap"].get("properties", None),
        top_list=config["umap"].get("top_list", None),
        n_neighbors=config["umap"].get("neighbors", None),
        min_dist=config["umap"].get("min_dist", None),
        init=config["umap"].get("init", None),
        metric=config["umap"].get("metric", None),
        norm=config["umap"].get("norm", None),
        seed=config["umap"].get("seed", None),
        xticks_names=config["xticks_names"],
        colors=config["colors"],
        plots_props=config["common"].get("plots", None),
    conda: "envs/analysis_packages.yml"
    script: "scripts/umap_scatter_plot.R"

rule volcano_human_scatter_plot_manual:
    input: result=join(config["results_folder"], get_steps_folder("diff_exp_human"))
    output: result=join(config["results_folder"], get_steps_folder("volcano_human_scatter_plot_manual"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_log=config["volcano_human_scatter_plot_manual"].get("diff_log", None),
        comparisons=config["volcano_human_scatter_plot_manual"].get("comparisons_human", None),
        parameters=config["volcano_human_scatter_plot_manual"].get("parameters", None),
        manual_plots_parameters=config["volcano_human_scatter_plot_manual"].get("manual_plots_parameters", None),
        colors=config["colors"],
        xticks_names=config["xticks_names"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/volcano_scatter_plot_manual.R"

rule volcano_human_scatter_plot:
    input: result=join(config["results_folder"], get_steps_folder("diff_exp_human"))
    output: result=join(config["results_folder"], get_steps_folder("volcano_human_scatter_plot"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_log=config["volcano_scatter_plot"].get("diff_log", None),
        comparisons_age=config["volcano_scatter_plot"].get("comparisons_human", None),
        parameters=config["volcano_scatter_plot"].get("parameters", None),
        colors=config["colors"],
        xticks_names=config["xticks_names"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/volcano_scatter_plot.R"

rule make_publication_folder:
   input:
   output: result=join(config["results_folder"], get_steps_folder("make_publication_folder"))
   params:
   conda: "envs/analysis_packages.yml"
   script: "scripts/extract_paper_figures.py"
