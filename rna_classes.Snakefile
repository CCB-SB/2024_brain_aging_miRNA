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

##### target rules #####
rule all:
    input: result_files


#------------------ plots and tables ------------------
rule mapping_stats_bar_plot:
    input: 
        mapping_info_rna_classes=config["mapping_info_rna_classes"],
        mapping_info_rna_classes_microglia=config["mapping_info_rna_classes_microglia"],
    output: result=join(config["results_folder"], get_steps_folder("mapping_stats_bar_plot"))
    params: 
        data_input=config["common"].get("input", None),
        plots_props=config["common"].get("plots", None),
        parameters_props=config["common"].get("parameters", None),
        grouping=config["mapping_stats_bar_plot"].get("grouping", None),
        image_size=config["mapping_stats_bar_plot"].get("image_size", None),
        xticks_names=config["xticks_names"],
        xticks_names_microglia=config["xticks_names_microglia"],
        colors=config["colors"],
        colors_microglia=config["colors_microglia"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/mapping_stats_bar_plot.R"

rule expressed_features_per_rna_class_bar_plot:
    input:
    output: result=join(config["results_folder"], get_steps_folder("expressed_features_per_rna_class_bar_plot"))
    params: 
        data_input=config["common"].get("input", None),
        grouping=config["expressed_features_per_rna_class_bar_plot"].get("grouping", None),
        min_detect_raw=config["expressed_features_per_rna_class_bar_plot"].get("min_detect_raw", None),
        detection_rate=config["expressed_features_per_rna_class_bar_plot"].get("detection_rate", None),
        colors=config["colors"],
        xticks_names=config["xticks_names"],
        plots_props=config["common"].get("plots", None),
        image_size=config["expressed_features_per_rna_class_bar_plot"].get("image_size", None),
        results_folder=config["results_folder"],
    conda: "envs/analysis_packages.yml"
    script: "scripts/expressed_features_per_rna_class_bar_plot.R"

# rule metadata_overview_tables:
#     input:
#     output: result=join(config["results_folder"], get_steps_folder("metadata_overview_tables"))
#     params: 
#         data_input=config["common"].get("input", None),
#         parameters_props=config["common"].get("parameters", None),
#         combinations=config["metadata_overview_tables"].get("combinations", None),
#         results_folder=config["results_folder"],
#     conda: "envs/analysis_packages.yml"
#     script: "scripts/metadata_overview_tables.R"

rule expressed_features_per_rna_class_line_plot:
    input: 
    output: result=join(config["results_folder"], get_steps_folder("expressed_features_per_rna_class_line_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        parameters=config["expressed_features_per_rna_class_line_plot"].get("parameters", None),
        manual_plots_parameters=config["expressed_features_per_rna_class_line_plot"].get("manual_plots_parameters", None),
        xticks_names=config["xticks_names"],
        colors=config["colors"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/expressed_features_per_rna_class_line_plot.R"

rule correlation_comp_age_per_rna_class_ridgeline_plot:
    input: 
    output: join(config["results_folder"], get_steps_folder("correlation_comp_age_per_rna_class_ridgeline_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        prop_attribute=config["correlation_comp_age_per_rna_class_ridgeline_plot"].get("props", None),
        corr_methods=config["correlation_comp_age_per_rna_class_ridgeline_plot"].get("method", None),
        corr_th=config["correlation_comp_age_per_rna_class_ridgeline_plot"].get("corr_th", None),
        p_val_adj=config["correlation_comp_age_per_rna_class_ridgeline_plot"].get("adjustment", None),
        adj_p_value=config["correlation_comp_age_per_rna_class_ridgeline_plot"].get("adj_p_value", None),
        plot_size=config["correlation_comp_age_per_rna_class_ridgeline_plot"].get("plot_size", None),
        xticks_names=config["xticks_names"],
        colors=config["colors"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/correlation_comp_age_per_rna_class_ridgeline_plot.R"

rule merge_expression_files_per_rna_classes:
    input: 
    output: result=join(config["results_folder"],get_steps_folder("merge_expression_files_per_rna_classes"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        results_folder=config["results_folder"],
    conda: "envs/analysis_packages.yml"
    script: "scripts/merge_expression_files_per_rna_classes.R"

rule umap_analysis:
    input: result=join(config["results_folder"], get_steps_folder("merge_expression_files_per_rna_classes"))
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
    
rule make_publication_folder:
   input:
   output: result=join(config["results_folder"], get_steps_folder("make_publication_folder"))
   params:
   conda: "envs/analysis_packages.yml"
   script: "scripts/extract_paper_figures.py"
