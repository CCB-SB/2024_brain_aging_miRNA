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

# # find which diff_exp_CA2 should be run
# we search for it in the modules
result_file_diff_exp = ""
for m in config['modules']:
    if m.startswith("diff_exp_CA2"):
        result_file_diff_exp = get_steps_folder(m)

##### target rules #####
rule all:
    input: result_files


#----------------- preprocessing stats -----------------
rule alignment_stats_bar_plot:
    input: mapping_info=config["mapping_info"]
    output: result=join(config["results_folder"], get_steps_folder("alignment_stats_bar_plot"))
    params: 
        data_input=config["common"].get("input", None),
        plots_props=config["common"].get("plots", None),
        parameters_props=config["common"].get("parameters", None),
        grouping=config["alignment_stats_bar_plot"].get("grouping", None),
        xticks_names=config["xticks_names"],
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/alignment_stats_bar_plot.R"

rule mirna_number_per_tissue_bar_plot:
    input:
    output: result=join(config["results_folder"], get_steps_folder("mirna_number_per_tissue_bar_plot"))
    params: 
        data_input=config["common"].get("input", None),
        tissue=config["mirna_number_per_tissue_bar_plot"].get("tissue", None),
        detection_rate_list=config["mirna_number_per_tissue_bar_plot"].get("detection_rate", None),
        min_detect_raw=config["mirna_number_per_tissue_bar_plot"].get("min_detect_raw", None),
        colors=config["colors"],
        xticks_names=config["xticks_names"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"],
    conda: "envs/analysis_packages.yml"
    script: "scripts/mirna_number_per_tissue_bar_plot.R"

#------------------ plots and tables ------------------
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
    
rule diff_exp_CA2_Experiment_diet_restriction:
    input:
    output: result=join(config["results_folder"], get_steps_folder("diff_exp_CA2_Experiment_diet_restriction"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        comparisons_supp_table=config["diff_exp_CA2_Experiment_diet_restriction"].get("comparisons", None),
        comparisons_volcano_brain_region=config["volcano_scatter_plot"].get("comparisons_brain_region_diet_restriction", None),
        log_list=config["diff_exp_CA2_Experiment_diet_restriction"].get("log_list", None),
        adjustment=config["diff_exp_CA2_Experiment_diet_restriction"]["adjustment"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/diff_exp.R"

rule diff_exp_CA2_Experiment_young_mouse_plasma_injection_type_retro_orbital:
    input:
    output: result=join(config["results_folder"], get_steps_folder("diff_exp_CA2_Experiment_young_mouse_plasma_injection_type_retro_orbital"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        comparisons_supp_table=config["diff_exp_CA2_Experiment_young_mouse_plasma_injection_type_retro_orbital"].get("comparisons", None),
        comparisons_volcano_brain_region=config["volcano_scatter_plot"].get("comparisons_brain_region_injection", None),
        log_list=config["diff_exp_CA2_Experiment_young_mouse_plasma_injection_type_retro_orbital"].get("log_list", None),
        adjustment=config["diff_exp_CA2_Experiment_young_mouse_plasma_injection_type_retro_orbital"]["adjustment"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/diff_exp.R"

rule volcano_CA2_Experiment_diet_restriction_scatter_plot:
    input: results_diff_exp=join(config["results_folder"], result_file_diff_exp)
    output: result=join(config["results_folder"], get_steps_folder("volcano_CA2_Experiment_diet_restriction_scatter_plot"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_log=config["volcano_scatter_plot"].get("diff_log", None),
        comparisons=config["diff_exp_CA2_Experiment_diet_restriction"].get("comparisons", None),
        comparisons_brain_region=config["volcano_scatter_plot"].get("comparisons_brain_region_diet_restriction", None),
        parameters=config["volcano_scatter_plot"].get("parameters", None),
        colors=config["colors"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/volcano_scatter_plot.R"

rule volcano_CA2_Experiment_young_mouse_plasma_injection_type_retro_orbital_scatter_plot:
    input: results_diff_exp=join(config["results_folder"], result_file_diff_exp)
    output: result=join(config["results_folder"], get_steps_folder("volcano_CA2_Experiment_young_mouse_plasma_injection_type_retro_orbital_scatter_plot"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_log=config["volcano_scatter_plot"].get("diff_log", None),
        comparisons=config["diff_exp_CA2_Experiment_young_mouse_plasma_injection_type_retro_orbital"].get("comparisons", None),
        comparisons_brain_region=config["volcano_scatter_plot"].get("comparisons_brain_region_injection", None),
        parameters=config["volcano_scatter_plot"].get("parameters", None),
        colors=config["colors"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/volcano_scatter_plot.R"

rule correlation_comp_age_heatmap_plot_table:
    input: feature_list=config["feature_list"]
    output: result_local=join(config["results_folder"], get_steps_folder("correlation_comp_age_heatmap_plot_table"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        method=config["correlation_comp_age_heatmap_plot_table"].get("method", None),
        file=config["correlation_comp_age_heatmap_plot_table"].get("file", None),
        props=config["correlation_comp_age_heatmap_plot_table"].get("props", None),
        adjustment=config["correlation_comp_age_heatmap_plot_table"].get("adjustment", None),
        corr_th=config["correlation_comp_age_heatmap_plot_table"].get("corr_th", None),
        colors=config["colors"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: 
        "scripts/correlation_comp_age_heatmap_plot_table.R"

rule pca_scatter_plot:
    input: mapping_info=config["mapping_info"]
    output: result=join(config["results_folder"], get_steps_folder("pca_scatter_plot"))
    params: 
        data_input=config["common"].get("input", None),
        parameters=config["pca_scatter_plot"].get("parameters", None),
        properties=config["pca_scatter_plot"].get("properties", None),
        plots_props=config["common"].get("plots", None),
        parameters_props=config["common"].get("parameters", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/pca_scatter_plot.R"

rule top_expressed_mirnas:
    input: 
    output: result=join(config["results_folder"], get_steps_folder("top_expressed_mirnas"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        results_folder=config["results_folder"],
        top_list=config["umap"].get("top_list", None),
    conda: "envs/analysis_packages.yml"
    script: "scripts/top_expressed_mirnas.R"

rule pvca_all_brain_regions_bar_plot:
    input: 
    output: result=join(config["results_folder"], get_steps_folder("pvca_all_brain_regions_bar_plot"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        props=config["pvca_all_brain_regions_bar_plot"].get("properties", None),
        parameters=config["pvca_all_brain_regions_bar_plot"].get("parameters", None),
        plots_props=config["common"].get("plots", None),
        xticks_names=config["xticks_names"],
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/pvca_bar_plot.R"

rule pvca_per_brain_regions_bar_plot:
    input: 
    output: result=join(config["results_folder"], get_steps_folder("pvca_per_brain_regions_bar_plot"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        props=config["pvca_per_brain_regions_bar_plot"].get("properties", None),
        parameters=config["pvca_per_brain_regions_bar_plot"].get("parameters", None),
        plot_parameters=config["pvca_per_brain_regions_bar_plot"].get("plot_parameters", None),
        plots_props=config["common"].get("plots", None),
        xticks_names=config["xticks_names"],
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/pvca_bar_plot.R"

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

rule mirna_expression_median_violin_plot:
   input: 
       result_local=join(config["results_folder"], get_steps_folder("correlation_comp_age_heatmap_plot_table")),
       results_diff_exp=join(config["results_folder"], result_file_diff_exp)
   output: result=join(config["results_folder"], get_steps_folder("mirna_expression_median_violin_plot"))
   params: 
       data_input=config["common"].get("input", None),
       files=config["mirna_expression_median_violin_plot"].get("files", None),
       features=config["mirna_expression_median_violin_plot"].get("features", None),
       comparisons=config["mirna_expression_median_violin_plot"].get("comparisons", None),
       parameters_props=config["mirna_expression_median_violin_plot"].get("parameters", None),
       manual_plots_parameters=config["mirna_expression_median_violin_plot"].get("manual_plots_parameters", None),
       plots_props=config["common"].get("plots", None),
       colors=config["colors"],
       xticks_names=config["xticks_names"],
       results_folder=config["results_folder"]
   conda: "envs/analysis_packages.yml"
   script: "scripts/mirna_expression_median_violin_plot.R"

rule target_genes_scatter_plot:
    input:
    output: result=join(config["results_folder"], get_steps_folder("target_genes_scatter_plot"))
    params:
        data_input=config["common"].get("input", None),
        feature=config["target_genes_scatter_plot"].get("feature", None),
        gene_list_folder_path=config["target_genes_scatter_plot"].get("gene_list_folder_path", None),
        split_prop=config["target_genes_scatter_plot"].get("split_prop", None),
        hline_th=config["target_genes_scatter_plot"].get("hline_th", None),
        xticks_names=config["xticks_names"],
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/target_genes_scatter_plot_CA2.R"

rule isomir_line_plot: 
    input: 
    output: result_local=join(config["results_folder"], get_steps_folder("isomir_line_plot"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        plot_group=config["isomir_line_plot"].get("plot_group", None),
        comp_group=config["isomir_line_plot"].get("comp_group", None),
        order_of_x_axis=config["isomir_line_plot"].get("order_of_x_axis", None),
        feature=config["isomir_line_plot"].get("feature", None),
        colours=config["isomir_line_plot"].get("colours", None),
        feature_col=config["isomir_line_plot"].get("feature_col", None),
        plot_layout=config["isomir_line_plot"].get("plot_layout", None),
        plot_size=config["isomir_line_plot"].get("plot_size", None),
        manual_plot=config["isomir_line_plot"].get("manual_plot", None),
        results_folder=config["results_folder"],
        xticks_names=config["xticks_names"],
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
    conda: "envs/analysis_packages.yml"
    script: 
        "scripts/isomir_line_plot.R"

rule isomir_heatmap_plot: 
    input: 
    output: result_local=join(config["results_folder"], get_steps_folder("isomir_heatmap_plot"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        plot_group=config["isomir_heatmap_plot"].get("plot_group", None),
        feature=config["isomir_heatmap_plot"].get("feature", None),
        top_list=config["isomir_heatmap_plot"].get("top_list", None),
        colours=config["isomir_heatmap_plot"].get("colours", None),
        feature_col=config["isomir_heatmap_plot"].get("feature_col", None),
        results_folder=config["results_folder"],
        xticks_names=config["xticks_names"],
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
    conda: "envs/analysis_packages.yml"
    script: 
        "scripts/isomir_heatmap_plot.R"
