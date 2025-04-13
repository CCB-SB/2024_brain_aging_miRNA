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

rule diff_exp:
    input:
    output: result=join(config["results_folder"], get_steps_folder("diff_exp"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        comparisons_supp_table=config["diff_exp"].get("comparisons", None),
        comparisons=config["volcano_scatter_plot"].get("comparisons", None),
        log_list=config["diff_exp"].get("log_list", None),
        adjustment=config["diff_exp"]["adjustment"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/diff_exp.R"

rule top_expressed_mirnas:
    input: 
    output: result=join(config["results_folder"],get_steps_folder("top_expressed_mirnas"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        results_folder=config["results_folder"],
        top_list=config["top_expressed_mirnas"].get("top_list", None),
    conda: "envs/analysis_packages.yml"
    script: "scripts/top_expressed_mirnas.R"

rule top_expressed_mirnas_heatmap_plot:
    input: result=join(config["results_folder"], get_steps_folder("top_expressed_mirnas"))
    output: result=join(config["results_folder"], get_steps_folder("top_expressed_mirnas_heatmap_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        top_list=config["top_expressed_mirnas_heatmap_plot"].get("top_list", None),
        sample_order=config["top_expressed_mirnas_heatmap_plot"].get("sample_order", None),
        properties=config["top_expressed_mirnas_heatmap_plot"].get("properties", None),
        xticks_names=config["xticks_names"],
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"],
    conda: "envs/analysis_packages.yml"
    script: "scripts/top_expressed_mirnas_heatmap_plot.R"

rule coefficient_of_variation:
    input:
    output: result=join(config["results_folder"], get_steps_folder("coefficient_of_variation"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        tissues=config["coefficient_of_variation"].get("tissues", None),
        top_list=config["coefficient_of_variation"].get("top_list", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/coefficient_of_variation.R"

rule hclust_expression_zscores_heatmap_plot:
    input: result=join(config["results_folder"], get_steps_folder("coefficient_of_variation"))
    output: result=join(config["results_folder"], get_steps_folder("hclust_expression_zscores_heatmap_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        top_list=config["hclust_expression_zscores_heatmap_plot"].get("top_list", None),
        sample_order=config["hclust_expression_zscores_heatmap_plot"].get("sample_order", None),
        properties=config["hclust_expression_zscores_heatmap_plot"].get("properties", None),
        diff_log=config["hclust_expression_zscores_heatmap_plot"].get("diff_log", None),
        xticks_names=config["xticks_names"],
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"],
    conda: "envs/analysis_packages.yml"
    script: "scripts/hclust_expression_zscores_heatmap_plot.R"

rule pca_scatter_plot:
    input: mapping_info=config["mapping_info"]
    output: result=join(config["results_folder"], get_steps_folder("pca_scatter_plot"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        parameters=config["pca_scatter_plot"].get("parameters", None),
        properties=config["pca_scatter_plot"].get("properties", None),
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/pca_scatter_plot.R"

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

rule mirna_expression_median_per_group_line_plot:
   input: result=join(config["results_folder"], get_steps_folder("diff_exp"))
   output: result=join(config["results_folder"], get_steps_folder("mirna_expression_median_per_group_line_plot"))
   params: 
       data_input=config["common"].get("input", None),
       diff_log=config["mirna_expression_median_per_group_line_plot"].get("diff_log", None),
       properties=config["mirna_expression_median_per_group_line_plot"].get("properties", None),
       features_black=config["mirna_expression_median_per_group_line_plot"].get("features_black", None),
       features_grey=config["mirna_expression_median_per_group_line_plot"].get("features_grey", None),
       plots_props=config["common"].get("plots", None),
       colors=config["colors"],
       xticks_names=config["xticks_names"],
       results_folder=config["results_folder"]
   conda: "envs/analysis_packages.yml"
   script: "scripts/mirna_expression_median_per_group_line_plot.R"

rule fc_versus_expr_scatter_plot:
    input: result=join(config["results_folder"], get_steps_folder("diff_exp"))
    output: result=join(config["results_folder"], get_steps_folder("fc_versus_expr_scatter_plot"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_exp_log=config["fc_versus_expr_scatter_plot"].get("diff_log", None),
        prop=config["fc_versus_expr_scatter_plot"].get("prop", None),
        comparisons=config["fc_versus_expr_scatter_plot"].get("comparisons", None),
        thresholds=config["fc_versus_expr_scatter_plot"].get("thresholds", None),
        xticks_names=config["xticks_names"],
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/fc_versus_expr_scatter_plot.R" 

rule volcano_scatter_plot:
    input: result=join(config["results_folder"], get_steps_folder("diff_exp"))
    output: result=join(config["results_folder"], get_steps_folder("volcano_scatter_plot"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_log=config["fc_versus_expr_scatter_plot"].get("diff_log", None),
        comparisons=config["volcano_scatter_plot"].get("comparisons", None),
        parameters=config["volcano_scatter_plot"].get("parameters", None),
        colors=config["colors"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/volcano_scatter_plot_microglia.R"
     
rule fc_versus_expr_scatter_plot_CA1:
    input: result=join(config["results_folder"], get_steps_folder("diff_exp_CA1"))
    output: result=join(config["results_folder"], get_steps_folder("fc_versus_expr_scatter_plot_CA1"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_exp_log=config["fc_versus_expr_scatter_plot"].get("diff_log", None),
        prop=config["fc_versus_expr_scatter_plot"].get("prop", None),
        comparisons_CA1=config["fc_versus_expr_scatter_plot"].get("comparisons_CA1", None),
        thresholds=config["fc_versus_expr_scatter_plot"].get("thresholds", None),
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/fc_versus_expr_scatter_plot.R" 
    
rule comparing_overlapping_feature_sets:
    input: result=join(config["results_folder"], get_steps_folder("mirna_number_per_tissue_bar_plot"))
    output: result=join(config["results_folder"], get_steps_folder("comparing_overlapping_feature_sets"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        sets=config["comparing_overlapping_feature_sets"].get("sets", None),
        colour_set=config["comparing_overlapping_feature_sets"].get("colour_set", None),
        plots_props=config["common"].get("plots", None),
        xticks_names=config["xticks_names"],
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/comparing_overlapping_feature_sets.R" 

rule make_publication_folder:
   input:
   output: result=join(config["results_folder"], get_steps_folder("make_publication_folder"))
   params:
   conda: "envs/analysis_packages.yml"
   script: "scripts/extract_paper_figures.py"

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
