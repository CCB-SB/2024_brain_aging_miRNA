from os.path import join

from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("5.4.0")


##### load config #####

configfile: "config.yaml"

def get_steps_folder(step_name):
    step_folder_path = f"{config['common'].get('input', None)['rna_class']}_{config['common'].get('input', None)['detection_rate']}/results_{config['common'].get('input', None)['data_sub_set_male']}_{config['common'].get('input', None)['data_sub_set_female']}/steps/{step_name}.done"
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

#-------------------- preprocessing --------------------
# no preproccing only one preprocessing for total dataset

#----------------- preprocessing stats -----------------


#------------------ plots and tables ------------------
rule diff_exp_analysis:
    input: 
        # result_local=join(config["results_folder"], get_steps_folder("aggregation_diff_exp_bar_plot")),
        #result_bar=join(config["results_folder"], get_steps_folder("aggregation_compact_diff_exp_bar_plot")),
        result_dots=join(config["results_folder"], get_steps_folder("global_sig_deregulated_mirnas_dot_plot")),
        result_heatmap=join(config["results_folder"], get_steps_folder("deregulated_comparison_tissue_heatmap"))
    output: join(config["results_folder"], get_steps_folder("diff_exp_analysis"))
    shell:
        "touch {output};"

rule deregulated_comparison_tissue_heatmap_male_female:
    input: 
    output: result=join(config["results_folder"], get_steps_folder("deregulated_comparison_tissue_heatmap_male_female"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_exp_log=config["deregulated_comparison_tissue_heatmap_male_female"].get("diff_log", None),
        thresholds=config["deregulated_comparison_tissue_heatmap_male_female"].get("thresholds", None),
        props=config["deregulated_comparison_tissue_heatmap_male_female"].get("props", None),
        xticks_names=config["xticks_names"],
        colors=config["colors"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: 
        "scripts/deregulated_comparison_tissue_heatmap_male_female.R"

#rule aggregation_diff_exp_bar_plot:
#    input: results_diff_exp=join(config["results_folder"], result_file_diff_exp)
#    output: result=join(config["results_folder"], get_steps_folder("aggregation_diff_exp_bar_plot"))
#    params:
#        data_input=config["common"].get("input", None),
#        parameters_props=config["common"].get("parameters", None),
#        diff_exp_log=config["aggregation_diff_exp_bar_plot"].get("diff_log", None),
#        thresholds=config["aggregation_diff_exp_bar_plot"].get("thresholds", None),
#        props=config["aggregation_diff_exp_bar_plot"].get("props", None),
#        plots_props=config["common"].get("plots", None),
#        colors=config["colors"],
#        results_folder=config["results_folder"]
#    conda: "envs/analysis_packages.yml"
#    script: "scripts/aggregation_diff_exp_bar_plot.R"

#rule aggregation_compact_diff_exp_bar_plot:
#    input: results_diff_exp=join(config["results_folder"], result_file_diff_exp)
#    output: result=join(config["results_folder"], get_steps_folder("aggregation_compact_diff_exp_bar_plot"))
#    params:
#        data_input=config["common"].get("input", None),
#        parameters_props=config["common"].get("parameters", None),
#        diff_exp_log=config["aggregation_diff_exp_bar_plot"].get("diff_log", None),
#        thresholds=config["aggregation_diff_exp_bar_plot"].get("thresholds", None),
#        props=config["aggregation_diff_exp_bar_plot"].get("props", None),
#        plots_props=config["common"].get("plots", None),
#        colors=config["colors"],
#        xticks_names=config["xticks_names"],
#        results_folder=config["results_folder"]
#    conda: "envs/analysis_packages.yml"
#    script: "scripts/aggregation_compact_diff_exp_bar_plot.R"

#rule global_sig_deregulated_mirnas_dot_plot:
#    input: result=join(config["results_folder"], get_steps_folder("aggregation_compact_diff_exp_bar_plot"))
#    output: result=join(config["results_folder"], get_steps_folder("global_sig_deregulated_mirnas_dot_plot"))
#    params:
#        data_input=config["common"].get("input", None),
#        parameters_props=config["common"].get("parameters", None),
#        diff_exp_log=config["aggregation_diff_exp_bar_plot"].get("diff_log", None),
#        thresholds=config["aggregation_diff_exp_bar_plot"].get("thresholds", None),
#        props=config["aggregation_diff_exp_bar_plot"].get("props", None),
#        plots_props=config["common"].get("plots", None),
#        colors=config["colors"],
#        results_folder=config["results_folder"]
#    conda: "envs/analysis_packages.yml"
#    script: "scripts/global_sig_deregulated_mirnas_dot_plot.R"

rule correlation_comp_age_heatmap_plot_table_male_female:
    input: feature_list=config["feature_list"]
    output: result_local=join(config["results_folder"], get_steps_folder("correlation_comp_age_heatmap_plot_table_male_female"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        method=config["correlation_comp_age_heatmap_plot_table_male_female"].get("method", None),
        props=config["correlation_comp_age_heatmap_plot_table_male_female"].get("props", None),
        adjustment=config["correlation_comp_age_heatmap_plot_table_male_female"].get("adjustment", None),
        corr_th=config["correlation_comp_age_heatmap_plot_table_male_female"].get("corr_th", None),
        xticks_names=config["xticks_names"],
        colors=config["colors"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: 
        "scripts/correlation_comp_age_heatmap_plot_table_male_female.R"

rule candidate_comparison_venn_plot_cluster_split_male_female:
    input: 
    output: result=join(config["results_folder"],get_steps_folder("candidate_comparison_venn_plot_cluster_split_male_female"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_exp_props=config["candidate_comparison_venn_plot_cluster_split_male_female"].get("diff_exp_props", None),
        corr_props=config["candidate_comparison_venn_plot_cluster_split_male_female"].get("corr_props", None),
        hclust_props=config["candidate_comparison_venn_plot_cluster_split_male_female"].get("hclust_props", None),
        props=config["candidate_comparison_venn_plot_cluster_split_male_female"].get("props", None),
        tissue_specific_fig_size=config["candidate_comparison_venn_plot_cluster_split_male_female"].get("tissue_specific", None),
        all_features_fig_size=config["candidate_comparison_venn_plot_cluster_split_male_female"].get("all_features", None),
        xticks_names=config["xticks_names"],
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/candidate_comparison_venn_plot_cluster_split_male_female.R"

rule candidate_comparison_scatter_plot_tissue_split_male_female:
    input: 
    output: result=join(config["results_folder"],get_steps_folder("candidate_comparison_scatter_plot_tissue_split_male_female"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_exp_props=config["candidate_comparison_scatter_plot_tissue_split_male_female"].get("diff_exp_props", None),
        corr_props=config["candidate_comparison_scatter_plot_tissue_split_male_female"].get("corr_props", None),
        hclust_props=config["candidate_comparison_scatter_plot_tissue_split_male_female"].get("hclust_props", None),
        props=config["candidate_comparison_scatter_plot_tissue_split_male_female"].get("props", None),
        xticks_names=config["xticks_names"],
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/candidate_comparison_scatter_plot_tissue_split_male_female.R"

rule mirna_expression_median_bar_box_plot:
    input: 
    output: result=join(config["results_folder"], get_steps_folder("mirna_expression_median_bar_box_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        files=config["mirna_expression_median_bar_box_plot"].get("files", None),
        features=config["mirna_expression_median_bar_box_plot"].get("features", None),
        parameters=config["mirna_expression_median_bar_box_plot"].get("parameters", None),
        xticks_names=config["xticks_names"],
        colors=config["colors"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/mirna_expression_median_bar_box_plot.R"

rule mirna_expression_over_time_box_bar_plot:
    input: 
        results_diff_exp=join(config["results_folder"], result_file_diff_exp),
        feature_list=config["feature_list"]
    output: result=join(config["results_folder"], get_steps_folder("mirna_expression_over_time_box_bar_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        files=config["mirna_expression_over_time_box_bar_plot"].get("files", None),
        features=config["mirna_expression_over_time_box_bar_plot"].get("features", None),
        parameters=config["mirna_expression_over_time_box_bar_plot"].get("parameters", None),
        manual_plots_parameters=config["mirna_expression_over_time_box_bar_plot"].get("manual_plots_parameters", None),
        corr_methods=config["correlation_comp_age_heatmap_plot_table_male_female"].get("method", None),
        xticks_names=config["xticks_names"],
        colors=config["colors"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/mirna_expression_over_time_box_bar_plot.R"

rule mirna_expression_over_time_box_bar_plot_male_female:
    input: 
        results_diff_exp=join(config["results_folder"], result_file_diff_exp),
        feature_list=config["feature_list"]
    output: result=join(config["results_folder"], get_steps_folder("mirna_expression_over_time_box_bar_plot_male_female"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        features=config["mirna_expression_over_time_box_bar_plot"].get("features", None),
        parameters=config["mirna_expression_over_time_box_bar_plot"].get("parameters", None),
        manual_plots_parameters=config["mirna_expression_over_time_box_bar_plot"].get("manual_plots_parameters", None),
        xticks_names=config["xticks_names"],
        colors=config["colors"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/mirna_expression_over_time_box_bar_plot_male_female.R"

rule mirna_expression_per_time_over_brain_region_box_bar_plot:
    input: results_diff_exp=join(config["results_folder"], result_file_diff_exp)
    output: result=join(config["results_folder"], get_steps_folder("mirna_expression_per_time_over_brain_region_box_bar_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        files=config["mirna_expression_per_time_over_brain_region_box_bar_plot"].get("files", None),
        features=config["mirna_expression_per_time_over_brain_region_box_bar_plot"].get("features", None),
        parameters=config["mirna_expression_per_time_over_brain_region_box_bar_plot"].get("parameters", None),
        xticks_names=config["xticks_names"],
        colors=config["colors"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/mirna_expression_per_time_over_brain_region_box_bar_plot.R"

#rule coefficient_of_variation:
#    input:
#    output: result=join(config["results_folder"], get_steps_folder("coefficient_of_variation"))
#    params:
#        data_input=config["common"].get("input", None),
#        parameters_props=config["common"].get("parameters", None),
#        tissues=config["coefficient_of_variation"].get("tissues", None),
#        top_list=config["coefficient_of_variation"].get("top_list", None),
#        results_folder=config["results_folder"]
#    conda: "envs/analysis_packages.yml"
#    script: "scripts/coefficient_of_variation.R"

#rule tissue_specificity_index:
#    input:
#    output: result=join(config["results_folder"], get_steps_folder("tissue_specificity_index"))
#    params:
#        data_input=config["common"].get("input", None),
#        parameters_props=config["common"].get("parameters", None),        
#        tissues=config["tissue_specificity_index"].get("tissues", None),
#        top_list=config["tissue_specificity_index"].get("top_list", None),
#        results_folder=config["results_folder"]
#    conda: "envs/analysis_packages.yml"
#    script: "scripts/tissue_specificity_index.R"

rule hclust_expression_brain_region_for_mirna_list_zscores_heatmap_plot_male_female:
    input: 
    output: result=join(config["results_folder"], get_steps_folder("hclust_expression_brain_region_for_mirna_list_zscores_heatmap_plot_male_female"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        prop=config["hclust_expression_brain_region_for_mirna_list_zscores_heatmap_plot_male_female"].get("prop", None),
        diff_exp_log=config["hclust_expression_brain_region_for_mirna_list_zscores_heatmap_plot_male_female"].get("diff_log", None),
        zscore_th=config["hclust_expression_brain_region_for_mirna_list_zscores_heatmap_plot_male_female"].get("zscore_th", None),
        method=config["hclust_expression_brain_region_for_mirna_list_zscores_heatmap_plot_male_female"].get("method", None),
        top_list=config["hclust_expression_brain_region_for_mirna_list_zscores_heatmap_plot_male_female"].get("top_list", None),
        params=config["hclust_expression_brain_region_for_mirna_list_zscores_heatmap_plot_male_female"].get("hclust_params", None),
        xticks_names=config["xticks_names"],
        plot_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/hclust_expression_brain_region_for_mirna_list_zscores_heatmap_plot_male_female.R"

rule hclust_expression_brain_region_for_mirna_list_zscores_venn_plot_male_female:
   input: 
   output: result=join(config["results_folder"], get_steps_folder("hclust_expression_brain_region_for_mirna_list_zscores_venn_plot_male_female"))
   params:
       data_input=config["common"].get("input", None),
       parameters_props=config["common"].get("parameters", None),
       prop=config["hclust_expression_brain_region_for_mirna_list_zscores_venn_plot_male_female"].get("prop", None),
       method=config["hclust_expression_brain_region_for_mirna_list_zscores_venn_plot_male_female"].get("method", None),
       top_list=config["hclust_expression_brain_region_for_mirna_list_zscores_venn_plot_male_female"].get("top_list", None),
       zscore_th=config["hclust_expression_brain_region_for_mirna_list_zscores_venn_plot_male_female"].get("zscore_th", None),
       xticks_names=config["xticks_names"],
       plot_props=config["common"].get("plots", None),
       colors=config["colors"],
       results_folder=config["results_folder"]
   conda: "envs/analysis_packages.yml"
   script: "scripts/hclust_expression_brain_region_for_mirna_list_zscores_venn_plot_male_female.R"
     
rule make_publication_folder:
   input:
   output: result=join(config["results_folder"], get_steps_folder("make_publication_folder"))
   params:
   conda: "envs/analysis_packages.yml"
   script: "scripts/extract_paper_figures.py"
