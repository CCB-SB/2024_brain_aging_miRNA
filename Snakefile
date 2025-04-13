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


#----------------- preprocessing stats -----------------
rule alignment_stats_bar_plot:
    input: mapping_info=config["mapping_info"]
    output: result=join(config["results_folder"], get_steps_folder("alignment_stats_bar_plot"))
    params: 
        data_input=config["common"].get("input", None),
        plots_props=config["common"].get("plots", None),
        parameters_props=config["common"].get("parameters", None),
        grouping=config["alignment_stats_bar_plot"].get("grouping", None),
        number_of_reads=config["alignment_stats_bar_plot"].get("number_of_reads", None),
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

rule diff_exp_CA1:
    input:
    output: result=join(config["results_folder"], get_steps_folder("diff_exp_CA1"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        comparisons_supp_table=config["diff_exp_CA1"].get("comparisons", None),
        comparisons_volcano_age=config["volcano_scatter_plot"].get("comparisons_age", None),
        comparisons_volcano_brain_region=config["volcano_scatter_plot"].get("comparisons_brain_region", None),
        log_list=config["diff_exp_CA1"].get("log_list", None),
        adjustment=config["diff_exp_CA1"]["adjustment"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/diff_exp.R"

rule diff_exp_CA1_age_3m_12m_15m:
    input:
    output: result=join(config["results_folder"], get_steps_folder("diff_exp_CA1_age_3m_12m_15m"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        comparisons_heatmap_fc=config["diff_exp_CA1_age_3m_12m_15m"].get("comparisons", None),
        comparisons_fc_scatter=config["fc_versus_expr_scatter_plot"].get("comparisons_CA1_age_3m_12m_15m", None),
        comparisons_volcano_brain_region=config["volcano_scatter_plot"].get("comparisons_brain_region", None),
        log_list=config["diff_exp_CA1"].get("log_list", None),
        adjustment=config["diff_exp_CA1"]["adjustment"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/diff_exp.R"

rule diff_exp_CA1_without_age_26m_28m:
    input:
    output: result=join(config["results_folder"], get_steps_folder("diff_exp_CA1_without_age_26m_28m"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        comparisons_volcano_sex=config["volcano_scatter_plot"].get("comparisons_sex", None),
        comparisons_sex_for_br_age=config["volcano_scatter_plot"].get("comparisons_sex_for_br_age", None),
        comparisons_volcano_brain_region=config["volcano_scatter_plot"].get("comparisons_brain_region", None),
        comparisons_volcano_age=config["volcano_scatter_plot"].get("comparisons_age_without_age_26m_28m", None),
        comparisons_fc_vs_median_expr=config["fc_versus_expr_scatter_plot"].get("comparisons_CA1_without_age_26m_28m", None),
        log_list=config["diff_exp_CA1"].get("log_list", None),
        adjustment=config["diff_exp_CA1"]["adjustment"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/diff_exp.R"

rule diff_exp_CA1_age_3m_12m_15m_TMS_brain_age_3m_12m_15m:
    input:
    output: result=join(config["results_folder"], get_steps_folder("diff_exp_CA1_age_3m_12m_15m_TMS_brain_age_3m_12m_15m"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        comparisons_heatmap_fc=config["diff_exp_CA1_age_3m_12m_15m"].get("comparisons", None),
        comparisons_volcano_brain_region_tms=config["volcano_scatter_plot"].get("comparisons_brain_region_tms", None),
        log_list=config["diff_exp_CA1"].get("log_list", None),
        adjustment=config["diff_exp_CA1"]["adjustment"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/diff_exp.R"

#rule diff_exp_values_table:
#    input: 
#        result=join(config["results_folder"], get_steps_folder("diff_exp")),
#        annot=config["annot"],
#        feature_list=config["feature_list"]
#    output: result=join(config["results_folder"], get_steps_folder("diff_exp_values_table"))
#    params: 
#        parameters_props=config["common"].get("parameters", None),
#        diff_exp_logs=config["diff_exp_values_table"].get("diff_exp_logs", None),
#        comparisons=config["diff_exp_values_table"].get("comparisons", None),
#        results_folder=config["results_folder"]
#    conda: "envs/analysis_packages.yml"
#    script: "scripts/diff_exp_values_table.R"

rule diff_exp_analysis:
    input: 
        # result_local=join(config["results_folder"], get_steps_folder("aggregation_diff_exp_bar_plot")),
        result_bar=join(config["results_folder"], get_steps_folder("aggregation_compact_diff_exp_bar_plot")),
        result_dots=join(config["results_folder"], get_steps_folder("global_sig_deregulated_mirnas_dot_plot")),
        result_heatmap=join(config["results_folder"], get_steps_folder("deregulated_comparison_tissue_heatmap"))
    output: join(config["results_folder"], get_steps_folder("diff_exp_analysis"))
    shell:
        "touch {output};"

rule deregulated_comparison_tissue_heatmap:
    input: results_diff_exp=join(config["results_folder"], result_file_diff_exp)
    output: result=join(config["results_folder"], get_steps_folder("deregulated_comparison_tissue_heatmap"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_exp_log=config["aggregation_diff_exp_bar_plot"].get("diff_log", None),
        thresholds=config["aggregation_diff_exp_bar_plot"].get("thresholds", None),
        props=config["aggregation_diff_exp_bar_plot"].get("props", None),
        xticks_names=config["xticks_names"],
        colors=config["colors"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: 
        "scripts/deregulated_comparison_tissue_heatmap.R"

rule aggregation_diff_exp_bar_plot:
    input: results_diff_exp=join(config["results_folder"], result_file_diff_exp)
    output: result=join(config["results_folder"], get_steps_folder("aggregation_diff_exp_bar_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_exp_log=config["aggregation_diff_exp_bar_plot"].get("diff_log", None),
        thresholds=config["aggregation_diff_exp_bar_plot"].get("thresholds", None),
        props=config["aggregation_diff_exp_bar_plot"].get("props", None),
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/aggregation_diff_exp_bar_plot.R"

rule aggregation_compact_diff_exp_bar_plot:
    input: results_diff_exp=join(config["results_folder"], result_file_diff_exp)
    output: result=join(config["results_folder"], get_steps_folder("aggregation_compact_diff_exp_bar_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_exp_log=config["aggregation_diff_exp_bar_plot"].get("diff_log", None),
        thresholds=config["aggregation_diff_exp_bar_plot"].get("thresholds", None),
        props=config["aggregation_diff_exp_bar_plot"].get("props", None),
        plots_props=config["common"].get("plots", None),
        xticks_names=config["xticks_names"],
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/aggregation_compact_diff_exp_bar_plot.R"

rule global_sig_deregulated_mirnas_dot_plot:
    input: result=join(config["results_folder"], get_steps_folder("aggregation_compact_diff_exp_bar_plot"))
    output: result=join(config["results_folder"], get_steps_folder("global_sig_deregulated_mirnas_dot_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_exp_log=config["aggregation_diff_exp_bar_plot"].get("diff_log", None),
        thresholds=config["aggregation_diff_exp_bar_plot"].get("thresholds", None),
        props=config["aggregation_diff_exp_bar_plot"].get("props", None),
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/global_sig_deregulated_mirnas_dot_plot.R"

rule deregulated_comparison_sex_tissue_heatmap:
    input: results_diff_exp=join(config["results_folder"], result_file_diff_exp)
    output: result=join(config["results_folder"], get_steps_folder("deregulated_comparison_sex_tissue_heatmap"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_exp_log=config["deregulated_comparison_sex_tissue_heatmap"].get("diff_log", None),
        thresholds=config["deregulated_comparison_sex_tissue_heatmap"].get("thresholds", None),
        props=config["deregulated_comparison_sex_tissue_heatmap"].get("props", None),
        xticks_names=config["xticks_names"],
        colors=config["colors"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: 
        "scripts/deregulated_comparison_sex_tissue_heatmap.R"

rule correlation:
    input: 
        result_table=join(config["results_folder"], get_steps_folder("correlation_comp_age_heatmap_plot_table")),
        # result_heatmap=join(config["results_folder"], get_steps_folder("correlation_comp_age_heatmap")),
        result_dots=join(config["results_folder"], get_steps_folder("correlation_comp_age_dot_plot")),
        # result=join(config["results_folder"], get_steps_folder("aggregation_correlation_bar_plot")),
        result_aggr=join(config["results_folder"], get_steps_folder("aggregation_compact_correlation_bar_plot"))
    output: join(config["results_folder"], get_steps_folder("correlation"))
    shell:
        "touch {output};"

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
        xticks_names=config["xticks_names"],
        colors=config["colors"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: 
        "scripts/correlation_comp_age_heatmap_plot_table.R"

rule correlation_comp_age_heatmap:
    input: result_local=join(config["results_folder"], get_steps_folder("correlation_comp_age_heatmap_plot_table"))
    output: result=join(config["results_folder"], get_steps_folder("correlation_comp_age_heatmap"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        thresholds=config["aggregation_correlation_bar_plot"].get("thresholds", None),
        method=config["correlation_comp_age_heatmap_plot_table"].get("method", None),
        props=config["correlation_comp_age_heatmap_plot_table"].get("props", None),
        adjustment=config["correlation_comp_age_heatmap_plot_table"].get("adjustment", None),
        colors=config["colors"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: 
        "scripts/correlation_comp_age_heatmap.R"

rule correlation_comp_age_dot_plot:
    input: result_local=join(config["results_folder"], get_steps_folder("aggregation_compact_correlation_bar_plot"))
    output: result=join(config["results_folder"], get_steps_folder("correlation_comp_age_dot_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        thresholds=config["aggregation_correlation_bar_plot"].get("thresholds", None),
        props=config["aggregation_correlation_bar_plot"].get("props", None),
        method=config["aggregation_correlation_bar_plot"].get("method", None),
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: 
        "scripts/correlation_comp_age_dot_plot.R"

rule aggregation_correlation_bar_plot:
    input: result=join(config["results_folder"], get_steps_folder("global_significant_correlated_mirna_with_age"))
    output: result=join(config["results_folder"], get_steps_folder("aggregation_correlation_bar_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        thresholds=config["aggregation_correlation_bar_plot"].get("thresholds", None),
        props=config["aggregation_correlation_bar_plot"].get("props", None),
        method=config["aggregation_correlation_bar_plot"].get("method", None),
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/aggregation_correlation_bar_plot.R"

rule aggregation_compact_correlation_bar_plot:
    input: result=join(config["results_folder"], get_steps_folder("correlation_comp_age_heatmap_plot_table"))
    output: result=join(config["results_folder"], get_steps_folder("aggregation_compact_correlation_bar_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        thresholds=config["aggregation_correlation_bar_plot"].get("thresholds", None),
        props=config["aggregation_correlation_bar_plot"].get("props", None),
        method=config["aggregation_correlation_bar_plot"].get("method", None),
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        xticks_names=config["xticks_names"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/aggregation_compact_correlation_bar_plot.R"

rule sig_tissue_specific_comparison:
    input: 
        result_diff_exp=join(config["results_folder"], get_steps_folder("aggregation_compact_diff_exp_bar_plot")),
        result_corr=join(config["results_folder"], get_steps_folder("aggregation_compact_correlation_bar_plot"))
    output: result=join(config["results_folder"],get_steps_folder("sig_tissue_specific_comparison"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_exp_props=config["sig_tissue_specific_comparison"].get("diff_exp_props", None),
        corr_props=config["sig_tissue_specific_comparison"].get("corr_props", None),
        props=config["sig_tissue_specific_comparison"].get("props", None),
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/sig_tissue_specific_comparison.R"

rule correlation_age_table:
    input:
    output: result_local=join(config["results_folder"], get_steps_folder("correlation_age_table"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        method=config["correlation_age_table"].get("method", None),
        adjustment=config["correlation_comp_age_heatmap_plot_table"].get("adjustment", None),
        prop=config["correlation_age_table"].get("prop", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: 
        "scripts/correlation_age_table.R"

rule mfuzz_ts_analysis_for_sig_tissue_specific_miRNAs_line_plot:
    input: result=join(config["results_folder"], get_steps_folder("sig_tissue_specific_comparison"))
    output: result=join(config["results_folder"], get_steps_folder("mfuzz_ts_analysis_for_sig_tissue_specific_miRNAs_line_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        cluster_props=config["mfuzz_ts_analysis_for_sig_tissue_specific_miRNAs_line_plot"].get("cluster_props", None),
        plot_category=config["mfuzz_ts_analysis_for_sig_tissue_specific_miRNAs_line_plot"].get("plot_category", None),
        diff_exp_props=config["mfuzz_ts_analysis_for_sig_tissue_specific_miRNAs_line_plot"].get("diff_exp_props", None),
        corr_props=config["mfuzz_ts_analysis_for_sig_tissue_specific_miRNAs_line_plot"].get("corr_props", None),
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/mfuzz_ts_analysis_for_sig_tissue_specific_miRNAs_line_plot.R"

rule sig_tissue_specific_miRNAs_line_plot:
    input: result=join(config["results_folder"], get_steps_folder("sig_tissue_specific_comparison"))
    output: result=join(config["results_folder"], get_steps_folder("sig_tissue_specific_miRNAs_line_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        plot_category=config["sig_tissue_specific_miRNAs_line_plot"].get("plot_category", None),
        diff_exp_props=config["sig_tissue_specific_miRNAs_line_plot"].get("diff_exp_props", None),
        corr_props=config["sig_tissue_specific_miRNAs_line_plot"].get("corr_props", None),
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/sig_tissue_specific_miRNAs_line_plot.R"

rule candidate_comparison_venn_plot:
    input: 
        result_diff_exp=join(config["results_folder"], get_steps_folder("aggregation_compact_diff_exp_bar_plot")),
        result_corr=join(config["results_folder"], get_steps_folder("aggregation_compact_correlation_bar_plot"))
    output: result=join(config["results_folder"],get_steps_folder("candidate_comparison_venn_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_exp_props=config["candidate_comparison_venn_plot"].get("diff_exp_props", None),
        corr_props=config["candidate_comparison_venn_plot"].get("corr_props", None),
        hclust_props=config["candidate_comparison_venn_plot"].get("hclust_props", None),
        props=config["candidate_comparison_venn_plot"].get("props", None),
        tissue_specific_fig_size=config["candidate_comparison_venn_plot"].get("tissue_specific", None),
        all_features_fig_size=config["candidate_comparison_venn_plot"].get("all_features", None),
        xticks_names=config["xticks_names"],
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/candidate_comparison_venn_plot.R"

rule candidate_comparison_venn_plot_cluster_split:
    input: 
        result_diff_exp=join(config["results_folder"], get_steps_folder("aggregation_compact_diff_exp_bar_plot")),
        result_corr=join(config["results_folder"], get_steps_folder("aggregation_compact_correlation_bar_plot"))
    output: result=join(config["results_folder"],get_steps_folder("candidate_comparison_venn_plot_cluster_split"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_exp_props=config["candidate_comparison_venn_plot_cluster_split"].get("diff_exp_props", None),
        corr_props=config["candidate_comparison_venn_plot_cluster_split"].get("corr_props", None),
        hclust_props=config["candidate_comparison_venn_plot_cluster_split"].get("hclust_props", None),
        props=config["candidate_comparison_venn_plot_cluster_split"].get("props", None),
        plot_layout=config["candidate_comparison_venn_plot_cluster_split"].get("plot_layout", None),
        xticks_names=config["xticks_names"],
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/candidate_comparison_venn_plot_cluster_split.R"


rule candidate_comparison_venn_plot_tissue_split:
    input: 
        result_diff_exp=join(config["results_folder"], get_steps_folder("aggregation_compact_diff_exp_bar_plot")),
        result_corr=join(config["results_folder"], get_steps_folder("aggregation_compact_correlation_bar_plot"))
    output: result=join(config["results_folder"],get_steps_folder("candidate_comparison_venn_plot_tissue_split"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_exp_props=config["candidate_comparison_venn_plot_tissue_split"].get("diff_exp_props", None),
        corr_props=config["candidate_comparison_venn_plot_tissue_split"].get("corr_props", None),
        hclust_props=config["candidate_comparison_venn_plot_tissue_split"].get("hclust_props", None),
        props=config["candidate_comparison_venn_plot_tissue_split"].get("props", None),
        plot_layout=config["candidate_comparison_venn_plot_tissue_split"].get("plot_layout", None),
        xticks_names=config["xticks_names"],
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/candidate_comparison_venn_plot_tissue_split.R"

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
    input: results_diff_exp=join(config["results_folder"], result_file_diff_exp)
    output: result=join(config["results_folder"], get_steps_folder("mirna_expression_over_time_box_bar_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        files=config["mirna_expression_over_time_box_bar_plot"].get("files", None),
        features=config["mirna_expression_over_time_box_bar_plot"].get("features", None),
        parameters=config["mirna_expression_over_time_box_bar_plot"].get("parameters", None),
        manual_plots_parameters=config["mirna_expression_over_time_box_bar_plot"].get("manual_plots_parameters", None),
        corr_methods=config["correlation_comp_age_heatmap_plot_table"].get("method", None),
        xticks_names=config["xticks_names"],
        colors=config["colors"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/mirna_expression_over_time_box_bar_plot.R"

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

rule tissue_specificity_index:
    input:
    output: result=join(config["results_folder"], get_steps_folder("tissue_specificity_index"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        tissues=config["tissue_specificity_index"].get("tissues", None),
        top_list=config["tissue_specificity_index"].get("top_list", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/tissue_specificity_index.R"

rule hclust_expression_brain_region_for_mirna_list_zscores_heatmap_plot:
    input: 
        requ=join(config["results_folder"], get_steps_folder("coefficient_of_variation")),
        result=join(config["results_folder"], get_steps_folder("tissue_specificity_index")),
        results_diff_exp=join(config["results_folder"], result_file_diff_exp)
    output: result=join(config["results_folder"], get_steps_folder("hclust_expression_brain_region_for_mirna_list_zscores_heatmap_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        prop=config["hclust_expression_brain_region_for_mirna_list_zscores_heatmap_plot"].get("prop", None),
        diff_exp_log=config["hclust_expression_brain_region_for_mirna_list_zscores_heatmap_plot"].get("diff_log", None),
        zscore_th=config["hclust_expression_brain_region_for_mirna_list_zscores_heatmap_plot"].get("zscore_th", None),
        method=config["hclust_expression_brain_region_for_mirna_list_zscores_heatmap_plot"].get("method", None),
        top_list=config["hclust_expression_brain_region_for_mirna_list_zscores_heatmap_plot"].get("top_list", None),
        params=config["hclust_expression_brain_region_for_mirna_list_zscores_heatmap_plot"].get("hclust_params", None),
        plot_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/hclust_expression_brain_region_for_mirna_list_zscores_heatmap_plot.R"

rule hclust_expression_brain_region_for_mirna_list_zscores_graph_venn_plot:
    input: requ_zscores=join(config["results_folder"], get_steps_folder("hclust_expression_brain_region_for_mirna_list_zscores_heatmap_plot"))
    output: result=join(config["results_folder"], get_steps_folder("hclust_expression_brain_region_for_mirna_list_zscores_graph_venn_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        prop=config["hclust_expression_brain_region_for_mirna_list_zscores_graph_venn_plot"].get("prop", None),
        method=config["hclust_expression_brain_region_for_mirna_list_zscores_graph_venn_plot"].get("method", None),
        top_list=config["hclust_expression_brain_region_for_mirna_list_zscores_graph_venn_plot"].get("top_list", None),
        diff_exp_log=config["hclust_expression_brain_region_for_mirna_list_zscores_graph_venn_plot"].get("diff_log", None),
        zscore_th=config["hclust_expression_brain_region_for_mirna_list_zscores_graph_venn_plot"].get("zscore_th", None),
        feature_keep_rate=config["hclust_expression_brain_region_for_mirna_list_zscores_graph_venn_plot"].get("feature_keep_rate", None),
        intersection_th=config["hclust_expression_brain_region_for_mirna_list_zscores_graph_venn_plot"].get("intersection_th", None),
        colors=config["colors"],
        plot_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/hclust_expression_brain_region_for_mirna_list_zscores_graph_venn_plot.R"

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
        xticks_names=config["xticks_names"],
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/fc_versus_expr_scatter_plot.R" 

rule fc_versus_expr_scatter_plot_CA1_age_3m_12m_15m:
    input: result=join(config["results_folder"], get_steps_folder("diff_exp_CA1_age_3m_12m_15m"))
    output: result=join(config["results_folder"], get_steps_folder("fc_versus_expr_scatter_plot_CA1_age_3m_12m_15m"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_exp_log=config["fc_versus_expr_scatter_plot"].get("diff_log", None),
        prop=config["fc_versus_expr_scatter_plot"].get("prop", None),
        comparisons_CA1_age_3m_12m_15m=config["fc_versus_expr_scatter_plot"].get("comparisons_CA1_age_3m_12m_15m", None),
        thresholds=config["fc_versus_expr_scatter_plot"].get("thresholds", None),
        xticks_names=config["xticks_names"],
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/fc_versus_expr_scatter_plot.R" 

rule fc_versus_expr_scatter_plot_CA1_without_age_26m_28m:
    input: result=join(config["results_folder"], get_steps_folder("diff_exp_CA1_without_age_26m_28m"))
    output: result=join(config["results_folder"], get_steps_folder("fc_versus_expr_scatter_plot_CA1_without_age_26m_28m"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_exp_log=config["fc_versus_expr_scatter_plot"].get("diff_log", None),
        prop=config["fc_versus_expr_scatter_plot"].get("prop", None),
        comparisons_CA1_without_age_26m_28m=config["fc_versus_expr_scatter_plot"].get("comparisons_CA1_without_age_26m_28m", None),
        thresholds=config["fc_versus_expr_scatter_plot"].get("thresholds", None),
        xticks_names=config["xticks_names"],
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/fc_versus_expr_scatter_plot.R" 

rule volcano_CA1_scatter_plot:
    input: result=join(config["results_folder"], get_steps_folder("diff_exp_CA1"))
    output: result=join(config["results_folder"], get_steps_folder("volcano_CA1_scatter_plot"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_log=config["fc_versus_expr_scatter_plot"].get("diff_log", None),
        comparisons_age=config["volcano_scatter_plot"].get("comparisons_age", None),
        comparisons_brain_region=config["volcano_scatter_plot"].get("comparisons_brain_region", None),
        parameters=config["volcano_scatter_plot"].get("parameters", None),
        colors=config["colors"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/volcano_scatter_plot.R"

rule volcano_CA1_age_3m_12m_15m_scatter_plot:
    input: result=join(config["results_folder"], get_steps_folder("diff_exp_CA1_age_3m_12m_15m"))
    output: result=join(config["results_folder"], get_steps_folder("volcano_CA1_age_3m_12m_15m_scatter_plot"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_log=config["fc_versus_expr_scatter_plot"].get("diff_log", None),
        comparisons_brain_region=config["volcano_scatter_plot"].get("comparisons_brain_region", None),
        comparisons_age=config["fc_versus_expr_scatter_plot"].get("comparisons_CA1_age_3m_12m_15m", None),
        parameters=config["volcano_scatter_plot"].get("parameters", None),
        colors=config["colors"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/volcano_scatter_plot.R"

rule volcano_CA1_without_age_26m_28m_scatter_plot:
    input: result=join(config["results_folder"], get_steps_folder("diff_exp_CA1_without_age_26m_28m"))
    output: result=join(config["results_folder"], get_steps_folder("volcano_CA1_without_age_26m_28m_scatter_plot"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_log=config["volcano_scatter_plot"].get("diff_log", None),
        comparisons_sex=config["volcano_scatter_plot"].get("comparisons_sex", None),
        omparisons_brain_region=config["volcano_scatter_plot"].get("comparisons_brain_region", None),
        comparisons_age=config["volcano_scatter_plot"].get("comparisons_age_without_age_26m_28m", None),
        parameters=config["volcano_scatter_plot"].get("parameters", None),
        manual_plots_parameters=config["volcano_scatter_plot"].get("manual_plots_parameters", None),
        colors=config["colors"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/volcano_scatter_plot.R"

rule volcano_CA1_without_age_26m_28m_scatter_plot_manual:
    input: result=join(config["results_folder"], get_steps_folder("diff_exp_CA1_without_age_26m_28m"))
    output: result=join(config["results_folder"], get_steps_folder("volcano_CA1_without_age_26m_28m_scatter_plot_manual"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_log=config["volcano_scatter_plot_manual"].get("diff_log", None),
        comparisons=config["volcano_scatter_plot_manual"].get("comparisons", None),
        parameters=config["volcano_scatter_plot_manual"].get("parameters", None),
        manual_plots_parameters=config["volcano_scatter_plot_manual"].get("manual_plots_parameters", None),
        colors=config["colors"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/volcano_scatter_plot_manual.R"

rule volcano_CA1_age_3m_12m_15m_TMS_brain_age_3m_12m_15m_scatter_plot:
    input: result=join(config["results_folder"], get_steps_folder("diff_exp_CA1_age_3m_12m_15m_TMS_brain_age_3m_12m_15m"))
    output: result=join(config["results_folder"], get_steps_folder("volcano_CA1_age_3m_12m_15m_TMS_brain_age_3m_12m_15m_scatter_plot"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        diff_log=config["fc_versus_expr_scatter_plot"].get("diff_log", None),
        comparisons_brain_region_tms=config["volcano_scatter_plot"].get("comparisons_brain_region_tms", None),
        parameters=config["volcano_scatter_plot"].get("parameters", None),
        colors=config["colors"],
        plots_props=config["common"].get("plots", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/volcano_scatter_plot.R"

#rule pathway_analysis_results_heatmap_plots:
#    input: featuer_list=config["feature_list"]
#    output: result=join(config["results_folder"], get_steps_folder("pathway_analysis_results_heatmap_plots"))
#    params:
#        parameters_props=config["common"].get("parameters", None),
#        filter_params=config["pathway_analysis_results_heatmap_plots"].get("filter_params", None),
#        plot_props=config["common"].get("plots", None),
#        results_folder=config["results_folder"]
#    conda: "envs/mouse_cold_p.yml"
#    script: "scripts/pathway_analysis_results_heatmap_plots.R"

#rule halushka_violin_plots:
#    input: 
#    output: result=join(config["results_folder"], get_steps_folder("halushka_violin_plots"))
#    params: 
#        parameters_props=config["common"].get("parameters", None),
#        plots_props=config["common"].get("plots", None),
#        results_folder=config["results_folder"]
#    conda: "envs/mouse_cold_p.yml"
#    script: "scripts/halushka_violin_plots.R"

# rule mfuzz_ts_clustering:
#     input: result=join(config["results_folder"], get_steps_folder("diff_exp_CA1"))
#     output: result=join(config["results_folder"], get_steps_folder("mfuzz_ts_clustering"))
#     params:
#         data_input=config["common"].get("input", None),
#         parameters_props=config["common"].get("parameters", None),
#         cluster_props=config["mfuzz_ts_clustering"].get("cluster_props", None),
#         plot_category=config["mfuzz_ts_clustering"].get("plot_category", None),
#         diff_log=config["mfuzz_ts_clustering"].get("diff_log", None),
#         plots_props=config["common"].get("plots", None),
#         colors=config["colors"],
#          xticks_names=config["xticks_names"],
#         results_folder=config["results_folder"]
#     conda: "envs/analysis_packages.yml"
#     script: "scripts/mfuzz_ts_clustering.R"

rule mfuzz_ts_analysis_overview_line_plot:
    input: 
        result_diff_exp=join(config["results_folder"], get_steps_folder("diff_exp_CA1")),
        #result_clustering=join(config["results_folder"], get_steps_folder("mfuzz_ts_clustering")),
        feature_list=config["feature_list"],
        cluster_results_folder_path=config["cluster_results_folder_path"]
    output: result=join(config["results_folder"], get_steps_folder("mfuzz_ts_analysis_overview_line_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        cluster_props=config["mfuzz_ts_clustering"].get("cluster_props", None),
        plot_category=config["mfuzz_ts_clustering"].get("plot_category", None),
        diff_log=config["mfuzz_ts_clustering"].get("diff_log", None),
        thresholds=config["mfuzz_ts_analysis_overview_line_plot"].get("thresholds", None),
        xticks_names=config["xticks_names"],
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/mfuzz_ts_analysis_overview_line_plot.R"

rule mfuzz_ts_analysis_membership_box_plot:
    input: 
        result_diff_exp=join(config["results_folder"], get_steps_folder("diff_exp_CA1")),
        #result_clustering=join(config["results_folder"], get_steps_folder("mfuzz_ts_clustering")),
        feature_list=config["feature_list"],
        cluster_results_folder_path=config["cluster_results_folder_path"]
    output: result=join(config["results_folder"], get_steps_folder("mfuzz_ts_analysis_membership_box_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        cluster_props=config["mfuzz_ts_clustering"].get("cluster_props", None),
        thresholds=config["mfuzz_ts_analysis_overview_line_plot"].get("thresholds", None),
        interesting_cluster=config["mfuzz_ts_analysis_membership_box_plot"].get("interesting_cluster", None),
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/mfuzz_ts_analysis_membership_box_plot.R"

rule mfuzz_ts_analysis_tissue_specific_clusters_line_plot:
    input: 
        result_diff_exp=join(config["results_folder"], get_steps_folder("diff_exp_CA1")),
        result_clustering=join(config["results_folder"], get_steps_folder("mfuzz_ts_analysis_overview_line_plot")),
        feature_list=config["feature_list"],
        cluster_results_folder_path=config["cluster_results_folder_path"]
    output: result=join(config["results_folder"], get_steps_folder("mfuzz_ts_analysis_tissue_specific_clusters_line_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        cluster_props=config["mfuzz_ts_clustering"].get("cluster_props", None),
        plot_category=config["mfuzz_ts_clustering"].get("plot_category", None),
        diff_log=config["mfuzz_ts_clustering"].get("diff_log", None),
        thresholds=config["mfuzz_ts_analysis_overview_line_plot"].get("thresholds", None),
        plots_props=config["common"].get("plots", None),
        xticks_names=config["xticks_names"],
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/mfuzz_ts_analysis_tissue_specific_clusters_line_plot.R"

rule mfuzz_ts_analysis_scatter_plot:
    input: 
        result_diff_exp=join(config["results_folder"], get_steps_folder("diff_exp_CA1")),
        result_clustering=join(config["results_folder"], get_steps_folder("mfuzz_ts_analysis_overview_line_plot")),
        feature_list=config["feature_list"],
        cluster_results_folder_path=config["cluster_results_folder_path"]
    output: result=join(config["results_folder"], get_steps_folder("mfuzz_ts_analysis_scatter_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        cluster_props=config["mfuzz_ts_clustering"].get("cluster_props", None),
        plot_category=config["mfuzz_ts_clustering"].get("plot_category", None),
        thresholds=config["mfuzz_ts_analysis_overview_line_plot"].get("thresholds", None),
        xticks_names=config["xticks_names"],
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/mfuzz_ts_analysis_scatter_plot.R"

rule mfuzz_ts_analysis_single_line_plot:
    input: 
        result_diff_exp=join(config["results_folder"], get_steps_folder("diff_exp_CA1")),
        #result_clustering=join(config["results_folder"], get_steps_folder("mfuzz_ts_clustering")),
        cluster_results_folder_path=config["cluster_results_folder_path"]
    output: result=join(config["results_folder"], get_steps_folder("mfuzz_ts_analysis_single_line_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        cluster_props=config["mfuzz_ts_clustering"].get("cluster_props", None),
        plot_category=config["mfuzz_ts_clustering"].get("plot_category", None),
        diff_log=config["mfuzz_ts_clustering"].get("diff_log", None),
        thresholds=config["mfuzz_ts_analysis_overview_line_plot"].get("thresholds", None),
        interesting_cluster=config["mfuzz_ts_analysis_single_line_plot"].get("interesting_cluster", None),
        xticks_names=config["xticks_names"],
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/mfuzz_ts_analysis_single_line_plot.R"

rule mfuzz_ts_analysis_age_development_clusters:
    input: 
        result_diff_exp=join(config["results_folder"], get_steps_folder("diff_exp_CA1")),
        #result_clustering=join(config["results_folder"], get_steps_folder("mfuzz_ts_clustering")),
        cluster_results_folder_path=config["cluster_results_folder_path"]
    output: result=join(config["results_folder"], get_steps_folder("mfuzz_ts_analysis_age_development_clusters"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        cluster_props=config["mfuzz_ts_clustering"].get("cluster_props", None),
        plot_category=config["mfuzz_ts_clustering"].get("plot_category", None),
        diff_log=config["mfuzz_ts_clustering"].get("diff_log", None),
        thresholds=config["mfuzz_ts_analysis_overview_line_plot"].get("thresholds", None),
        cluster_groups=config["mfuzz_ts_analysis_age_development_clusters"].get("cluster_groups", None),
        xticks_names=config["xticks_names"],
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/mfuzz_ts_analysis_age_development_clusters.R"

rule fc_brain_region_vs_tms_brain_upset_plot:
    input:
    output: result=join(config["results_folder"], get_steps_folder("fc_brain_region_vs_tms_brain_upset_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        params=config["fc_brain_region_vs_tms_brain_upset_plot"].get("parameters", None),
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/fc_brain_region_vs_tms_brain_upset_plot.R"

rule pathway_analysis_bubble_plot:
    input:
    output: result=join(config["results_folder"], get_steps_folder("pathway_analysis_bubble_plot"))
    params:
        data_input=config["common"].get("input", None),
        sig_lvl=config["pathway_analysis_bubble_plot"].get("sig_lvl", None),
        method=config["pathway_analysis_bubble_plot"].get("method", None),
        comparisons=config["pathway_analysis_bubble_plot"].get("comparisons", None),
        xticks_names=config["xticks_names"],
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/pathway_analysis_bubble_plot.R"

rule pathway_analysis_bubble_plot_age_correlated:
    input:
    output: result=join(config["results_folder"], get_steps_folder("pathway_analysis_bubble_plot_age_correlated"))
    params:
        data_input=config["common"].get("input", None),
        sig_lvl=config["pathway_analysis_bubble_plot_age_correlated"].get("sig_lvl", None),
        method=config["pathway_analysis_bubble_plot_age_correlated"].get("method", None),
        comparisons=config["pathway_analysis_bubble_plot_age_correlated"].get("comparisons", None),
        xticks_names=config["xticks_names"],
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/pathway_analysis_bubble_plot_age_correlated.R"

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
    script: "scripts/target_genes_scatter_plot_CA1.R"

rule preprocessing_mRNA_data:
    input:
    output: result=join(config["results_folder"], get_steps_folder("preprocessing_mRNA_data"))
    params:
        data_input=config["common"].get("input", None),
    conda: "envs/analysis_packages.yml"
    script: "scripts/preprocessing_mRNA_data.R"

rule correlation_mirna_mrna:
    input: result=join(config["results_folder"], get_steps_folder("preprocessing_mRNA_data"))
    output: result_local=join(config["results_folder"], get_steps_folder("correlation_mirna_mrna"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        method=config["correlation_mirna_mrna"].get("method", None),
        adjustment=config["correlation_mirna_mrna"].get("adjustment", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: 
        "scripts/correlation_mirna_mrna.R"

rule correlation_mirna_mrna_heatmap_plot:
    input: 
        #targetscan=config["targetscan_mirna_targets"],
        #result_local=join(config["results_folder"], get_steps_folder("correlation_mirna_mrna")),
        result=join(config["results_folder"], get_steps_folder("preprocessing_mRNA_data"))
    output: result_local=join(config["results_folder"], get_steps_folder("correlation_mirna_mrna_heatmap_plot"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        method=config["correlation_mirna_mrna"].get("method", None),
        adjustment=config["correlation_mirna_mrna"].get("adjustment", None),
        corr_th=config["correlation_mirna_mrna"].get("corr_th", None),
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: 
        "scripts/correlation_mirna_mrna_heatmap_plot.R"

rule correlation_mirna_target_genes_cas:
    input: result=join(config["results_folder"], get_steps_folder("preprocessing_mRNA_data"))
    output: result_local=join(config["results_folder"], get_steps_folder("correlation_mirna_target_genes_cas"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        feature=config["correlation_mirna_target_genes_cas"].get("feature", None),
        gene_list_cas_folder_path=config["correlation_mirna_target_genes_cas"].get("gene_list_cas_folder_path", None),
        gene_list_target_genes_folder_path=config["correlation_mirna_target_genes_cas"].get("gene_list_target_genes_folder_path", None),
        split_prop=config["correlation_mirna_target_genes_cas"].get("split_prop", None),
        corr_params=config["correlation_mirna_target_genes_cas"].get("corr_params", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: 
        "scripts/correlation_mirna_target_genes_cas.R"

rule correlation_mirna_target_genes_cas_upset_plot:
    input: 
        result=join(config["results_folder"], get_steps_folder("preprocessing_mRNA_data")),
        result_clustering=join(config["results_folder"], get_steps_folder("correlation_mirna_target_genes_cas"))
    output: result_local=join(config["results_folder"], get_steps_folder("correlation_mirna_target_genes_cas_upset_plot"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        feature=config["correlation_mirna_target_genes_cas"].get("feature", None),
        gene_list_cas_folder_path=config["correlation_mirna_target_genes_cas"].get("gene_list_cas_folder_path", None),
        gene_list_target_genes_folder_path=config["correlation_mirna_target_genes_cas"].get("gene_list_target_genes_folder_path", None),
        split_prop=config["correlation_mirna_target_genes_cas"].get("split_prop", None),
        corr_params=config["correlation_mirna_target_genes_cas"].get("corr_params", None),
        colour_genes=config["correlation_mirna_target_genes_cas"].get("colour_genes", None),
        plots_props=config["common"].get("plots", None),
        xticks_names=config["xticks_names"],
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: 
        "scripts/correlation_mirna_target_genes_cas_upset_plot.R"

rule correlation_mirna_target_genes_cas_scatter_plot:
    input:
        result=join(config["results_folder"], get_steps_folder("preprocessing_mRNA_data")),
        result_cas=join(config["results_folder"], get_steps_folder("correlation_mirna_target_genes_cas"))
    output: result_local=join(config["results_folder"], get_steps_folder("correlation_mirna_target_genes_cas_scatter_plot"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        feature=config["correlation_mirna_target_genes_cas"].get("feature", None),
        gene_list_cas_folder_path=config["correlation_mirna_target_genes_cas"].get("gene_list_cas_folder_path", None),
        gene_list_target_genes_folder_path=config["correlation_mirna_target_genes_cas"].get("gene_list_target_genes_folder_path", None),
        save_folder=config["correlation_mirna_target_genes_cas"].get("save_folder", None),
        split_prop_scatter=config["correlation_mirna_target_genes_cas"].get("split_prop_scatter", None),
        corr_params=config["correlation_mirna_target_genes_cas"].get("corr_params", None),
        manual_plots_parameters=config["correlation_mirna_target_genes_cas"].get("manual_plots_parameters", None),
        plots_props=config["common"].get("plots", None),
        xticks_names=config["xticks_names"],
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/correlation_mirna_target_genes_cas_scatter_plot.R"

rule correlation_mirna_target_genes_martin:
    input: result=join(config["results_folder"], get_steps_folder("preprocessing_mRNA_data"))
    output: result_local=join(config["results_folder"], get_steps_folder("correlation_mirna_target_genes_martin"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        feature=config["correlation_mirna_target_genes_martin"].get("feature", None),
        gene_list_cas_folder_path=config["correlation_mirna_target_genes_martin"].get("gene_list_cas_folder_path", None),
        gene_list_target_genes=config["correlation_mirna_target_genes_martin"].get("gene_list_target_genes", None),
        split_prop=config["correlation_mirna_target_genes_martin"].get("split_prop", None),
        corr_params=config["correlation_mirna_target_genes_martin"].get("corr_params", None),
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/correlation_mirna_target_genes_martin.R"
    
rule correlation_mirna_target_genes_martin_scatter_plot:
    input:
        result=join(config["results_folder"], get_steps_folder("preprocessing_mRNA_data")),
        result_martin=join(config["results_folder"], get_steps_folder("correlation_mirna_target_genes_martin")),
    output: result_local=join(config["results_folder"], get_steps_folder("correlation_mirna_target_genes_martin_scatter_plot"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        feature=config["correlation_mirna_target_genes_cas"].get("feature", None),
        gene_list_cas_folder_path=config["correlation_mirna_target_genes_martin"].get("gene_list_cas_folder_path", None),
        gene_list_target_genes_folder_path=config["correlation_mirna_target_genes_martin"].get("gene_list_target_genes_folder_path", None),
        save_folder=config["correlation_mirna_target_genes_martin"].get("save_folder", None),
        split_prop_scatter=config["correlation_mirna_target_genes_martin"].get("split_prop_scatter", None),
        corr_params=config["correlation_mirna_target_genes_martin"].get("corr_params", None),
        manual_plots_parameters=config["correlation_mirna_target_genes_martin"].get("manual_plots_parameters", None),
        plots_props=config["common"].get("plots", None),
        xticks_names=config["xticks_names"],
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/correlation_mirna_target_genes_cas_scatter_plot.R"

rule correlation_mirna_target_genes_mtor_table_scatter_plot:
    input: result=join(config["results_folder"], get_steps_folder("preprocessing_mRNA_data"))
    output: result_local=join(config["results_folder"], get_steps_folder("correlation_mirna_target_genes_mtor_table_scatter_plot"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        feature=config["correlation_mirna_target_genes_mtor_table_scatter_plot"].get("feature", None),
        gene_list_cas_folder_path=config["correlation_mirna_target_genes_mtor_table_scatter_plot"].get("gene_list_cas_folder_path", None),
        gene_list_target_genes=config["correlation_mirna_target_genes_mtor_table_scatter_plot"].get("gene_list_target_genes", None),
        gene_list_target_genes_folder_path=config["correlation_mirna_target_genes_mtor_table_scatter_plot"].get("gene_list_target_genes_folder_path", None),
        gene_list_mtor=config["correlation_mirna_target_genes_mtor_table_scatter_plot"].get("gene_list_mtor", None),
        split_prop=config["correlation_mirna_target_genes_mtor_table_scatter_plot"].get("split_prop", None),
        corr_prop=config["correlation_mirna_target_genes_mtor_table_scatter_plot"].get("corr_prop", None),
        corr_params=config["correlation_mirna_target_genes_mtor_table_scatter_plot"].get("corr_params", None),
        manual_plots_parameters=config["correlation_mirna_target_genes_mtor_table_scatter_plot"].get("manual_plots_parameters", None),
        plots_props=config["common"].get("plots", None),
        xticks_names=config["xticks_names"],
        colors=config["colors"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/correlation_mirna_target_genes_mtor_table_scatter_plot.R"

rule external_bar_plot:
    input: 
    output: result_local=join(config["results_folder"], get_steps_folder("external_bar_plot"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        raw_data_path=config["external_bar_plot"].get("raw_data_path", None),
        id_col=config["external_bar_plot"].get("id_col", None),
        plot_list=config["external_bar_plot"].get("plot_list", None),
        #rna=config["external_bar_plot"].get("rna", None),
        #plot_titles=config["external_bar_plot"].get("plot_titles", None),
        columns_mean=config["external_bar_plot"].get("columns_mean", None),
        columns_fc=config["external_bar_plot"].get("columns_fc", None),
        results_folder=config["external_bar_plot"].get("results_folder", None),
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
    conda: "envs/analysis_packages.yml"
    script: 
        "scripts/external_bar_plot.R"

rule external_protein_plot:
    input: 
    output: result_local=join(config["results_folder"], get_steps_folder("external_protein_plot"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        raw_data_path=config["external_protein_plot"].get("raw_data_path", None),
        age_colour=config["external_protein_plot"].get("age_colour", None),
        results_folder=config["external_protein_plot"].get("results_folder", None),
        intensity_col=config["external_protein_plot"].get("intensity_col", None),
        plot_parameters=config["external_protein_plot"].get("plot_parameters", None),
        plots_props=config["common"].get("plots", None),
        xticks_names=config["xticks_names"],
        colors=config["colors"],
    conda: "envs/analysis_packages.yml"
    script: 
        "scripts/external_protein_plot.R"

rule overlap_sig_correlated_paper_bar_plot:
    input: result=join(config["results_folder"], get_steps_folder("correlation_comp_age_heatmap_plot_table"))
    output: result=join(config["results_folder"], get_steps_folder("overlap_sig_correlated_paper_bar_plot"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        raw_data_path=config["overlap_sig_correlated_paper_bar_plot"].get("raw_data_path", None),
        thresholds=config["overlap_sig_correlated_paper_bar_plot"].get("thresholds", None),
        props=config["overlap_sig_correlated_paper_bar_plot"].get("props", None),
        method=config["overlap_sig_correlated_paper_bar_plot"].get("method", None),
        plot_size=config["overlap_sig_correlated_paper_bar_plot"].get("plot_size", None),
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        xticks_names=config["xticks_names"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/overlap_sig_correlated_paper_bar_plot.R"

rule mirna_neighborhood:
    input: result=join(config["results_folder"], get_steps_folder("correlation_comp_age_heatmap_plot_table"))
    output: result=join(config["results_folder"], get_steps_folder("mirna_neighborhood"))
    params:
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        thresholds=config["mirna_neighborhood"].get("thresholds", None),
        method=config["mirna_neighborhood"].get("method", None),
        props=config["mirna_neighborhood"].get("props", None),
        plot_layout=config["mirna_neighborhood"].get("plot_layout", None),
        plot_size=config["mirna_neighborhood"].get("plot_size", None),
        plots_props=config["common"].get("plots", None),
        colors=config["colors"],
        xticks_names=config["xticks_names"],
        results_folder=config["results_folder"]
    conda: "envs/analysis_packages.yml"
    script: "scripts/mirna_neighborhood.R"

rule isomir_line_plot: 
    input: 
    output: result_local=join(config["results_folder"], get_steps_folder("isomir_line_plot"))
    params: 
        data_input=config["common"].get("input", None),
        parameters_props=config["common"].get("parameters", None),
        plot_group=config["isomir_line_plot"].get("plot_group", None),
        comp_group=config["isomir_line_plot"].get("comp_group", None),
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
        "scripts/isomir_line_plot_numeric.R"

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

rule make_publication_folder:
   input:
   output: result=join(config["results_folder"], get_steps_folder("make_publication_folder"))
   params:
   conda: "envs/analysis_packages.yml"
   script: "scripts/extract_paper_figures.py"
