#!/usr/bin/env python3

"""
Extract the relevant figures for the publication.
In addition, rename them to the appropriate figure name.
Finally, remove the source figure folder.

Usage:
Run this script from the folder where the "results" folder is located!
"""

from pathlib import Path
import shutil

# import file_mapping.py
# from .file_mapping import file_mapping

prefix = "results/results_rpmm/"
file_mapping = {
    "fig_1b": prefix + "ncRNA_10p/results_CA1_aging/figures/corr_plots/age/ridgeline_plot/corr_method=spearman_corr_th=0.5_features_with_age_per_brain_region_6x9",
    "fig_1c": prefix + "tRNA_10p/results_CA1_aging_male_CA1_aging_without_age=26m_28m_female/figures/expr/box_plots/time_single_corr_method=spearman_poly_deg=2/tRNA-Glu-TTC-1-1_cor_th_plx_tRNA-Glu-TTC-1-1_cor_th_plx",
    "fig_1d": prefix + "tRNA_10p/results_CA1_aging/figures/umap_plots/scatter_plots/brain_region/all/std/umap_neigh25_mdist0.25_euclidean_initrandom_seed42_norm_std_6x9",
    "fig_1e": prefix + "tRNA_10p/results_CA1_aging/figures/pvca_plots/bar_plots/all_brain_regions_9x6",
    "fig_1f": prefix + "miRNA_10p/results_CA1_aging/figures/umap_plots/scatter_plots/brain_region/all/std/umap_neigh10_mdist0.25_euclidean_initrandom_seed42_norm_std_6x9",
    "fig_1g": prefix + "miRNA_10p/results_CA1_aging/figures/pvca_plots/bar_plots/all_brain_regions_9x6",

    "fig_2a": prefix + "miRNA_10p/results_CA1_aging_age=3m_12m_15m_male_CA1_aging_age=3m_12m_15m_female/figures/cv_clustering/heatmap_complex/50_zscore_th=0.5_num_clusters=4_removed_zero_show_legend=TRUE_row=miRNA_portrait_mod_15x12",
    "fig_2b": prefix + "miRNA_10p/results_CA1_aging_age=3m_12m_15m_male_CA1_aging_age=3m_12m_15m_female/figures/cv_clustering/venn_plots/50_zscore_th=0.5_clustered_by_cv_vertical",
    "fig_2c": prefix + "miRNA_10p/results_CA1_aging_without_age=26m_28m_without_pon/figures/pvca_plots/scatter_plots/combined_brain_region_results_age_sex_6x6",
    "fig_2d": prefix + "miRNA_10p/results_CA1_aging_without_age=26m_28m_without_pon/figures/volcano/scatter_plots/merged/sex__Male_vs_Female_brain_region=cor_cp_th_olf.ttest.adj_2x2_6x6",
    "fig_2e": prefix + "miRNA_10p/results_CA1_aging_without_age=26m_28m_without_pon/figures/pathway/bubble_plot/gsea_sex__Male_vs_Female_brain_region=cor_cp_th",

    "fig_3a": prefix + "miRNA_10p/results_CA1_aging_without_pon_male_CA1_aging_without_age=26m_28m_without_pon_female/figures/expr/box_plots/time_single_corr_method=spearman_poly_deg=2/mmu-miR-9-5p_med_cor_olf_mmu-miR-9-5p_med_cor_olf",
    "fig_3b_male": prefix + "miRNA_10p/results_CA1_aging_without_pon_male/figures/corr_plots/bar_plots/corr_method=spearman_th=0.5_sig_aggregation_for_2_in_brain_region_complete_9x3_flip_left",
    "fig_3b_female": prefix + "miRNA_10p/results_CA1_aging_without_age=26m_28m_without_pon_female/figures/corr_plots/bar_plots/corr_method=spearman_th=0.5_sig_aggregation_for_2_in_brain_region_complete_9x3_flip_right",
    "fig_3c_male": prefix + "miRNA_10p/results_CA1_aging_without_pon_male/figures/corr_plots/heatmap_complex/global_miRNAs/corr_method=spearman_global_miRNA_over_age_corr_th=0.5_for_2_brain_region_sig_dots",
    "fig_3c_female": prefix + "miRNA_10p/results_CA1_aging_without_age=26m_28m_without_pon_female/figures/corr_plots/heatmap_complex/global_miRNAs/corr_method=spearman_global_miRNA_over_age_corr_th=0.5_for_2_brain_region_sig_dots",
    "fig_3d": prefix + "miRNA_10p/results_CA1_aging_without_pon_male_CA1_aging_without_age=26m_28m_without_pon_female/figures/diff_exp/heatmap_complex/num_fc_th_1.5_per_brain_region/diff_exp_global_de_reg_age_per_miRNA_adj",
    "fig_3e_male": prefix + "miRNA_10p/results_CA1_aging_without_pon_male/figures/diff_exp/heatmap_complex/global_miRNAs/diff_exp_global_sig_de_reg_miRNA_at_least_in_1_age_comps_for_2_brain_region",
    "fig_3e_female": prefix + "miRNA_10p/results_CA1_aging_without_age=26m_28m_without_pon_female/figures/diff_exp/heatmap_complex/global_miRNAs/diff_exp_global_sig_de_reg_miRNA_at_least_in_1_age_comps_for_2_brain_region",
    "fig_3f_male": prefix + "miRNA_10p/results_CA1_aging_without_pon_male_CA1_aging_without_age=26m_28m_without_pon_female/figures/candidate_comparison_brain_region_split/scatter_plots/candidate_comp_sig_diff_exp_both_miRNA_at_least_in_1_age_comps_for_2_brain_region_sig_corr_method=spearman_both_for_2_brain_region_miRNA_corr_th=0.5_hclust_method=cv_zscore_th=0.5_top=50_miRNAs_log10_male",
    "fig_3f_female": prefix + "miRNA_10p/results_CA1_aging_without_pon_male_CA1_aging_without_age=26m_28m_without_pon_female/figures/candidate_comparison_brain_region_split/scatter_plots/candidate_comp_sig_diff_exp_both_miRNA_at_least_in_1_age_comps_for_2_brain_region_sig_corr_method=spearman_both_for_2_brain_region_miRNA_corr_th=0.5_hclust_method=cv_zscore_th=0.5_top=50_miRNAs_log10_female",
    "fig_3g": prefix + "miRNA_10p/results_ROSMAP/figures/deregulated_miRNAs_human_age_bin_3_groups/deregulation_fc_th=1.5_per_age_bin_3_groups__92_vs_71_msex=1_age_bin_3_groups__92_vs_71_msex=0_12x3.5",
    "fig_3h": prefix + "miRNA_10p/results_ROSMAP/figures/volcano/scatter_plots/merged/age_bin_3_groups__92_vs_71_msex=1_0.cohend_estimate_1x2_6x3",

    "fig_4a": prefix + "miRNA_10p/results_CA1_aging/figures/diff_exp/heatmap_complex/num_fc_th_1.5_per_brain_region/diff_exp_global_de_reg_age_per_miRNA_adj",
    "fig_4b": prefix + "miRNA_10p/results_CA1_aging/figures/candidate_comparison_brain_region_split/upset_plots/candidate_comp_sig_diff_exp_both_miRNA_at_least_in_1_age_comps_for_2_brain_region_sig_corr_method=spearman_both_for_2_brain_region_miRNA_corr_th=0.5_intersection_th=5_6x9",
    "fig_4c": prefix + "miRNA_10p/results_CA1_aging/figures/mfuzz_clustering/k=53/0.15/scatter_plots/cluster_overview_tissue_occ_th=0.3_miRNA_occ_th=4_9x6",
    "fig_4d": prefix + "miRNA_10p/results_CA1_aging/figures/mfuzz_clustering/k=53/0.15/line_plots/tissue_specific_clusters_more_0.3_from_one_brain_region_center_lines_colouring_adj",
    "fig_4e_cl_41": prefix + "miRNA_10p/results_CA1_aging/figures/mfuzz_clustering/k=53/0.15/line_plots/single/cluster_41_expression_6x6",
    "fig_4e_cl_28": prefix + "miRNA_10p/results_CA1_aging/figures/mfuzz_clustering/k=53/0.15/line_plots/single/cluster_28_expression_6x6",
    "fig_4e_cl_38": prefix + "miRNA_10p/results_CA1_aging/figures/mfuzz_clustering/k=53/0.15/line_plots/single/cluster_38_expression_6x6",
    "fig_4e_cl_24": prefix + "miRNA_10p/results_CA1_aging/figures/mfuzz_clustering/k=53/0.15/line_plots/single/cluster_24_expression_6x6",
    "fig_4f": prefix + "miRNA_10p/results_CA1_aging/figures/mfuzz_clustering/k=53/0.15/scatter_plots/cluster_analysis_tissue_occ_th=0.3_miRNA_occ_th=4_9x6_label",
    "fig_4g": prefix + "miRNA_10p/results_CA1_aging/figures/mfuzz_clustering/k=53/0.15/line_plots/age_development/mmu-miR-155-5p_center_lines_for_clusters_11_25_33_28_23_51_22_legend_side=right_9x6",
    "fig_4h": prefix + "miRNA_10p/results_CA1_aging/figures/expr/box_plots/time_single_corr_method=spearman_poly_deg=2/mmu-miR-155-5p/mmu-miR-155-5p_vis_svz_cp_cer_cc_pon_med_18x4",

    "fig_5a": "results/results_external/paper_10.1523_JNEUROSCI.3582-16.2017/figures/bar_plots/miR-155-5p",
    "fig_5b": prefix + "miRNA_10p/results_microglia/figures/set_overlap/venn_plots/overlap_young_old_paper_label",
    "fig_5c": prefix + "miRNA_10p/results_microglia/figures/top_expressed/heatmap_complex/top_25_log10_miRNAs",
    "fig_5d": prefix + "miRNA_10p/results_microglia/figures/cv_clustering/heatmap_complex/top_25_miRNAs",
    "fig_5e": prefix + "miRNA_10p/results_microglia/figures/expr/line_plots/prop_comparison_for_3_miRNAs_log10_quad",
    "fig_5f": prefix + "miRNA_10p/results_microglia/figures/volcano/scatter_plots/single/age__old_vs_young.ttest.raw.labels_quad",
    "fig_5g": prefix + "iso_10p/results_CA1_aging/figures/iso/line_plots/mmu-miR-155-5p_vis_svz_cer_cc_pon",
    "fig_5h": prefix + "iso_10p/results_CA1_aging/figures/iso/heatmap_complex/zscore=row/top25_isomiRs_zscores=row_mmu-miR-155-5p_9x6",
    "fig_5i": prefix + "miRNA_10p/results_CA1_aging/figures/mmu-miR-155-5p_target_genes_corr/upset_plots/target_genes_all/spearman/upset_mmu-miR-155-5p_corr_m=spearman_corr_th=0.3_intersection_th=1_12x9",
    "fig_5j": prefix + "miRNA_10p/results_CA1_aging/figures/mmu-miR-155-5p_target_genes_corr/scatter_plots/target_genes_all/spearman/scatter_plot_PTPRJ_vs_mmu-miR-155-5p_corr_m=spearman_corr_th=0.3_coloured_by_brain_region_ARNTL_CCNH_ZFP322A_PTPRJ_9x9",


    "supp_fig_1a": prefix + "miRNA_10p/results_CA1_aging/figures/mapping_statistics/bar_plots/alignment/mapping_stats_grouping_brain_region_age_9x6",    
    "supp_fig_1b": prefix + "ncRNA_10p/results_CA1_aging/figures/mapping_statistics/bar_plots/mapping/grouped_by_brain_region_sample_6x9",   
    "supp_fig_1c": prefix + "ncRNA_10p/results_CA1_aging/figures/umap_plots/scatter_plots/brain_region/all/std/umap_neigh25_mdist0.25_euclidean_initrandom_seed42_norm_std_6x6",
    "supp_fig_1d": prefix + "ncRNA_10p/results_CA1_aging/figures/umap_plots/scatter_plots/sex/all/std/umap_neigh25_mdist0.25_euclidean_initrandom_seed42_norm_std_6x6", 
    "supp_fig_1e": prefix + "ncRNA_10p/results_CA1_aging/figures/umap_plots/scatter_plots/age/all/std/umap_neigh25_mdist0.25_euclidean_initrandom_seed42_norm_std_6x6",   
    "supp_fig_1f": prefix + "ncRNA_10p/results_CA1_aging/figures/expressed_per_brain_region_per_rna_class/bar_plots/grouped_by_brain_region_6x6",   
    "supp_fig_1g": prefix + "ncRNA_10p/results_CA1_aging/figures/expressed_per_brain_region_per_rna_class/line_plots/poly_deg=3/grouped_by_cc_cer_cor_cp_ent_hi_hi2_hy_med_olf_plx_pon_svz_th_vis_age_12x12",   

    "supp_fig_2a": prefix + "tRNA_10p/results_CA1_aging/figures/umap_plots/scatter_plots/sex/all/std/umap_neigh25_mdist0.25_euclidean_initrandom_seed42_norm_std_6x6",
    "supp_fig_2b": prefix + "tRNA_10p/results_CA1_aging/figures/umap_plots/scatter_plots/age/all/std/umap_neigh25_mdist0.25_euclidean_initrandom_seed42_norm_std_6x6",
    "supp_fig_2c": prefix + "miRNA_10p/results_CA1_aging/figures/umap_plots/scatter_plots/sex/all/std/umap_neigh10_mdist0.25_euclidean_initrandom_seed42_norm_std_6x6",
    "supp_fig_2d": prefix + "miRNA_10p/results_CA1_aging/figures/umap_plots/scatter_plots/age/all/std/umap_neigh10_mdist0.25_euclidean_initrandom_seed42_norm_std_6x6",
    "supp_fig_2e": prefix + "miRNA_10p/results_CA1_aging/figures/expressed_per_brain_region/upset_plots/number_of_expressed_miRNAs_per_brain_region_th=5_counts_detection_rate=10_intersection_th=5_6x6",
    "supp_fig_2f": prefix + "miRNA_10p/results_CA1_aging_age=3m_12m_15m/figures/umap_plots/scatter_plots/brain_region/all/std/umap_neigh5_mdist0.25_euclidean_initrandom_seed42_norm_std_6x6",
    "supp_fig_2g": prefix + "miRNA_10p/results_CA1_aging_age=3m_12m_15m/figures/umap_plots/scatter_plots/sex/all/std/umap_neigh5_mdist0.25_euclidean_initrandom_seed42_norm_std_6x6",
    "supp_fig_2h": prefix + "miRNA_10p/results_CA1_aging_age=3m_12m_15m/figures/umap_plots/scatter_plots/age/all/std/umap_neigh5_mdist0.25_euclidean_initrandom_seed42_norm_std_6x6",
    "supp_fig_2i": prefix + "miRNA_10p/results_CA1_aging_age=3m_12m_15m/figures/cv_clustering/heatmap_complex/50_zscore_th=0.5_num_clusters=4_removed_zero_rows_rownames_pos=left_show_legend=TRUE_row=miRNA_portrait_mod",
    "supp_fig_2j": prefix + "miRNA_10p/results_CA1_aging/figures/expr/box_plots/brain_region/mmu-miR-9-5p",

    "supp_fig_3a": prefix + "miRNA_10p/results_CA1_aging_without_pon_male_CA1_aging_without_age=26m_28m_without_pon_female/figures/expr/box_plots/time_single_corr_method=spearman_poly_deg=2/mmu-miR-9-5p_cc_cer_cp_ent_hi_hi2_hy_plx_svz_th_vis_mmu-miR-9-5p_cc_cer_cp_ent_hi_hi2_hy_plx_svz_th_vis",
    "supp_fig_3b": prefix + "miRNA_10p/results_CA1_aging_without_pon_male_CA1_aging_without_age=26m_28m_without_pon_female/figures/corr_plots/heatmap_complex/all_miRNA_with_age/corr_method=spearman_miRNA_over_age_per_brain_region_filtered_corr_th=0.5",
    "supp_fig_3c": prefix + "miRNA_10p/results_CA1_aging_without_pon_male_CA1_aging_without_age=26m_28m_without_pon_female/figures/diff_exp/heatmap_complex/num_fc_th_1.5_per_brain_region/diff_exp_global_sig_de_reg_age_per_miRNA_adj",
    
    "supp_fig_4a_male": prefix + "miRNA_10p/results_CA1_aging_without_pon_male/figures/diff_exp/bar_plots/diff_exp_sig_aggregation_at_least_in_1_age_comps_for_2_brain_region_complete_9x3_flip_left",
    "supp_fig_4a_female": prefix + "miRNA_10p/results_CA1_aging_without_age=26m_28m_without_pon_female/figures/diff_exp/bar_plots/diff_exp_sig_aggregation_at_least_in_1_age_comps_for_2_brain_region_complete_9x3_flip_right",
    "supp_fig_4b_male": prefix + "miRNA_10p/results_CA1_aging_without_pon_male/figures/candidate_comparison_brain_region_split/upset_plots/candidate_comp_sig_diff_exp_both_miRNA_at_least_in_1_age_comps_for_2_brain_region_sig_corr_method=spearman_both_for_2_brain_region_miRNA_corr_th=0.5_intersection_th=1_6x9",
    "supp_fig_4b_female": prefix + "miRNA_10p/results_CA1_aging_without_age=26m_28m_without_pon_female/figures/candidate_comparison_brain_region_split/upset_plots/candidate_comp_sig_diff_exp_both_miRNA_at_least_in_1_age_comps_for_2_brain_region_sig_corr_method=spearman_both_for_2_brain_region_miRNA_corr_th=0.5_intersection_th=1_6x9",
    "supp_fig_4c": prefix + "miRNA_10p/results_CA1_aging/figures/corr_plots/heatmap_complex/all_miRNA_with_age/corr_method=spearman_miRNA_over_age_per_brain_region_filtered_corr_th=0.5",
    "supp_fig_4d": prefix + "miRNA_10p/results_CA1_aging/figures/corr_plots/bar_plots/corr_method=spearman_th=0.5_sig_aggregation_for_2_in_brain_region_complete_9x3_flip_right",
    "supp_fig_4e": prefix + "miRNA_10p/results_CA1_aging/figures/diff_exp/heatmap_complex/num_fc_th_1.5_per_brain_region/diff_exp_global_sig_de_reg_age_per_miRNA_adj",
    "supp_fig_4f": prefix + "miRNA_10p/results_CA1_aging/figures/diff_exp/bar_plots/diff_exp_sig_aggregation_at_least_in_1_age_comps_for_2_brain_region_complete_9x3_flip_right",
    "supp_fig_4g": prefix + "miRNA_10p/results_CA1_aging/figures/fc_vs_expr/scatter_plots/logfc_vs_expr_coloured=brain_region",
    "supp_fig_4h": prefix + "miRNA_10p/results_CA1_aging/figures/corr_plots/overlap_bar_plots/corr_method=spearman_sig_corr_th=0.5_grouped_by_brain_region_age_6x6",
    "supp_fig_4i": prefix + "miRNA_10p/results_CA1_aging/figures/mirna_neighborhood/heatmap_plots/count_heatmap_corr_method=spearman_corr_th=0.5_miRNAs_with_age_12x4.5",

    "supp_fig_5a": prefix + "miRNA_10p/results_CA1_aging/figures/mirna_neighborhood/manhattan_plots/manhattan_corr_method=spearman_corr_th=0.5_miRNAs_with_age_18x28",
    
    "supp_fig_6a": prefix + "miRNA_10p/results_CA1_aging/figures/pathway/bubble_plot_age_correlated/gsea_num=50_age_corr_brain_region=cc_cer_cor_cp_ent_hi_hi2_hy_med_olf_plx_pon_svz_th_vis_18x28",

    "supp_fig_7a": prefix + "miRNA_10p/results_CA1_aging/figures/mfuzz_clustering/k=53/0.15/line_plots/cluster_expression_adj",

    "supp_fig_8a": prefix + "miRNA_10p/results_CA1_aging/figures/mfuzz_clustering/k=53/0.15/line_plots/age_development/mmu-miR-146a-5p_center_lines_for_clusters_11_22_25_28_33_41_51_legend_side=right_9x4",
    "supp_fig_8b": prefix + "miRNA_10p/results_CA1_aging/figures/expr/box_plots/time_corr_method=spearman_poly_deg=2/mmu-miR-146a-5p",
    "supp_fig_8c": prefix + "miRNA_10p/results_CA1_aging/figures/mfuzz_clustering/k=53/0.15/line_plots/age_development/mmu-miR-5100_center_lines_for_clusters_4_18_22_23_25_30_41_51_legend_side=right_9x4",
    "supp_fig_8d": prefix + "miRNA_10p/results_CA1_aging/figures/expr/box_plots/time_single_corr_method=spearman_poly_deg=2/mmu-miR-155-5p/mmu-miR-155-5p_hi2_hy_th_cor_ent_hi_olf_plx_18x4",

    "supp_fig_9a": "results/results_external/paper_10.1523_JNEUROSCI.3582-16.2017/figures/bar_plots/miR-146a-5p",
    "supp_fig_9b": prefix + "miRNA_10p/results_microglia/figures/mapping_statistics/bar_plots/alignment/mapping_stats_6x6",
    "supp_fig_9c": prefix + "ncRNA_10p/results_microglia/figures/mapping_statistics/bar_plots/mapping/grouped_by_age_sample_6x6",
    "supp_fig_9d": prefix + "miRNA_10p/results_CA2_diet_restriction/figures/mapping_statistics/bar_plots/alignment/mapping_stats_grouping_brain_region_treatment_6x6",
    "supp_fig_9e": prefix + "miRNA_10p/results_CA2_young_mouse_plasma/figures/mapping_statistics/bar_plots/alignment/mapping_stats_grouping_brain_region_treatment_6x6",
    "supp_fig_9f": prefix + "ncRNA_10p/results_CA2_diet_restriction/figures/mapping_statistics/bar_plots/mapping/grouped_by_brain_region_sample_6x4.5",
    "supp_fig_9g": prefix + "ncRNA_10p/results_CA2_young_mouse_plasma/figures/mapping_statistics/bar_plots/mapping/grouped_by_brain_region_sample_6x4.5",
    "supp_fig_9h": prefix + "miRNA_10p/results_CA2_diet_restriction/figures/volcano/scatter_plots/single/6x6/treatment__treatment_vs_control.ttest.adj.labels",
    "supp_fig_9i": prefix + "miRNA_10p/results_CA2_diet_restriction/figures/expr/violin_plots/median_brain_region_single/mmu-miR-155-5p/cer_cor_cp_ent_hi_hi2_hy_med_olf_plx_pon_svz_th_vis",
    "supp_fig_9j": prefix + "miRNA_10p/results_CA2_young_mouse_plasma/figures/volcano/scatter_plots/single/6x6/treatment__old mice with young plasma_vs_PBS.ttest.adj.labels",
    "supp_fig_9k": prefix + "miRNA_10p/results_CA2_young_mouse_plasma/figures/expr/violin_plots/median_brain_region_single/mmu-miR-155-5p/cc_cer_cor_cp_hi_med_olf_plx_pon_th_vis",

    "supp_fig_10a": prefix + "iso_10p/results_CA1_aging/figures/iso/line_plots/mmu-miR-155-5p_cor_cp_ent_hi_hi2_hy_med_olf_plx_th",
    "supp_fig_10b": prefix + "miRNA_10p/results_CA1_aging/figures/mmu-miR-155-5p_target_genes_martin_corr/scatter_plots/target_genes_all/spearman/scatter_plot_ADAM23_vs_mmu-miR-155-5p_corr_m=spearman_corr_th=0.3_coloured_by_brain_region_NRCAM_NSG2_PCSK5_REPS2_ADAM23_8x10.5",
    "supp_fig_10c": prefix + "miRNA_10p/results_CA1_aging/figures/mmu-miR-155-5p_target_genes_mtor/corr_method=spearman_mmu-miR-155-5p_per_brain_region_corr_th=0.3_tissue_flipped_10x14",
    "supp_fig_10d": "results/results_external/web_server_aging_b6_proteomics/figures/box_plots/protein_intensity_plot_per_sex_over_age_Mef2a_hippocampus_4x8",
    }

# relevant file types
file_types = ["png", "svg"]
#file_types = ["png", , "svg", "csv", "xlsx"]

# snakemake output folder
source_figures_folder = Path("results/")
# check if it is really a folder
if source_figures_folder.is_dir():
    print(f"{source_figures_folder} exists and is a folder.")
else:
    print(f"{source_figures_folder} is not a folder ... aborting!")
    raise SystemExit("Bye!")

# create output folder if non-existant
pub_figures_folder = Path("publication/")
pub_figures_folder.mkdir(parents=True, exist_ok=True)

# loop over all figures
for pub_figure, source_paths in file_mapping.items():
    # if the dictionary item is not a list (multiple sub-plots), then make one with 1 item
    if not (isinstance(source_paths, list) or isinstance(source_paths, tuple)):
        source_paths = [source_paths]
    for i, source_path in enumerate(source_paths):
        source_path = Path(source_path)
        for file_type in file_types:
            # build full filename with ending and check if it exists
            source_full_path = Path(f"{source_path}.{file_type}")
            if not source_full_path.is_file():
                print(f"{pub_figure}: {source_full_path} does not exist ... skipping!")
                continue
            # only add number if len(source_paths) > 1
            if len(source_paths) > 1:
                target_full_path = pub_figures_folder / Path(f"{pub_figure}_{i+1:.0f}.{file_type}")
            else:
                target_full_path = pub_figures_folder / Path(f"{pub_figure}.{file_type}")
            shutil.copyfile(source_full_path, target_full_path)
            #print(f"{pub_figure} creating {target_full_path} from {source_full_path}")

# remove figures folder
#print(f"removing {source_figures_folder}!")
#shutil.rmtree(source_figures_folder)

