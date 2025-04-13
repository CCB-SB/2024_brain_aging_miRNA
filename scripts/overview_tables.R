suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(openxlsx))


# ----- create overview table for all data sets including the number of samples after sequencing and filtering and number of features after mapping and filtering -----
# load data
expr_filtered = c()
expr = c()

# CA1
expr_filtered[["CA1"]] = fread("../../data_CA1/without_mid/sncrna-pipeline-input_CA1/miRNA_filtered_quantification_rpmmm_norm.detection_rate_10p_per_brain_region.csv", sep='\t', colClasses=c(ID="character")) 
expr[["CA1"]] = fread("../../data_CA1/without_mid/sncrna-pipeline-input_CA1/miRNA_complete_quantification_rpmm_norm.csv", sep='\t', colClasses=c(ID="character")) 

# CA2
expr_filtered[["CA2_diet_restriction"]] = fread("../../data_CA2/without_mid/sncrna-pipeline-input_CA2_Experiment=diet restriction_quantity_control/miRNA_filtered_quantification_rpmm_norm.detection_rate_10p_per_brain_region.csv", sep='\t', colClasses=c(ID="character")) 
expr[["CA2_diet_restriction"]] = fread("../../data_CA2/without_mid/sncrna-pipeline-input_CA2_Experiment=diet restriction/miRNA_complete_quantification_rpmm_norm.csv", sep='\t', colClasses=c(ID="character")) 

expr_filtered[["CA2_injection"]] = fread("../../data_CA2/without_mid/sncrna-pipeline-input_CA2_Experiment=young mouse plasma_injection type=retro-orbital_quantity_control/miRNA_filtered_quantification_rpmm_norm.detection_rate_10p_per_brain_region.csv", sep='\t', colClasses=c(ID="character")) 
expr[["CA2_injection"]] = fread("../../data_CA2/without_mid/sncrna-pipeline-input_CA2_Experiment=young mouse plasma_injection type=retro-orbital/miRNA_complete_quantification_rpmm_norm.csv", sep='\t', colClasses=c(ID="character")) 

# microglia
expr_filtered[["microglia"]] = fread("../../data_microglia/sncrna-pipeline-input_microglia/miRNA_filtered_quantification_rpmm_norm.detection_rate_10p_per_age.csv", sep='\t', colClasses=c(ID="character")) 
expr[["microglia"]] = fread("../../data_microglia/sncrna-pipeline-input_microglia/miRNA_complete_quantification_rpmm_norm.csv", sep='\t', colClasses=c(ID="character")) 

info_table = c()
for (data_set in names(expr_filtered)) {
  info_table = rbind(info_table, c(data_set, dim(expr[[data_set]])[2] - 1, dim(expr_filtered[[data_set]])[2] - 1, dim(expr[[data_set]])[1], dim(expr_filtered[[data_set]])[1]))
}

colnames(info_table) = c("data_set", "sequenced_samples", "filtered_samples", "mapped_miRNAs", "filtered_miRNAs")

# save table
fwrite(info_table, "results_rpmm/samples_filtering_alignment_against_mmu10_over_2e6_reads/without_mid/miRNA_10p/overview_table.csv", sep = "\t", row.names = FALSE)
write.xlsx(info_table, "results_rpmm/samples_filtering_alignment_against_mmu10_over_2e6_reads/without_mid/miRNA_10p/overview_table.xlsx", colNames = TRUE, rowNames = FALSE, append = FALSE)


# ----- create overview table for all data sets including the samples we used for analysis, their matadata information and the information resulting from the alignment and mapping
# load data

annot = c()
alignment_info = c() # alignment against the mouse genome
mapping_info = c() # mapping of the aligned reads against the RNA types especially miRNAs


# CA1
annot[["CA1"]] = fread("../../data_CA1/without_mid/sncrna-pipeline-input_CA1/annotation_CA1_filtered.csv", sep='\t', colClasses=c(ID="character")) 
annot[["CA1"]]$V1 = c()
annot[["CA1"]]$sampleID = c()
annot[["CA1"]][["Sample ID"]] = c()
annot[["CA1"]]$ID = c()

alignment_info[["CA1"]] = fread("../../data_CA1_CA2/mapping_infos/mapping_vs_mm10_mis1.csv", sep='\t', colClasses=c(ID="character")) 
mask = (alignment_info[["CA1"]]$fastq_name == "BA10_CA1_cp_26_1")
alignment_info[["CA1"]] = alignment_info[["CA1"]][!mask,]
tmp = (alignment_info[["CA1"]]$fastq_name %in% annot[["CA1"]]$fastq_name)
alignment_info[["CA1"]] = alignment_info[["CA1"]][tmp,]

mapping_info[["CA1"]] = fread("../../data_CA1_CA2/mapping_infos/overall_composition.detailed_mod.csv")
mask = mapping_info[["CA1"]]$Sample == "BA10_CA1_cp_26_1"
mapping_info[["CA1"]] = mapping_info[["CA1"]][!mask,]
tmp = (mapping_info[["CA1"]]$Sample %in% annot[["CA1"]]$fastq_name)
mapping_info[["CA1"]] = mapping_info[["CA1"]][tmp,]

#rna_composition_perc = fread("../../data_CA1_CA2/mapping_infos/overall_composition.detailed_mod.mapped_perc.csv")
#overall_composition = fread("../../data_CA1_CA2/mapping_infos/overall_composition.csv")
#overall_composition_rownames = overall_composition$V1
#overall_composition$V1 = c()
#overall_composition_colnames = colnames(overall_composition)
#overall_composition = transpose(overall_composition)
#colnames(overall_composition) = overall_composition_rownames
#rownames(overall_composition) = overall_composition_colnames


# CA2
annot[["CA2_diet_restriction"]] = fread("../../data_CA2/without_mid/sncrna-pipeline-input_CA2_Experiment=diet restriction_quantity_control/annotation_CA2_Experiment=diet restriction_quantity_control_filtered.csv", sep='\t', colClasses=c(ID="character")) 
annot[["CA2_diet_restriction"]]$V1  = c()
annot[["CA2_diet_restriction"]]$sampleID = c()
annot[["CA2_diet_restriction"]][["Sample ID"]] = c()
annot[["CA2_diet_restriction"]]$ID = c()

alignment_info[["CA2_diet_restriction"]] = fread("../../data_CA1_CA2/mapping_infos/mapping_vs_mm10_mis1.csv", sep='\t', colClasses=c(ID="character")) 
tmp = (alignment_info[["CA2_diet_restriction"]]$fastq_name %in% annot[["CA2_diet_restriction"]]$fastq_name)
alignment_info[["CA2_diet_restriction"]] = alignment_info[["CA2_diet_restriction"]][tmp,]

mapping_info[["CA2_diet_restriction"]] = fread("../../data_CA1_CA2/mapping_infos/overall_composition.detailed_mod.csv")
tmp = (mapping_info[["CA2_diet_restriction"]]$Sample %in% annot[["CA2_diet_restriction"]]$fastq_name)
mapping_info[["CA2_diet_restriction"]] = mapping_info[["CA2_diet_restriction"]][tmp,]

annot[["CA2_injection"]] = fread("../../data_CA2/without_mid/sncrna-pipeline-input_CA2_Experiment=young mouse plasma_injection type=retro-orbital_quantity_control/annotation_CA2_Experiment=young mouse plasma_injection type=retro-orbital_quantity_control_filtered.csv", sep='\t', colClasses=c(ID="character")) 
annot[["CA2_injection"]]$V1 = c()
annot[["CA2_injection"]]$sampleID = c()
annot[["CA2_injection"]][["Sample ID"]] = c()
annot[["CA2_injection"]]$ID = c()

alignment_info[["CA2_injection"]] = fread("../../data_CA1_CA2/mapping_infos/mapping_vs_mm10_mis1.csv", sep='\t', colClasses=c(ID="character")) 
tmp = (alignment_info[["CA2_injection"]]$fastq_name %in% annot[["CA2_injection"]]$fastq_name)
alignment_info[["CA2_injection"]] = alignment_info[["CA2_injection"]][tmp,]

mapping_info[["CA2_injection"]] = fread("../../data_CA1_CA2/mapping_infos/overall_composition.detailed_mod.csv")
tmp = (mapping_info[["CA2_injection"]]$Sample %in% annot[["CA2_injection"]]$fastq_name)
mapping_info[["CA2_injection"]] = mapping_info[["CA2_injection"]][tmp,]

# microglia
annot[["microglia"]] = fread("../../data_microglia/sncrna-pipeline-input_microglia/annotation_microglia_filtered.csv", sep='\t', colClasses=c(ID="character")) 
annot[["microglia"]]$Name = c()
annot[["microglia"]]$ID = c()
annot[["microglia"]]$group = c()

alignment_info[["microglia"]] = fread("../../data_microglia/mapping_infos/mapping_vs_mm10_mis1.csv", sep='\t', colClasses=c(ID="character")) 
tmp = (alignment_info[["microglia"]]$fastq_name %in% annot[["microglia"]]$fastq_name)
alignment_info[["microglia"]] = alignment_info[["microglia"]][tmp,]

mapping_info[["microglia"]] = fread("../../data_microglia/mapping_infos/overall_composition.detailed_mod.csv")
tmp = (mapping_info[["microglia"]]$Sample %in% annot[["microglia"]]$fastq_name)
mapping_info[["microglia"]] = mapping_info[["microglia"]][tmp,]

info_table = c()
for (data_set in names(annot)) {
  keep_cols = c("fastq_name", "reads_processed", "reads_aligned")
  keep_cols_names = c("fastq_name", "reads_sequenced", "reads_aligned_against_mmu")
  tmp = (colnames(alignment_info[[data_set]]) %in% keep_cols)
  filtered_alignment_info = alignment_info[[data_set]][,..tmp]
  filtered_alignment_info = filtered_alignment_info[,..keep_cols]
  colnames(filtered_alignment_info) = keep_cols_names
  annot[[data_set]] = merge(annot[[data_set]], filtered_alignment_info, by = "fastq_name")
  
  keep_cols = c("fastq_name", "total assigned reads", "miRNA")
  keep_cols_names = c("fastq_name", "aligned_without_suppressed_m=100", "reads_mapped_against_miRNA")
  mapping_info[[data_set]]$fastq_name = mapping_info[[data_set]]$Sample 
  tmp = (colnames(mapping_info[[data_set]]) %in% keep_cols)
  filtered_mapping_info = mapping_info[[data_set]][,..tmp]
  filtered_mapping_info = filtered_mapping_info[,..keep_cols]
  colnames(filtered_mapping_info) = keep_cols_names
  annot[[data_set]] = merge(annot[[data_set]], filtered_mapping_info, by = "fastq_name")
  
  annot[[data_set]]$data_set = data_set
  
  info_table = rbind(info_table, annot[[data_set]], fill=TRUE)
  #info_table = rbind(info_table, c(data_set, dim(expr[[data_set]])[2] - 1, dim(expr_filtered[[data_set]])[2] - 1, dim(expr[[data_set]])[1], dim(expr_filtered[[data_set]])[1]))
}

mask = info_table[, c("data_set", setdiff(names(info_table), "data_set"))]
info_table = info_table[,..mask]

# save table
fwrite(info_table, "results_rpmm/samples_filtering_alignment_against_mmu10_over_2e6_reads/without_mid/miRNA_10p/metadata_table.csv", sep = "\t", row.names = FALSE)
write.xlsx(info_table, "results_rpmm/samples_filtering_alignment_against_mmu10_over_2e6_reads/without_mid/miRNA_10p/metadata_table.xlsx", colNames = TRUE, rowNames = FALSE, append = FALSE)


## ----- create overview table for all data sets including the number of matadata information 
## load data
#metadata_overview = c()
#
## CA1
#metadata_overview[["CA1_sex"]] = fread("../../data_CA1/without_mid/sncrna-pipeline-input_CA1/sample_overview/overview_brain_region_sex_SampleID_total_828_SampleIDs.csv", sep='\t', colClasses=c(ID="character")) 
#metadata_overview[["CA1_age"]] = fread("../../data_CA1/without_mid/sncrna-pipeline-input_CA1/sample_overview/overview_brain_region_age_SampleID_total_828_SampleIDs.csv", sep='\t', colClasses=c(ID="character")) 
#metadata_overview[["CA1_age"]]$all = c()
#
#metadata_overview[["CA1_age_male"]] = fread("../../data_CA1/without_mid/sncrna-pipeline-input_CA1_male/sample_overview/overview_brain_region_age_SampleID_total_479_SampleIDs.csv", sep='\t', colClasses=c(ID="character")) 
#metadata_overview[["CA1_age_male"]]$all = c()
#metadata_overview[["CA1_age_female"]] = fread("../../data_CA1/without_mid/sncrna-pipeline-input_CA1_without_age=26m_28m_female/sample_overview/overview_brain_region_age_SampleID_total_349_SampleIDs.csv", sep='\t', colClasses=c(ID="character")) 
#metadata_overview[["CA1_age_female"]]$all = c()
#
## CA2
#metadata_overview[["CA2_diet_restriction_sex"]] = fread("../../data_CA2/without_mid/sncrna-pipeline-input_CA2_Experiment=diet restriction_quantity_control/sample_overview/overview_brain_region_sex_SampleID_total_113_SampleIDs.csv", sep='\t', colClasses=c(ID="character")) 
#metadata_overview[["CA2_diet_restriction_experiment"]] = fread("../../data_CA2/without_mid/sncrna-pipeline-input_CA2_Experiment=diet restriction_quantity_control/sample_overview/overview_brain_region_treatment_SampleID_total_113_SampleIDs.csv", sep='\t', colClasses=c(ID="character")) 
#
#metadata_overview[["CA2_injection_sex"]] = fread("../../data_CA2/without_mid/sncrna-pipeline-input_CA2_Experiment=young mouse plasma_injection type=retro-orbital_quantity_control/sample_overview/overview_brain_region_sex_SampleID_total_68_SampleIDs.csv", sep='\t', colClasses=c(ID="character")) 
#metadata_overview[["CA2_injection_experiment"]] = fread("../../data_CA2/without_mid/sncrna-pipeline-input_CA2_Experiment=young mouse plasma_injection type=retro-orbital_quantity_control/sample_overview/overview_brain_region_treatment_SampleID_total_68_SampleIDs.csv", sep='\t', colClasses=c(ID="character")) 
#
## microglia
#metadata_overview[["microglia_age"]] = fread("../../data_microglia/sncrna-pipeline-input_microglia/sample_overview/overview_age_group_SampleID_total_8_SampleIDs.csv", sep='\t', colClasses=c(ID="character")) 


## ----- load diff exp table for all, male and female samples separate and save only relevant columns
## load data
# all
diff_exp = fread("results_rpmm/samples_filtering_alignment_against_mmu10_over_2e6_reads/without_mid/miRNA_10p/results_CA1/matrices/diff_exp/diff_exp_log2.csv", sep='\t', colClasses=c(ID="character")) 

diff_exp_filtered = data.table(feature=diff_exp[["miRNA"]])
colnames_list = c("miRNA")                   
for (br in unique(annot$CA1$brain_region)) {
  prop = "brain_region"
  sub_string = sprintf("_%s=%s", prop, br)
  comp_prop = "age"
  
  col_part1 = sprintf("num_%s__3%s", comp_prop, sub_string)
  col_part2 = sprintf("median_%s__3%s", comp_prop, sub_string)
  plot_df_part1 = data.table(
    num_control=diff_exp[[col_part1]], 
    median=diff_exp[[col_part2]]
  )
  diff_exp_filtered = cbind(diff_exp_filtered, plot_df_part1)
  
  colnames_list = append(colnames_list, col_part1)
  colnames_list = append(colnames_list, col_part2)
  
  for (age in c(12, 15, 18, 21, 26, 28)) {
    
    plot_df_part2 = data.table(num_case=diff_exp[[sprintf("num_%s__%s%s", comp_prop, age, sub_string)]],
                               median=diff_exp[[sprintf("median_%s__%s%s", comp_prop, age, sub_string)]],
                               ttest_adj_pval=diff_exp[[sprintf("ttest_adjp_%s__%s_vs_3%s", comp_prop, age, sub_string)]],
                               fc=diff_exp[[sprintf("fc_%s__%s_vs_3%s", comp_prop, age, sub_string)]]
    )
    diff_exp_filtered = cbind(diff_exp_filtered, plot_df_part2)
    
    colnames_list = append(colnames_list, sprintf("num_%s__%s%s", comp_prop, age, sub_string))
    colnames_list = append(colnames_list, sprintf("median_%s__%s%s", comp_prop, age, sub_string))
    colnames_list = append(colnames_list, sprintf("ttest_adjp_%s__%s_vs_3%s", comp_prop, age, sub_string)) 
    colnames_list = append(colnames_list, sprintf("fc_%s__%s_vs_3%s", comp_prop, age, sub_string))
  }
}
colnames(diff_exp_filtered) = colnames_list

# save table
fwrite(diff_exp_filtered, "results_rpmm/samples_filtering_alignment_against_mmu10_over_2e6_reads/without_mid/miRNA_10p/results_CA1/matrices/diff_exp/diff_exp_log2_filtered.csv", sep = "\t", row.names = FALSE)
write.xlsx(diff_exp_filtered, "results_rpmm/samples_filtering_alignment_against_mmu10_over_2e6_reads/without_mid/miRNA_10p/results_CA1/matrices/diff_exp/diff_exp_log2_filtered.xlsx", colNames = TRUE, rowNames = FALSE, append = FALSE)

# male
diff_exp_male = fread("results_rpmm/samples_filtering_alignment_against_mmu10_over_2e6_reads/without_mid/miRNA_10p/results_CA1_without_pon_male/matrices/diff_exp/diff_exp_log2.csv", sep='\t', colClasses=c(ID="character")) 

diff_exp_filtered_male = data.table(feature=diff_exp_male[["miRNA"]])
colnames_list = c("miRNA")                   
for (br in unique(annot$CA1$brain_region)) {
  if (br == "pon") {
    next
  }
  prop = "brain_region"
  sub_string = sprintf("_%s=%s", prop, br)
  comp_prop = "age"
  
  col_part1 = sprintf("num_%s__3%s", comp_prop, sub_string)
  col_part2 = sprintf("median_%s__3%s", comp_prop, sub_string)
  plot_df_part1 = data.table(
                       num_control=diff_exp_male[[col_part1]], 
                       median=diff_exp[[col_part2]]
                       )
  diff_exp_filtered_male = cbind(diff_exp_filtered_male, plot_df_part1)

  colnames_list = append(colnames_list, col_part1)
  colnames_list = append(colnames_list, col_part2)
  
  for (age in c(12, 15, 18, 21, 26, 28)) {
    
    plot_df_part2 = data.table(num_case=diff_exp[[sprintf("num_%s__%s%s", comp_prop, age, sub_string)]],
                               median=diff_exp_male[[sprintf("mean_%s__%s%s", comp_prop, age, sub_string)]],
                               ttest_adj_pval=diff_exp_male[[sprintf("ttest_adjp_%s__%s_vs_3%s", comp_prop, age, sub_string)]],
                               fc=diff_exp_male[[sprintf("fc_%s__%s_vs_3%s", comp_prop, age, sub_string)]]
                              )
    diff_exp_filtered_male = cbind(diff_exp_filtered_male, plot_df_part2)
    
    colnames_list = append(colnames_list, sprintf("num_%s__%s%s", comp_prop, age, sub_string))
    colnames_list = append(colnames_list, sprintf("mean_%s__%s%s", comp_prop, age, sub_string))
    colnames_list = append(colnames_list, sprintf("ttest_adjp_%s__%s_vs_3%s", comp_prop, age, sub_string)) 
    colnames_list = append(colnames_list, sprintf("fc_%s__%s_vs_3%s", comp_prop, age, sub_string))
  }
}
colnames(diff_exp_filtered_male) = colnames_list

# save table
fwrite(diff_exp_filtered_male, "results_rpmm/samples_filtering_alignment_against_mmu10_over_2e6_reads/without_mid/miRNA_10p/results_CA1_without_pon_male/matrices/diff_exp/diff_exp_log2_filtered.csv", sep = "\t", row.names = FALSE)
write.xlsx(diff_exp_filtered_male, "results_rpmm/samples_filtering_alignment_against_mmu10_over_2e6_reads/without_mid/miRNA_10p/results_CA1_without_pon_male/matrices/diff_exp/diff_exp_log2_filtered.xlsx", colNames = TRUE, rowNames = FALSE, append = FALSE)

# female
# load data
diff_exp_female = fread("results_rpmm/samples_filtering_alignment_against_mmu10_over_2e6_reads/without_mid/miRNA_10p/results_CA1_without_age=26m_28m_female/matrices/diff_exp/diff_exp_log2.csv", sep='\t', colClasses=c(ID="character")) 

diff_exp_filtered_female = data.table(feature=diff_exp_female[["miRNA"]])
colnames_list = c("miRNA")                   
for (br in unique(annot$CA1$brain_region)) {
  if (br == "pon") {
    next
  }
  prop = "brain_region"
  sub_string = sprintf("_%s=%s", prop, br)
  comp_prop = "age"
  
  col_part1 = sprintf("num_%s__3%s", comp_prop, sub_string)
  col_part2 = sprintf("median_%s__3%s", comp_prop, sub_string)
  plot_df_part1 = data.table(
    num_control=diff_exp_female[[col_part1]],
    median=diff_exp[[col_part2]]
  )
  diff_exp_filtered_female = cbind(diff_exp_filtered_female, plot_df_part1)
  
  colnames_list = append(colnames_list, col_part1)
  colnames_list = append(colnames_list, col_part2)
  
  for (age in c(12, 15, 18, 21)) {
    
    plot_df_part2 = data.table(num_case=diff_exp[[sprintf("num_%s__%s%s", comp_prop, age, sub_string)]],
                               median=diff_exp_female[[sprintf("median_%s__%s%s", comp_prop, age, sub_string)]],
                               ttest_adj_pval=diff_exp_female[[sprintf("ttest_adjp_%s__%s_vs_3%s", comp_prop, age, sub_string)]],
                               fc=diff_exp_female[[sprintf("fc_%s__%s_vs_3%s", comp_prop, age, sub_string)]]
    )
    diff_exp_filtered_female = cbind(diff_exp_filtered_female, plot_df_part2)
    
    colnames_list = append(colnames_list, sprintf("num_%s__%s%s", comp_prop, age, sub_string))
    colnames_list = append(colnames_list, sprintf("median_%s__%s%s", comp_prop, age, sub_string))
    colnames_list = append(colnames_list, sprintf("ttest_adjp_%s__%s_vs_3%s", comp_prop, age, sub_string)) 
    colnames_list = append(colnames_list, sprintf("fc_%s__%s_vs_3%s", comp_prop, age, sub_string))
  }
}
colnames(diff_exp_filtered_female) = colnames_list

# save table
fwrite(diff_exp_filtered_female, "results_rpmm/samples_filtering_alignment_against_mmu10_over_2e6_reads/without_mid/miRNA_10p/results_CA1_without_age=26m_28m_without_pon_female/matrices/diff_exp/diff_exp_log2_filtered.csv", sep = "\t", row.names = FALSE)
write.xlsx(diff_exp_filtered_female, "results_rpmm/samples_filtering_alignment_against_mmu10_over_2e6_reads/without_mid/miRNA_10p/results_CA1_without_age=26m_28m_without_pon_female/matrices/diff_exp/diff_exp_log2_filtered.xlsx", colNames = TRUE, rowNames = FALSE, append = FALSE)


# write correlation form mmu-miR-155-5p and its target genes in one excel

