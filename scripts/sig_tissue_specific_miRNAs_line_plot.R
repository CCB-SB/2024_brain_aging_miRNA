suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
#suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(ggrepel))

#snakemake = readRDS("snakemake_sig_tissue_specific_miRNAs_line_plot.rds")

data.table::setDTthreads(1)
set.seed(snakemake@params$parameters_props$set_seed)

print("sig_tissue_specific_miRNAs_line_plot")


#---------------------------------- Functions ----------------------------------
reashape_features = function(data, rna_class) {

  diff_exp_clean = data$diff_exp[data$diff_exp != ""]
  corr_clean = data$corr[data$corr != ""]
  
  merged_list = unique(c(diff_exp_clean, corr_clean))
  result = data.frame("rna_class_column" = merged_list, "method"="tmp")
  colnames(result)[1] = rna_class

  for (rna in merged_list) {
    if ((rna %in% diff_exp_clean == TRUE) & (rna %in% corr_clean == TRUE)) {
      result[result[[rna_class]] == rna,]$method = "both"
    } else if ((rna %in% diff_exp_clean == TRUE) & (rna %in% corr_clean == FALSE)) {
      # If only one value is empty, add the non-empty value to the result list
      result[result[[rna_class]] == rna,]$method = "diff_exp"
    } else {
      result[result[[rna_class]] == rna,]$method = "corr"
    }
  }
  
  return(result)
}

line_plot = function(median_table, feature_distribution_reshaped, rna_class, data_input, output_folder_path_fig, cluster_range, threads, plot_props){
  
  # Combine data for ggplot
  median_table_zscored_melted = melt(median_table)#, id.vars = "row_names")
  
  #data_plot = merge(median_table_zscored_melted, acore_dt, by="Var1")
  data_plot = merge(median_table_zscored_melted, feature_distribution_reshaped, by = "feature", all = TRUE)
  
  #colnames(data_plot) = gsub(paste(rna_class), "feature", colnames(data_plot))
  
  #data_plot$cluster = factor(data_plot$method)
  method_colours = c("diff_exp" = "#d62728", "both" = "#9467bd", "corr" = "#1f77b4")
  
  p = ggplot(data_plot, aes(x = variable, y = value, group = feature)) + 
    geom_line(aes(color=method)) +
    geom_text_repel(aes(label = feature), data = data_plot[data_plot$variable == sort(unique(data_plot$variable))[1],], color = "black", size = 3) +
    ##scale_linetype_manual(values=c("diff_exp" = "solid", "both" = "dashed", "corr" = "dotted")) + 
    scale_colour_manual(values=method_colours) +
    #scale_x_continuous(breaks = sort(unique(data_plot$variable))) +#, labels = 1:10) 
    xlab("") + 
    ylab(sprintf("standardised expression\nvalues (%s)", data_input$norm)) +
    theme_classic() +
    theme(legend.position="bottom", 
          legend.margin = margin(t = -15, r = 0, b = 0, l = 0),
          text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          axis.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          axis.title = element_text(family = plot_props$font_family, size = plot_props$font_size),
          plot.title = element_text(family = plot_props$font_family, size = plot_props$font_size_header),
          legend.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          legend.title = element_text(family = plot_props$font_family, size = plot_props$font_size)
    ) +
    guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 2))
  p
  
  # save bar plot
  ggsave(sprintf("%s.png", output_folder_path_fig), p, dpi=plot_props$dpi, width=plot_props$image_width, height=plot_props$image_height, units = plot_props$image_units)
  ggsave(sprintf("%s.svg", output_folder_path_fig), p, width=plot_props$image_width, height=plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
}


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
plot_category = snakemake@params$plot_category
diff_exp_props = snakemake@params$diff_exp_props
corr_props = snakemake@params$corr_props
colours = snakemake@params$colors
plot_props = snakemake@params$plots_props
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_sig_tissue_specific_miRNAs_line_plot.rds")

expr = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
diff_exp = fread(sprintf("%s/%s_%s/results_%s/matrices/diff_exp/diff_exp%s.csv",results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, diff_exp_props$diff_log), sep='\t') 

expr_rownames = expr[[data_input$rna_class]]
expr[[data_input$rna_class]] = c()

folder_path_tab = sprintf("%s/%s_%s/results_%s/matrices/sig_tissue_specific_comp", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)

output_folder_path_median_fig = sprintf("%s/%s_%s/results_%s/figures/sig_tissue_specific_comp/line_plots/median", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_path_median_fig, recursive=TRUE)
output_folder_path_mean_fig = sprintf("%s/%s_%s/results_%s/figures/sig_tissue_specific_comp/line_plots/mean", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_path_mean_fig, recursive=TRUE)


#------------------------------------ Script ----------------------------------- 
direction = list(c("both", "both"), c("up", "pos"), c("down", "neg"))

for (tissue in unique(annot[[plot_category[1]]])){
  for (num_dereg_combinations in diff_exp_props$num_dereg_combinations_th) {
    for (num_dereg_tissues in diff_exp_props$num_dereg_tissues_th) {
      for (m in corr_props$method) {
        for (th in corr_props$corr_th) {
          for (num_corr_tissues in corr_props$num_corr_tissues_th) {
            for (dir in direction) {
              file_name = sprintf("sig_%s_specific_comparison_diff_exp_%s_%s_at_least_in_%s_%s_comps_for_%s_%s_corr_method=%s_%s_for_%s_%s_%s_corr_th=%s", tissue, dir[1], data_input$rna_class, num_dereg_combinations, plot_category[2], num_dereg_tissues, plot_category[1], m, dir[2], num_corr_tissues, plot_category[1], data_input$rna_class, th)
              matching_files = list.files(sprintf("%s/venn_plots", folder_path_tab), pattern = sprintf("%s.csv", file_name))
              
              print(file_name)
              
              if (length(matching_files) > 0) {
                feature_distribution = fread(sprintf("%s/venn_plots/%s.csv", folder_path_tab, file_name), sep='\t')
                feature_list = unique(c(feature_distribution$diff_exp, feature_distribution$corr)[nzchar(c(feature_distribution$diff_exp, feature_distribution$corr))])
                
                # filter expr with feature_list
                diff_exp_filtered = diff_exp[diff_exp[[data_input$rna_class]] %in% feature_list,]
                
                median_table = data.table(feature=diff_exp_filtered[[data_input$rna_class]])
                #median_table_log2 = data.table(feature=diff_exp_filtered[[data_input$rna_class]])
                mean_table = data.table(feature=diff_exp_filtered[[data_input$rna_class]])
                for (tp in paste(sort(unique(annot[[plot_category[2]]])))) {
                  median_table[[tp]] = 2^(as.numeric(diff_exp_filtered[[sprintf("median_%s__%s_%s=%s", plot_category[2], tp, plot_category[1], tissue)]]))
                  #median_table_log2[[tp]] = as.numeric(diff_exp_filtered[[sprintf("median_%s__%s_%s=%s", plot_category[2], tp, plot_category[1], tissue)]])
                  mean_table[[tp]] = 2^(as.numeric(diff_exp_filtered[[sprintf("mean_%s__%s_%s=%s", plot_category[2], tp, plot_category[1], tissue)]]))
                }
                
                if (dim(median_table)[1] < 4) {
                  print(sprintf("Not enough features for %s: %s", plot_category[1], tissue))
                  next
                } else if ((unique(is.na(feature_distribution$diff_exp)) == TRUE) || (unique(is.na(feature_distribution$corr)) == TRUE)) {
                  if (unique(is.na(feature_distribution$diff_exp)) == TRUE) {
                    diff_exp_or_corr = "diff_exp"
                  } else {
                    diff_exp_or_corr = "corr"
                  }
                  print(sprintf("Not enough features for %s: %s and method: %s", plot_category[1], tissue, diff_exp_or_corr))
                  next
                }
                
                median_table_df = data.frame(median_table, check.names = FALSE)
                #median_table_log2_df = data.frame(median_table_log2, check.names = FALSE)
                mean_table_df = data.frame(mean_table, check.names = FALSE)
                
                # rearrange features
                feature_distribution_reshaped = reashape_features(feature_distribution, "feature")

                file_name_output = sprintf("sig_%s_specific_comparison_diff_exp_%s_%s_at_least_in_%s_%s_comps_for_%s_%s_corr_method=%s_%s_for_%s_%s_%s_corr_th=%s", tissue, dir[1], data_input$rna_class, num_dereg_combinations, plot_category[2], num_dereg_tissues, plot_category[1], m, dir[2], num_corr_tissues, plot_category[1], data_input$rna_class, th)
                output_folder_path_median_fig_tissue = sprintf("%s/%s", output_folder_path_median_fig, file_name_output)
                line_plot(median_table_df, feature_distribution_reshaped, data_input$rna_class, data_input, output_folder_path_median_fig_tissue, cluster_range, threads, plot_props)
              
                #file_name_output = sprintf("sig_%s_specific_comparison_diff_exp_%s_%s_at_least_in_%s_%s_comps_for_%s_%s_corr_method=%s_%s_for_%s_%s_%s_corr_th=%s_log2", tissue, dir[1], data_input$rna_class, num_dereg_combinations, plot_category[2], num_dereg_tissues, plot_category[1], m, dir[2], num_corr_tissues, plot_category[1], data_input$rna_class, th)
                #output_folder_path_median_fig_tissue = sprintf("%s/%s", output_folder_path_median_fig, file_name_output)
                #line_plot(median_table_log2_df, feature_distribution_reshaped, data_input$rna_class, data_input, output_folder_path_median_fig_tissue, cluster_range, threads, plot_props)
            
                file_name_output = sprintf("sig_%s_specific_comparison_diff_exp_%s_%s_at_least_in_%s_%s_comps_for_%s_%s_corr_method=%s_%s_for_%s_%s_%s_corr_th=%s", tissue, dir[1], data_input$rna_class, num_dereg_combinations, plot_category[2], num_dereg_tissues, plot_category[1], m, dir[2], num_corr_tissues, plot_category[1], data_input$rna_class, th)
                output_folder_path_mean_fig_tissue = sprintf("%s/%s", output_folder_path_mean_fig, file_name_output)
                line_plot(mean_table_df, feature_distribution_reshaped, data_input$rna_class, data_input, output_folder_path_mean_fig_tissue, cluster_range, threads, plot_props)
                
              } else {
                print(sprintf("No features for %s: %s", plot_category[1], tissue))
                next
              }
            }
          }
        }
      }
    }
  }  
}
  

#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])

