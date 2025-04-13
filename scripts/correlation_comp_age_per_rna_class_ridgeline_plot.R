suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(openxlsx))
#library(corrplot)
suppressPackageStartupMessages(library(Hmisc))

#snakemake = readRDS("snakemake_correlation_comp_age_per_rna_class_ridgeline_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("correlation_comp_age_per_rna_class_ridgeline_plot")


#---------------------------------- Functions ----------------------------------


#---------------------------------- Load data ---------------------------------- 
data_input = snakemake@params$data_input
props = snakemake@params$prop_attribute
method = snakemake@params$corr_methods
corr_ths = snakemake@params$corr_th
adjustment = snakemake@params$p_val_adj
adj_p_value = snakemake@params$adj_p_value
ID = data_input$identifier_column
xticks_names = snakemake@params$xticks_names
plots_props = snakemake@params$plots_props
colors = snakemake@params$colors
plot_size = snakemake@params$plot_size
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_correlation_comp_age_per_rna_class_ridgeline_plot.rds")
#stop()

annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 


#------------------------------------ Script ----------------------------------- 
for (m in method) {
  for (corr_th in corr_ths) {
      
    #for (prop in props) {
    split_group = props[1]
    corr_variable = props[2]
    
    output_folder_fig = sprintf("%s/%s_%s/results_%s/figures/corr_plots/age/ridgeline_plot", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
    dir.create(output_folder_fig, recursive=TRUE)
    output_folder_fig_2 = sprintf("%s/%s_%s/results_%s/figures/corr_plots/age/histo_plot", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
    dir.create(output_folder_fig_2, recursive=TRUE)
    output_folder_tab = sprintf("%s/%s_%s/results_%s/matrices/corr_plots/age/ridgeline_plot", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
    dir.create(output_folder_tab, recursive=TRUE)
    
    corr_values = c()
    corr_values_per_split_group = c()
    aggregated_corr = c()
    for (rna_cl in data_input$rna_classes) {
      print(rna_cl)
      expr = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, rna_cl, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
      fixed_cols = c(colnames(expr)[[1]], annot[[data_input$identifier_column]])
      expr = expr[, ..fixed_cols]
      #print(unique(annot[[data_input$identifier_column]] == colnames(expr)[2:dim(expr)[[2]]]))

      if (rna_cl == "miRNA") {
        feature_col_name = rna_cl
      } else {
        feature_col_name = "RNA"
      }
        
      correlation = c()
      corr_pvalues_adj = c()
      pos_corr = c()
      neg_corr = c()
      sig_pos_corr = c()
      sig_neg_corr = c()
      for(group in unique(annot[[split_group]])){
        groups_ID = annot[annot[[split_group]] == group,][[ID]]
        tmp = colnames(expr) %in% groups_ID
        expr_sort = data.frame(feature=expr[[feature_col_name]], expr[,..tmp], check.names = FALSE)
        
        tmp = annot[[ID]] %in% groups_ID
        corr_var = annot[tmp,][[corr_variable]]
        
        if (!is.numeric(corr_var)) {
          unique_strings = unique(corr_var)
          string_to_number = setNames(seq_along(unique_strings), unique_strings)
          
          numeric_corr_var = string_to_number[corr_var]
        } else {
          numeric_corr_var = corr_var
        }
        
        if (length(corr_var) <= 4) {
          print(sprintf("For group %s, the number of samples is %s (must at least be 5)", group, length(corr_var)))
          print("skip")
        } else {
          correlation[[group]] = as.data.frame(cor(t(expr_sort[,2:dim(expr_sort)[2]]), numeric_corr_var, method = m))
          # if all expr for every timepoint for a sample is 0 for a feature
          # then the correlation is NA
          # we then force the correlation to be 0
          na_cor = is.na(correlation[[group]])
          # uncomment this to show the corresponding lines in the expr_matrix
          # print(expr_sort[,2:dim(expr_sort)[2]][na_cor, ])
          # replace NA with 0
          correlation[[group]][na_cor, ] = 0
          #rownames(correlation[[group]]) = expr_sort[,1]
          
          # get pvalues from rcorr
          # $P gets p-values
          pvalues_tmp = as.data.frame(rcorr(t(expr_sort[,2:dim(expr_sort)[2]]), numeric_corr_var, type=m)$P)
          # select y because x=matrix, y=corr_var
          # and remove last value because it is correlation of corr_var with corr_var = NA
          corr_pvalues = head(pvalues_tmp$y, -1)
          # adjustment
          corr_pvalues_adj[[group]] = as.data.frame(p.adjust(corr_pvalues, method=adjustment))
          # remove all pvalues == NA and replace them with 1 (then it is not significant)
          na_pvalue = is.na(corr_pvalues_adj[[group]])
          corr_pvalues_adj[[group]][na_pvalue] = 1
          #rownames(corr_pvalues_adj[[group]]) = expr_sort[,1]
          # check if anything is still NA
          #print(sprintf("NA in correlation: %s", any(is.na(correlation[[group]]))))
          #print(sprintf("NA in pvalues: %s", any(is.na(corr_pvalues_adj[[group]]))))
          
          pos_corr[[group]] = length(correlation[[group]][correlation[[group]] >= corr_th])
          neg_corr[[group]] = length(correlation[[group]][correlation[[group]] <= -corr_th])
          sig_pos_corr[[group]] = length(correlation[[group]][((correlation[[group]] >= corr_th) & (corr_pvalues_adj[[group]] < adj_p_value))])
          sig_neg_corr[[group]] = length(correlation[[group]][((correlation[[group]] <= -corr_th) & (corr_pvalues_adj[[group]] < adj_p_value))])
        }
      }
    
      correlation_df = as.data.frame(correlation)
      colnames(correlation_df) = unique(annot[[split_group]])
      
      corr_pvalues_df = as.data.frame(corr_pvalues_adj)
      colnames(corr_pvalues_df) = colnames(correlation_df)
      
      # save tables
      correlation_df = data.frame(feature = expr[[feature_col_name]], correlation_df)
      colnames(correlation_df) = gsub("feature", paste(rna_cl), colnames(correlation_df))
      fwrite(correlation_df, sprintf("%s/corr_method=%s_%ss_with_%s_per_%s.csv", output_folder_tab, m, rna_cl, corr_variable, split_group), sep = "\t", row.names = TRUE)
      write.xlsx(correlation_df, sprintf("%s/corr_method=%s_%ss_with_%s_per_%s.xlsx", output_folder_tab, m, rna_cl, corr_variable, split_group), colNames = TRUE, rowNames = TRUE, append = FALSE)
      
      corr_pvalues_df = data.frame(feature = expr[[feature_col_name]], corr_pvalues_df)
      colnames(corr_pvalues_df) = gsub("feature", paste(rna_cl), colnames(corr_pvalues_df))
      fwrite(corr_pvalues_df, sprintf("%s/corr_method=%s_adj_pvalues_%ss_with_%s_per_%s.csv", output_folder_tab, m, rna_cl, corr_variable, split_group), sep = "\t", row.names = TRUE)
      write.xlsx(corr_pvalues_df, sprintf("%s/corr_method=%s_adj_pvalues_%ss_with_%s_per_%s.xlsx", output_folder_tab, m, rna_cl, corr_variable, split_group), colNames = TRUE, rowNames = TRUE, append = FALSE)
     
      corr_values_per_split_group[[rna_cl]] = correlation_df[,2:ncol(correlation_df)]
      
      corr_values[[rna_cl]] = unname(unlist(correlation_df[,2:ncol(correlation_df)]))
      
      aggregated_corr_rna_class = data.frame(group = unique(annot[[split_group]]), num_pos_corr = unname(unlist(pos_corr)), num_neg_corr = unname(unlist(neg_corr)), num_sig_pos_corr = unname(unlist(sig_pos_corr)), num_sig_neg_corr = unname(unlist(sig_neg_corr)))
      tmp = 2:ncol(aggregated_corr_rna_class)
      aggregated_corr_rna_class = rbind(aggregated_corr_rna_class, c("total", colSums(aggregated_corr_rna_class[, tmp])))
      aggregated_corr_rna_class$rna_class = rna_cl
      aggregated_corr[[rna_cl]] = aggregated_corr_rna_class
    }
    
    corr_values_melted = melt(corr_values)
    colnames(corr_values_melted) = c("value", "variable")
    
    mask = (names(unlist(colors$rna_classes)) %in% unique(corr_values_melted$variable))
    colours_rna_class_filtered = unlist(colors$rna_classes)[mask]
    
    median_values = aggregate(value ~ variable, data = corr_values_melted, FUN = median, na.rm = TRUE)
    median_values = median_values[order(-median_values$value), ]  # Sort in descending order
    
    corr_values_melted$variable_fac = factor(corr_values_melted$variable, levels = median_values$variable)
    
    #corr_values_melted$variable_fac = factor(corr_values_melted$variable, levels = rev(names(colors$rna_classes)))
    
    if (m == "spearman") {
      x_axis_title = sprintf("Spearman corr. (%s) per %s", corr_variable, xticks_names$categories[[split_group]])
    } else {
      x_axis_title = sprintf("Pearson corr. (%s) per %s", corr_variable, xticks_names$categories[[split_group]])
    }
    
    ridgeline_plot = ggplot(corr_values_melted, aes(x = value, y = variable_fac, fill = variable_fac)) +
      #geom_density_ridges(scale = 1, alpha = 0.7, size = 0.3, show.legend = FALSE) +  # Ridgeline plot
      geom_density_ridges_gradient(quantile_lines = TRUE, quantiles = 2, size = 0.15, show.legend = FALSE) + 
      geom_vline(xintercept = c(-corr_th, 0, corr_th), linetype = "dashed", color = "grey", size = 0.3) +  # Threshold lines
      scale_y_discrete(expand = expansion(mult = c(0, 0))) +  # Space at the bottom for x-axis
      #scale_x_continuous(limits = c(-1,1), breaks = c(-1, -0.5, 0, 0.5, 1)) +  # Space at the bottom for x-axis
      scale_fill_manual(values = colours_rna_class_filtered) +
      xlab(x_axis_title) +
      theme_classic() +
      theme(
        #legend.position = "none", 
        #legend.title = element_blank(), #element_text(family = plots_props$font_family, size = plots_props$font_size),
        #legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size, margin = margin(0, 0, 0, 0.3, "pt")),
        text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        #axis.text = element_text(angle = 0, family = plots_props$font_family, size = plots_props$font_size),
        #axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
        axis.text.y = element_text(family = plots_props$font_family, size = plots_props$font_size, vjust = -0.5, margin = margin(r = -10)),  # Labels for matrix names
        axis.title.y = element_blank(),  # Remove y-axis title
        axis.text.x = element_text(family = plots_props$font_family, size = plots_props$font_size, margin = margin(t = -2)),
        axis.title.x = element_text(family = plots_props$font_family, size = plots_props$font_size), #, hjust = 0
        axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        #strip.text.y.left = element_text(family = plots_props$font_family, angle = 0, size=plots_props$font_size),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        #panel.spacing = unit(0, "pt"),
        #plot.background = element_rect(fill = "transparent"),
      ) +
      coord_cartesian(clip = "off")  # Prevents cutting off labels
    #ridgeline_plot
  
    
    # Step 1: Compute the median for sorting
    median_values = as.data.table(corr_values_melted)[, .(median_value = median(value)), by = variable]
    
    # Step 2: Order by median and update factor levels
    sorted_levels = median_values[order(median_values$median_value)]$variable
    corr_values_melted$variable_fac = factor(corr_values_melted$variable, levels = sorted_levels)
    
    histo_plot = ggplot(corr_values_melted, aes(x = value, fill = variable_fac)) +
      geom_histogram(bins = 30, alpha = 0.7, position = "identity") +  # Adjust bins if needed
      #geom_vline(data = median_values, aes(xintercept = median_value, group = variable), color = "black", size = 0.7) + # Meidan line in each histo
      facet_wrap(~variable_fac, scales = "free_y", ncol = 1, strip.position = "left") +  # Labels on left
      geom_vline(xintercept = c(-corr_th, 0, corr_th), linetype = "dashed", color = "grey", size = 0.3) +  # Threshold lines
      scale_fill_manual(values = colours_rna_class_filtered) +
      xlab(x_axis_title) +
      theme_classic() +
      theme(
        legend.position = "none",  # Hide legend if not needed
        strip.text.y.left = element_text(family = plots_props$font_family, size=plots_props$font_size, angle = 0, hjust = 1),  # Move facet labels to the left
        strip.placement = "outside",  # Keep labels outside the histograms
        panel.spacing = unit(0.5, "lines"),  # Adjust space between histograms
        #legend.title = element_blank(), #element_text(family = plots_props$font_family, size = plots_props$font_size),
        #legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size, margin = margin(0, 0, 0, 0.3, "pt")),
        text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        #axis.text = element_text(angle = 0, family = plots_props$font_family, size = plots_props$font_size),
        #axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
        axis.text.y = element_text(family = plots_props$font_family, size = plots_props$font_size, vjust = -0.5, margin = margin(r = -10)),  # Labels for matrix names
        axis.title.y = element_blank(),  # Remove y-axis title
        axis.text.x = element_text(family = plots_props$font_family, size = plots_props$font_size, margin = margin(t = -2)),
        axis.title.x = element_text(family = plots_props$font_family, size = plots_props$font_size, hjust = 1),  # Remove x-axis title
        axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        #plot.background = element_rect(fill = "transparent"),
      ) +
      coord_cartesian(clip = "off")  # Prevents cutting off labels
    histo_plot
    
    # save
    width = plot_size$image_width
    height = plot_size$image_height
    ggsave(sprintf("%s/corr_method=%s_corr_th=%s_features_with_%s_per_%s_%sx%s.png", output_folder_fig, m, corr_th, corr_variable, split_group, height, width), ridgeline_plot, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
    ggsave(sprintf("%s/corr_method=%s_corr_th=%s_features_with_%s_per_%s_%sx%s.svg", output_folder_fig, m, corr_th, corr_variable, split_group, height, width), ridgeline_plot, width = width, height = height, unit = plots_props$image_units, dpi = plots_props$dpi)

    height = 2*plot_size$image_height
    ggsave(sprintf("%s/corr_method=%s_corr_th=%s_features_with_%s_per_%s_%sx%s.png", output_folder_fig_2, m, corr_th, corr_variable, split_group, height, width), histo_plot, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
    ggsave(sprintf("%s/corr_method=%s_corr_th=%s_features_with_%s_per_%s_%sx%s.svg", output_folder_fig_2, m, corr_th, corr_variable, split_group, height, width), histo_plot, width = width, height = height, unit = plots_props$image_units, dpi = plots_props$dpi)
    
    aggregated_corr_combined = do.call(rbind, aggregated_corr)
    fwrite(aggregated_corr_combined, sprintf("%s/number_of_corr_corr_method=%s_corr_th=%s_%ss_with_%s_per_%s.csv", output_folder_tab, m, corr_th, rna_cl, corr_variable, split_group), sep = "\t", row.names = FALSE)
    write.xlsx(aggregated_corr_combined, sprintf("%s/number_of_corr_corr_method=%s_corr_th=%s_%ss_with_%s_per_%s.xlsx", output_folder_tab, m, corr_th, rna_cl, corr_variable, split_group), colNames = TRUE, rowNames = FALSE, append = FALSE)
  }
  print(sprintf("Done: %s", m))
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
