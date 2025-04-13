suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(latex2exp))
suppressPackageStartupMessages(library(openxlsx))


#snakemake = readRDS("snakemake_external_bar_plot.rds")

set.seed(snakemake@params$parameters_porps$set_seed)

print("external_bar_plot") 


#---------------------------------- Load data ---------------------------------- 
raw_data_path = snakemake@params$raw_data_path
id_col = snakemake@params$id_col
#rna = snakemake@params$rna
plot_list = snakemake@params$plot_list
mean_ct_col_mapping = unlist(snakemake@params$columns_mean)
fc_col_mapping = unlist(snakemake@params$columns_fc)
#plot_titles = snakemake@params$plot_titles
xticks_names = snakemake@params$xticks_names
plots_props = snakemake@params$plots_props
colours = snakemake@params$colors
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_external_bar_plot.rds")

df = fread(raw_data_path, sep=',', header=T)

bar_colour =  "#335F70"


# --------------------------------- Script -------------------------------------
# output folder
output_folder_bar = sprintf("%s/figures/bar_plots", results_folder)
dir.create(output_folder_bar, recursive=TRUE)
  
for (i in 1:length(plot_list)) {
      plot_info = plot_list[[i]]

      rna = plot_info$rna
      plot_titles = plot_info$plot_titles

      # means
      # get cols
      tmp = (colnames(df) %in% c(id_col, names(mean_ct_col_mapping)))
      df_mean = df[,..tmp]  
      # get feature
      df_mean_feature = df_mean[df_mean$`Our name` == rna,]
      df_mean_feature$`Our name` = c()
      # exchange colnames by names for the plot
      colnames(df_mean_feature) = unname(mean_ct_col_mapping[colnames(df_mean_feature)])

      df_mean_feature_melt = melt(df_mean_feature)
      colnames(df_mean_feature_melt) = c("variable", "mean_value")
      df_mean_feature_melt$mean_value_power_2 = 2^df_mean_feature_melt$mean_value
      df_mean_feature_melt$mean_value_log10 = log10(df_mean_feature_melt$mean_value)

      # fc
      # get cols
      tmp = (colnames(df) %in% c(id_col, names(fc_col_mapping)))
      df_fc = df[,..tmp]  
      # get feature
      df_fc_feature = df_fc[df_fc$`Our name` == rna,]
      df_fc_feature$`Our name` = c()
      # exchange colnames by names for the plot
      colnames(df_fc_feature) = unname(fc_col_mapping[colnames(df_fc_feature)])

      plot_df = melt(df_fc_feature)
      colnames(plot_df) = c("variable", "fc_value")

      # add the property which the fc was calculated to as 1
      plot_df = rbind(plot_df, data.frame(variable = setdiff(unname(mean_ct_col_mapping), unname(fc_col_mapping)), fc_value = log2(1)))
      
      plot_df$fc_value_power_2 = 2^(plot_df$fc_value)

      plot_df = merge(plot_df, df_mean_feature_melt, by = "variable")

      # reproduce the ordering of the bars
      plot_df$variable_factored = factor(plot_df$variable, levels = rev(unname(mean_ct_col_mapping)))

      # bar plot
      p = ggplot(plot_df, aes(x=variable_factored, y=fc_value_power_2)) + #, fill=tissue 
      geom_bar(stat="identity", fill=bar_colour, width=0.5) +
      #scale_fill_manual(values=colours[[tissue]]) +
      scale_y_continuous(expand = c(0,0)) +
      xlab(plot_titles$x_lab) +
      ylab(plot_titles$y_lab) +
      ggtitle(plot_titles$plot_title)  +
      coord_flip() +
      theme_classic() +
      theme(legend.position="none",
            text = element_text(family = plots_props$font_family, size = plots_props$font_size),
            axis.text = element_text(angle = 0, family = plots_props$font_family, size = plots_props$font_size),
            axis.text.x = element_text(family = plots_props$font_family, size = plots_props$font_size), #, angle = 90, vjust = 0.5, hjust=1 
            axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
            plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header),
            legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
            plot.margin = unit(c(0.1,0.5,0.1,-0.3), "cm"), #t, r, b, l
      )

      file_name = rna
      plot_width = 6
      plot_height = 6
      ggsave(file = sprintf("%s/%s.svg", output_folder_bar, file_name), plot = p, width = plot_width, height = plot_height, unit = plots_props$image_units, dpi = plots_props$dpi)
      ggsave(file = sprintf("%s/%s.png", output_folder_bar, file_name), plot = p, width = plot_width, height = plot_height, unit = plots_props$image_units, dpi = plots_props$dpi)     
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
