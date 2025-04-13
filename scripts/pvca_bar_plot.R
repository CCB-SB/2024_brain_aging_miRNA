suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(plotly))
suppressMessages(library(lme4))
suppressMessages(library(openxlsx))
suppressMessages(library(reshape2))
suppressMessages(library(ggrepel))

#snakemake = readRDS("snakemake_pvca_bar_plot.rds")

snakemake@source("./custom_pvca.R")

data.table::setDTthreads(1)
set.seed(snakemake@params$parameters_props$set_seed)

print("pvca_bar_plot")


#---------------------------------- Functions ----------------------------------
pvca_calculation_and_plotting = function(tbl, annot, props, params, xticks_names, plots_props, output_path_fig, output_path_tab) {
  prop_mapping = make.names(props)
  names(prop_mapping) = props
  
  # keep only expression
  expr = tbl[, sapply(tbl, is.numeric), with=F]
  
  # remove constants
  expr = expr[apply(expr, 1, var) != 0,]
  
  # check if this is logarithmized data
  if((is.null(params$force) || !params.force) && max(expr) >= 100) {
    stop("It seems that your data is not logarithmized. PVCA expects logarithmized data! You can force execution by setting params.force=True")
  }
  
  # Make factors
  annot <- annot[,names(prop_mapping) := lapply(.SD, as.factor), .SDcols = names(prop_mapping)]
  
  # Prepare for PVCA
  #expr.mat <- as.matrix(expr)
  pheno_dat <- data.frame(annot)
  rownames(pheno_dat) <- pheno_dat[[data_input$identifier_column]]
  pheno_dat <- pheno_dat[match(colnames(expr), rownames(pheno_dat)),]

  # Run PVCA
  pvcaObj <- pvcaBatchAssess(expr, pheno_dat, batch.factors=unname(prop_mapping), threshold=params$min_var, threads=1, skip.unique=T)

  # Prepare visualization
  pvcaObj.df <- data.frame(variance=round(pvcaObj$dat, 3)*100, variable=gsub(".", " ", sub("resid", "Residual", pvcaObj$label, fixed=T), fixed=T), variance_label=paste0(as.character(round(pvcaObj$dat , 3)*100), "%"))
  if(!params$keep_zeros){
    pvcaObj.df <- pvcaObj.df[pvcaObj.df$variance != 0,]
  }
  other_var <- sum(pvcaObj.df[pvcaObj.df$variance < 100*params$min_cutoff,]$variance)
  if (params$min_cutoff != 0) {
    pvcaObj.df <- rbind(pvcaObj.df, data.frame(variance=other_var, variable="Other", variance_label=sub(",", ".", paste0(as.character(other_var), "%"))), stringsAsFactors=FALSE)
  }
  pvcaObj.df <- pvcaObj.df[pvcaObj.df$variance >= 100*params$min_cutoff,]
  pvcaObj.df <- pvcaObj.df[order(-pvcaObj.df$variance),]
  pvcaObj.df$variance <- pmax(pvcaObj.df$variance, 1)
  pvcaObj.df$variable <- factor(pvcaObj.df$variable, levels=rev(pvcaObj.df$variable))
  
  # remove possible _ in names with xticksnames from the config file
  variable_name = c()
  for (var in pvcaObj.df$variable) {
    if (var != "Residual") {
      if (grepl(":", var)) {
        var_contain = strsplit(var, ":")[[1]]
        var_contain_list = c()
        for (var_c in var_contain) {
          tmp = strsplit(xticks_names$categories[[sprintf("%s_capital", var_c)]], " \\(")[[1]][1]
          var_contain_list = append(var_contain_list, tmp)
        }
        var_contain_list = paste(var_contain_list, collapse = ":")
        variable_name = append(variable_name, var_contain_list)
      } else {
        tmp = strsplit(xticks_names$categories[[sprintf("%s_capital", var)]], " \\(")[[1]][1]
        variable_name = append(variable_name, tmp) 
      }
    } else {
      variable_name = append(variable_name, "Residual")
    }
  }
  
  pvcaObj.df$variable_name = variable_name
  pvcaObj.df$variable_name <- factor(pvcaObj.df$variable_name, levels=rev(pvcaObj.df$variable_name))
  
  # Create plot
  pvca_plot = ggplot(pvcaObj.df, aes(x = variance, y = variable_name)) + #aes(variable, variance)
    geom_bar(stat = "identity", width = 0.75, fill = "#D5A021") +
    #geom_point(shape = 21, fill="lightgrey", color="black", size = 1, stroke = 1) +
    scale_x_log10(limits = c(1, 900), expand = c(0, 0)) + #scale_y_log10(limits = c(1, 900), expand = c(0, 0))
    geom_text(aes(label=variance_label), nudge_x=log10(3.6), size = plots_props$font_size/11.04*3.88, angle = 0, nudge_y=0) + #geom_text(aes(label=variance_label), nudge_x=0, nudge_y=log10(3.6), size=params$plots_props$font_size/11.04*3.88, angle = 90)
    labs(x="Observed variance (%)", y="") + #labs(x="", y="Observed variance (%, log10)")
    theme_classic() +
    theme(aspect.ratio = (plots_props$image_width / plots_props$image_width), 
          axis.text.x = element_text(size = plots_props$font_size), #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=params$plots_props$font_size)
          legend.position="none", 
          text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
          plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header),
          legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size)
    )
  
  # save figures
  ggsave(sprintf("%s.png", output_path_fig), pvca_plot, dpi=plots_props$dpi, width = ((2 * plots_props$image_width) / 3), height=plots_props$image_height, units = plots_props$image_units)
  ggsave(sprintf("%s.svg", output_path_fig), pvca_plot, width = ((2 * plots_props$image_width) / 3), height=plots_props$image_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
  
  ggsave(sprintf("%s_9x6.png", output_path_fig), pvca_plot, dpi=plots_props$dpi, width = plots_props$image_width, height=plots_props$image_height, units = plots_props$image_units)
  ggsave(sprintf("%s_9x6.svg", output_path_fig), pvca_plot, width = plots_props$image_width, height=plots_props$image_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
  
  # Basic piechart
  #pie_chart = ggplot(pvcaObj.df, aes(x = variance, y = variable_name)) +
  #  geom_bar(stat="identity", width=1, color="#D5A021") +
  #  coord_polar("y", start=0) +
  #  scale_x_log10(limits = c(1, 900), expand = c(0, 0)) +
  #  theme_void() + 
  #  theme(legend.position="none") +
  #  geom_text(aes(label=variance_label), size = plots_props$font_size/11.04*3.88) +
  #  labs(x="Observed variance (%)", y="") +
  #  scale_fill_brewer(palette="Set1")
  #  theme_classic() +
  #  theme(aspect.ratio = (plots_props$image_width / plots_props$image_width), 
  #        axis.text.x = element_text(size = plots_props$font_size), #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=params$plots_props$font_size)
  #        legend.position="none", 
  #        text = element_text(family = plots_props$font_family, size = plots_props$font_size),
  #        axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
  #        axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
  #        plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header),
  #        legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
  #        legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size)
  #  )
   
  #pie_chart = ggplot(pvcaObj.df, aes(x = "", y = variance, fill = variable_name)) +
  #  geom_bar(stat = "identity", width = 1, color = "white") +
  #  coord_polar("y", start = 0) +
  #  theme_void() + 
  #  theme(legend.title = element_blank()) +
  #  labs(title = "Proportion of Variance Explained (PVCA)") +
  #  geom_text(aes(label = variance_label), position = position_stack(vjust = 0.5), size = 3.5) +
  #  scale_fill_brewer(palette = "Set3")
  
  #height = plots_props$image_height
  #width = plots_props$image_width
  #ggsave(sprintf("%s_pie_chart_%sx%s.png", output_path_fig, height, width), pie_chart, dpi=plots_props$dpi, width=width , height=height, units = plots_props$image_units)
  #ggsave(sprintf("%s_pie_chart_%sx%s.svg", output_path_fig, height, width), pie_chart, width=width, height=height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
    
  #pvca_plot_ply = ggplotly(pvca_plot) %>% toWebGL()
  #saveRDS(pvca_plot_ply, file=sprintf("%s/pvca.rds", output_path_fig))
  
  # save pvca values to table
  pvcaObj.df$variable_name = c()
  if (output_path_tab != "none") {
    fwrite(pvcaObj.df, sprintf("%s.csv", output_path_tab), sep = "\t", row.names = TRUE)
    write.xlsx(pvcaObj.df, sprintf("%s.xlsx", output_path_tab), colNames = TRUE, rowNames = TRUE, append = FALSE)
  }
  
  return(pvcaObj.df)
}

pvca_per_tissue_plotting = function(pvca_data_table, props, params, plot_params, xticks_names, plots_props, colors, output_path_fig) {
  feature_x = props[1] 
  feature_y = props[2]
  
  dot_size_variable = sprintf("%s.%s", feature_x, feature_y)
  tmp_x = strsplit(xticks_names$categories[[sprintf("%s_capital", feature_x)]], " \\(")[[1]][1]
  tmp_y = strsplit(xticks_names$categories[[sprintf("%s_capital", feature_y)]], " \\(")[[1]][1]
  legend_title_point_size = sprintf("%s:%s", tmp_x, tmp_y)
  dot_labels = unlist(xticks_names[[params$tissue]])[as.character(pvca_data_table[[params$tissue]])]
  
  scatter_pvca_plot = ggplot(pvca_data_table, aes_string(x=feature_x, y=feature_y)) + 
    geom_point(aes_string(size=dot_size_variable, colour = factor(pvca_data_table[[params$colour]]))) +
    scale_color_manual(values=colors) +
    geom_text_repel(label=dot_labels, size = (1/2.82) * plots_props$font_size, family = plots_props$font_family) + #pvca_data_table[[params$tissue]]
    labs(x=sprintf("%s (observed variance (%%))", xticks_names$categories[[sprintf("%s_capital", feature_x)]]), y=sprintf("%s (observed variance (%%))", xticks_names$categories[[sprintf("%s_capital", feature_y)]])) +
    guides(color = FALSE) +
    theme_classic() +
    ylim(min(pvca_data_table[[feature_y]]), max(pvca_data_table[[feature_y]]) + 1) +
    theme(#aspect.ratio = (plots_props$image_width / plots_props$image_width), 
          legend.position="bottom", 
          text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
          plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header),
          legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
          legend.margin=margin(t = -8), legend.key.size = unit(0.3, "cm"), legend.spacing.y = unit(0.4, "cm")
    ) +
    guides(size=guide_legend(legend_title_point_size))
  
  if (length(plot_params) > 0){
    scatter_pvca_plot = scatter_pvca_plot + geom_vline(xintercept = plot_params$vline, linetype="dashed", color = "grey") +
                                            geom_hline(yintercept = plot_params$hline, linetype="dashed", color = "grey")
  }
  
  # save figures
  ggsave(sprintf("%s_%sx%s.png", output_path_fig, plots_props$image_width, plots_props$image_height), scatter_pvca_plot, dpi=plots_props$dpi, width = plots_props$image_width, height=plots_props$image_height, units = plots_props$image_units)
  ggsave(sprintf("%s_%sx%s.svg", output_path_fig, plots_props$image_width, plots_props$image_height), scatter_pvca_plot, width = plots_props$image_width, height=plots_props$image_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
  
  ggsave(sprintf("%s_%sx%s.png", output_path_fig, plots_props$image_height, plots_props$image_height), scatter_pvca_plot, dpi=plots_props$dpi, width = plots_props$image_height, height=plots_props$image_height, units = plots_props$image_units)
  ggsave(sprintf("%s_%sx%s.svg", output_path_fig, plots_props$image_height, plots_props$image_height), scatter_pvca_plot, width = plots_props$image_height, height=plots_props$image_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
}


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
props = snakemake@params$props
params = snakemake@params$parameters
plot_params = snakemake@params$plot_parameters
plots_props = snakemake@params$plots_props
xticks_names = snakemake@params$xticks_names
colors = unlist(snakemake@params$colors[[params$colour]])
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_pvca_bar_plot.rds")

# load data files
tbl = fread(sprintf("%s_%s/%s_%s_quantification_%s_log%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, params$log, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 

# creating output folder
output_folder_fig_bar = sprintf("%s/%s_%s/results_%s/figures/pvca_plots/bar_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_fig_bar, recursive=TRUE)
output_folder_tab = sprintf("%s/%s_%s/results_%s/matrices/pvca", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_tab, recursive=TRUE)


#------------------------------------ Script ----------------------------------- 
if (params$tissue == "all_brain_regions") {
  # create output paths
  output_path_fig = sprintf("%s/%s", output_folder_fig_bar, params$tissue)
  output_path_tab = sprintf("%s/%s", output_folder_tab, params$tissue)
  
  pvcaObj.df = pvca_calculation_and_plotting(tbl, annot, props, params, xticks_names, plots_props, output_path_fig, output_path_tab)
  
  pvcaObj.df$variance_label = c()
  pvcaObj.df_reordered = transpose(pvcaObj.df, make.names = "variable")
  pvca_data_table = data.frame(pvcaObj.df_reordered)
  
  #print(sprintf("Done: %s", params$tissues))
  
} else {
  output_folder_fig_scatter = sprintf("%s/%s_%s/results_%s/figures/pvca_plots/scatter_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
  dir.create(output_folder_fig_scatter, recursive=TRUE)
  
  pvca_data_table = c()
  for (ti in unique(annot[[params$tissue]])) {
    # filter expr for single tissue samples
    tissue_ids = annot[annot[[params$colour]] == ti,][[data_input$identifier_column]]
    tmp = colnames(tbl) %in% tissue_ids
    tbl_tissue = tbl[,..tmp]
        
    #annot_tissue = annot[annot[[params$colour]] == ti,]
    #annot_tissue$mouse_ID = droplevels(annot_tissue$mouse_ID)
    #annot_tissue$mouse_ID = factor(annot_tissue$mouse_ID)
    
    output_path_fig = sprintf("%s/%s", output_folder_fig_bar, ti)
    pvcaObj.df = pvca_calculation_and_plotting(tbl_tissue, annot, props, params, xticks_names, plots_props, output_path_fig, "none")
    
    pvcaObj.df$variance_label = c()
    pvcaObj.df_reordered = transpose(pvcaObj.df, make.names = "variable")
    pvcaObj.df_reordered_df = data.frame(placeholder = ti, pvcaObj.df_reordered) 
    names(pvcaObj.df_reordered_df)[names(pvcaObj.df_reordered_df) == "placeholder"] = paste0(params$tissue) 
    
    pvca_data_table = rbind(pvca_data_table, pvcaObj.df_reordered_df)
    
    print(sprintf("Done: %s", ti))
  }
  
  # save pvca values for all tissues 
  fwrite(pvca_data_table, sprintf("%s/combined_%s_results.csv", output_folder_tab, params$colour), sep = "\t", row.names = TRUE)
  write.xlsx(pvca_data_table, sprintf("%s/combined_%s_results.xlsx", output_folder_tab, params$colour), colNames = TRUE, rowNames = TRUE, append = FALSE)
  
  output_path_fig = sprintf("%s/combined_%s_results_%s", output_folder_fig_scatter, params$colour, paste(props, collapse = "_"))
  pvca_per_tissue_plotting(pvca_data_table, props, params, plot_params, xticks_names, plots_props, colors, output_path_fig)
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])

