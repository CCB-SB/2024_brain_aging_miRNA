suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))


#snakemake = readRDS("snakemake_expression_analysis.rds")

data.table::setDTthreads(1)
set.seed(snakemake@params$parameters_props$set_seed)

print("expression_analysis")


################################## Functions ###################################

expr_histo_density_plot= function(expr, th, colour_list, save_path, dpi, image_width, image_height){
  expr_melt = melt(expr[,-"miRNA"])
  expr_melt$group = str_split_fixed(expr_melt$variable, "_", 3)[,2]
  expr_melt$value_log10 = log10(expr_melt$value)
  
  p = ggplot(expr_melt, aes(x = value_log10, color = group)) +
    scale_color_manual(values = colour_list) +
    geom_vline(xintercept = log10(th), color = highlight_colour, size = 1) +
    xlab("expression value (log10)") +
    ylab("occurence") +
    ggtitle("Density") +
    #labs(color = "samples") +
    facet_grid(variable ~ .) +
    geom_density(alpha = .2) +
    theme_classic()
  #coord_cartesian(xlim = c(0,1000))
  dir.create(save_path, recursive=TRUE)
  ggsave(sprintf("%s/density_log10_th=%s.png", save_path, th), p, dpi = dpi, width = image_width, height = ncol(expr)*2, units = "cm")
  
  p = ggplot(expr_melt, aes(x = value_log10, color = group)) +
    geom_histogram(position = "stack", binwidth = 1) +
    geom_vline(xintercept = log10(th), color = highlight_colour, size = 1) +
    scale_fill_manual(values = colour_list) +
    xlab("expression value (log10)") +
    ylab("occurence") +
    ggtitle("Histogram") +
    #labs(color = "samples") +
    facet_grid(variable ~ .) +
    #coord_cartesian(xlim = c(0,25)) +
    theme_classic() + 
    theme(legend.position = "none")
  ggsave(sprintf("%s/histo_log10_th=%s.png", save_path, th), p, dpi = dpi, width = image_width, height = ncol(expr)*2, units = "cm")
}

sample_fewest_reads = function(raw){
  samples = colnames(raw[,-"miRNA"])
  reads = c()
  for (sample in samples){
    reads = append(reads, colSums(raw[,..sample]))
  }
  sample_min = append(names(which.min(reads)), min(reads))
  
  return(sample_min)
}

expr_histo_density_plot_raw_samples_fewest_reads = function(raw, sample_min, th, save_path, dpi, image_width, image_height){
  background_noise_th_raw = th * as.numeric(sample_min[2]) / 1e6
  
  tmp = colnames(raw) %in% sample_min[1]
  raw_melt = melt(raw[, ..tmp][, -"miRNA"])
  p = ggplot(raw_melt, aes(x = value, color = variable)) +
    geom_density(alpha = .2) +
    #scale_color_manual(values = colour_list) +
    geom_vline(xintercept = log10(background_noise_th_raw), color = highlight_colour, size = 0.5) +
    xlab("expression value (log10)") +
    ylab("occurence") +
    ggtitle("Density") +
    #labs(color = "samples") +
    #facet_grid(variable ~ .) +
    theme_classic() + 
    theme(legend.position = "bottom")
  #coord_cartesian(xlim = c(0,1000))
  dir.create(save_path, recursive=TRUE)
  ggsave(sprintf("%s/density_raw_sample_fewest_reads_log10_th=%s.png", save_path, background_noise_th_raw), p, dpi = dpi, width = image_width, height = image_height, units = "cm")
  
  p = ggplot(raw_melt, aes(x = value, color = variable)) +
    geom_histogram(position = "identity", binwidth = 1) +
    geom_vline(xintercept = background_noise_th_raw, color = highlight_colour, size = 0.5) +
    #scale_fill_manual(values = colour_list) +
    xlab("expression value (log10)") +
    ylab("occurence") +
    ggtitle("Histogramm") +
    #labs(color = "samples") +
    #facet_grid(variable ~ .) +
    #geom_density(alpha = .2) +
    coord_cartesian(xlim = c(0,25)) +
    theme_classic() +
    theme(legend.position = "bottom")
  ggsave(sprintf("%s/histo_raw_sample_fewest_reads_log10_th=%s.png", save_path, th), p, dpi = dpi, width = image_width, height = image_height, units = "cm")
}

log_expr = function(expr, basis, save_path){
  #expr_log2 = log2(expr[,-"miRNA"])
  expr_log = log2(expr[,-"miRNA"] + 1) / log2(basis)
  #expr_log[expr_log == -Inf] <- 0
  #expr_log[expr_log == -NaN] <- 0
  expr_log$miRNA = as.character(expr$miRNA)
  #expr_log$miRNA = paste(expr_log$miRNA, "", sep="")
  setcolorder(expr_log, c("miRNA", colnames(expr_log[,-"miRNA"])))
  dir.create(save_path, recursive=TRUE)
  fwrite(expr_log, sprintf("results/miRNA_filtered_quantification_rpmm_norm_log%s.csv", save_path, basis), row.names = FALSE, sep='\t')
}

box_plot_expr_thresholded = function(expr, th_list, save_path, dpi, image_height, image_width){
  plot_df = melt(expr)
  plot_df$th = rep(0, each = dim(plot_df)[1])
  for (th in th_list) {
    expr_bin = ifelse(expr >= log10(th), 1, 0)
    expr_thresholded = expr * expr_bin
    expr_thresholded_melted = melt(expr_thresholded)
    expr_thresholded_melted$th = rep(sprintf("log10(%s)", th), each = dim(expr_thresholded_melted)[1])
    plot_df = rbind(plot_df, expr_thresholded_melted)
  }
  plot_df[plot_df$value == 0]$value <- NA
  
  p <- ggplot(plot_df, aes(x=variable, y=value, group=variable)) + 
    geom_boxplot(aes(fill=variable)) +
    scale_fill_manual(values = sample_colors) +
    xlab("samples") +
    ylab("expression values (rpmm, log10)") +
    #theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    facet_grid(. ~ th)
  dir.create(save_path, recursive=TRUE)
  ggsave(sprintf("%s/box_plot_log10_without_0.png", save_path), p, dpi = dpi, width = image_width, height = image_height, units = "cm")
}

background_noise_filtering = function(expr, th, detection_rate){
  expr_bin = ifelse(expr[,-"miRNA"] >= th, 1, 0)
  expr_filtered_detect = expr[,-"miRNA"] * expr_bin
  #expr_filtered$miRNA = expr$miRNA
  
  if (length(unique(rowMeans(expr_bin) >= detection_rate/100)) == 2) {
    expr_filtered = expr_filtered_detect[rowMeans(expr_bin) >= detection_rate/100,]
  } else if ((length(unique(rowMeans(expr_bin) >= detection_rate/100) == 1)) & (unique(rowMeans(expr_bin) >= detection_rate/100) == TRUE)) {
    expr_filtered = expr_filtered_detect[rowMeans(expr_bin) >= detection_rate/100,]
  }
  
  return(expr_filtered)
}

miRNAs_per_sample = function(expr){
  samples = colnames(expr[,-"miRNA"])
  sample_names = str_split_fixed(samples, "_", 2)[,2]
  sample_mapping = data.frame(sample_names = sample_names, samples=samples)
  list_miRNAs = c()
  for (sample in samples) {
    tmp = expr[,..sample]
    tmp_1 = tmp[,..sample] > 0
    miRNA_names = expr$miRNA
    # tmp_2 = sample_mapping[sample_mapping$samples == sample,]$sample_names
    
    
    list_miRNAs[[sample]] = miRNA_names[tmp_1]
  }
  
  return(list_miRNAs)
}

num_miRNAs_per_sample = function(list_miRNAs){
  samples = names(list_miRNAs)
  num_miRNAs = c()
  for (sample in samples) {
    num_miRNAs = append(num_miRNAs, length(list_miRNAs[sample][[1]]))
  }
  num_miRNAs_df = data.frame(samples = samples, num_miRNAs = num_miRNAs)
  
  return(num_miRNAs_df)
}

bar_plot_num_miRNAs = function(annot, num_miRNAs, th, colour_list, save_path, dpi, image_height, image_width, background_noise_th){
  tissue_list = c()
  for (samp in num_miRNAs$samples) {
    tissue_list = append(tissue_list, annot[annot$ID == samp,][[tissue]])
  }
  num_miRNAs$group = tissue_list
  
  dir.create(save_path, recursive=TRUE)
  p = ggplot(data=num_miRNAs, aes(x=samples, y=num_miRNAs, fill=group)) + 
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values=colour_list) + 
    geom_text(aes(label=num_miRNAs), position=position_dodge(width=0), vjust = 0.5, hjust = 1.1, angle = 90, size = 1) +
    xlab("") +
    ylab(sprintf("number of miRNAs\n(expression >= %srpmm)", th)) +
    theme_classic() + 
    theme(legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 3))
  
  ggsave(sprintf("%s/expressed_miRNA_per_sample_th=%srpmm.png", save_path, background_noise_th), p, dpi = dpi, height = image_height, width = image_width, unit="cm")
}

upset_plot_miRNAs_per_sample = function(miRNA_lists, th, save_path, dpi, image_width, image_height){
  plot_table = make_comb_mat(miRNA_lists)
  plot_table_thresholded = plot_table[comb_size(plot_table) >= th]
  
  cs = comb_size(plot_table_thresholded)
  nc = ncol(plot_table_thresholded)
  #plot_table_thresholded = plot_table_thresholded[1:(dim(plot_table_thresholded)[1]-1),]
  
  dir.create(save_path, recursive=TRUE)
  png(sprintf("%s/threshold=%s.png", save_path, th), res=dpi, width=image_width, height=image_height)
  ht = ComplexHeatmap::UpSet(plot_table_thresholded, comb_order = order(comb_size(plot_table_thresholded), decreasing = TRUE),
                             top_annotation = HeatmapAnnotation("Number of\n intersecting\n miRNAs\n" = anno_barplot(comb_size(plot_table_thresholded), ylim = c(0, max(comb_size(plot_table_thresholded))*1.1),border = FALSE, gp = gpar(fill = "black"), 
                                                                                                                     labels_gp = gpar(col = "black", fontsize = 10), 
                                                                                                                     height = unit(4, "cm")), annotation_name_side = "left", annotation_name_rot = 90),
                             #top_annotation = upset_top_annotation(plot_table_thresholded, add_numbers = TRUE),
                             left_annotation = rowAnnotation("Number of miRNAs\n in the samples" = anno_barplot(-set_size(plot_table_thresholded), baseline = -0, axis_param = list(at = c(0, -200, -400, -600, -800), labels = c(0, 200, 400, 600,800), labels_rot = 0), 
                                                                                                                border = FALSE, gp = gpar(fill = sample_colors, col = sample_colors), width = unit(4, "cm"))),
                             right_annotation = NULL)
  
  #decorate_annotation("Intersection\nsize", {grid.text(cs[co], x = 1:nc, y = unit(cs[co], "native") + unit(1, "mm"), gp = gpar(fontsize = 5), just = "bottom", default.units = "native")})
  ht = draw(ht); co = column_order(ht)
  dev.off()
}


#---------------------------------- Load data ---------------------------------- 
data_input = snakemake@params$data_input
tissue = snakemake@params$tissue
background_noise_th_list = snakemake@params$th_rpmmm
xticks_names = snakemake@params$xticks_names
plots_props = snakemake@params$plots_props
colours = snakemake@params$colors
results_folder = snakemake@params$results_folder

expr = fread(sprintf("%s/data_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
raw_count = fread(sprintf("%s_%s/%s_%s_quantification_raw.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering), sep='\t', header=T)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 

# save rdata
#saveRDS(snakemake, file = "snakemake_expression_analysis.rds")

output_folder_fig = sprintf("%s/%s_%s/results_%s/figures/expressed_per_%s/analysis", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, tissue)
dir.create(output_folder_fig, recursive=TRUE)
output_folder_tab = sprintf("%s/%s_%s/results_%s/matrices/expressed_per_%s/analysis", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, tissue)
dir.create(output_folder_tab, recursive=TRUE)


################################## Parameters ##################################

# mathematical parameters
intersecting_th = 5
box_plot_th_list = rep(5, dim(expr)[1])

highlight_colour = "#E61549"


#################################### Script ####################################
for (ti in unique(annot[[tissue]])){
  print(ti)
  
  save_path_plots = sprintf("%s/%s", output_folder_fig, ti)
  dir.create(sprintf(save_path_plots), recursive=TRUE)
  save_path_tables = sprintf("%s/%s", output_folder_tab, ti)
  dir.create(sprintf(save_path_tables), recursive=TRUE)
  
  samples = annot[annot[[tissue]] == ti,][[data_input$identifier_column]]
  tmp = colnames(expr) %in% c("miRNA", samples)
  expr_ti = raw_count[,..tmp]
  tmp = colnames(raw_count) %in% c("miRNA", samples)
  raw_count_ti = raw_count[,..tmp]
  background_noise_th = background_noise_th_list[ti]
  
  # to check if background_noise_th is set appropriated
  # density and histogram plot of expression values per sample
  expr_histo_density_plot(expr_ti, unlist(background_noise_th), unlist(colours[[tissue]]), sprintf("%s/background_noise_correction_plots", save_path_plots), plots_props$dpi, plots_props$image_width, plots_props$image_height)
  
  # samples with the lowest number of reads analysis of the selected threshold 
  sample_min = sample_fewest_reads(raw_count_ti)
  
  # density and histogram plot of raw expression values of the sample with the lowest reads 
  expr_histo_density_plot_raw_samples_fewest_reads(raw_count_ti, sample_min, unlist(background_noise_th), sprintf("%s/background_noise_correction_plots", save_path_plots), plots_props$dpi, plots_props$image_width, plots_props$image_height)
  
  # compare different background noise thresholds with log10(expression)
  # expr_log10 = log_expr(expr_filtered, 10, FALSE)
  # box_plot_expr_thresholded(expr_log10, box_plot_th_list, sprintf("%s/box_plots", save_path_plots), dpi, image_height, image_width)
  
  # removes background noise by thresholding rpmm normalised data with a threshold
  expr_filtered = background_noise_filtering(expr_ti, background_noise_th)
  
  # log data for snc_analysis-pipeline
  # rpmm normalised data
  # log_expr(expr, 2, sprintf("%s", save_path_tables))
  # filtered rpmm normalised data
  # log_expr(expr_filtered, 2, sprintf("%s", save_path_tables))
  
  # generates a list of lists containing for every sample which miRNA is expressed (>0, due to filtering >= 40)
  miRNA_lists = miRNAs_per_sample(expr_filtered)
  
  # calculates how many miRNAs per sample are expressed
  num_miRNAs = num_miRNAs_per_sample(miRNA_lists)
  dir.create(sprintf("%s/num_miRNAs", save_path_tables), recursive=TRUE)
  fwrite(num_miRNAs, sprintf("%s/num_miRNAs/th_expressed_miRNAs_th=%srpmm.csv", sprintf("%s", save_path_tables), background_noise_th))
  
  # Bar plot for expressed miRNAs 
  bar_plot_num_miRNAs(annot, num_miRNAs, background_noise_th, unlist(colours[[tissue]]), sprintf("%s/bar_plots/th", save_path_plots), plots_props$dpi, plots_props$image_height, plots_props$image_width, background_noise_th)
  
  # Upset plot for the expressed miRNAs per sample
  # upset_plot_miRNAs_per_sample(miRNA_lists, intersecting_th, sprintf("%s/upset_plots", save_path_plots), dpi, image_width, image_height)
  
  
  
  # un-thresholded 
  # generates a list of lists containing for every sample which miRNA is expressed (>0, due to filtering >= 40)
  miRNA_lists = miRNAs_per_sample(expr_ti)
  
  # calculates how many miRNAs per sample are expressed
  num_miRNAs = num_miRNAs_per_sample(miRNA_lists)
  dir.create(sprintf("%s/num_miRNAs", save_path_tables), recursive=TRUE)
  fwrite(num_miRNAs, sprintf("%s/num_miRNAs/expressed_miRNAs_in_spike-in.csv", save_path_tables))
  
  # Bar plot for expressed miRNAs 
  bar_plot_num_miRNAs(annot, num_miRNAs, background_noise_th, unlist(colours[[tissue]]), sprintf("%s/bar_plots/expressed", save_path_plots), plots_props$dpi, plots_props$image_height, plots_props$image_width, background_noise_th)
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
