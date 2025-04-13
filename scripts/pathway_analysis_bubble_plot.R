suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))

#snakemake = readRDS("snakemake_pathway_analysis_bubble_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("pathway_analysis_bubble_plot")


#---------------------------------- Functions ----------------------------------


#---------------------------------- Load data ---------------------------------- 
data_input = snakemake@params$data_input
sig_lvl = snakemake@params$sig_lvl
methods = snakemake@params$method
comps = snakemake@params$comparisons
plots_props = snakemake@params$plots_props
xticks_names = snakemake@params$xticks_names
colors = snakemake@params$colors
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_pathway_analysis_bubble_plot.rds")

output_folder_fig = sprintf("%s/%s_%s/results_%s/figures/pathway/bubble_plot", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, data_input$rna_class)
dir.create(output_folder_fig, recursive=TRUE)
input_folder_tab = sprintf("%s/%s_%s/results_%s/matrices/mieaa", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)


#------------------------------------ Script ----------------------------------- 
for (m in methods) {
  tissues = c()
  pathway_results = c()
  for(i in 1:length(comps)){
    c = comps[[i]]
    prop = names(comps)[i]
    c$g1 = c$g[1]
    c$g2 = c$g[2]
    
    sub_string = ifelse(!is.null(c$subset), sprintf("_%s=%s", c$subset[1], c$subset[2]), "")
    tissues = append(tissues, c$subset[2])
    
    pathway_table = fread(sprintf("data/mieaa_CA1_without_old_sex/gsea/sex/%s_vs_%s%s_all/GO_Annotations_indirect_mature.txt", c$g1, c$g2, sub_string)) #Male_vs_Female_brain_region=cc_all
  
    colnames(pathway_table) = gsub("#Name", "Pathway", colnames(pathway_table))
    colnames(pathway_table) = gsub("P-value", "pvalue", colnames(pathway_table))
    
    tmp = (pathway_table$pvalue < sig_lvl)
    pathway_table_sig = pathway_table[tmp, ]
    
    pathway_table_sig$Regulation_direction = gsub("0", "depleted", pathway_table_sig$Regulation_direction)
    pathway_table_sig$Regulation_direction = gsub("1", "enriched", pathway_table_sig$Regulation_direction)
    
    pathway_table_sig$Pathway_short = str_split_fixed(pathway_table_sig$Pathway, " GO", 2)[,1]

    pathway_table_sig$tissue = rep(c$subset[2], length(pathway_table_sig$Regulation_direction))
    
    pathway_results = rbind(pathway_results, pathway_table_sig)
  }
  
  pathway_results$Regulation_direction_sig = gsub("depleted", -1, pathway_results$Regulation_direction)
  pathway_results$Regulation_direction_sig = gsub("enriched", 1, pathway_results$Regulation_direction_sig)
  
  pathway_results$pvalue_neglog10 = -log10(as.numeric(pathway_results$pvalue))
  
  pathway_results$sort_col = -log10(as.numeric(pathway_results$pvalue)) * as.numeric(pathway_results$Regulation_direction_sig)
  order_idx = order(pathway_results$sort_col, decreasing = TRUE)
  pathway_results_sorted = pathway_results[order_idx,]
  
  if (dim(pathway_results_sorted[pathway_results_sorted$sort_col > 0,])[1] >= 10 ) {
    first_rows = head(pathway_results_sorted, 10)
  } else {
    first_rows = pathway_results_sorted[pathway_results_sorted$sort_col > 0,]
  }
  if (dim(pathway_results_sorted[pathway_results_sorted$sort_col < 0,])[1] >= 10 ) {
    last_rows = tail(pathway_results_sorted, 10)
  } else {
    last_rows = pathway_results_sorted[pathway_results_sorted$sort_col < 0,]
  }
  
  pathway_results_sorted_filtered = rbind(first_rows, last_rows)
  
  color_list = unlist(colors[["direction"]])[c("enriched", "depleted")]
  
  bubble_plot = ggplot(pathway_results_sorted_filtered, aes(x = tissue, y = reorder(Pathway_short, pvalue_neglog10))) +
    geom_point(aes(size = pvalue_neglog10, color = Regulation_direction)) +
    scale_size_continuous(range=c(0.5,2.5), limits = -log10(c(0.05, 0.001)), labels = c(0.05, 0.01, 0.005, 0.001), breaks = -log10(c(0.05, 0.01, 0.005, 0.001))) +
    scale_color_manual(values = color_list) +
    xlab("") +
    ylab("") +
    theme_classic() +
    theme(legend.position = "bottom", 
          legend.box = "vertical",
          #plot.background = element_rect(fill="transparent"),
          #panel.background = element_rect(fill="transparent"),
          text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          axis.text.y = element_text(family = plots_props$font_family, size = 4), 
          axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header), 
          legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size, face = "bold"),
          plot.margin = margin(0.1,0.1,0,-0.3, "cm"),
          legend.spacing.y = unit(0, "cm"),
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(-0.8,0,0,-4, "cm"),
          legend.spacing.x = unit(0, "cm")
          ) +
    guides(
      color = guide_legend(title = "", nrow = 1, title.position = "left"),  # Color legend in one row
      size = guide_legend(title = "P-value", nrow = 1, title.position = "top")   # Shape legend in one row
    )
  bubble_plot
  ggsave(sprintf("%s/%s_%s__%s_vs_%s_%s=%s.png", output_folder_fig, m, prop, c$g1, c$g2, c$subset[1], paste(tissues, collapse= "_")), bubble_plot, dpi=plots_props$dpi, width= 2/3 * plots_props$image_width, height=plots_props$image_height, units = plots_props$image_units)
  ggsave(sprintf("%s/%s_%s__%s_vs_%s_%s=%s.svg", output_folder_fig, m, prop, c$g1, c$g2, c$subset[1], paste(tissues, collapse= "_")), bubble_plot, width= 2/3 * plots_props$image_width, height=plots_props$image_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])

