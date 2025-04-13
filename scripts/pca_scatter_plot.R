suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(cowplot))


#snakemake = readRDS("snakemake_pca_scatter_plot.rds")

source("./scripts/iwanthue.R")

set.seed(snakemake@params$parameters_porps$set_seed)

print("pca_scatter_plot") 


#---------------------------------- Functions ----------------------------------
plot_scatter = function(plot_df, annot, imp, props, prop_mapping, colour_para, major_prop, colors, plots_props, legend_title, output_path, output_filename) {
  # set font-family for all plots
  #theme_set(theme_cowplot(font_family = plots_props$font_family))
  
  if(is.vector(props) && length(props) == 2){
    p = ggplot(plot_df, aes_string(x = "D1", y = "D2", color = prop_mapping[major_prop], shape = prop_mapping[props[2]], text = "ID"))
  } else {
    if (props == "reads_aligned") {
      p = ggplot(plot_df, aes_string(x = "D1", y = "D2", color = colour_para, text = "ID"))
    } else {
      p = ggplot(plot_df, aes_string(x = "D1", y = "D2", color = prop_mapping[major_prop], text = "ID"))
    }
  }
  
  #p = p + geom_point(size = 1)
  # for phd thesis change point size
  p = p + geom_point(size = 1.5)
  
  # if you want the pca plot cutted to a specific area
  #p = p + xlim(-15, 18) + ylim(-10, 15)
  
  if (length(colors) == 1 && colors == "viridis_c") {
    p = p + scale_colour_viridis_c(name=major_prop)
  } else if (length(colors) == 1 && colors == "viridis_d") {
    p = p + scale_colour_viridis_d(name=major_prop)
  } else if (length(props) == 1 && colour_para == "reads_aligned") {
    p = p + scale_colour_viridis_c(name="reads_aligned")
  } else if (length(props) == 1 && colour_para == "reads_aligned_log10") {
    p = p + scale_colour_viridis_c(name="reads_aligned_log10")
  } else {
    
    p = p + scale_color_manual(values = colors)
  }
  
  if(is.vector(props) && length(props) == 2){
    p = p + #scale_shape(name=props[2]) 
    scale_shape_manual(values=c(16, 17, 15, 3, 8, 7, 10))
  }
  
  if(nrow(plot_df) < 25){
    p = p + geom_text_repel(aes(label = plot_df$ID), show.legend = FALSE)
  }
  
  if (major_prop != FALSE) {
    if(!is.numeric(annot[[major_prop]])) {
      p = p + theme_classic() + theme(legend.position="bottom")
    }
  } else {
    p = p + theme_classic() + theme(legend.position="bottom", legend.key.width = unit(plots_props$image_width/10, plots_props$image_units)) 
  }
  
  p = p + xlab(sprintf("PC1 (%.2f%%)", imp["Proportion of Variance","PC1"]*100)) +
    ylab(sprintf("PC2 (%.2f%%)", imp["Proportion of Variance", "PC2"]*100))
  
  p = p + theme(text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
                axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
                axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size), 
                plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header), 
                legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
                legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size, face = "bold"),
                legend.margin=margin(t = -8), legend.key.size = unit(0.3, "cm"), legend.spacing.y = unit(0.4, "cm"))
  
  needed_rows = 1
  if(is.vector(props) && length(props) == 2){
    total_len = sum(nchar(unique(as.character(plot_df[[prop_mapping[major_prop]]])))) + sum(nchar(unique(as.character(plot_df[[prop_mapping[props[2]]]]))))
    if(nlevels(plot_df[[prop_mapping[major_prop]]]) + nlevels(plot_df[[prop_mapping[props[2]]]]) > 4 || total_len > 50){
      needed_rows = ceiling(total_len / 50)
      #p = p + theme(legend.box = "horizontal")
      p = p + guides(col=FALSE, fill=guide_legend(nrow=needed_rows, title=legend_title), shape=guide_legend(nrow=needed_rows, title=legend_title))
    }
  } else if (props != "reads_aligned") {
    if(length(colors) != 1 || colors == "viridis_d") {
      total_len = sum(nchar(unique(as.character(plot_df[[prop_mapping[major_prop]]]))), na.rm=T)
      if(nlevels(plot_df[[prop_mapping[major_prop]]]) >= 4 || total_len > 50){
        needed_rows = ceiling(total_len / 30)
        print(sprintf("Needed rows: %s", needed_rows))
        p = p + guides(col=guide_legend(nrow=needed_rows, title=legend_title))
      }
    } 
    
    # for phd thesis remove legend
    p = p + guides(col = "none")
    
  } else {
    p = p + guides(col = guide_colourbar(title = legend_title)) 
  }
  
  hight = plots_props$image_height
  width = plots_props$image_width
  ggsave(file = sprintf("%s_%sx%s.svg", output_path, hight, width), plot = p, width = width, height = hight, unit = plots_props$image_units, dpi = plots_props$dpi)
  ggsave(file = sprintf("%s_%sx%s.png", output_path, hight, width), plot = p, width = width, height = hight, unit = plots_props$image_units, dpi = plots_props$dpi)
  
  hight = 9
  width = 9
  ggsave(file = sprintf("%s_%sx%s.svg", output_path, hight, width), plot = p, width = width, height = hight, unit = plots_props$image_units, dpi = plots_props$dpi)
  ggsave(file = sprintf("%s_%sx%s.png", output_path, hight, width), plot = p, width = width, height = hight, unit = plots_props$image_units, dpi = plots_props$dpi)

  print(sprintf("Done: %s", output_filename))
}

#---------------------------------- Load data ---------------------------------- 
data_input = snakemake@params$data_input
properties = snakemake@params$properties
params = snakemake@params$parameters
plots_props = snakemake@params$plots_props
snakemake_colors = snakemake@params$colors
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_pca_scatter_plot.rds")

expr = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 

mapping_info = fread(snakemake@input$mapping_info, sep='\t')
# for some reason, the colum Sample if named fastq_names sometimes
if (!("Sample" %in% colnames(mapping_info)) && ("fastq_names" %in% colnames(mapping_info))) {
  mapping_info$Sample = mapping_info$fastq_names
}
tmp = str_split_fixed(mapping_info$Sample, "_", 5)
mapping_info$ID = paste(tmp[,2], tmp[,3], tmp[,4], sep = "_")

output_folder_fig = sprintf("%s/%s_%s/results_%s/figures/pca_plots/scatter_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_fig, recursive=TRUE)


# --------------------------------- Script -------------------------------------
for (props in properties) {
  # remove constants
  expr_filtered = expr[apply(expr[,-1], 1, var) != 0,]
  #print(sprintf("Dim original: %s, Dim filtered: %s", dim(expr)[1], dim(expr_filtered)[1]))

  comp = prcomp(t(expr_filtered[,-1]), center = T, scale=T)
  imp = summary(comp)$importance
  comp_df = as.data.frame(comp$x)
  plot_df = data.frame(D1=comp_df$PC1, D2=comp_df$PC2, ID=colnames(expr_filtered[,-1]), stringsAsFactors = F)
  
  prop_mapping = make.names(props)
  names(prop_mapping) = props
  
  if (length(props) > 1){
    for(a in props){
      plot_df[[prop_mapping[a]]] = as.factor(unlist(lapply(plot_df$ID, function(n) annot[annot[[data_input$identifier_column]] == n, a, with=F][1])))
    }
    major_prop = props[1]
    
    # create legend title
    if (props[2] == "Time") {
      legend_title = "Timepoint"
    } else {
      legend_title = props[2]
    }
    
    # create output path
    output_filename = sprintf("%s_%s", props[1], props[2])
    
  } else {
    if (props == "reads_aligned") {
      # unlist here
      plot_df[[prop_mapping[props]]] = unlist(lapply(plot_df$ID, function(n) mapping_info[(mapping_info$ID == n), props, with=F][1]))
      #print(mapping_info)
      plot_df$reads_aligned_log10 = log10(plot_df$reads_aligned + 1)
      major_prop = FALSE
      plot_df[[prop_mapping[props]]] = as.numeric(plot_df[[prop_mapping[props]]])
    } else {
      plot_df[[prop_mapping[props]]] = as.factor(unlist(lapply(plot_df$ID, function(n) annot[annot[[data_input$identifier_column]] == n, props, with=F][1])))
      major_prop = props
    }
    
    legend_title = props
    
    # create output path
    output_filename = sprintf("%s", props)
  }
  
  # set output file name and create folder
  output_path = sprintf("%s/%s", output_folder_fig, output_filename)
  
  if(!is.null(snakemake_colors)){
    if (major_prop != FALSE) {
      colors = unlist(snakemake_colors[[major_prop]])
      if(length(colors) == 1 && colors == "auto") {
        colors = unname(iwanthue(nlevels(plot_df[[prop_mapping[major_prop]]])))
      } else if(length(colors) == 1 && colors == "viridis_c"){
        # make continuous instead of discrete
        plot_df[[prop_mapping[major_prop]]] = unlist(lapply(plot_df$ID, function(n) annot[annot[[data_input$identifier_column]] == n, props, with=F][1]))
      } else if(length(colors) > 1 && length(colors) < length(unique(prop_mapping[major_prop]))) {
        stop(sprintf("Not enough colors provided for the property %s! Expected %d, got %d (%s)", major_prop, length(unique(prop_mapping[major_prop])), length(color), paste(names(colors), sep=', ')))
      } 
    } 
  } else {
    colors = unname(iwanthue(nlevels(plot_df[[prop_mapping[major_prop]]])))
  }
  
  if (length(props) == 1  && props == "reads_aligned") {
    plot_scatter(plot_df, annot, imp, props, prop_mapping, "reads_aligned", major_prop, colors, plots_props, legend_title, output_path, output_filename)
    plot_scatter(plot_df, annot, imp, props, prop_mapping, "reads_aligned_log10", major_prop, colors, plots_props, sprintf("%s_log10", legend_title), sprintf("%s_log10", output_path), sprintf("%s_log10", output_filename))
  } else {
    plot_scatter(plot_df, annot, imp, props, prop_mapping, FALSE, major_prop, colors, plots_props, legend_title, output_path, output_filename)
  }
  
  #max_list = c()
  #for (feature in colnames(comp$rotation)) {
  #  max_list[[feature]] = sort(comp$rotation[,feature], decreasing = TRUE)[1:5] 
  #  #argmax_list[[feature]] = names(sort(comp$rotation[,feature], decreasing = TRUE)[1:5]) 
  #}
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])



