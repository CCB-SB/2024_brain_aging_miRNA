suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(viridisLite))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(openxlsx))

suppressPackageStartupMessages(library(RhpcBLASctl))

suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))

#snakemake = readRDS("snakemake_expressed_features_per_rna_class_bar_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("expressed_features_per_rna_class_bar_plot")


#---------------------------------- Functions ----------------------------------
feature_filtering_for_grouping = function(raw_count_ti, min_detect_raw, detection_rate){
  raw_count_ti_bin = ifelse(raw_count_ti >= min_detect_raw, 1, 0)
  if (length(unique(rowMeans(raw_count_ti_bin) >= detection_rate/100)) == 2) {
    raw_count_ti_filtered = raw_count_ti[rowMeans(raw_count_ti_bin) >= detection_rate/100,]
  } else if ((length(unique(rowMeans(raw_count_ti_bin) >= detection_rate/100)) == 1) & (unique(rowMeans(raw_count_ti_bin) >= detection_rate/100) == TRUE)) {
    raw_count_ti_filtered = raw_count_ti[rowMeans(raw_count_ti_bin) >= detection_rate/100,]
  }
  return(raw_count_ti_filtered)
}


averaged_number_of_reads = function(raw_count_ti_filtered){
  if (nrow(raw_count_ti_filtered) == 0) {
    return(0)
  } else {
    raw_count_ti_filtered_row_sums = rowSums(raw_count_ti_filtered)
    raw_count_ti_filtered_row_sums_mean = mean(raw_count_ti_filtered_row_sums)
    #samples = colnames(raw_count_ti_filtered)
    #sample_names = str_split_fixed(samples, "_", 2)[,2]
    #sample_mapping = data.frame(sample_names = sample_names, samples=samples)
    #list_miRNAs = c()
    #for (sample in samples) {
    #  tmp = expr[,..sample]
    #  tmp_1 = tmp[,..sample] > 0
    #  miRNA_names = expr$miRNA
    #  # tmp_2 = sample_mapping[sample_mapping$samples == sample,]$sample_names
    #  
    #  
    #  list_miRNAs[[sample]] = miRNA_names[tmp_1]
    #}
    
    return(raw_count_ti_filtered_row_sums_mean)
  }
}

modify_facet_appearance <- function(plot = NULL, strip.background.x.fill = NULL, strip.background.y.fill = NULL, strip.background.x.col = NULL,
                                    strip.background.y.col = NULL, strip.text.x.col = NULL, strip.text.y.col = NULL){
  
  if(is.null(plot)){stop("A ggplot (gg class) needs to be provided!")}
  
  # Generate gtable object to modify the facet strips:
  g <- ggplot_gtable(ggplot_build(plot))
  
  # Get the locations of the right and top facets in g:
  stripy <- which(grepl('strip-r|strip-l', g$layout$name)) # account for when strip positions are switched r-l and/or t-b in facet_grid(switch = )
  stripx <- which(grepl('strip-t|strip-b', g$layout$name))
  
  # Check that the provided value arrays have the same length as strips the plot has:
  lx <- c(length(strip.background.x.fill), length(strip.background.x.col), length(strip.text.x.col))
  if(!all(lx==length(stripx) | lx==0)){stop("x The provided vectors with values need to have the same length and the number of facets in the plot!")}
  ly <- c(length(strip.background.y.fill), length(strip.background.y.col), length(strip.text.y.col))
  if(!all(ly==length(stripy) | ly==0)){stop("y The provided vectors with values need to have the same length and the number of facets in the plot!")}
  
  # Change the strips on the y axis:
  for (i in seq_along(stripy)){ # if no strips in the right, the loop will not be executed as seq_along(stripy) will be integer(0)
    
    # Change strip fill and (border) colour :
    j1 <- which(grepl('strip.background.y', g$grobs[[stripy[i]]]$grobs[[1]]$childrenOrder))
    if(!is.null(strip.background.y.fill[i])){g$grobs[[stripy[i]]]$grobs[[1]]$children[[j1]]$gp$fill <- strip.background.y.fill[i]} # fill
    if(!is.null(strip.background.y.col[i])){g$grobs[[stripy[i]]]$grobs[[1]]$children[[j1]]$gp$col <- strip.background.y.col[i]} # border colour
    
    # Change color of text:
    j2 <- which(grepl('strip.text.y', g$grobs[[stripy[i]]]$grobs[[1]]$childrenOrder))
    if(!is.null(strip.text.y.col[i])){g$grobs[[stripy[i]]]$grobs[[1]]$children[[j2]]$children[[1]]$gp$col <- strip.text.y.col[i]}
    
  }
  
  # Same but for the x axis:
  for (i in seq_along(stripx)){
    
    # Change strip fill and (border) colour :
    j1 <- which(grepl('strip.background.x', g$grobs[[stripx[i]]]$grobs[[1]]$childrenOrder))
    if(!is.null(strip.background.x.fill[i])){g$grobs[[stripx[i]]]$grobs[[1]]$children[[j1]]$gp$fill <- strip.background.x.fill[i]} # fill
    if(!is.null(strip.background.x.col[i])){g$grobs[[stripx[i]]]$grobs[[1]]$children[[j1]]$gp$col <- strip.background.x.col[i]} # border colour
    
    # Change color of text:
    j2 <- which(grepl('strip.text.x', g$grobs[[stripx[i]]]$grobs[[1]]$childrenOrder))
    if(!is.null(strip.text.x.col[i])){g$grobs[[stripx[i]]]$grobs[[1]]$children[[j2]]$children[[1]]$gp$col <- strip.text.x.col[i]}
    
  }
  
  return(g)
}
# Note that it returns a gtable object. This can be ploted with plot() or grid::draw.grid().
# patchwork can handle the addition of such gtable to a layout with other ggplot objects,
# but be sure to use patchwork::wrap_ggplot_grob(g) for proper alignment of plots!
# See: https://patchwork.data-imaginist.com/reference/wrap_ggplot_grob.html


# ---------------------- Load data ----------------------
data_input = snakemake@params$data_input
plots_props = snakemake@params$plots_props
results_folder = snakemake@params$results_folder
grouping = snakemake@params$grouping
min_detect_raw = snakemake@params$min_detect_raw
detection_rate = snakemake@params$detection_rate
image_size = snakemake@params$image_size
xticks_names = snakemake@params$xticks_names
colors = snakemake@params$colors

# save rdata
#saveRDS(snakemake, file = sprintf("snakemake_expressed_features_per_rna_class_bar_plot.rds"))
#stop()

annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
# rna_composition_complete = fread(snakemake@input$mapping_info_rna_classes, sep='\t')

output_folder_fig = sprintf("%s/%s_%s/results_%s/figures/expressed_per_brain_region_per_rna_class/bar_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_fig, recursive=TRUE)
output_folder_tab = sprintf("%s/%s_%s/results_%s/matrices/expressed_per_brain_region_per_rna_class/bar_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_tab, recursive=TRUE)

# ----------------------- Script -----------------------
number_of_reads = c()
for (rna_class in data_input$rna_classes) {
  print(rna_class)
  # load data
  parts = strsplit(data_input$data_sub_set, "_")[[1]]
  suffix = paste(parts[-1], collapse = "_")
  #raw_count = fread(sprintf("%s_%s/%s_%s_quantification_raw.csv", data_input$raw_data_folder, data_input$data_sub_set, rna_class, data_input$feature_filtering), sep='\t', header=T)
  raw_count = fread(sprintf("%s_%s/%s_%s_expression_raw_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, rna_class, data_input$feature_filtering ,suffix), sep='\t', header=T)
  
  number_of_reads_ti = c()
  for (ti in unique(annot[[grouping]])) {
    #print(ti)
    mask = annot[annot[[grouping]] == ti,][[data_input$identifier_column]]
    raw_count_ti = raw_count[, ..mask]
    raw_count_ti_filtered = feature_filtering_for_grouping(raw_count_ti, min_detect_raw, detection_rate)
    number_of_reads_ti[[ti]] = averaged_number_of_reads(raw_count_ti_filtered)
  }
  
  number_of_reads[[rna_class]] = number_of_reads_ti
}
  
rna_comp_df = as.data.table(melt(number_of_reads))
colnames(rna_comp_df) = c("value", "grouping", "variable")

rna_comp_df$value_per = rna_comp_df$value
for (ti in unique(annot[[grouping]])){
  total = sum(rna_comp_df[rna_comp_df$grouping == ti,]$value)
  rna_comp_df[rna_comp_df$grouping == ti,]$value_per = (rna_comp_df[rna_comp_df$grouping == ti,]$value * 100) / total
}

rna_comp_df[, variable:=factor(variable, levels=rna_comp_df[, mean(value_per), by="variable"][order(V1, decreasing = T)]$variable)]

fwrite(rna_comp_df, sprintf("%s/grouped_by_%s.csv", output_folder_tab, grouping), sep = "\t", row.names = TRUE)
write.xlsx(rna_comp_df, sprintf("%s/grouped_by_%s.xlsx", output_folder_tab, grouping), colNames = TRUE, rowNames = TRUE, append = FALSE)

# sum all types < 0.5 mean expression in category "Other"
# other_rna_classes = rna_comp_df[,mean(value_per), by="variable"][V1 < 0.5]$variable
# rna_comp_df[variable %in% other_rna_classes, variable:="other"]

# filter for wanted calsses
#observed_classes = names(colors$rna_classes)
#observed_classes = observed_classes[observed_classes != "mRNA"]
#rna_comp_df[!(variable %in% observed_classes), variable:="others"]
##rna_comp_df = rna_comp_df[, list(value_per=sum(value_per)), by=c("Sample", "variable", "value_per", "Protocol")]

# to order the colours in the same way like the protocols in the plot
#tmp = str_split_fixed(rna_comp_df$grouping, "_", 2)
#rna_comp_df$plot_names_grouping = unname(xticks_names[[grouping]][tmp[,1]])
rna_comp_df$plot_names_grouping = unname(xticks_names[[grouping]][rna_comp_df$grouping])
rna_comp_df$plot_names_grouping = factor(rna_comp_df$plot_names_grouping, levels=unname(xticks_names[[grouping]]))

##cellline_colors_ordered = cellline_colors[unique(rna_comp_df$Protocol)]

mask = (names(unlist(colors$rna_classes)) %in% unique(rna_comp_df$variable))
colours_rna_class_filtered = unlist(colors$rna_classes)[mask]

colours_grouping = colors[[grouping]]
names(colours_grouping) = unname(xticks_names[[grouping]][names(colours_grouping)])
mask = (names(unlist(colours_grouping)) %in% unique(rna_comp_df$plot_names_grouping))
colours_grouping_filtered = unlist(colours_grouping)[mask]

colours_grouping_text = colors[[sprintf("%s_text", grouping)]]
names(colours_grouping_text) = unname(xticks_names[[grouping]][names(colours_grouping_text)])
mask = (names(unlist(colours_grouping_text)) %in% unique(rna_comp_df$plot_names_grouping))
colours_grouping_text_filtered = unlist(colours_grouping_text)[mask]

#print(rna_comp_df[rna_comp_df$variable == "miRNA",])
#print(mean(rna_comp_df[rna_comp_df$variable == "miRNA",]$value_per))

#if (length(grouping) == 1) {
factors_x_levels = rna_comp_df[variable == "miRNA"][order(-value_per)]$grouping
rna_comp_df$plot_names_grouping_fac = factor(rna_comp_df$plot_names_grouping, levels = unlist(unname(xticks_names[[grouping]][factors_x_levels])))

colours_grouping_filtered_sorted = colours_grouping_filtered[unlist(unname(xticks_names[[grouping]][factors_x_levels]))]
colours_grouping_text_filtered_sorted = colours_grouping_text_filtered[unlist(unname(xticks_names[[grouping]][factors_x_levels]))]

#} else {
#  if (subgrouping != "sample") {
#    factors_x_levels = rna_comp_df[variable == "miRNA"][order(value_per)]$grouping
#  } else {
#    rna_comp_df$grouping = sprintf("%s_%s", rna_comp_df$plot_names_grouping, rna_comp_df$Sample)
#    factors_x_levels = rna_comp_df[variable == "miRNA"][order(value_per)]$grouping
#  }
#}

if (image_size$image_width >= 9) {
  ncol_legend = 6
  left_margin_legend = -10
  legend_spacing_x = 0.5
  legend_spacing_key_text = -13
  plot_margin_left = 5.5
  x_axis_title_left = 0.5
} else if (image_size$image_width >= 6) {
  ncol_legend = 4
  left_margin_legend = -47
  legend_spacing_x = 0.3
  legend_spacing_key_text = -8
  plot_margin_left = -10
  x_axis_title_left = 0.6
} else {
  ncol_legend = 4
  left_margin_legend = -47
  legend_spacing_x = 0.3
  legend_spacing_key_text = -8
  plot_margin_left = -8
  x_axis_title_left = 0.5
}

# Rows sorted by colour mapping
rna_comp_df$variable_fac = factor(rna_comp_df$variable, levels = names(colours_rna_class_filtered))

compo_grouped = ggplot(rna_comp_df, aes(x=grouping, #factor(grouping, levels=factors_x_levels)
                                        y=value_per, fill=variable_fac)) + 
  geom_col() +
  facet_grid(rows=vars(plot_names_grouping_fac), scales = "free_y", space = "free", switch = "y") +
  ##theme_cowplot(11) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1, scale=1, suffix=""), expand = expansion(mult=c(0, 0.05))) +
  coord_flip() + 
  xlab("") +
  ylab("Expressed count composition (%)") +
  scale_fill_manual(values=colours_rna_class_filtered,
                    labels=function(x) str_trim(gsub("_", " ", gsub("ensembl", "", x))), 
                    name="Class",
                    guide=guide_legend(ncol=ncol_legend, title.position = "top")
  ) +
  theme_classic() +
  theme(legend.position = "bottom", 
        text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header),
        legend.title = element_blank(), #element_text(family = plots_props$font_family, size = plots_props$font_size),
        legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size, margin = margin(0, 0, 0, legend_spacing_key_text, "pt")),
        axis.text = element_text(angle = 0, family = plots_props$font_family, size = plots_props$font_size),
        axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(),
        strip.text.y.left = element_text(family = plots_props$font_family, angle = 0, size=plots_props$font_size),
        axis.title.x = element_text(family = plots_props$font_family, size = plots_props$font_size, hjust = x_axis_title_left),
        panel.spacing = unit(0, "pt"),
        plot.margin = margin(5.5, 10, 5.5, plot_margin_left, "pt"),
        legend.margin = margin(-15, 0, 0, left_margin_legend, "pt"), legend.key.size = unit(0.3, "cm"), legend.spacing.y = unit(0.4, "cm"), legend.spacing.x = unit(legend_spacing_x, "cm"),
  )

#compo_grouped

compo_grouped_mod = modify_facet_appearance(compo_grouped,
                                            strip.background.y.fill = colours_grouping_filtered_sorted,
                                            strip.background.y.col = colours_grouping_filtered_sorted, #rep("white", length(cellline_colors)),
                                            strip.text.y.col = colours_grouping_text_filtered_sorted #rep("black", length(cellline_colors))
                                            )

width = image_size$image_width
height = image_size$image_height
ggsave(sprintf("%s/grouped_by_%s_%sx%s.png", output_folder_fig, grouping, height, width), compo_grouped_mod, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
ggsave(sprintf("%s/grouped_by_%s_%sx%s.svg", output_folder_fig, grouping, height, width), compo_grouped_mod, width = width, height = height, unit = plots_props$image_units, dpi = plots_props$dpi)


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
