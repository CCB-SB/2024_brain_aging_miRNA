suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(cowplot))
#suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(viridisLite))
#suppressPackageStartupMessages(library(gghalves))
#suppressPackageStartupMessages(library(ggridges))
#suppressPackageStartupMessages(library(pbapply))
#suppressPackageStartupMessages(library(extrafont))
#suppressPackageStartupMessages(library(vioplot))
#suppressPackageStartupMessages(library(ggpointdensity))
#suppressPackageStartupMessages(library(ggbeeswarm))
#suppressPackageStartupMessages(library(uwot))
#suppressPackageStartupMessages(library(UpSetR))
#suppressPackageStartupMessages(library(aricode))
## for Vennerable
## conda install -c conda-forge r-devtools
## conda install -c bioconda bioconductor-rbgl
## suppressPackageStartupMessages(library(devtools))
## install_github("js229/Vennerable")
#suppressPackageStartupMessages(library(Vennerable))
#suppressPackageStartupMessages(library(ComplexUpset))
suppressPackageStartupMessages(library(stringr))
#suppressPackageStartupMessages(library(tidyverse))
#suppressPackageStartupMessages(library(gtable))
#suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(svglite))

#suppressPackageStartupMessages(library(rbioapi))
suppressPackageStartupMessages(library(RhpcBLASctl))

suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))

#snakemake = readRDS("snakemake_mapping_stats_bar_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("mapping_stats_bar_plot")


#---------------------------------- Functions ----------------------------------
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
image_size = snakemake@params$image_size

# save rdata
#saveRDS(snakemake, file = sprintf("snakemake_mapping_stats_bar_plot.rds"))

metadata = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 

if (data_input$data_sub_set == "microglia") {
  
  xticks_names = snakemake@params$xticks_names_microglia
  colors = snakemake@params$colors_microglia
  
  rna_composition_complete = fread(snakemake@input$mapping_info_rna_classes_microglia, sep='\t')
  
} else {
  
  xticks_names = snakemake@params$xticks_names
  colors = snakemake@params$colors
  
  rna_composition_complete = fread(snakemake@input$mapping_info_rna_classes, sep='\t')
  
}
# after loading of rna_composition_complete
# read lnc_classes_mapping
lnc_mapping = fread("data/lncRNA_mapping.tsv", sep = "\t", check.names=F)
# get the class names with a 1
lnc_classes = lnc_mapping[lnc_mapping$lncRNA == 1, ]$class
# get the columns from rna_composition_complete which belong to lnc
lnc_cols = colnames(rna_composition_complete)[colnames(rna_composition_complete) %in% lnc_classes]
# sum up to percentages
lnc_perc = rowSums(rna_composition_complete[, ..lnc_cols])
# remove columns that contribute to lncRNA
rna_composition_complete[, lnc_cols] = c()
# add lncRNA column
rna_composition_complete$lncRNA = lnc_perc

output_folder = sprintf("%s/%s_%s/results_%s/figures/mapping_statistics/bar_plots/mapping", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder, recursive=TRUE)


# ----------------------- Script -----------------------
# #tobias
# blas_p = blas_get_num_procs()
# blas_set_num_threads(1)
# omp_p = omp_get_num_procs()
# omp_set_num_threads(1)

#setnames(rna_composition, "miRNA_primary_transcript", "pre-miRNA")

for (group in grouping) {
  print(group)
  group = unlist(unname(group))
  # one samples was two times sequenced due to problems in the lab
  # we removed the sample for further analysis for which less reads could be aligned
  # in the quality control files this sample is still included, therefore we have to remove it again
  #mask = !rna_composition_complete$Sample == "BA10_CA1_cp_26_1"
  #rna_composition = rna_composition_complete[mask,]
  # filter the rna_compostion table by the samples included in metadata table
  mask = (rna_composition_complete$Sample %in% metadata$fastq_name)
  rna_composition = rna_composition_complete[mask,]
  if (length(group) == 1){
    
    split_group = group
    tmp = metadata[match(rna_composition$Sample, metadata$fastq_name)][[split_group]]
    rna_composition[, split_group:=tmp]
    
    # Select numeric columns only
    numeric_cols = sapply(rna_composition, is.numeric)  # Logical vector identifying numeric columns
    numeric_cols["split_group"] = TRUE  # Keep split_group for grouping
    
    # Apply aggregation only to numeric columns
    rna_composition_merged = aggregate(rna_composition[, ..numeric_cols], by = list(rna_composition$split_group), FUN = mean)
    
    # remove old split_group column
    rna_composition_merged$split_group = c() 
    
    # Rename Group.1 column (generated by aggregate)
    colnames(rna_composition_merged)[1] = "split_group"
    
  } else if (length(group) == 2) {
    
    split_group = group[1]
    subsplit_group = group[2]
    
    if (subsplit_group != "sample") {
      tmp = metadata[match(rna_composition$Sample, metadata$fastq_name)][[split_group]]
      rna_composition[, split_group:=tmp]
      tmp = metadata[match(rna_composition$Sample, metadata$fastq_name)][[subsplit_group]]
      rna_composition[, subsplit_group:=tmp]
      
      rna_composition$split_group_merged = sprintf("%s_%s", rna_composition$split_group, rna_composition$subsplit_group)
      
      # Select numeric columns only
      numeric_cols = sapply(rna_composition, is.numeric)  # Logical vector identifying numeric columns
      numeric_cols["split_group_merged"] = TRUE  # Keep split_group_merged for grouping
      
      # Apply aggregation only to numeric columns
      rna_composition_merged = aggregate(rna_composition[, ..numeric_cols], by = list(rna_composition$split_group_merged), FUN = mean)
      
      # remove old split_group_merged column
      rna_composition_merged$split_group_merged = c() 
      rna_composition_merged$subsplit_group = c() 
      
      # Rename Group.1 column (generated by aggregate)
      colnames(rna_composition_merged)[1] = "split_group"
      
      #tmp = str_split_fixed(rna_composition_merged$split_group_merged, "_", 2)
      #rna_composition_merged$split_group = tmp[,1]
      #rna_composition_merged$subsplit_group = tmp[,2]

    } else{
      tmp = metadata[match(rna_composition$Sample, metadata$fastq_name)][[split_group]]
      rna_composition[, split_group:=tmp]
      rna_composition_merged = rna_composition 
    }
  } else {
    print(sprintf("You entered to many split groups. Skip %s", group))
    next
  }
  rna_comp_df = as.data.table(melt(rna_composition_merged))
  ##rna_comp_df[, Group := factor(annot$`Sample Type`[match(Sample, annot$Name)], levels=type_levels)]
  #rna_comp_df[, split_group:=metadata$split_group[match(Sample, metadata$Sample)]]
  
  rna_comp_df[, variable:=factor(variable, levels=rna_comp_df[, mean(value), by="variable"][order(V1, decreasing = T)]$variable)]
  
  # only take samples that are present in annot
  ##rna_comp_df = rna_comp_df[rna_comp_df$Sample %in% metadata$Sample, ]
  
  # sum all types < 0.5 mean expression in category "Other"
  # other_rna_classes = rna_comp_df[,mean(value), by="variable"][V1 < 0.5]$variable
  # rna_comp_df[variable %in% other_rna_classes, variable:="other"]
  # filter for wanted calsses
  observed_classes = names(colors$rna_classes)
  observed_classes = observed_classes[observed_classes != "mRNA"]
  rna_comp_df[!(variable %in% observed_classes), variable:="others"]
  ##rna_comp_df = rna_comp_df[, list(value=sum(value)), by=c("Sample", "variable", "value", "Protocol")]
  
  # to order the colours in the same way like the protocols in the plot
  tmp = str_split_fixed(rna_comp_df$split_group, "_", 2)
  rna_comp_df$plot_names_split_group = unname(xticks_names[[split_group]][tmp[,1]])
  rna_comp_df$plot_names_split_group = factor(rna_comp_df$plot_names_split_group, levels=unname(xticks_names[[split_group]]))
  
  ##cellline_colors_ordered = cellline_colors[unique(rna_comp_df$Protocol)]
  
  mask = (names(unlist(colors$rna_classes)) %in% unique(rna_comp_df$variable))
  colours_rna_class_filtered = unlist(colors$rna_classes)[mask]
  
  colours_split_group = colors[[split_group]]
  names(colours_split_group) = unname(xticks_names[[split_group]][names(colours_split_group)])
  mask = (names(unlist(colours_split_group)) %in% unique(rna_comp_df$plot_names_split_group))
  colours_split_group_filtered = unlist(colours_split_group)[mask]
  
  colours_split_group_text = colors[[sprintf("%s_text", split_group)]]
  names(colours_split_group_text) = unname(xticks_names[[split_group]][names(colours_split_group_text)])
  mask = (names(unlist(colours_split_group_text)) %in% unique(rna_comp_df$plot_names_split_group))
  colours_split_group_text_filtered = unlist(colours_split_group_text)[mask]
  
  #print(rna_comp_df[rna_comp_df$variable == "miRNA",])
  #print(mean(rna_comp_df[rna_comp_df$variable == "miRNA",]$value))
  
  if (length(group) == 1) {
    factors_x_levels = rna_comp_df[variable == "miRNA"][order(value)]$split_group
  } else {
    if (subsplit_group != "sample") {
      factors_x_levels = rna_comp_df[variable == "miRNA"][order(value)]$split_group
    } else {
      rna_comp_df$split_group = sprintf("%s_%s", rna_comp_df$plot_names_split_group, rna_comp_df$Sample)
      factors_x_levels = rna_comp_df[variable == "miRNA"][order(value)]$split_group
    }
  }
  
  if (image_size$image_width >= 9) {
    ncol_legend = 6
    left_margin_legend = -10
    legend_spacing_x = 0.5
    legend_spacing_key_text = -13
    plot_margin_left = 5.5
  } else if (image_size$image_width > 6) {
    ncol_legend = 4
    left_margin_legend = -30
    legend_spacing_x = 0.3
    legend_spacing_key_text = 0
    plot_margin_left = -10
  } else {
    ncol_legend = 4
    left_margin_legend = -10
    legend_spacing_x = 0.5
    legend_spacing_key_text = -8
    plot_margin_left = -8.5
  }

  for (rna_cl in unique(rna_comp_df$variable)) {
    print(sprintf("Average of the %s: %s", rna_cl, mean(rna_comp_df[rna_comp_df$variable == rna_cl,]$value)))
  }

  rna_comp_df$variable = factor(rna_comp_df$variable, levels = names(colours_rna_class_filtered))
  
  compo_grouped = ggplot(rna_comp_df, aes(x=factor(split_group, levels=factors_x_levels),
                                          y=value, fill=variable)) + 
    geom_col() +
    facet_grid(rows=vars(plot_names_split_group), scales = "free_y", space = "free", switch = "y") +
    ##theme_cowplot(11) + 
    scale_y_continuous(labels = scales::percent_format(accuracy = 1, scale=1, suffix=""), expand = expansion(mult=c(0, 0.05))) +
    coord_flip() + 
    xlab("") +
    ylab("Read composition (%)") +
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
          #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4),
          panel.spacing = unit(0, "pt"),
          plot.margin = margin(5.5, 10, 5.5, plot_margin_left, "pt"),
          legend.margin = margin(-15, 0, 0, left_margin_legend, "pt"), legend.key.size = unit(0.3, "cm"), legend.spacing.y = unit(0.4, "cm"), legend.spacing.x = unit(legend_spacing_x, "cm"),
    )
  
  compo_grouped

  compo_grouped_mod = modify_facet_appearance(compo_grouped,
                                              strip.background.y.fill = colours_split_group_filtered,
                                              strip.background.y.col = colours_split_group_filtered, #rep("white", length(cellline_colors)),
                                              strip.text.y.col = colours_split_group_text_filtered #rep("black", length(cellline_colors))
                                              )
  
  width = image_size$image_width
  height = image_size$image_height
  ggsave(sprintf("%s/grouped_by_%s_%sx%s.png", output_folder, paste(group, collapse = "_"), height, width), compo_grouped_mod, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
  ggsave(sprintf("%s/grouped_by_%s_%sx%s.svg", output_folder, paste(group, collapse = "_"), height, width), compo_grouped_mod, width = width, height = height, unit = plots_props$image_units, dpi = plots_props$dpi)

  # compo_grouped = compo_grouped + theme(strip.text.y.left = element_text(family = plots_props$font_family, angle = 0, size=plots_props$font_size))
  # num_protocols = length(unique(rna_comp_df$plot_names_split_group))
  # compo_grouped_mod = modify_facet_appearance(compo_grouped,
  #                                             strip.background.y.fill = rep("grey", num_protocols),
  #                                             strip.background.y.col = rep("grey", num_protocols),
  #                                             strip.text.y.col =rep("black", num_protocols)
  #                                             )
  # 
  # ggsave(sprintf("%s/grouped_by_%s_%sx%s_2.png", output_folder, paste(group, collapse = "_"), height, width), compo_grouped_mod, dpi = plots_props$dpi, width = width, height = height, units = plots_props$image_units)
  # ggsave(sprintf("%s/grouped_by_%s_%sx%s_2.svg", output_folder, paste(group, collapse = "_"), height, width), compo_grouped_mod, width = width, height = height, unit = plots_props$image_units, dpi = plots_props$dpi)
}


#idx = c("male", "female")
#
#p1 = expand_grid(id=idx, id2=unique(rna_comp_df$Group2), x=1:100) %>%
#  mutate(y=rnorm(n=n())) %>%
#  ggplot(rna_comp_df, aes(x=factor(Sample, levels=rna_comp_df[variable == "miRNA"][order(value)]$Sample),
#                                        y=value, fill=variable)) +
#  geom_col() +
#  facet_wrap(~id + id2, nrow = 4, ncol=8) +
#  scale_y_continuous(labels = scales::percent_format(accuracy = 1, scale=1), expand = expansion(mult=c(0, 0.05))) +
#  coord_flip() +
#  scale_fill_manual(values=rna_class_colors, 
#                    labels=function(x) str_trim(gsub("_", " ", gsub("ensembl", "", x))), 
#                    name="Class",
#                    guide=guide_legend(ncol=6, title.position = "top")) +
#  xlab("") +
#  ylab("Read composition") +
#  theme_classic() +
#  theme(legend.position = "bottom", 
#        text = element_text(family = plots_props$font_family, size = plots_props$font_size),
#        plot.title = element_blank(), #element_text(family = plots_props$font_family, size = plots_props$font_size_header),
#        legend.title = element_blank(), #element_text(family = plots_props$font_family, size = plots_props$font_size),
#        legend.text=element_text(family = plots_props$font_family, size = 4),
#        axis.text = element_text(angle = 0, family = plots_props$font_family, size = plots_props$font_size),
#        axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
#        axis.text.y = element_blank(),
#        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
#        strip.text.y.left = element_text(angle = 0, size=9),
#        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4),
#        #panel.spacing = unit(0, "pt"),
#        #legend.margin = margin(-25, 0, 0, 0, "pt") #legend.key.size = unit(0.3, "cm"), legend.spacing.y = unit(0.4, "cm")
#        )
#
#g <- ggplot_gtable(ggplot_build(p1))
#
#stript <- grep("strip", g$layout$name)
#
#grid_cols <- sort(unique(g$layout[stript,]$l))
#t_vals <- rep(sort(unique(g$layout[stript,]$t)), each = length(grid_cols)/2)
#l_vals <- rep(grid_cols[seq_along(grid_cols) %% 2 == 1], length = length(t_vals))
#r_vals <- rep(grid_cols[seq_along(grid_cols) %% 2 == 0], length = length(t_vals))
#labs   <- levels(as.factor(p1$data$id))
#
#for(i in seq_along(labs))
#{
#  filler <- rectGrob(y = 0.7, height = 0.6, gp = gpar(fill = "gray80", col = NA))
#  tg    <- textGrob(label = labs[i], y = 0.75, gp = gpar(cex = 0.8))
#  g     <- gtable_add_grob(g, filler, t = t_vals[i], l = l_vals[i], r = r_vals[i], 
#                           name = paste0("filler", i))
#  g     <- gtable_add_grob(g, tg, t = t_vals[i], l = l_vals[i], r = r_vals[i], 
#                           name = paste0("textlab", i))
#}
#
#grid.newpage()
#grid.draw(g)
#
#ggsave("../results/report/compo_grouped_complex.png", g, dpi=600, width=9, height=6, units = "cm")
#ggsave("../results/report/compo_grouped_complex.svg", g, width = 9 height = 6, unit = "cm", dpi = 600)
#
# mir_expr = fread("results/quantification/overview/std_quantification_rpmmm_norm.csv")
# mir_expr[, miRNA:=paste(precursor, miRNA, sep='|')]
# mir_expr[, precursor:=NULL]
# 
# detected_mirs = data.frame(num_mirnas=colSums(mir_expr[, 2:ncol(mir_expr)] > 0))
# detected_mirs$Group = metadata$Group[match(rownames(detected_mirs), metadata$Sample)]
# detected_mirs = setDT(detected_mirs)
# 
# detected_mirs_p = ggplot(detected_mirs, aes(x=factor(Group, levels=detected_mirs[, mean(num_mirnas), by="Group"][order(V1)]$Group), y=num_mirnas, color=Group, fill=Group)) + 
#   #geom_quasirandom(alpha = 1, width = 0.2) +
#   geom_half_point() +
#   geom_half_boxplot(outlier.colour = NA) +
#   #theme_cowplot(11) +
#   scale_y_continuous(breaks=scales::pretty_breaks(10)) +
#   scale_fill_manual(values=cellline_colors) + 
#   scale_color_manual(values=cellline_colors) +
#   coord_flip() +
#   xlab("") + 
#   ylab("Detected miRNAs") + 
#   theme_classic() +
#   theme(legend.position = "none", 
#         text = element_text(family = "Arial", size = 9),
#         plot.title = element_blank(), #element_text(family = "Arial", size = 12),
#         legend.title = element_blank(), #element_text(family = "Arial", size = 9),
#         legend.text=element_text(family = "Arial", size = 4),
#         axis.text = element_text(angle = 0, family = "Arial", size = 9),
#         axis.title = element_text(family = "Arial", size = 9),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(), 
#         axis.line.y = element_blank(),
#         strip.text.y.left = element_text(family = "Arial", angle = 0, size=9),
#         #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4),
#         panel.spacing = unit(0, "pt"),
#         legend.margin = margin(-15, 0, 0, 0, "pt"), legend.key.size = unit(0.3, "cm"), legend.spacing.y = unit(0.4, "cm")
#         )
# 
# ggsave(sprintf("results/report/%s/detected_mirs_th=0rpmmm%s_new_tissue_order.png", folder_name, file_name), detected_mirs_p, dpi=600, width=9, height=6, units = "cm")
# ggsave(sprintf("results/report/%s/detected_mirs_th=0rpmmm%s_new_tissue_order.svg", folder_name, file_name), detected_mirs_p, width = 9, height = 6, unit = "cm", dpi = 600)


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
