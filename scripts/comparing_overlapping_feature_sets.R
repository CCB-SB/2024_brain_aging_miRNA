suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(viridisLite))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ggnetwork))
suppressPackageStartupMessages(library(ggVennDiagram))
suppressPackageStartupMessages(library(eulerr))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))


#snakemake = readRDS("snakemake_comparing_overlapping_feature_sets.rds")

data.table::setDTthreads(1)
set.seed(snakemake@params$parameters_props$set_seed)

print("hclust_comparing_overlapping_feature_sets")


#---------------------------------- Functions ----------------------------------
# write the tissue name on the side of a venn
#tissue_name_side = function(plot, tissue_group, tissue, xticks_names) {
#  plot_with_label = plot + annotation_custom(
#    grob = textGrob(sprintf("%s", xticks_names[[tissue_group]][tissue]), rot = 270, vjust = 0.5, hjust = 0.5),
#    xmin = 1000, xmax = 1000, ymin = -Inf, ymax = Inf
#  )
#  return(plot_with_label)
#}

create_venn_plot = function(venn_input, colour_set, colour_list_text, data_input, xticks_names, plot_props, output_folder_venn, file_name) {
  # manually calculate the stuff ggVennDiagramm should do automatically
  # this way we cat color the edges manually
  venn = Venn(venn_input)
  features_cluster_manual_venn = process_data(venn)

  items = features_cluster_manual_venn@region$item
  item_combined = c()
  for (item in items) {
    item_combined = append(item_combined, do.call(paste, c(list(item), collapse = "\n")))
  }
  # remove mmu-miR
  item_combined = gsub("mmu-miR-", "", item_combined)
  features_cluster_manual_venn@region$item_combined = item_combined
  
  # we need to specify all colours for the overlaps
  colour_list_with_overlap = unlist(colour_set)[features_cluster_manual_venn@region$name]
  #for (origin in features_cluster_manual_venn@region$name) {
  #  if (origin %in% names(colour_list)) {
  #    colour_list_with_overlap[[origin]] = colour_list[[origin]]
  #  } else {
  #    colour_list_with_overlap[[origin]] = "red"
  #  }
  #}

  # plot with label list
  p_labels = ggVennDiagram(venn_input, label_alpha=0,
                           set_size = plot_props$font_size_header * 25.4/72, #set_color = colour_list,
                           label="none", label_color = "black", label_size = plot_props$font_size * 25.4/72, 
                           edge_size = 0.5,
  ) +
  # text with boxes
  # geom_sf_label(aes(label = item_combined), data = venn_region(features_cluster_manual_venn), alpha=0.5) + 
  # no boxes
  geom_sf_text(aes(label = item_combined, color = name), data = venn_region(features_cluster_manual_venn), size=6 * 25.4/72) + 
  geom_sf(aes(color = name), data = venn_setedge(features_cluster_manual_venn), show.legend = F) +
  # expand the x axis
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  # geom_sf(aes(fill = name), data = venn_region(features_cluster_manual_venn), show.legend = F) +
  # scale_fill_manual(values =  alpha(colour_list, .2)) +
  scale_color_manual(values =  alpha(colour_list_with_overlap, .2)) +
  #scale_fill_manual(values = c("grey","grey","grey","grey")) +
  # scale_fill_gradient(low = "#D5FBFF", high = "#0092BF")
  #scale_fill_gradient(low = "white", high = "white", guide = "none") +
  #scale_fill_gradient(low = "#FFFFC6", high = "#4B86B9", name = sprintf("Number of %ss\n(top %s regarding %s, th=%szscore)", data_input$rna_class, i, m, zscore_th), guide="colourbar") +
  #scale_fill_distiller(palette = "Greens", name = sprintf("Number of %ss\n(top %s regarding %s, th=%szscore)", data_input$rna_class, i, m, zscore_th), guide="colourbar") +
  #scale_fill_viridis_c() + 
  #ggtitle(sprintf("Number of %ss\n(top %s regarding %s, th=%szscore)", data_input$rna_class, i, m, zscore_th)) +
  scale_y_continuous(expand = expansion(mult = .1)) +
  theme(text = element_text(size = plot_props$font_size, family = plot_props$font_family),
        #plot.title = element_text(family = plot_props$font_family, size = plot_props$font_size),
        legend.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
        legend.title = element_text(family = plot_props$font_family, size = plot_props$font_size, face = "bold"),
        legend.position = "none",# legend.box = "horizontal",
        plot.margin = unit(c(-1,-1,-1,-1), "cm"),
        #panel.spacing = margin(t = 0, r = 0, b = -3, l = 0, unit = "cm"),
  ) #+ 
  #guides(fill = guide_colourbar(title.position = "top", title.hjust = 0), size = guide_legend(title.position = "top", title.hjust = 0))
  
  # change internal fill aes to name
  p_labels$layers[[1]]$mapping = aes(fill = name)
  # set the area color (fill)
  p_labels = p_labels + 
                      scale_fill_manual(values = colour_list_with_overlap) +
                      theme(legend.position = "none")
  p_labels = p_labels + 
                      scale_color_manual(values = c("black","black","black","black", "black","black","black", "black","black","black")) +
                      theme(legend.position = "none")
  #p_labels
  plot_height = plot_props$image_height
  plot_width = 2/3 * plot_props$image_width
  ggsave(sprintf("%s/%s_label.png", output_folder_venn, file_name), p_labels, dpi=plot_props$dpi, width=plot_width, height=plot_height, units = plot_props$image_units)
  ggsave(sprintf("%s/%s_label.svg", output_folder_venn, file_name), p_labels, width=plot_width, height=plot_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
  
  
  # plot with counts
  p = ggVennDiagram(venn_input, label_alpha=1,
                    set_size = plot_props$font_size_header * 25.4/72, #set_color = colour_list_filtered,
                    label = "count", label_color = "black", label_size = 6 * 25.4/72, label_geom = "text",
                    edge_size = 0.5,
  ) +
  # text with boxes
  # geom_sf_label(aes(label = item_combined), data = venn_region(features_cluster_manual_venn), alpha=0.5) + 
  # no boxes
  #geom_sf_text(aes(label = item_combined), data = venn_region(features_cluster_manual_venn), size=6 * 25.4/72) + 
  geom_sf(aes(color = name), data = venn_setedge(features_cluster_manual_venn), show.legend = F) +
  # expand the x axis
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  # geom_sf(aes(fill = name), data = venn_region(feature_cluster_manual_venn), show.legend = F) +
  scale_color_manual(values =  alpha(colour_list_with_overlap, .2)) +
  # scale_fill_manual(values = c("grey","grey","grey","grey")) +
  # scale_fill_gradient(low = "#D5FBFF", high = "#0092BF")
  #scale_fill_gradient(low = "white", high = "white", guide = "none") +
  #scale_fill_gradient(low = "#FFFFC6", high = "#4B86B9", name = sprintf("Number of %ss\n(top %s regarding %s, th=%szscore)", data_input$rna_class, i, m, zscore_th), guide="colourbar") +
  #scale_fill_distiller(palette = "Greens", name = sprintf("Number of %ss\n(top %s regarding %s, th=%szscore)", data_input$rna_class, i, m, zscore_th), guide="colourbar") +
  #scale_fill_viridis_c() + 
  #ggtitle(sprintf("Number of %ss\n(top %s regarding %s, th=%szscore)", data_input$rna_class, i, m, zscore_th)) +
  #ggtitle(plot_title) +
  scale_y_continuous(expand = expansion(mult = .1)) +
  theme(text = element_text(size = plot_props$font_size, family = plot_props$font_family),
        #plot.title = element_text(family = plot_props$font_family, size = plot_props$font_size),
        legend.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
        legend.title = element_text(family = plot_props$font_family, size = plot_props$font_size, face = "bold"),
        legend.position = "none",# legend.box = "horizontal",
        plot.margin = unit(c(-1,-1,-1,-1), "cm"),
        #panel.spacing = margin(t = 0, r = 0, b = -3, l = 0, unit = "cm"),
  ) #+ 
  #guides(fill = guide_colourbar(title.position = "top", title.hjust = 0), size = guide_legend(title.position = "top", title.hjust = 0))
  

  # change internal fill aes to name
  p$layers[[1]]$mapping = aes(fill = name)
  # set the color of the areas (fill)
  p = p + 
        scale_fill_manual(values = colour_list_with_overlap) +
        theme(legend.position = "none")
  
  # set the colors of the text (does not work)
  p$layers[[4]]$mapping = aes(label = count, color = name)
  #p = p +
  #    scale_color_manual(values = colour_list_text_filtered) +
  #    theme(legend.position = "none")
  
  # set the colors of the outlines (color)
  p = p + 
        scale_color_manual(values = c("black","black","black","black","black","black","black","black","black","black")) +
        theme(legend.position = "none")
  #p

  plot_height = plot_props$image_height
  plot_width = 2/3 * plot_props$image_width
  ggsave(sprintf("%s/%s.png", output_folder_venn, file_name), p, dpi=plot_props$dpi, width=plot_width, height=plot_height, units = plot_props$image_units)
  ggsave(sprintf("%s/%s.svg", output_folder_venn, file_name), p, width=plot_width, height=plot_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
}

create_euler_plot = function(hclust, features_cluster, i, zscore_th, keep_rate, cluster_colours_list, plot_props, output_folder_euler) {
  plot_df = euler(features_cluster)
  
  labels = plot_df$original 
  labels[labels == 0] = NA
  
  p = plot(plot_df, 
           quantities = labels,
           fill = cluster_colours_list,
           labels = list(font = plot_props$font_size, family = plot_props$font_family)
           )
  
  #png(sprintf("%s/%s_zscore_th=%s_num_clusters=%s_18x6.png", output_folder_euler, i, zscore_th, length(unique(hclust$cluster))), width = 2 * plot_props$image_width, height = plot_props$image_height, units = plot_props$image_units, res = plot_props$dpi)
  #print(p)
  #dev.off()
  #svglite(sprintf("%s/%s_zscore_th=%s_num_clusters=%s_18x6.svg", output_folder_euler, i, zscore_th, length(unique(hclust$cluster))), width = 2 * plot_props$image_width / 2.54, height = plot_props$image_height*1.25 / 2.54)
  #print(p)
  #dev.off()
  png(sprintf("%s/%s_zscore_th=%s_num_clusters=%s_%s_keep_rate=%s.png", output_folder_euler, i, zscore_th, length(unique(hclust$cluster)), data_input$rna_class, keep_rate), width = plot_props$image_width, height = plot_props$image_height, units = plot_props$image_units, res = plot_props$dpi)
  print(p)
  dev.off()
  svglite(sprintf("%s/%s_zscore_th=%s_num_clusters=%s_%s_keep_rate=%s.svg", output_folder_euler, i, zscore_th, length(unique(hclust$cluster)), data_input$rna_class, keep_rate), width = plot_props$image_width / 2.54, height = plot_props$image_height*1.25 / 2.54)
  print(p)
  dev.off()
}


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
sets = snakemake@params$sets
colour_set = snakemake@params$colour_set
xticks_names = snakemake@params$xticks_names
plot_props = snakemake@params$plots_props
colours = snakemake@params$colors
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_comparing_overlapping_feature_sets.rds")
 

#------------------------------------ Script ----------------------------------- 
output_folder_venn = sprintf("%s/%s_%s/results_%s/figures/set_overlap/venn_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_venn, recursive=TRUE)

set_list = c()
colour_list = c()
colour_text_list = c()
#i = 1
for (set in sets) {
  set_list[[set$name]] = fread(set$file_path, sep='\t', header=F)
  #set_list[[set$name]]$V2 = set$name
  #if (length(set$colour) != 1) {
  #  colour_list = append(colour_list, unname(colours[[set$colour[2]]][names(sets[i])]))
  #} else {
  #  colour_list = append(colour_list, set$colour)
  #}
  colour_text_list = append(colour_text_list, "black")
  #i = i + 1
}
#colour_list = unlist(colour_list)
#names(colour_list) = names(set_list)
names(colour_text_list) =  names(set_list)

# transform all_feature_list into list of
venn_input = list()
for (origin in names(set_list)) {
  venn_input[[origin]] = set_list[[origin]]$V1
}

file_name = sprintf("overlap_%s", paste(names(sets), collapse = "_"))
create_venn_plot(venn_input, colour_set, colour_text_list, data_input, xticks_names, plot_props, output_folder_venn, file_name)
            

#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])

