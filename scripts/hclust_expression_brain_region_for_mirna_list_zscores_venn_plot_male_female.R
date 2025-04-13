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
suppressPackageStartupMessages(library(openxlsx))

# save rdata
#saveRDS(snakemake, file = "snakemake_hclust_expression_brain_region_for_mirna_list_zscores_venn_plot_male_female.rds")
#stop()

#snakemake = readRDS("snakemake_hclust_expression_brain_region_for_mirna_list_zscores_venn_plot_male_female.rds")

data.table::setDTthreads(1)
set.seed(snakemake@params$parameters_props$set_seed)

print("hclust_expression_brain_region_for_mirna_list_zscores_venn_plot_male_female")


#---------------------------------- Functions ----------------------------------
scale = function (x, rows, columns) {
  if(rows){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    x = (x - m)/s
  }
  if(columns){
    x = t(x)
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    x = (x - m)/s
    x = t(x)
  }
  return(x)
}

# write the tissue name on the side of a venn
tissue_name_side = function(plot, tissue_group, tissue, xticks_names) {
  plot_with_label = plot + annotation_custom(
    grob = textGrob(sprintf("%s", xticks_names[[tissue_group]][tissue]), rot = 270, vjust = 0.5, hjust = 0.5),
    xmin = 1000, xmax = 1000, ymin = -Inf, ymax = Inf
  )
  return(plot_with_label)
}

create_venn_plot = function(membership_table, plot_title, plot_layout, plot_index, colour_list_filtered, colour_list_text_filtered) {
  
  # manually calculate the stuff ggVennDiagramm should do automatically
  # this way we cat color the edges manually
  venn = Venn(membership_table)
  features_cluster_manual_venn = process_data(venn)
   
  #print(features_cluster_manual_venn)
  items = features_cluster_manual_venn@region$item
  item_combined = c()
  for (item in items) {
    item_combined = append(item_combined, do.call(paste, c(list(item), collapse = "\n")))
  }
  # remove mmu-miR
  item_combined = gsub("mmu-miR-", "", item_combined)
  features_cluster_manual_venn@region$item_combined = item_combined

  # plot with label list
  if ((plot_layout == "normal") | ((plot_layout == "vertical") & (plot_index != 1))) {
    p_labels = ggVennDiagram(membership_table, label_alpha=0,
                             set_size = 0, #plot_props$font_size_header * 25.4/72,# set_color = colour_list_filtered,
                             label="none", label_color = "black", label_size = plot_props$font_size * 25.4/72, 
                             edge_size = 0.5,
    )
  } else {
    p_labels = ggVennDiagram(membership_table, label_alpha=0,
                             set_size = 0, #plot_props$font_size * 25.4/72,# set_color = colour_list_filtered,
                             label="none", label_color = "black", label_size = plot_props$font_size * 25.4/72, 
                             edge_size = 0.5,
    )
  }
  p_labels +
    # text with boxes
    # geom_sf_label(aes(label = item_combined), data = venn_region(features_cluster_manual_venn), alpha=0.5) + 
    # no boxes
    geom_sf_text(aes(label = item_combined, color = name), data = venn_region(features_cluster_manual_venn), size=6 * 25.4/72) + 
    geom_sf(aes(color = name), data = venn_setedge(features_cluster_manual_venn), show.legend = F) +
    # expand the x axis
    scale_x_continuous(expand = expansion(mult = 0.2)) +
    # geom_sf(aes(fill = name), data = venn_region(features_cluster_manual_venn), show.legend = F) +
    # scale_fill_manual(values =  alpha(colour_list_filtered, .2)) +
    # scale_fill_manual(values = c("grey","grey","grey","grey")) +
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
  
  if (plot_layout == "normal") {
    p_labels = p_labels + 
      ggtitle(plot_title) +
      theme(plot.title = element_text(family = plot_props$font_family, size = plot_props$font_size))
  }
  
  # change internal fill aes to name
  p_labels$layers[[1]]$mapping = aes(fill = name)
  # set the area color (fill)
  p_labels = p_labels + 
                      scale_fill_manual(values = colour_list_filtered) +
                      theme(legend.position = "none")
  p_labels = p_labels + 
                      scale_color_manual(values = c("black","black","black","black")) +
                      theme(legend.position = "none")
  #p_labels
  
  
  # plot with counts
  if ((plot_layout == "normal") | ((plot_layout == "vertical") & (plot_index != 1))) {
    p = ggVennDiagram(membership_table, label_alpha=1,
                      set_size = 0, #plot_props$font_size_header * 25.4/72, #set_color = colour_list_filtered,
                      label = "count", label_color = "black", label_size = 6 * 25.4/72, label_geom = "text",
                      edge_size = 0.5,
    )
  } else {
    p = ggVennDiagram(membership_table, label_alpha=1,
                      set_size = 0, #plot_props$font_size * 25.4/72, #set_color = colour_list_filtered,
                      label = "count", label_color = "black", label_size = 6 * 25.4/72, label_geom = "text",
                      edge_size = 0.5,
    )
  }

  p +
  # text with boxes
  # geom_sf_label(aes(label = item_combined), data = venn_region(features_cluster_manual_venn), alpha=0.5) + 
  # no boxes
  #geom_sf_text(aes(label = item_combined), data = venn_region(features_cluster_manual_venn), size=6 * 25.4/72) + 
  geom_sf(aes(color = name), data = venn_setedge(features_cluster_manual_venn), show.legend = F) +
  # expand the x axis
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  # geom_sf(aes(fill = name), data = venn_region(feature_cluster_manual_venn), show.legend = F) +
  # scale_fill_manual(values =  alpha(colour_list_filtered, .2)) +
  # scale_fill_manual(values = c("grey","grey","grey","grey")) +
  # scale_fill_gradient(low = "#D5FBFF", high = "#0092BF")
  scale_fill_gradient(low = "white", high = "white", guide = "none") +
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
  
  if (plot_layout == "normal") {
    p = p + 
      ggtitle(plot_title) +
      theme(plot.title = element_text(family = plot_props$font_family, size = plot_props$font_size))
  } else {
    p = p +
      theme(
        plot.margin = unit(c(0,0,0,0), "cm"),
        #plot.margin = unit(c(0.3,0,-0.3,0), "cm"),
        #panel.spacing = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
      ) #+ 
  }
  
  # change internal fill aes to name
  p$layers[[1]]$mapping = aes(fill = name)
  # set the color of the areas (fill)
  p = p + 
        scale_fill_manual(values = colour_list_filtered) +
        theme(legend.position = "none")
  
  # set the colors of the text (does not work)
  p$layers[[4]]$mapping = aes(label = count, color = name)
  #p = p +
  #    scale_color_manual(values = colour_list_text_filtered) +
  #    theme(legend.position = "none")
  
  # set the colors of the outlines (color)
  p = p + 
        scale_color_manual(values = c("black","black","black","black")) +
        theme(legend.position = "none")
  #p
  
  return(list("plot" = p, "plot_labels" = p_labels))
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
method = snakemake@params$method
top_list = snakemake@params$top_list
zscore_th_list = snakemake@params$zscore_th
prop = snakemake@params$prop
xticks_names = snakemake@params$xticks_names
plot_props = snakemake@params$plot_props
colours = snakemake@params$colors
results_folder = snakemake@params$results_folder

expr = c()
annot = c()

expr$Male = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set_male, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot$Male = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set_male, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 

expr$Female = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set_female, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot$Female = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set_female, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 


#------------------------------------ Script ----------------------------------- 
method_split = data.frame(stri_split_fixed(method, "_"))[1,]
for (m in method_split) {
  #output_folder_graph = sprintf("%s/%s_%s/results_%s/figures/%s_clustering/graph_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m)
  #dir.create(output_folder_graph, recursive=TRUE)
  output_folder_venn = sprintf("%s/%s_%s/results_%s_%s/figures/%s_clustering/venn_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set_male, data_input$data_sub_set_female, m)
  dir.create(output_folder_venn, recursive=TRUE)
  output_folder_tab = sprintf("%s/%s_%s/results_%s_%s/matrices/%s_clustering/venn_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set_male, data_input$data_sub_set_female, m)
  dir.create(output_folder_tab, recursive=TRUE)
  output_folder_upset = sprintf("%s/%s_%s/results_%s_%s/figures/%s_clustering/upset_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set_male, data_input$data_sub_set_female, m)
  dir.create(output_folder_upset, recursive=TRUE)
  #output_folder_tab = sprintf("%s/%s_%s/results_%s/matrices/%s_clustering/graph_plot_info", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m)
  #dir.create(output_folder_tab, recursive=TRUE)
  
  for (i in top_list){
    print(i)
    #feature_list = c()
    #feature_list_tmp = fread(sprintf("%s/%s_%s/results_%s/matrices/%s/%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set_male, m, i))
    #feature_list$male = feature_list_tmp[[data_input$rna_class]]
    #feature_list_tmp = fread(sprintf("%s/%s_%s/results_%s/matrices/%s/%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set_female, m, i))
    #feature_list$female = feature_list_tmp[[data_input$rna_class]]
    #print(feature_list)
    for (zscore_th in zscore_th_list) {
      
      file_to_load = c()
      file_to_load$Male = sprintf("%s/%s_%s/results_%s/matrices/%s_clustering/%s_zscore_th=%s_clustered_by_%s_removed_zero_rows_plot_df.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set_male, m, i, zscore_th, m)
      file_to_load$Female = sprintf("%s/%s_%s/results_%s/matrices/%s_clustering/%s_zscore_th=%s_clustered_by_%s_removed_zero_rows_plot_df.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set_female, m, i, zscore_th, m)
      
      membership_table = list()
      for (sex_key in c("Male", "Female")) {
      # check if the file exists, if not. There were no rows left after removing all zeros in "hclust_expression_brain_region_for_miRNA_list_heatmap_plot"
        #file_to_load = sprintf("%s/%s_%s/results_%s/matrices/%s_clustering/%s_zscore_th=%s_clustered_by_%s_removed_zero_rows.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set_male, m, i, zscore_th, m)
        
        if (file.exists(file_to_load[[sex_key]])) {
          hclust = fread(file_to_load[[sex_key]])
        } else {
          print(sprintf("skipping %s because it does not exist!", file_to_load))
          next
        }
        #print(hclust)
        for (ti in unique(annot[[sex_key]][[prop]])) {
          # Ensure the list structure exists
          #if (!is.list(membership_table[[ti]])) {
          #  membership_table[[ti]] = list()
          #}
          membership_table[[ti]][[sex_key]] = hclust[hclust[[ti]] == 1,]$V1
        }
      }
            
      colour_list = c(unlist(colours[["sex"]]), unlist(colours[[prop]]))
      colour_list_text = c(unlist(colours[["sex_text"]]), unlist(colours[[sprintf("%s_text", prop)]]))
      venn_plot = list()
      venn_plot_labels = list()
      venn_plot_vertical = list()
      venn_plot_vertical_name = list()
      venn_plot_vertical_labels = list()
      venn_plot_vertical_labels_name = list()
      #upset_plot = list()
      #upset_plot_labels = list()
      suppl_tab = c()
      plot_index = 1
      for (ti in names(colours[[prop]])) {
      #for (ti in unique(annot$Male[[prop]])) {
        if (ti %in% unique(annot$Male[[prop]])) {
          colour_list_filtered = colour_list[c("Male", "Female", paste(ti))]
          colour_list_text_filtered = colour_list_text[c("Male", "Female", paste(ti))]
          names(colour_list_filtered) = c("Male", "Female", "Male..Female")
          names(colour_list_text_filtered) = c("Male", "Female", "Male..Female")
          # venn plot
          plot_title = xticks_names[[prop]][ti]
          layout = "normal"
          result = create_venn_plot(membership_table[[ti]], plot_title, layout, plot_index, colour_list_filtered, colour_list_text_filtered)
          venn_plot[[ti]] = result$plot
          venn_plot_labels[[ti]] = result$plot_labels
          
          layout = "vertical"
          result = create_venn_plot(membership_table[[ti]], plot_title, layout, plot_index, colour_list_filtered, colour_list_text_filtered)
  
          # write tissue name to one side of the plot
          #venn_plot_vertical[[ti]] = tissue_name_side(result$plot, prop, ti, xticks_names)
          #venn_plot_vertical_labels[[ti]] = tissue_name_side(result$plot_labels, prop, ti, xticks_names)
          
          venn_plot_vertical[[ti]] = result$plot
          venn_plot_vertical_labels[[ti]] = result$plot_labels
          
          venn_plot_vertical_name[[ti]] = textGrob(sprintf("%s", xticks_names[[prop]][ti]), just = "center", gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family)) #rot = 270
          venn_plot_vertical_labels_name[[ti]] = textGrob(sprintf("%s", xticks_names[[prop]][ti]), just = "center", gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family)) #rot = 270
          
          #vertical_plot_title = annotation_custom(
          #  grob = textGrob(sprintf("%s", xticks_names[[prop]][ti]), rot = 270, vjust = 0.5, hjust = 0.5),
          #  xmin = 1000, xmax = 1000, ymin = -Inf, ymax = Inf
          #)
          
          # upset plot
          #result = create_upset_plot(membership_table[[ti]], colour_list_filtered)
          #upset_plot[[ti]] = result$plot
          #upset_plot_labels[[ti]] = result$plot_labels
          
          # euler plot
          #create_euler_plot(hclust, features_cluster, i, zscore_th, keep_rate, cluster_colours_list, plot_props, output_folder_euler)
          
          plot_index = plot_index + 1
          
          # creat suppl table for publication containg the venn plot details
          both = intersect(membership_table[[ti]]$Male, membership_table[[ti]]$Female)
          only_male = setdiff(membership_table[[ti]]$Male, both)
          only_female = setdiff(membership_table[[ti]]$Female, both)
          suppl_tab = rbind(suppl_tab, data.frame(xticks_names[[prop]][[ti]], paste(only_male, collapse = ", "), paste(only_female, collapse = ", "), paste(both, collapse = ", ")))
        }
      }
      
      suppl_tab_df = data.frame(suppl_tab)
      colnames(suppl_tab_df) = c("brain_region", "only_male", "only_female", "both")
      
      # save mat
      file_name_tab = sprintf("%s/%s_zscore_th=%s_clustered_by_%s", output_folder_tab, i, zscore_th, m)
      write.table(suppl_tab_df, sprintf("%s.csv", file_name_tab), sep='\t', row.names = F)
      #write.xlsx(suppl_tab_df, sprintf("%s.xlsx", file_name_tab), colNames = TRUE, rowNames = False, append = FALSE)
      
      # 5x3 plot
      venn_plots = arrangeGrob(grobs=venn_plot, ncol=3, nrow=5)
      venn_plots_labels = arrangeGrob(grobs=venn_plot_labels, ncol=3, nrow=5)
      
      # vertical plot
      top_text = "    Male   Female"
      venn_plots_vertical_venns = arrangeGrob(grobs=venn_plot_vertical, ncol=1, nrow=15, top = textGrob(top_text, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family)))
      venn_plots_vertical_labels_venns = arrangeGrob(grobs=venn_plot_vertical_labels, ncol=1, nrow=15, top = textGrob(top_text, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family)))
      venn_plots_vertical_names = arrangeGrob(grobs=append(list(nullGrob()), venn_plot_vertical_name), ncol=1, nrow=16, heights=c(0.6, rep(1, 15)))
      venn_plots_vertical_labels_names = arrangeGrob(grobs=append(list(nullGrob()), venn_plot_vertical_labels_name), ncol=1, nrow=16, heights=c(0.6, rep(1, 15)))

      #figure_final_vertical = ggarrange(venn_plots_vertical_venns, venn_plots_vertical_names, widths=c(1.0, 0.5))
      figure_final_vertical = arrangeGrob(venn_plots_vertical_venns, venn_plots_vertical_names, widths=c(0.7, 0.5)) #widths=c(1, 0.25)
      figure_final_vertical_labels = arrangeGrob(venn_plots_vertical_labels_venns, venn_plots_vertical_labels_names, widths=c(0.7, 0.5)) #widths=c(1, 0.25)
                                                        
      #upset_plots = arrangeGrob(grobs=upset_plot, ncol=3, nrow=5, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family)), left = textGrob(y_axis_title, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family), rot = 90))
      #upset_plots_labels = arrangeGrob(grobs=upset_plot_labels, ncol=3, nrow=5, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family)), left = textGrob(y_axis_title, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family), rot = 90))
      
      # save venn plot
      ggsave(sprintf("%s/%s_zscore_th=%s_clustered_by_%s.png", output_folder_venn, i, zscore_th, m), venn_plots, dpi=plot_props$dpi, width=2/3*plot_props$image_width, height=plot_props$image_height, units = plot_props$image_units)
      ggsave(sprintf("%s/%s_zscore_th=%s_clustered_by_%s.svg", output_folder_venn, i, zscore_th, m), venn_plots, width=2/3*plot_props$image_width, height=plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
      
      ggsave(sprintf("%s/%s_zscore_th=%s_clustered_by_%s_labels.png", output_folder_venn, i, zscore_th, m), venn_plots_labels, dpi=plot_props$dpi, width=2/3*plot_props$image_width, height=plot_props$image_height, units = plot_props$image_units)
      ggsave(sprintf("%s/%s_zscore_th=%s_clustered_by_%s_labels.svg", output_folder_venn, i, zscore_th, m), venn_plots_labels, width=2/3*plot_props$image_width, height=plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
      
      # save vertical plot
      ggsave(sprintf("%s/%s_zscore_th=%s_clustered_by_%s_vertical.png", output_folder_venn, i, zscore_th, m), figure_final_vertical, dpi=plot_props$dpi, width=1/3*plot_props$image_width, height=2*plot_props$image_height, units = plot_props$image_units)
      ggsave(sprintf("%s/%s_zscore_th=%s_clustered_by_%s_vertical.svg", output_folder_venn, i, zscore_th, m), figure_final_vertical, width=1/3*plot_props$image_width, height=2*plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
      
      ggsave(sprintf("%s/%s_zscore_th=%s_clustered_by_%s_vertical_labels.png", output_folder_venn, i, zscore_th, m), figure_final_vertical_labels, dpi=plot_props$dpi, width=1/3*plot_props$image_width, height=2*plot_props$image_height, units = plot_props$image_units)
      ggsave(sprintf("%s/%s_zscore_th=%s_clustered_by_%s_vertical_labels.svg", output_folder_venn, i, zscore_th, m), figure_final_vertical_labels, width=1/3*plot_props$image_width, height=2*plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
      
      
      # save upset plot
      #ggsave(sprintf("%s/%s_zscore_th=%s_clustered_by_%s.png", output_folder_upset, i, zscore_th, m), upset_plot, dpi=plot_props$dpi, width=2/3*plot_props$image_width, height=plot_props$image_height, units = plot_props$image_units)
      #ggsave(sprintf("%s/%s_zscore_th=%s_clustered_by_%s.svg", output_folder_upset, i, zscore_th, m), upset_plot, width=2/3*plot_props$image_width, height=plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
      
      #ggsave(sprintf("%s/%s_zscore_th=%s_clustered_by_%s_labels.png", output_folder_upset, i, zscore_th, m), upset_plots_labels, dpi=plot_props$dpi, width=2/3*plot_props$image_width, height=plot_props$image_height, units = plot_props$image_units)
      #ggsave(sprintf("%s/%s_zscore_th=%s_clustered_by_%s_labels.svg", output_folder_upset, i, zscore_th, m), upset_plots_labels, width=2/3*plot_props$image_width, height=plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
      
      rm(venn_plot)
      rm(venn_plot)
      rm(venn_plot_labels)
      rm(venn_plot_labels)
      rm(venn_plots)
      rm(venn_plots)
      rm(venn_plots_labels)
      rm(venn_plots_labels)
      
      rm(venn_plot_vertical)
      rm(venn_plot_vertical)
      rm(venn_plot_vertical_labels)
      rm(venn_plot_vertical_labels)
      rm(venn_plots_vertical)
      rm(venn_plots_vertical)
      rm(venn_plots_vertical_labels)
      rm(venn_plots_vertical_labels)
      
      rm(upset_plot)
      rm(upset_plot)
      rm(upset_plot_labels)
      rm(upset_plot_labels)
      rm(upset_plots)
      rm(upset_plots)
      rm(upset_plots_labels)
      rm(upset_plots_labels)
      
      gc()
    }
  }
}

#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])

