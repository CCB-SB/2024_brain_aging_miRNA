suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(viridisLite))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ggnetwork))
suppressPackageStartupMessages(library(ggVennDiagram))
suppressPackageStartupMessages(library(eulerr))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))


#snakemake = readRDS("snakemake_sig_tissue_specific_comparison.rds")

data.table::setDTthreads(1)
set.seed(snakemake@params$parameters_props$set_seed)

print("sig_tissue_specific_comparison")


#---------------------------------- Functions ----------------------------------
replace_null_with_empty_vector = function(input) {
  # because of some strange behaviour, the venn function thinks c() != NULL != character(0)
  # it seems that it needs the last
  if (length(input) == 0) {
    return(character(0))
  } else {
    return(input)
  }
}

reduce_landscape_spacing = function() {
  return(theme(plot.margin=unit(c(t=0, r=-0.25, b=0, l=-0.25), "cm")))
}

# write the tissue name on the left side of a venn
tissue_name_side = function(plot, tissue) {
  plot_with_label = plot + annotation_custom(
    grob = textGrob(sprintf("%s", tissue), rot = 90, vjust = 0.5, hjust = 0.5),
    xmin = 50, xmax = -Inf, ymin = -Inf, ymax = Inf
    )
  return(plot_with_label)
}

# write the tissue name on the top side of a venn
tissue_name_top = function(plot, tissue) {
  plot_with_label = plot + annotation_custom(
    grob = textGrob(sprintf("%s", tissue), rot = 0, vjust = 0.5, hjust = 0.5),
    xmin = -Inf, xmax = Inf, ymin = 800, ymax = Inf
  )
  return(plot_with_label)
}

venn_plot = function(diff_exp_ts_rnas, corr_ts_rnas, plot_props, file_name, file_name_tab, show_set_name = TRUE) {
  plot_df = list(diff_exp=diff_exp_ts_rnas, corr=corr_ts_rnas)
  # manually calculate the stuff ggVennDiagramm should do automatically
  # this way we can color the edges manually
  venn = Venn(plot_df)
  features_cluster_manual_venn = process_data(venn)

  items = features_cluster_manual_venn@region$item
  item_combined = c()
  for (item in items) {
    item_combined = append(item_combined, do.call(paste, c(list(item), collapse = "\n")))
  }
  
  # remove mmu-miR
  item_combined = gsub("mmu-miR-", "", item_combined)
  features_cluster_manual_venn@region$item_combined = item_combined

  
  cluster_colours_list = c("diff_exp"="grey", "corr"="red")
  # use colorRampPalette to create function that interpolates colors 
  colfunc <- colorRampPalette(cluster_colours_list)
  # call function and create vector of 15 colors
  col <- colfunc(15)
  
  cmap = colorRamp2(c(-20, 0, 20), hcl_palette = "YlOrRd")
  
  if (show_set_name) {
    set_size = plot_props$font_size_header /3
  } else {
    set_size = 0
  }
  p_labels = ggVennDiagram(plot_df,
                    set_size = set_size,
                    label="none", label_color = "black", label_size = plot_props$font_size /3,
                    edge_size = 0.5,
                    ) +
    # text with boxes
    # geom_sf_label(aes(label = item_combined), data = venn_region(features_cluster_manual_venn), alpha=0.5) + 
    # no boxes
    geom_sf_text(aes(label = item_combined), data = venn_region(features_cluster_manual_venn), size=plot_props$font_size /3) + 
    geom_sf(aes(color = name), data = venn_setedge(features_cluster_manual_venn), show.legend = F) +
    #scale_color_manual(values = cluster_colours_list) +
    # expand the x axis
    scale_x_continuous(expand = expansion(mult = 0.2)) +
    scale_y_continuous(expand = expansion(mult = 0.3)) +
    #scale_fill_distiller(palette = "YlOrRd", direction = -1, values=c(0, 5, 10)) +
    # geom_sf(aes(fill = name), data = venn_region(features_cluster_manual_venn), show.legend = F) +
    # scale_fill_manual(values =  alpha(cluster_colours_list, .2)) +
    #scale_fill_manual(values = c("grey","grey", "grey", "grey")) +
    # scale_fill_gradient(low = "#D5FBFF", high = "#0092BF") +
    scale_fill_gradient(low = "#F7D748", high = "#d62728") +
    scale_color_manual(values = c("black", "black", "black", "black")) + 
    #scale_fill_viridis_c() + 
    #scale_y_continuous(expand = expansion(mult = .1)) +
    theme(text = element_text(size = plot_props$font_size, family = plot_props$font_family),
          plot.title = element_text(family = plot_props$font_family, size = plot_props$font_size_header),
          legend.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          legend.title = element_text(family = plot_props$font_family, size = plot_props$font_size, face = "bold"),
          legend.position = "none",# legend.box = "horizontal",
          #plot.margin = unit(c(-0.25,0,0,0), "cm"),
          #panel.spacing = margin(t = 0, r = 0, b = -3, l = 0, unit = "cm"),
    ) #+ 
    #guides(fill = guide_colourbar(title.position = "top", title.hjust = 0), size = guide_legend(title.position = "top", title.hjust = 0))
  #p
  
  ggsave(sprintf("%s_labels.png", file_name), p_labels, dpi=plot_props$dpi, width=plot_props$image_width, height=plot_props$image_height, units = plot_props$image_units)
  ggsave(sprintf("%s_labels.svg", file_name), p_labels, width=plot_props$image_width, height=plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)

  if (show_set_name) {
    set_size = plot_props$font_size_header /3
  } else {
    set_size = 0
  }
  p = ggVennDiagram(plot_df, 
                    set_size = set_size,
                    label="count", label_color = "black", label_size = plot_props$font_size /3,
                    edge_size = 0.5,
  ) +
    # text with boxes
    # geom_sf_label(aes(label = item_combined), data = venn_region(features_cluster_manual_venn), alpha=0.5) + 
    # no boxes
    #geom_sf_text(aes(label = item_combined), data = venn_region(features_cluster_manual_venn), size=6 * 25.4/72) + 
    geom_sf(aes(color = name), data = venn_setedge(features_cluster_manual_venn), show.legend = F) +
    #scale_color_manual(values = cluster_colours_list) +
    # expand the x axis
    scale_x_continuous(expand = expansion(mult = 0.2)) +
    scale_y_continuous(expand = expansion(mult = 0.3)) +
    # geom_sf(aes(fill = name), data = venn_region(features_cluster_manual_venn), show.legend = F) +
    # scale_fill_manual(values =  alpha(cluster_colours_list, .2)) +
    # scale_fill_manual(values = c("grey","grey","grey","grey")) +
    # scale_fill_gradient(low = "#D5FBFF", high = "#0092BF") +
    scale_fill_gradient(low = "#F7D748", high = "#d62728") +
    scale_color_manual(values = c("black", "black", "black", "black")) + 
    #scale_fill_gradient(low = "white", high = "white", guide = "none") + 
    #scale_fill_viridis_c() + 
    theme(text = element_text(size = plot_props$font_size, family = plot_props$font_family),
          plot.title = element_text(family = plot_props$font_family, size = plot_props$font_size_header),
          legend.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          legend.title = element_text(family = plot_props$font_family, size = plot_props$font_size, face = "bold"),
          legend.position = "none",# legend.box = "horizontal",
          #plot.margin = unit(c(-0.25,0,0,0), "cm"),
          #panel.spacing = margin(t = 0, r = 0, b = -3, l = 0, unit = "cm"),
    ) #+ 
  #guides(fill = guide_colourbar(title.position = "top", title.hjust = 0), size = guide_legend(title.position = "top", title.hjust = 0))
  #p

  ggsave(sprintf("%s.png", file_name), p, dpi=plot_props$dpi, width=plot_props$image_width, height=plot_props$image_height, units = plot_props$image_units)
  ggsave(sprintf("%s.svg", file_name), p, width=plot_props$image_width, height=plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
  
  # pad both lists
  max_length = max(length(plot_df$diff_exp), length(plot_df$corr))
  if (max_length == 0) {
    # make a least one Element if both are empty
    max_length = 1
  }
  list_diff_exp_padded = c(plot_df$diff_exp, rep(NA, max_length - length(plot_df$diff_exp)))
  list_corr_padded = c(plot_df$corr, rep(NA, max_length - length(plot_df$corr)))
  
  plot_df_padded = data.frame(diff_exp=list_diff_exp_padded, corr=list_corr_padded)

  fwrite(plot_df_padded, sprintf("%s.csv", file_name_tab), sep = "\t", row.names = FALSE)
  write.xlsx(plot_df_padded, sprintf("%s.xlsx", file_name_tab), colNames = TRUE, rowNames = FALSE, append = FALSE)
  
  return(list(p_labels = p_labels, p = p))
}

euler_plot = function(diff_exp_ts_rnas, corr_ts_rnas, plot_props, file_name) {
  plot_df = euler(featurs_cluster)
  
  labels = plot_df$original 
  labels[labels == 0] = NA
  
  p = plot(plot_df, 
           quantities = labels,
           fill = cluster_colours_list,
           labels = list(font = plot_props$font_size, family = plot_props$font_family)
           )
  
  png(sprintf("%s.png", file_name), width = plot_props$image_width, height = plot_props$image_height, units = plot_props$image_units, res = plot_props$dpi)
  print(p)
  dev.off()
  svglite(sprintf("%s.svg", file_name), width = plot_props$image_width / 2.54, height = plot_props$image_height*1.25 / 2.54)
  print(p)
  dev.off()
}


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
diff_exp_props = snakemake@params$diff_exp_props
corr_props = snakemake@params$corr_props
props = snakemake@params$props
colours = snakemake@params$colours
plot_props = snakemake@params$plots_props
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_sig_tissue_specific_comparison.rds")

output_folder_venn = sprintf("%s/%s_%s/results_%s/figures/sig_tissue_specific_comp/venn_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_venn, recursive=TRUE)
output_folder_venn_tab = sprintf("%s/%s_%s/results_%s/matrices/sig_tissue_specific_comp/venn_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_venn_tab, recursive=TRUE)
#output_folder_euler = sprintf("%s/%s_%s/results_%s/figures/sig_tissue_specific_comp/euler_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
#dir.create(output_folder_euler, recursive=TRUE)
#output_folder_tab = sprintf("%s/%s_%s/results_%s/matrices/sig_tissue_specific_comp", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
#dir.create(output_folder_tab, recursive=TRUE)


#------------------------------------ Script ----------------------------------- 
for (prop in props) {
  
  tissue = prop[1]
  time = prop[2]
  for (num_dereg_combinations in diff_exp_props$num_dereg_combinations_th) {
    for (num_dereg_tissues in diff_exp_props$num_dereg_tissues_th) {
      for (m in corr_props$method) {
        for (th in corr_props$corr_th) {
          for (num_corr_tissues in corr_props$num_corr_tissues_th) {
            # diff exp
            file_to_load_up = sprintf("%s/%s_%s/results_%s/matrices/aggregation_diff_exp/diff_exp_sig_up_%s_at_least_in_%s_%s_comps_for_%s_%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue)
            file_to_load_down = sprintf("%s/%s_%s/results_%s/matrices/aggregation_diff_exp/diff_exp_sig_down_%s_at_least_in_%s_%s_comps_for_%s_%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue)
            
            print(file_to_load_up)
            
            diff_exp_list_up = fread(file_to_load_up)
            diff_exp_list_down = fread(file_to_load_down)
              
            #diff_exp_list_up$direction = "up"
            #diff_exp_list_down$direction = "down"
            
            #diff_exp_list = rbind(diff_exp_list_up, diff_exp_list_down, fill=TRUE)
            #diff_exp_list[is.na(diff_exp_list)] = 0
    
            diff_exp_list_up_rownames = diff_exp_list_up$V1
            diff_exp_list_up$V1 = c()
            tmp = (rowSums(diff_exp_list_up) == 1)
            diff_exp_list_up_ts = diff_exp_list_up[tmp,]
            diff_exp_list_up_rownames_ts = diff_exp_list_up_rownames[tmp]
            
            diff_exp_list_down_rownames = diff_exp_list_down$V1
            diff_exp_list_down$V1 = c()
            tmp = (rowSums(diff_exp_list_down) == 1)
            diff_exp_list_down_ts = diff_exp_list_down[tmp,]
            diff_exp_list_down_rownames_ts = diff_exp_list_down_rownames[tmp]
            
            #diff_exp_dir = diff_exp_list$dir
            #diff_exp_list$dir = c()
            #diff_exp_dir_ts = diff_exp_dir[tmp]
    
            
            # corr
            file_to_load_pos = sprintf("%s/%s_%s/results_%s/matrices/aggregation_corr/corr_method=%s_sig_pos_for_%s_%s_%s_corr_th=%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m, num_corr_tissues, tissue, data_input$rna_class, th)
            file_to_load_neg = sprintf("%s/%s_%s/results_%s/matrices/aggregation_corr/corr_method=%s_sig_neg_for_%s_%s_%s_corr_th=%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m, num_corr_tissues, tissue, data_input$rna_class, th)
            
            corr_list_pos = fread(file_to_load_pos)
            corr_list_neg = fread(file_to_load_neg)
            
            #if ((nrow(corr_list_pos) == 0) && (nrow(corr_list_neg) == 0)) {
            #  print("no pos and no neg correlated")
            #  next
            #}
            
            #corr_list_pos$direction = "pos"
            #corr_list_neg$direction = "neg"
            
            #corr_list = rbind(corr_list_pos, corr_list_neg, fill=TRUE)
            #corr_list[is.na(corr_list)] = 0
            
            corr_list_pos_rownames = corr_list_pos$V1
            corr_list_pos$V1 = c()
            tmp = (rowSums(corr_list_pos) == 1)
            corr_list_pos_ts = corr_list_pos[tmp,]
            corr_list_pos_rownames_ts = corr_list_pos_rownames[tmp]
            
            corr_list_neg_rownames = corr_list_neg$V1
            corr_list_neg$V1 = c()
            tmp = (rowSums(corr_list_neg) == 1)
            corr_list_neg_ts = corr_list_neg[tmp,]
            corr_list_neg_rownames_ts = corr_list_neg_rownames[tmp]
            
            #corr_dir = corr_list_dir$dir
            #corr_list$dir = c()
            #corr_dir_ts = corr_dir[tmp]
  
            # merge          
            tissue_list_up_pos = union(colnames(diff_exp_list_up_ts), colnames(corr_list_pos_ts))
            tissue_list_down_neg = union(colnames(diff_exp_list_down_ts), colnames(corr_list_neg_ts))
            tissue_list_merged = unique(c(tissue_list_up_pos, tissue_list_down_neg))
            
            
            # venn plot
            # both
            # only if not empty
            if (length(tissue_list_merged) > 0) {
              p_both_labels = list()
              p_both = list()
              p_both_landscape_labels = list()
              p_both_landscape = list()
              i = 1
              for (ti in tissue_list_merged) {
                diff_exp_list_up_ts_rnas = diff_exp_list_up_rownames_ts[diff_exp_list_up_ts[[ti]] == 1] 
                corr_list_pos_ts_rnas = corr_list_pos_rownames_ts[corr_list_pos_ts[[ti]] == 1]
                
                diff_exp_list_down_ts_rnas = diff_exp_list_down_rownames_ts[diff_exp_list_down_ts[[ti]] == 1] 
                corr_list_neg_ts_rnas = corr_list_neg_rownames_ts[corr_list_neg_ts[[ti]] == 1]
                
                if (length(intersect(diff_exp_list_up_ts_rnas, corr_list_neg_ts_rnas)) > 0) {
                  print(sprintf("!Same %s in sig. up regulated and sig negative correlated in tissue %s for given parameters!", data_input$rna_class, ti))
                  opt = options(show.error.messages = FALSE)
                  on.exit(options(opt))
                  stop()
                } else if (length(intersect(diff_exp_list_down_ts_rnas, corr_list_pos_ts_rnas)) > 0) {
                  print(sprintf("!Same %s in sig. down regulated and sig postive correlated in tissue %s for given parameters!", data_input$rna_class, ti))
                  opt = options(show.error.messages = FALSE)
                  on.exit(options(opt))
                  stop()
                }
                
                diff_exp_list_ts_rnas = unique(c(diff_exp_list_up_ts_rnas, diff_exp_list_down_ts_rnas))
                corr_list_ts_rnas = unique(c(corr_list_pos_ts_rnas, corr_list_neg_ts_rnas))
                
                # remove NULL
                diff_exp_list_ts_rnas = replace_null_with_empty_vector(diff_exp_list_ts_rnas)
                corr_list_ts_rnas = replace_null_with_empty_vector(corr_list_ts_rnas)
                
                # if both sets are empty
                if ((length(diff_exp_list_ts_rnas) == 0) && (length(corr_list_ts_rnas) == 0)) {
                  print(sprintf("no sig corr or dregulated %ss", data_input$rna_class))
                  next
                }
                
                file_name = sprintf("sig_%s_specific_comparison_diff_exp_both_%s_at_least_in_%s_%s_comps_for_%s_%s_corr_method=%s_both_for_%s_%s_%s_corr_th=%s", ti, data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th)
                file_path_venn_plot = sprintf("%s/%s", output_folder_venn, file_name)
                file_path_venn_tab = sprintf("%s/%s", output_folder_venn_tab, file_name)
  
                # first, we run the plot function to fill without set name for the variable p_output
                # then, we plot it again with show_set_name = TRUE to save the right plot
                # and fill the variable p_output_with_set_name
                p_output = venn_plot(diff_exp_list_ts_rnas, corr_list_ts_rnas, plot_props, file_path_venn_plot, file_path_venn_tab, show_set_name = FALSE)
                p_output_with_set_name = venn_plot(diff_exp_list_ts_rnas, corr_list_ts_rnas, plot_props, file_path_venn_plot, file_path_venn_tab, show_set_name = TRUE)
  
                #p_both_labels[[i]] = p_output$p_labels + theme(plot.title = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) + labs(title = sprintf("%s", ti))
                #p_both[[i]] = p_output$p + theme(plot.title = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) + labs(title = sprintf("%s", ti))
                
                # write the tissue names on the left side (vertical)
                if (i == 1) {
                  # first plot with set names
                  p_both_labels[[i]] = tissue_name_side(p_output_with_set_name$p_labels, ti)
                  p_both[[i]] = tissue_name_side(p_output_with_set_name$p, ti)
                } else {
                  # the rest without
                  p_both_labels[[i]] = tissue_name_side(p_output$p_labels, ti)
                  p_both[[i]] = tissue_name_side(p_output$p, ti)
                }
                # on top (landscape)
                # never print set name
                p_both_landscape_labels[[i]] = tissue_name_top(p_output$p_labels, ti) + reduce_landscape_spacing()
                p_both_landscape[[i]] = tissue_name_top(p_output$p, ti) + reduce_landscape_spacing()
                
                # euler plot
                #file_name_euler = sprintf("%s/sig_tissue_specific_comparison_diff_exp_%s_%s_at_least_in_%s_%s_comps_for_%s_%s_corr_method=%s_%s_for_%s_%s_%s_corr_th=%s", output_folder_euler, diff_exp_dir, data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, corr_dir, num_corr_tissues, tissue, data_input$rna_class, th)
                #euler_plot(diff_exp_ts_rnas, corr_ts_rnas, plot_props, file_name_euler)
                
                i = i + 1
              }
              
              # merged plots
              # both
              number_of_plots = length(tissue_list_merged)
              file_name_venn_merged = sprintf("sig_specific_comparison_diff_exp_both_%s_at_least_in_%s_%s_comps_for_%s_%s_corr_method=%s_both_for_%s_%s_%s_corr_th=%s", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th)
              file_path_venn_merged = sprintf("%s/%s", output_folder_venn, file_name_venn_merged)
              
              # vertical
              # creating concatinated plot
              all_p_labels = arrangeGrob(grobs=p_both_labels, ncol=1, nrow=number_of_plots, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
              all_p = arrangeGrob(grobs=p_both, ncol=1, nrow=number_of_plots, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
              # save plots
              ggsave(sprintf("%s_vertical_labels.png", file_path_venn_merged), all_p_labels, dpi=plot_props$dpi, width=0.5*plot_props$image_width, height=0.5*number_of_plots*plot_props$image_height, units = plot_props$image_units)
              ggsave(sprintf("%s_vertical_labels.svg", file_path_venn_merged), all_p_labels, width=0.5*plot_props$image_width, height=0.5*number_of_plots*plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
              
              ggsave(sprintf("%s_vertical.png", file_path_venn_merged), all_p, dpi=plot_props$dpi, width=0.5*plot_props$image_width, height=0.5*number_of_plots*plot_props$image_height, units = plot_props$image_units)
              ggsave(sprintf("%s_vertical.svg", file_path_venn_merged), all_p, width=0.5*plot_props$image_width, height=0.5*number_of_plots*plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
            
              # horizontal
              # creating concatinated plot
              all_p_labels = arrangeGrob(grobs=p_both_landscape_labels, ncol=number_of_plots, nrow=1, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
              all_p = arrangeGrob(grobs=p_both_landscape, ncol=number_of_plots, nrow=1, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
              
              # save plots
              ggsave(sprintf("%s_landscape_labels.png", file_path_venn_merged), all_p_labels, dpi=plot_props$dpi, width=2*plot_props$image_width, height=(2 * plot_props$image_height / 3), units = plot_props$image_units)
              ggsave(sprintf("%s_landscape_labels.svg", file_path_venn_merged), all_p_labels, width=2*plot_props$image_width, height=(2 * plot_props$image_height / 3), units =plot_props$image_units, limitsize = FALSE, scale = 1)
              
              ggsave(sprintf("%s_landscape.png", file_path_venn_merged), all_p, dpi=plot_props$dpi, width=2*plot_props$image_width, height=(2 * plot_props$image_height/ 3), units = plot_props$image_units)
              ggsave(sprintf("%s_landscape.svg", file_path_venn_merged), all_p, width=2*plot_props$image_width, height=(2 * plot_props$image_height / 3), units =plot_props$image_units, limitsize = FALSE, scale = 1)
            }
            
            # up and pos
            # only if not empty
            if (length(tissue_list_up_pos) > 0) {
              p_up_pos_labels = list()
              p_up_pos = list()
              p_up_pos_landscape_labels = list()
              p_up_pos_landscape = list()
              i = 1
              for (ti in tissue_list_up_pos) {
                diff_exp_list_up_ts_rnas = diff_exp_list_up_rownames_ts[diff_exp_list_up_ts[[ti]] == 1] 
                corr_list_pos_ts_rnas = corr_list_pos_rownames_ts[corr_list_pos_ts[[ti]] == 1]
                
                # remove NULL
                diff_exp_list_up_ts_rnas = replace_null_with_empty_vector(diff_exp_list_up_ts_rnas)
                corr_list_pos_ts_rnas = replace_null_with_empty_vector(corr_list_pos_ts_rnas)
                
                # if both sets are empty
                if ((length(diff_exp_list_up_ts_rnas) == 0) && (length(corr_list_pos_ts_rnas) == 0)) {
                  print("no sig up and no sig pos")
                  next
                }
                
                file_name_up_pos = sprintf("sig_%s_specific_comparison_diff_exp_up_%s_at_least_in_%s_%s_comps_for_%s_%s_corr_method=%s_pos_for_%s_%s_%s_corr_th=%s", ti, data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th)
                file_path_up_pos_venn_plot = sprintf("%s/%s", output_folder_venn, file_name_up_pos)
                file_path_up_pos_venn_tab = sprintf("%s/%s", output_folder_venn_tab, file_name_up_pos)
                
                p_output = venn_plot(diff_exp_list_up_ts_rnas, corr_list_pos_ts_rnas, plot_props, file_path_venn_plot, file_path_venn_tab, show_set_name = FALSE)
                p_output_with_set_name = venn_plot(diff_exp_list_up_ts_rnas, corr_list_pos_ts_rnas, plot_props, file_path_venn_plot, file_path_venn_tab, show_set_name = TRUE)
                
                # write the tissue names on the left side (vertical)
                if (i == 1) {
                  # first plot with set names
                  p_up_pos_labels[[i]] = tissue_name_side(p_output_with_set_name$p_labels, ti)
                  p_up_pos[[i]] = tissue_name_side(p_output_with_set_name$p, ti)
                } else {
                  # the rest without
                  p_up_pos_labels[[i]] = tissue_name_side(p_output$p_labels, ti)
                  p_up_pos[[i]] = tissue_name_side(p_output$p, ti)
                }
                # on top (landscape)
                # never print set name
                p_up_pos_landscape_labels[[i]] = tissue_name_top(p_output$p_labels, ti) + reduce_landscape_spacing()
                p_up_pos_landscape[[i]] = tissue_name_top(p_output$p, ti) + reduce_landscape_spacing()
                
                #p_output = venn_plot(diff_exp_list_up_ts_rnas, corr_list_pos_ts_rnas, plot_props, file_path_up_pos_venn_plot, file_path_up_pos_venn_tab)
                #p_up_pos_labels[[i]] = p_output$p_labels
                #p_up_pos[[i]] = p_output$p
                
                i = i + 1
              }
              
              # merged plots
              number_of_plots = length(tissue_list_up_pos)
              file_name_venn_merged = sprintf("sig_specific_comparison_diff_exp_up_%s_at_least_in_%s_%s_comps_for_%s_%s_corr_method=%s_pos_for_%s_%s_%s_corr_th=%s", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th)
              file_path_venn_merged = sprintf("%s/%s", output_folder_venn, file_name_venn_merged)
              
              # vertical
              # creating concatinated plot
              all_p_labels = arrangeGrob(grobs=p_up_pos_labels, ncol=1, nrow=number_of_plots, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
              all_p = arrangeGrob(grobs=p_up_pos, ncol=1, nrow=number_of_plots, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
              
              # save plots
              ggsave(sprintf("%s_vertical_labels.png", file_path_venn_merged), all_p_labels, dpi=plot_props$dpi, width=0.5*plot_props$image_width, height=0.5*number_of_plots*plot_props$image_height, units = plot_props$image_units)
              ggsave(sprintf("%s_vertical_labels.svg", file_path_venn_merged), all_p_labels, width=0.5*plot_props$image_width, height=0.5*number_of_plots*plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
              
              ggsave(sprintf("%s_vertical.png", file_path_venn_merged), all_p, dpi=plot_props$dpi, width=0.5*plot_props$image_width, height=0.5*number_of_plots*plot_props$image_height, units = plot_props$image_units)
              ggsave(sprintf("%s_vertical.svg", file_path_venn_merged), all_p, width=0.5*plot_props$image_width, height=0.5*number_of_plots*plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
              
              # horizontal
              # creating concatinated plot
              all_p_labels = arrangeGrob(grobs=p_up_pos_landscape_labels, ncol=number_of_plots, nrow=1, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
              all_p = arrangeGrob(grobs=p_up_pos_landscape, ncol=number_of_plots, nrow=1, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
              
              # save plots
              ggsave(sprintf("%s_landscape_labels.png", file_path_venn_merged), all_p_labels, dpi=plot_props$dpi, width=2*plot_props$image_width, height=(2 * plot_props$image_height /3), units = plot_props$image_units)
              ggsave(sprintf("%s_landscape_labels.svg", file_path_venn_merged), all_p_labels, width=2*plot_props$image_width, height=(2 * plot_props$image_height / 3), units =plot_props$image_units, limitsize = FALSE, scale = 1)
              
              ggsave(sprintf("%s_landscape.png", file_path_venn_merged), all_p, dpi=plot_props$dpi, width=2*plot_props$image_width, height=(2 * plot_props$image_height / 3), units = plot_props$image_units)
              ggsave(sprintf("%s_landscape.svg", file_path_venn_merged), all_p, width=2*plot_props$image_width, height=(2 * plot_props$image_height / 3), units =plot_props$image_units, limitsize = FALSE, scale = 1)
            }
            
            # down and neg
            # only if not empty
            if (length(tissue_list_down_neg) > 0) {
              p_down_neg_labels = list()
              p_down_neg = list()
              p_down_neg_landscape_labels = list()
              p_down_neg_landscape = list()
              i = 1
              for (ti in tissue_list_down_neg) {
                diff_exp_list_down_ts_rnas = diff_exp_list_down_rownames_ts[diff_exp_list_down_ts[[ti]] == 1] 
                corr_list_neg_ts_rnas = corr_list_neg_rownames_ts[corr_list_neg_ts[[ti]] == 1]
                
                # remove NULL
                diff_exp_list_down_ts_rnas = replace_null_with_empty_vector(diff_exp_list_down_ts_rnas)
                corr_list_neg_ts_rnas = replace_null_with_empty_vector(corr_list_neg_ts_rnas)
                
                # if both sets are empty
                if ((length(diff_exp_list_down_ts_rnas) == 0) && (length(corr_list_neg_ts_rnas) == 0)) {
                  print("no sig down and no sig neg")
                  next
                }
                
                file_name_down_neg = sprintf("sig_%s_specific_comparison_diff_exp_down_%s_at_least_in_%s_%s_comps_for_%s_%s_corr_method=%s_neg_for_%s_%s_%s_corr_th=%s", ti, data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th)
                file_path_down_neg_venn_plot = sprintf("%s/%s", output_folder_venn, file_name_down_neg)
                file_path_down_neg_venn_tab = sprintf("%s/%s", output_folder_venn_tab, file_name_down_neg)
                
                p_output = venn_plot(diff_exp_list_down_ts_rnas, corr_list_neg_ts_rnas, plot_props, file_path_venn_plot, file_path_venn_tab, show_set_name = FALSE)
                p_output_with_set_name = venn_plot(diff_exp_list_down_ts_rnas, corr_list_neg_ts_rnas, plot_props, file_path_venn_plot, file_path_venn_tab, show_set_name = TRUE)
                
                # write the tissue names on the left side (vertical)
                if (i == 1) {
                  # first plot with set names
                  p_down_neg_labels[[i]] = tissue_name_side(p_output_with_set_name$p_labels, ti)
                  p_down_neg[[i]] = tissue_name_side(p_output_with_set_name$p, ti)
                } else {
                  # the rest without
                  p_down_neg_labels[[i]] = tissue_name_side(p_output$p_labels, ti)
                  p_down_neg[[i]] = tissue_name_side(p_output$p, ti)
                }
                # on top (landscape)
                # never print set name
                p_down_neg_landscape_labels[[i]] = tissue_name_top(p_output$p_labels, ti) + reduce_landscape_spacing()
                p_down_neg_landscape[[i]] = tissue_name_top(p_output$p, ti) + reduce_landscape_spacing()
                
                #p_output = venn_plot(diff_exp_list_down_ts_rnas, corr_list_neg_ts_rnas, plot_props, file_path_down_neg_venn_plot, file_path_down_neg_venn_tab)
                #p_down_neg_labels[[i]] = p_output$p_labels
                #p_down_neg[[i]] = p_output$p
                
                i = i + 1
              }
              # merged plots
              # down and neg
              number_of_plots = length(tissue_list_down_neg)
              file_name_venn_merged = sprintf("sig_specific_comparison_diff_exp_down_%s_at_least_in_%s_%s_comps_for_%s_%s_corr_method=%s_neg_for_%s_%s_%s_corr_th=%s", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th)
              file_path_venn_merged = sprintf("%s/%s", output_folder_venn, file_name_venn_merged)
              
              # vertical
              # creating concatinated plot
              all_p_labels = arrangeGrob(grobs=p_down_neg_labels, ncol=1, nrow=number_of_plots, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
              all_p = arrangeGrob(grobs=p_down_neg, ncol=1, nrow=number_of_plots, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
              
              # save plots
              ggsave(sprintf("%s_vertical_labels.png", file_path_venn_merged), all_p_labels, dpi=plot_props$dpi, width=0.5*plot_props$image_width, height=0.5*number_of_plots*plot_props$image_height, units = plot_props$image_units)
              ggsave(sprintf("%s_vertical_labels.svg", file_path_venn_merged), all_p_labels, width=0.5*plot_props$image_width, height=0.5*number_of_plots*plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
              
              ggsave(sprintf("%s_vertical.png", file_path_venn_merged), all_p, dpi=plot_props$dpi, width=0.5*plot_props$image_width, height=0.5*number_of_plots*plot_props$image_height, units = plot_props$image_units)
              ggsave(sprintf("%s_vertical.svg", file_path_venn_merged), all_p, width=0.5*plot_props$image_width, height=0.5*number_of_plots*plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
            
              # horizontal
              # creating concatinated plot
              all_p_labels = arrangeGrob(grobs=p_down_neg_landscape_labels, ncol=number_of_plots, nrow=1, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
              all_p = arrangeGrob(grobs=p_down_neg_landscape, ncol=number_of_plots, nrow=1, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
              
              # save plots
              ggsave(sprintf("%s_landscape_labels.png", file_path_venn_merged), all_p_labels, dpi=plot_props$dpi, width=2*plot_props$image_width, height=(2 * plot_props$image_height / 3), units = plot_props$image_units)
              ggsave(sprintf("%s_landscape_labels.svg", file_path_venn_merged), all_p_labels, width=2*plot_props$image_width, height=(2 * plot_props$image_height / 3), units =plot_props$image_units, limitsize = FALSE, scale = 1)
              
              ggsave(sprintf("%s_landscape.png", file_path_venn_merged), all_p, dpi=plot_props$dpi, width=2*plot_props$image_width, height=(2 * plot_props$image_height / 3), units = plot_props$image_units)
              ggsave(sprintf("%s_landscape.svg", file_path_venn_merged), all_p, width=2*plot_props$image_width, height=(2 * plot_props$image_height / 3), units =plot_props$image_units, limitsize = FALSE, scale = 1)
            }
          }
        }
      }
    }
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
