suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))

#snakemake = readRDS("snakemake_isomir_line_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("isomir_line_plot")


#---------------------------------- Functions ----------------------------------
rep_along = function(x, values) {
  rep(values, length.out = length(x))
}


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
plot_group = snakemake@params$plot_group
comp_group = snakemake@params$comp_group
feature = snakemake@params$feature
feature_colours = snakemake@params$colours
feature_col = snakemake@params$feature_col
manual_plot_list = snakemake@params$manual_plot
plots_props = snakemake@params$plots_props
xticks_names = snakemake@params$xticks_names
colours = snakemake@params$colours
colors = snakemake@params$colors
plot_layout = snakemake@params$plot_layout
plot_size = snakemake@params$plot_size
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_isomir_line_plot.rds")
#stop()

expr = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 


#------------------------------------ Script ----------------------------------- 
for (i in 1:length(manual_plot_list)) {
  manual_plot = manual_plot_list[[i]]
  
  # create output folders
  output_folder = sprintf("%s/%s_%s/results_%s/figures/%s/line_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, data_input$rna_class)
  dir.create(output_folder, recursive=TRUE)
  
  expr_feature = expr[expr[[feature_col]] == feature,]
  
  print(sprintf("There are %s different iso types for %s", dim(expr_feature)[1] - 1, feature ))
  
  if (length(expr_feature$iso_type) != length(unique(expr_feature$iso_type))) {
    for (duplicate in expr_feature$iso_type[!is.na(expr_feature$iso_type)]) {
      if (sum(expr_feature$iso_type[!is.na(expr_feature$iso_type)] == duplicate) > 1) {
        print(sprintf("There are more than one occurances of the same iso type (%s) for %s present", duplicate, feature))
        #print(sprintf("Skip %s", type))
      }
    }
  }
  
  expr_feature_plot_group_medians = list()
  for (group_1 in unique(annot[[plot_group]])) {
    group_1 = as.character(group_1)
    group_ID = annot[annot[[plot_group]] == group_1,][[data_input$identifier_column]]
    expr_feature_plot_group = expr_feature[, ..group_ID]
    
    expr_feature_plot_group_comp_medians = c()
    for (group_2 in unique(annot[[comp_group]])) {
      group_2 = as.character(group_2)
      group_ID = annot[(annot[[comp_group]] == group_2) & (annot[[plot_group]] == group_1),][[data_input$identifier_column]]
      expr_feature_plot_group_comp = expr_feature_plot_group[,..group_ID]
      expr_feature_plot_group_comp_medians = rbind(expr_feature_plot_group_comp_medians, data.frame(iso_type=expr_feature$iso_type, median=rowMedians(as.matrix(expr_feature_plot_group_comp)), comp=group_2))
    }
    expr_feature_plot_group_medians[[group_1]] = expr_feature_plot_group_comp_medians
  }
  
  maxi = c()
  for (group_1 in manual_plot$props) {
    maxi[[group_1]] = max(expr_feature_plot_group_medians[[group_1]]$median)
  }
  y_manual_height = max(unlist(maxi))
    
  line_plots = list()
  y_maxi = c()
  for (group_1 in unique(annot[[plot_group]])) {
    plot_df = expr_feature_plot_group_medians[[group_1]]
    plot_df$comp = as.numeric(plot_df$comp)
    #plot_df[[params$plot_group]] = str_split_fixed(plot_df$variable, "_", 3)[,2]
    
    plot_df = plot_df[!is.na(plot_df$iso_type),]
    
    plot_df$colour = "rest"
    plot_df[plot_df$iso_type == "0F_0T",]$colour = "canonical"
    
    plot_title = sprintf("%s", unlist(xticks_names[[plot_group]][[group_1]]))
    
    font_size_correction = 2.845
    
    x_min = min(unique(annot[[comp_group]]))
    x_max = max(unique(annot[[comp_group]])) + 2
    
    y_maxi[[group_1]] = max(expr_feature_plot_group_medians[[group_1]]$median)
    
    # create bar plot
    line_plots[[group_1]] = ggplot(plot_df, aes(x = comp, y = median, group=iso_type, color=colour)) + 
                            geom_line() +
                            ggtitle(plot_title) +
                            xlab("") +
                            ylab("") +
                            scale_x_continuous(limits = c(x_min, x_max), breaks=c(3, 10, 20, 30)) +
                            #scale_y_continuous(limits = c(0, 1.3*max(plot_df$median)), expand = c(0,0)) +
                            scale_color_manual(values = unlist(colours)) +
                            theme_classic() +
                            theme(#plot.margin = unit(c(0.1,-1.9,0.25,-1.75), plots_props$image_units), #top, right, bottom, left
                                  plot.margin = unit(c(0,-6,0,-6), plots_props$image_units), #top, right, bottom, left
                                  #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9), 
                                  plot.background = element_rect(fill='transparent', color=NA),
                                  legend.position="none", 
                                  text = element_text(family = plots_props$font_family, size = plots_props$font_size),
                                  axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
                                  axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
                                  plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header),
                                  legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
                                  legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size)
                            )
    
  }
   
  y_height = max(unlist(y_maxi))
  
  x_axis_title = xticks_names$categories[[sprintf("%s_capital", comp_group)]]
  y_axis_title = sprintf("Expr (%s)", gsub("_norm", "" ,data_input$norm))
    
  # creating concatinated plot
  plot_rows = plot_layout[1]
  plot_cols = plot_layout[2]
  
  # determine which figures are not on the left side
  left_col = c()
  for (row in 1:plot_rows - 1) {
    left_col = append(left_col, plot_cols*row + 1)
  }
  # determine which figures are not on the left side
  right_col = c()
  for (row in 1:plot_rows) {
    right_col = append(right_col, plot_cols*row)
  }
  
  line_plots_mod_margins = line_plots
  for (i in 1:length(line_plots_mod_margins)) {
    line_plots_mod_margins[[i]] = line_plots_mod_margins[[i]] + 
      scale_y_continuous(limits = c(0, y_height), expand = c(0,0))
    
    if ((!(i %in% left_col)) & (!(i %in% right_col))) {
      line_plots_mod_margins[[i]] = line_plots_mod_margins[[i]] + 
        theme(plot.margin = margin(0, -4, -7, -4, "pt")
        )
    }
    else if (i %in% left_col) {
      line_plots_mod_margins[[i]] = line_plots_mod_margins[[i]] + 
        theme(plot.margin = margin(0, -4, -7, -4, "pt")
        )
    } else if (i %in% right_col) {
      line_plots_mod_margins[[i]] = line_plots_mod_margins[[i]] + 
        theme(plot.margin = margin(0, 4, -7, -12, "pt")
        )
    }
    if (!(i %in% left_col)) {
      line_plots_mod_margins[[i]] = line_plots_mod_margins[[i]] + 
        scale_y_continuous(labels = function(breaks) {rep_along(breaks, "  ")}, limits = c(0, y_height), expand = c(0,0)) +
        theme(axis.ticks.y = element_line(color = "transparent"))
    }
  }
  
  all_line_plots = arrangeGrob(grobs=line_plots_mod_margins, ncol=plot_cols, nrow=plot_rows, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90), right = textGrob(feature, gp=gpar(fontsize=plots_props$font_size_header, fontfamily=plots_props$font_family), rot = 270))
  
  # save bar plot
  plot_width = plot_size[1]
  plot_height = plot_size[2]
  ggsave(sprintf("%s/%s.png", output_folder, feature), all_line_plots, dpi=plots_props$dpi, width=plot_width, height=plot_height, units = plots_props$image_units)
  ggsave(sprintf("%s/%s.svg", output_folder, feature), all_line_plots, width=plot_width, height=plot_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
  
  # manual plot
  # sort the figures given the order from the config
  for (group_1 in manual_plot$prop) {
    line_plots[[group_1]] = line_plots[[group_1]] + ylim(0, y_manual_height)
  }
    
  tmp = names(line_plots) %in% manual_plot$prop
  line_plots_manual = line_plots[tmp][manual_plot$prop]
  
  plot_rows = manual_plot$plot_layout[1]
  plot_cols = manual_plot$plot_layout[2]
  
  # determine which figures are not on the left side
  left_col = c()
  for (row in 1:plot_rows - 1) {
    left_col = append(left_col, plot_cols*row + 1)
  }
  # determine which figures are not on the left side
  right_col = c()
  for (row in 1:plot_rows) {
    right_col = append(right_col, plot_cols*row)
  }
  
  line_plots_manual_mod_margins = line_plots_manual
  for (i in 1:length(line_plots_manual_mod_margins)) {
    line_plots_manual_mod_margins[[i]] = line_plots_manual_mod_margins[[i]] + 
      scale_y_continuous(limits = c(0, y_manual_height), expand = c(0,0))
    
    if ((!(i %in% left_col)) & (!(i %in% right_col))) {
      line_plots_manual_mod_margins[[i]] = line_plots_manual_mod_margins[[i]] + 
        theme(plot.margin = margin(10, -2, -7, -12, "pt")
        )
    }
    else if (i %in% left_col) {
      line_plots_manual_mod_margins[[i]] = line_plots_manual_mod_margins[[i]] + 
        theme(plot.margin = margin(10, -6, -7, -8, "pt")
        )
    } else if (i %in% right_col) {
      line_plots_manual_mod_margins[[i]] = line_plots_manual_mod_margins[[i]] + 
        theme(plot.margin = margin(10, 2, -7, -16, "pt")
        )
    }
    if (!(i %in% left_col)) {
      line_plots_manual_mod_margins[[i]] = line_plots_manual_mod_margins[[i]] + 
        #scale_y_continuous(labels = function(breaks) {rep_along(breaks, "  ")}, limits = c(0, y_manual_height), expand = c(0,0)) +
        theme(axis.ticks.y = element_line(color = "transparent"),
              axis.text.y = element_text(color = "transparent")
              )
    }
  }
  
  all_line_plots_manual = arrangeGrob(grobs=line_plots_manual_mod_margins, ncol=plot_cols, nrow=plot_rows, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90), right = textGrob(feature, gp=gpar(fontsize=plots_props$font_size_header, fontfamily=plots_props$font_family), rot = 270))
  
  plot_width = manual_plot$plot_size[1]
  plot_height = manual_plot$plot_size[2]
  ggsave(sprintf("%s/%s_%s.png", output_folder, feature, paste(manual_plot$prop, collapse = "_")), all_line_plots_manual, dpi=plots_props$dpi, width=plot_width, height=plot_height, units = plots_props$image_units)
  ggsave(sprintf("%s/%s_%s.svg", output_folder, feature, paste(manual_plot$prop, collapse = "_")), all_line_plots_manual, width=plot_width, height=plot_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
  
  #if (manual_plots_parameters$plot_nrows == 1) {
  #  p_bar_manual_adj = list()
  #  p_box_manual_adj = list()
  #  for (ti in manual_plots_parameters$tissues) {
  #    #if (ti != manual_plots_parameters$tissues[1]){
  #    if (ti == manual_plots_parameters$tissues[length(manual_plots_parameters$tissues)]) {
  #      print(ti)
  #      p_bar_manual_adj[[ti]] = p_bar_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,0.25,-0.5,-0.1), "cm"))
  #      p_box_manual_adj[[ti]] = p_box_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,0.25,-0.5,-0.1), "cm"))
  #    #} else if (ti == manual_plots_parameters$tissues[length(manual_plots_parameters$tissues)]) {
  #    } else if (ti != manual_plots_parameters$tissues[1]){
  #      p_bar_manual_adj[[ti]] = p_bar_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,-0.1,-0.5,-0.1), "cm"))
  #      p_box_manual_adj[[ti]] = p_box_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,-0.1,-0.5,-0.1), "cm"))
  #    } else {
  #      p_bar_manual_adj[[ti]] = p_bar_manual[[ti]] + theme(plot.margin = unit(c(0,-0.1,-0.5,-0.35), "cm"))
  #      p_box_manual_adj[[ti]] = p_box_manual[[ti]] + theme(plot.margin = unit(c(0,-0.1,-0.5,-0.35), "cm"))
  #    } 
  #    #p_bar_manual_adj[[ti]] = p_bar_manual_adj[[ti]] + theme(plot.margin = unit(c(0,-0.15,-0.5,-0.15), "cm"))
  #    #p_box_manual_adj[[ti]] = p_box_manual_adj[[ti]] + theme(plot.margin = unit(c(0,-0.15,-0.5,-0.15), "cm"))
  #  }
  #  
  #  all_p_bar_manual = arrangeGrob(grobs=p_bar_manual_adj, ncol=manual_plots_parameters$plot_ncols, nrow=manual_plots_parameters$plot_nrows, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title_small, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
  #  all_p_box_manual = arrangeGrob(grobs=p_box_manual_adj, ncol=manual_plots_parameters$plot_ncols, nrow=manual_plots_parameters$plot_nrows, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title_small, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
  #  
  #  width = 2 * plots_props$image_width
  #  height = (plots_props$image_height*1.5)/2
  #  ggsave(sprintf("%s/%s_%sx%s.png", output_folder_bar_single_rna, paste(manual_plots_parameters$tissues, collapse = "_"), width, height), all_p_bar_manual, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
  #  ggsave(sprintf("%s/%s_%sx%s.svg", output_folder_bar_single_rna, paste(manual_plots_parameters$tissues, collapse = "_"), width, height), all_p_bar_manual, width=width, height=height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
  #  
  #  ggsave(sprintf("%s/%s_%sx%s.png", output_folder_box_single_rna, paste(manual_plots_parameters$tissues, collapse = "_"), width, height), all_p_box_manual, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
  #  ggsave(sprintf("%s/%s_%sx%s.svg", output_folder_box_single_rna, paste(manual_plots_parameters$tissues, collapse = "_"), width, height), all_p_box_manual, width=width, height=height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
  
  rm(line_plots)
  rm(line_plots_manual)
  rm(line_plots_mod_margins)
  rm(line_plots_manual_mod_margins)
  
  rm(all_line_plots_manual)
  rm(all_line_plots)
  
  gc()
}
  
#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
