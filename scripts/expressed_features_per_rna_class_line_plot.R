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

#snakemake = readRDS("snakemake_expressed_features_per_rna_class_line_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("expressed_features_per_rna_class_line_plot")


#---------------------------------- Functions ----------------------------------
feature_filtering_for_grouping = function(raw_count_ti_tp, min_detect_raw, detection_rate){
  raw_count_bin = ifelse(raw_count_ti_tp >= min_detect_raw, 1, 0)
  if (length(unique(rowMeans(raw_count_bin) >= detection_rate/100)) == 2) {
    raw_count_ti_tp_filtered = raw_count_ti_tp[rowMeans(raw_count_bin) >= detection_rate/100,]
  } else if ((length(unique(rowMeans(raw_count_bin) >= detection_rate/100) == 1)) & (unique(rowMeans(raw_count_bin) >= detection_rate/100) == TRUE)) {
    raw_count_ti_tp_filtered = raw_count_ti_tp[rowMeans(raw_count_bin) >= detection_rate/100,]
  }

  return(raw_count_ti_tp_filtered)
}


averaged_number_of_reads = function(raw_count_ti_tp_filtered){
  if (nrow(raw_count_ti_tp_filtered) == 0) {
    return(0)
  } else {
    raw_count_ti_tp_filtered_row_sums = rowSums(raw_count_ti_tp_filtered)
    raw_count_ti_tp_filtered_row_sums_mean = mean(raw_count_ti_tp_filtered_row_sums)
    #samples = colnames(raw_count_ti_tp_filtered)
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
    
    return(raw_count_ti_tp_filtered_row_sums_mean)
  }
}

create_steps = function(start, end) {
  seq(from = ceiling(start / 10) * 10, to = ceiling(end / 10) * 10, by = 10)
}

rep_along = function(x, values) {
  rep(values, length.out = length(x))
}


# ---------------------- Load data ----------------------
data_input = snakemake@params$data_input
plots_props = snakemake@params$plots_props
results_folder = snakemake@params$results_folder
parameters = snakemake@params$parameters
manual_plots_parameters = snakemake@params$manual_plots_parameters
xticks_names = snakemake@params$xticks_names
colors = snakemake@params$colors

# save rdata
#saveRDS(snakemake, file = sprintf("snakemake_expressed_features_per_rna_class_line_plot.rds"))
#stop()

grouping = parameters$grouping
split_group = grouping[1]
sub_split_group = grouping[2]
min_detect_raw = parameters$min_detect_raw
detection_rate = parameters$detection_rate
curves = parameters$curve
image_size = parameters$plot_size
image_layout = parameters$plot_layout
split_group_manual = manual_plots_parameters$tissues
image_size_manual = manual_plots_parameters$plot_size
image_layout_manual = manual_plots_parameters$plot_layout

annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
# rna_composition_complete = fread(snakemake@input$mapping_info_rna_classes, sep='\t')

output_folder_tab = sprintf("%s/%s_%s/results_%s/matrices/expressed_per_brain_region_per_rna_class/line_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_tab, recursive=TRUE)

# ----------------------- Script -----------------------
number_of_reads = c()
for (rna_class in data_input$rna_classes) {
  print(rna_class)
  # load data
  parts = strsplit(data_input$data_sub_set, "_")[[1]]
  suffix = paste(parts[-1], collapse = "_")
  raw_count = fread(sprintf("%s_%s/%s_%s_expression_raw_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, rna_class, data_input$feature_filtering ,suffix), sep='\t', header=T)
  
  number_of_reads_ti_tp = c()
  for (ti in unique(annot[[split_group]])) {
    #print(ti)
    for (tp in unique(annot[[sub_split_group]])) {
      mask = annot[((annot[[split_group]] == ti) & (annot[[sub_split_group]] == tp)),][[data_input$identifier_column]]
      raw_count_ti_tp = raw_count[, ..mask]
      raw_count_ti_tp_filtered = feature_filtering_for_grouping(raw_count_ti_tp, min_detect_raw, detection_rate)
      number_of_reads_ti_tp[[sprintf("%s_%s", ti, tp)]] = averaged_number_of_reads(raw_count_ti_tp_filtered)
    }
  }

  number_of_reads[[rna_class]] = number_of_reads_ti_tp
}

rna_comp_df = as.data.table(melt(number_of_reads))
colnames(rna_comp_df) = c("value", "grouping", "variable")

rna_comp_df$value_per = rna_comp_df$value
for (ti in unique(annot[[split_group]])) {
  for (tp in unique(annot[[sub_split_group]])) {
    gr = sprintf("%s_%s", ti, tp)
    total = sum(rna_comp_df[rna_comp_df$grouping == gr,]$value)
    rna_comp_df[rna_comp_df$grouping == gr,]$value_per = (rna_comp_df[rna_comp_df$grouping == gr,]$value * 100) / total
  }
}  

write.table(number_of_reads, sprintf("%s/expressed_count_composition_per.csv", output_folder_tab), sep='\t', row.names = F)
write.xlsx(number_of_reads, sprintf("%s/expressed_count_composition_per.xlsx", output_folder_tab), colNames = TRUE, rowNames = FALSE, append = FALSE)

tmp = str_split_fixed(rna_comp_df$grouping, "_", 2)
rna_comp_df$split_group = tmp[,1]
rna_comp_df$sub_split_group = tmp[,2] 

#rna_comp_df[, variable:=factor(variable, levels=rna_comp_df[, mean(value_per), by="variable"][order(V1, decreasing = T)]$variable)]

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

#rna_comp_df$plot_names_grouping = unname(xticks_names[[grouping]][rna_comp_df$grouping])
#rna_comp_df$plot_names_grouping = factor(rna_comp_df$plot_names_grouping, levels=unname(xticks_names[[grouping]]))

##cellline_colors_ordered = cellline_colors[unique(rna_comp_df$Protocol)]

mask = (names(unlist(colors$rna_classes)) %in% unique(rna_comp_df$variable))
colours_rna_class_filtered = unlist(colors$rna_classes)[mask]

#colours_grouping = colors[[grouping]]
#names(colours_grouping) = unname(xticks_names[[grouping]][names(colours_grouping)])
#mask = (names(unlist(colours_grouping)) %in% unique(rna_comp_df$plot_names_grouping))
#colours_grouping_filtered = unlist(colours_grouping)[mask]

#colours_grouping_text = colors[[sprintf("%s_text", grouping)]]
#names(colours_grouping_text) = unname(xticks_names[[grouping]][names(colours_grouping_text)])
#mask = (names(unlist(colours_grouping_text)) %in% unique(rna_comp_df$plot_names_grouping))
#colours_grouping_text_filtered = unlist(colours_grouping_text)[mask]

#print(rna_comp_df[rna_comp_df$variable == "miRNA",])
#print(mean(rna_comp_df[rna_comp_df$variable == "miRNA",]$value_per))

#if (length(grouping) == 1) {
factors_x_levels = rna_comp_df[variable == "miRNA"][order(value_per)]$grouping
#} else {
#  if (subgrouping != "sample") {
#    factors_x_levels = rna_comp_df[variable == "miRNA"][order(value_per)]$grouping
#  } else {
#    rna_comp_df$grouping = sprintf("%s_%s", rna_comp_df$plot_names_grouping, rna_comp_df$Sample)
#    factors_x_levels = rna_comp_df[variable == "miRNA"][order(value_per)]$grouping
#  }
#}

rna_comp_df$variable = factor(rna_comp_df$variable, levels = names(colours_rna_class_filtered))
rna_comp_df$sub_split_group = as.numeric(rna_comp_df$sub_split_group)
#rna_comp_df$sub_split_group = factor(rna_comp_df$sub_split_group, levels = unique(sort(rna_comp_df$sub_split_group)))

#curve_elem = curves$c2
for (curve_elem in curves) {
  curve_type = curve_elem[1]
  if (curve_type == "poly") {
    poly_deg = as.numeric(curve_elem[2])
    folder_name = sprintf("%s_deg=%s", curve_type, poly_deg)
  } else {
    folder_name = curve_type
  }
  print(folder_name)
  output_folder_fig = sprintf("%s/%s_%s/results_%s/figures/expressed_per_brain_region_per_rna_class/line_plots/%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, folder_name)
  dir.create(output_folder_fig, recursive=TRUE)

  output_folder_fig_single = sprintf("%s/%s_%s/results_%s/figures/expressed_per_brain_region_per_rna_class/line_plots/%s/%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, folder_name, split_group)
  dir.create(output_folder_fig_single, recursive=TRUE)

  maxi = max(rna_comp_df$value_per) + 5
  mask = (rna_comp_df$split_group %in% split_group_manual)
  maxi_manual = max(rna_comp_df[mask,]$value_per) + 5
  
  figure_range = append(min(rna_comp_df$sub_split_group), create_steps(min(rna_comp_df$sub_split_group), max(rna_comp_df$sub_split_group)))
  
  compo_grouped_list = list()
  for (ti in unique(annot[[split_group]])) {
    #print(ti)
    rna_comp_df_ti = rna_comp_df[rna_comp_df$split_group == ti,]
    
    compo_grouped = ggplot(rna_comp_df_ti, aes(x = sub_split_group, 
                                                     y = value_per, 
                                                     color = variable)) +
      geom_point(size = 1) # Scatter plot

    if (curve_type == "poly") {
      compo_grouped = compo_grouped + geom_smooth(data = rna_comp_df_ti, mapping = aes(x = sub_split_group, y = value_per, group = variable), method = "lm", formula = y ~ poly(x, poly_deg), se = FALSE, size = 0.5) # Polynomial fit per variable
    } else {
      compo_grouped = compo_grouped + geom_smooth(data = rna_comp_df_ti, mapping = aes(x = sub_split_group, y = value_per, group = variable), method = "loess", se = FALSE, size = 0.5) # Local polynomial regression fitting
    }
    
    compo_grouped = compo_grouped + 
      scale_color_manual(values = colours_rna_class_filtered) +  # Use predefined colors
      scale_x_continuous(limits = range(figure_range), breaks = figure_range) +
      ggtitle(xticks_names[[split_group]][ti]) +
      theme_classic() +
      theme(legend.position = "none", 
            text = element_text(family = plots_props$font_family, size = plots_props$font_size),
            legend.title = element_blank(), #element_text(family = plots_props$font_family, size = plots_props$font_size),
            #legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size, margin = margin(0, 0, 0, 0.3, "pt")),
            axis.text = element_text(angle = 0, family = plots_props$font_family, size = plots_props$font_size),
            axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
            #axis.text.y = element_blank(),
            #axis.ticks.y = element_blank(), 
            #axis.line.y = element_blank(),
            #strip.text.y.left = element_text(family = plots_props$font_family, angle = 0, size=plots_props$font_size),
            #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4),
            #axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed
            panel.spacing = unit(0, "pt"),
            # legend.margin = margin(0, 0, 0, 0, "pt"), legend.key.size = unit(0.3, "cm"), legend.spacing.y = unit(0.4, "cm"), legend.spacing.x = unit(0.3, "cm"),
            plot.background = element_rect(fill = "transparent"),
      )
    #compo_grouped[[ti]] 
    compo_grouped_list[[ti]] = compo_grouped + 
                                    scale_y_continuous(limits = c(0, maxi), breaks = round(pretty(c(0, maxi), n = 4), -1), expand = c(0,0)) +
                                    #ylim(0, maxi) +
                                    labs(x = "", y = "") +
                                    theme(plot.margin = margin(0, 4, 0, 0, "pt"),
                                          plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header, margin = margin(b=0)),
                                    )
    
    compo_grouped_single = compo_grouped + 
                                    labs(x = xticks_names$categories[[sprintf("%s_capital", sub_split_group)]], y = "Expressed count\ncomposition (%)") +
                                    theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
                                          plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header),
                                    )

    width = 6
    height = 6
    ggsave(sprintf("%s/grouped_by_%s_%s_%sx%s.png", output_folder_fig_single, ti, sub_split_group, height, width), compo_grouped_single, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
    ggsave(sprintf("%s/grouped_by_%s_%s_%sx%s.svg", output_folder_fig_single, ti, sub_split_group, height, width), compo_grouped_single, width = width, height = height, unit = plots_props$image_units, dpi = plots_props$dpi)
  }
  
  # plot all split_group objects
  # sort figures accroding to colour map
  compo_grouped_sorted = compo_grouped_list[names(colors[[split_group]])]
  
  # determine which figures are not on the left side
  left_col = c()
  for (row in 1:image_layout$plot_nrows - 1) {
    left_col = append(left_col, image_layout$plot_ncols*row + 1)
  }
  # determine which figures are not on the left side
  right_col = c()
  for (row in 1:image_layout$plot_nrows) {
    right_col = append(right_col, image_layout$plot_ncols*row)
  }
  
  compo_grouped_sorted_mod_margins = compo_grouped_sorted
  for (i in 1:length(compo_grouped_sorted_mod_margins)) {
    if ((!(i %in% left_col)) & (!(i %in% right_col))) {
      compo_grouped_sorted_mod_margins[[i]] = compo_grouped_sorted_mod_margins[[i]] + 
        theme(plot.margin = margin(0, -4, -7, -4, "pt")
        )
    }
    else if (i %in% left_col) {
      compo_grouped_sorted_mod_margins[[i]] = compo_grouped_sorted_mod_margins[[i]] + 
        theme(plot.margin = margin(0, -4, -7, -5, "pt")
        )
    } else if (i %in% right_col) {
      compo_grouped_sorted_mod_margins[[i]] = compo_grouped_sorted_mod_margins[[i]] + 
        theme(plot.margin = margin(0, 4, -7, -4, "pt")
        )
    }
    if (!(i %in% left_col)) {
      compo_grouped_sorted_mod_margins[[i]] = compo_grouped_sorted_mod_margins[[i]] + 
        #scale_y_continuous(labels = function(breaks) {rep_along(breaks, "  ")}, limits = c(0, maxi), expand = c(0,0)) +
        theme(axis.ticks.y = element_line(color = "transparent"),
              axis.text.y = element_text(color = "transparent")
        )
    }
  }
  
  # get the fixed width from a plot from the left
  #panel_index = which(ggplotGrob(compo_grouped_sorted_mod_margins[[1]])$layout$name == "panel")
  #fixed_widths = ggplotGrob(compo_grouped_sorted_mod_margins[[1]])$widths#[panel_index]
  
  #compo_grouped_sorted_mod_margins_axis = compo_grouped_sorted_mod_margins
  #for (i in 1:length(compo_grouped_sorted_mod_margins_axis)) {
  #  if (!(i %in% left_col)) {
  #    compo_grouped_sorted_mod_margins_axis[[i]] = compo_grouped_sorted_mod_margins_axis[[i]] + 
  #      theme(axis.text.y = element_blank(),
  #            axis.ticks.y = element_blank(),
  #            )
  #    # apply the fixed width
  #    compo_grouped_sorted_mod_margins_axis[[i]] = ggplotGrob(compo_grouped_sorted_mod_margins_axis[[i]])
  #    compo_grouped_sorted_mod_margins_axis[[i]]$widths = fixed_widths
  #  } 
  #}
  
  # merge figures to one big
  all_compo_grouped_sorted = arrangeGrob(grobs=compo_grouped_sorted_mod_margins, ncol=image_layout$plot_ncols, nrow=image_layout$plot_nrows,
                                         #padding = unit(0, "cm"),
                                         bottom = textGrob(xticks_names$categories[[sprintf("%s_capital", sub_split_group)]], gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)),
                                         left = textGrob("Expressed count composition (%)", gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
  
  # save figure
  width = image_size$image_width
  height = image_size$image_height
  ggsave(sprintf("%s/grouped_by_%s_%sx%s.png", output_folder_fig, paste(grouping, collapse = "_"), height, width), all_compo_grouped_sorted, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
  ggsave(sprintf("%s/grouped_by_%s_%sx%s.svg", output_folder_fig, paste(grouping, collapse = "_"), height, width), all_compo_grouped_sorted, width = width, height = height, unit = plots_props$image_units, dpi = plots_props$dpi)
  
  
  # for manually plot containing only given split_group objects
  mask = (names(compo_grouped_list) %in% split_group_manual)
  compo_grouped_filtered = compo_grouped_list[mask]
  compo_grouped_filtered_sorted = compo_grouped_filtered[split_group_manual]
  
  # determine which figures are not on the left side
  left_col = c()
  for (row in 1:image_layout_manual$plot_nrows - 1) {
    left_col = append(left_col, image_layout_manual$plot_ncols*row + 1)
  }
  # determine which figures are not on the left side
  right_col = c()
  for (row in 1:image_layout_manual$plot_nrows) {
    right_col = append(right_col, image_layout_manual$plot_ncols*row)
  }
  
  compo_grouped_filtered_sorted_mod_margins = compo_grouped_filtered_sorted
  for (i in 1:length(compo_grouped_filtered_sorted_mod_margins)) {
    
    compo_grouped_filtered_sorted_mod_margins[[i]] = compo_grouped_filtered_sorted_mod_margins[[i]] + scale_y_continuous(limits = c(0, maxi_manual), breaks = round(pretty(c(0, maxi), n = 4), -1), expand = c(0,0))
    
    if ((!(i %in% left_col)) & (!(i %in% right_col))) {
      compo_grouped_filtered_sorted_mod_margins[[i]] = compo_grouped_filtered_sorted_mod_margins[[i]] + 
        theme(plot.margin = margin(2, 4, -10, -8, "pt")
        )
    }
    else if (i %in% left_col) {
      compo_grouped_filtered_sorted_mod_margins[[i]] = compo_grouped_filtered_sorted_mod_margins[[i]] + 
        theme(plot.margin = margin(2, 4, -10, -8, "pt")
        )
    } else if (i %in% right_col) {
      compo_grouped_filtered_sorted_mod_margins[[i]] = compo_grouped_filtered_sorted_mod_margins[[i]] + 
        theme(plot.margin = margin(2, 4, -10, -8, "pt")
        )
    }
    if (!(i %in% left_col)) {
      compo_grouped_filtered_sorted_mod_margins[[i]] = compo_grouped_filtered_sorted_mod_margins[[i]] + 
        #scale_y_continuous(labels = function(breaks) {rep_along(breaks, "  ")}, limits = c(0, maxi_manual), expand = c(0,0)) +
        theme(axis.ticks.y = element_line(color = "transparent"),
              axis.text.y = element_text(color = "transparent")
        )
    }
  }
  
  # merge figures to one big
  all_compo_grouped_filtered_sorted = arrangeGrob(grobs=compo_grouped_filtered_sorted_mod_margins, ncol=image_layout_manual$plot_ncols, nrow=image_layout_manual$plot_nrows,
                                                   bottom = textGrob(xticks_names$categories[[sprintf("%s_capital", sub_split_group)]], gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)),
                                                   left = textGrob("Expressed count composition (%)", gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
  
  # save figure
  width = image_size_manual$image_width
  height = image_size_manual$image_height
  ggsave(sprintf("%s/grouped_by_%s_%s_%sx%s.png", output_folder_fig, paste(split_group_manual, collapse = "_"), sub_split_group, height, width), all_compo_grouped_filtered_sorted, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
  ggsave(sprintf("%s/grouped_by_%s_%s_%sx%s.svg", output_folder_fig, paste(split_group_manual, collapse = "_"), sub_split_group, height, width), all_compo_grouped_filtered_sorted, width = width, height = height, unit = plots_props$image_units, dpi = plots_props$dpi)

  rm(compo_grouped)
  rm(compo_grouped_single)
  rm(compo_grouped_list)
  rm(compo_grouped_sorted)
  rm(compo_grouped_sorted_mod_margins)
  rm(all_compo_grouped_sorted)
  rm(compo_grouped_filtered)
  rm(compo_grouped_filtered_sorted)
  rm(all_compo_grouped_filtered_sorted)
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])

