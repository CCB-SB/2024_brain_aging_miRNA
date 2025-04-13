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
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(showtext))
suppressPackageStartupMessages(library(shadowtext))

suppressPackageStartupMessages(library(png))



font_add("fa-solid", regular = "fonts/Font Awesome 6 Free-Solid-900.otf")
font_add("fa-regular", regular = "fonts/Font Awesome 6 Free-Regular-400.otf")

#suppressPackageStartupMessages(library(extrafont))
#font_import(paths = c("fonts/"), prompt = FALSE)

#suppressPackageStartupMessages(library(systemfonts))
#suppressPackageStartupMessages(library(ggtext))

#register_font("fa-regular", "fonts/Font Awesome 6 Free-Regular-400.otf")
#register_font("fa-solid",   "fonts/Font Awesome 6 Free-Solid-900.otf")


#snakemake = readRDS("snakemake_candidate_comparison_scatter_plot_tissue_split_male_female.rds")

data.table::setDTthreads(1)
set.seed(snakemake@params$parameters_props$set_seed)

print("candidate_comparison_scatter_plot_tissue_split_male_female")


#---------------------------------- Functions ----------------------------------
save_and_reload_png = function(plot_list, plot_props) {
  # the plots with the awesome icons do not like arrangeGrob
  # we write each single plot to a temporary file and load it again from the PNG file
  # then we pack it into grobs and return a list
  out_list = list()
  keys = names(plot_list)
  for (key in keys) {
    # save plot
    ggsave("/tmp/tmp_image.png", plot_list[[key]], dpi=plot_props$dpi, width=plot_props$image_width, height=plot_props$image_height, units = plot_props$image_units)
    # load plot
    png_image = readPNG("/tmp/tmp_image.png")
    out_list[[key]] = rasterGrob(png_image)
  }
  return(out_list)
}


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
  return(theme(plot.margin=unit(c(t=-0.15, r=-0.25, b=-0.15, l=-0.25), "cm")))
}

# write the tissue name on the left side of a venn
tissue_name_side = function(plot, tissue_group, tissue, xticks_names) {
  plot_with_label = plot + annotation_custom(
    grob = textGrob(sprintf("%s", xticks_names[[tissue_group]][tissue]), rot = 90, vjust = 0.5, hjust = 0.5),
    xmin = 50, xmax = -Inf, ymin = -Inf, ymax = Inf
    )
  return(plot_with_label)
}

# write the tissue name on the top side of a venn
tissue_name_top = function(plot, tissue_group, tissue, xticks_names) {
  plot_with_label = plot + annotation_custom(
    grob = textGrob(sprintf("%s", xticks_names[[tissue_group]][tissue]), rot = 0, vjust = 0.5, hjust = 0.5),
    xmin = -Inf, xmax = Inf, ymin = 800, ymax = Inf
  )
  return(plot_with_label)
}

create_scatter_plot = function(num_tissue_spec_features_df, tissue, x_axis_label, log_toggle, maxi_male_female, colours, plot_props) {
  #print(log_toggle)
  #print(maxi_male_female)
  
  # create the stuff for the colored rectangles and their labels
  # we add 3.5, 2.5, 1.5 as position for the horizontal lines
  #levels_with_dividers = c("4.5", "IV", "3.5", "III", "2.5", "II", "1.5", "I", "0.5")
  #half_yticks = c("0.5", "1.5", "2.5", "3.5", "4.5")
  tissues_sorted = unlist(xticks_names[[tissue]][names(colours[[tissue]])])
  tissues_sorted$mid = c()
  tissues_sorted = unname(unlist(tissues_sorted))
  start_step = 1
  step_size = 2
  levels_with_dividers = c(step_size*length(tissues_sorted) + start_step, rbind(tissues_sorted, rev(seq(start_step, step_size*length(tissues_sorted), step_size))))
  half_yticks = seq(start_step, step_size*length(tissues_sorted) + start_step, step_size)
  
  # create plot df
  # add zero row if any tissue is missing
  missing_tissues = names(xticks_names[[tissue]])[!(names(xticks_names[[tissue]]) %in% rownames(num_tissue_spec_features_df))]
  num_tissue_spec_features_df[missing_tissues, ] = 0
  # add tissue column and meld
  num_tissue_spec_features_df$tissue = rownames(num_tissue_spec_features_df)
  rownames(num_tissue_spec_features_df) = c()
  plot_df = melt(num_tissue_spec_features_df, id.vars = "tissue")
  
  # shapes adjustment
  #shape_formats = c("circle-arrow-up", "circle-arrow-down", "circle-plus", "circle-minus", "square", "circle")
  shape_formats = c("\uf0aa", "\uf0ab", "\uf055", "\uf056", "\uf45c", "\uf111")
  shape_families = c("fa-solid", "fa-solid", "fa-solid", "fa-solid", "fa-solid", "fa-solid")
  shape_colors = c(colours$direction$up, colours$direction$down, colours$direction$pos, colours$direction$neg, "#DB3085", "#17becf")
  names(shape_formats) = c("Up reg.", "Down reg.", "Pos. corr.", "Neg. corr.", "Region spec.", "Overlap")
  names(shape_families) = names(shape_formats)
  names(shape_colors) = names(shape_formats)
  
  shape_formats_list = c()
  shape_famlies_list = c()
  tissue_colours_list = c()
  for (i in 1:dim(plot_df)[1]) {
    shape_formats_list = append(shape_formats_list, shape_formats[plot_df[i,]$variable])
    shape_famlies_list = append(shape_famlies_list, shape_families[plot_df[i,]$variable])
    tissue_colours_list = append(tissue_colours_list, colours[[tissue]][plot_df[i,]$tissue])
  }
  plot_df$shapes = shape_formats_list
  plot_df$font_family = shape_famlies_list
  plot_df$tissue_colour = tissue_colours_list
  
  names(tissue_colours_list) = unlist(xticks_names[[tissue]][plot_df$tissue])
  
  rects = data.frame(ystart = half_yticks[-length(half_yticks)], yend = half_yticks[-1], col = tissues_sorted)

  plot_df$tissue_factor = factor(unlist(xticks_names[[tissue]][plot_df$tissue]), levels = tissues_sorted)
  
  showtext_auto()
  showtext_opts(dpi = plot_props$dpi)
  
  #print(sort(unique(plot_df[plot_df$value != 0,]$tissue)))
  
  scatter_plot = ggplot(plot_df[plot_df$value != 0,], aes(x=value, y=tissue_factor)) +
    #geom_text(aes(label=shapes, family=font_family, color=variable), size=2, show.legend=FALSE, alpha = 1) +
    scale_color_manual(values=shape_colors) +
    xlab(x_axis_label) +
    ylab("") +
    xlim(0, maxi_male_female) +
    coord_flip() +
    # only show breaks for roman numbers
    scale_y_discrete(breaks=plot_df$tissue_factor, limits=levels_with_dividers, expand=c(0, 0)) +
    #ylim(half_yticks[c(1, length(half_yticks))]) + 
    theme_classic() +
    #coord_cartesian(ylim = c(1, 4)) +  # Limit y-axis to the range of factors
    theme(plot.background = element_rect(fill='transparent', color=NA),
          legend.position="bottom", 
          text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          axis.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          axis.title = element_text(family = plot_props$font_family, size = plot_props$font_size),
          axis.ticks.x = element_blank(),
          #axis.line.x = element_line(arrow = grid::arrow()),
          axis.line = element_blank(),
          axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0, family = plot_props$font_family, size = plot_props$font_size),
          #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9), 
          plot.title = element_text(family = plot_props$font_family, size = plot_props$font_size_header),
          legend.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          legend.title = element_text(family = plot_props$font_family, size = plot_props$font_size),
          plot.margin = unit(c(0.25,0.25,0.5,0.25), "cm")
    )
  
  
  #label_function = trans_new(name="foo", transform=function(x) 10^x, inverse=function(x) x)

  if (log_toggle) {
    # maxi = round(ceiling(maxi_male_female * 10) / 10, digits = 1) # a bit over the actual max value
    maxi = log10(ceiling(10^maxi_male_female / 10) * 10) # rounf to nearest "Zehner Stelle"
    maxi_lower_int = as.integer(floor(maxi_male_female)) # the last integer
    log_breaks = c(seq(0, maxi_lower_int, 1), maxi)
    log_labels = as.integer(c(0, round(10^log_breaks[-1], digits=0)))
    print(log_breaks)
    print(log_labels)
    scatter_plot = scatter_plot + scale_x_continuous(breaks=log_breaks, labels=log_labels, limits=c(0, max(log_breaks)))
  }

  # add horizontal lines
  for (hline_pos in half_yticks[-c(1, length(half_yticks))]) {
    # durchgehende Linie
    scatter_plot = scatter_plot + geom_hline(yintercept=hline_pos, size = 0.3, color = "white")
    # Anfang und Ende variabel
    #scatter_plot = scatter_plot + geom_segment(x=0, xend=Inf, y=hline_pos, yend=hline_pos)
    
  } 
  scatter_plot_normal = scatter_plot + 
    # included coloured background according to the tissue coloures
    geom_rect(data=rects, aes(x=NULL, y=NULL, xmin=-Inf, xmax=Inf, ymin=ystart, ymax=yend, fill=col), alpha = 0.4, show.legend=FALSE) + 
    scale_fill_manual(values=tissue_colours_list) +
    geom_text(aes(label=shapes, family=font_family, color=variable), size=3, show.legend=FALSE, alpha = 0.75, position = position_jitter(width = 0, height = 0.25, seed = snakemake@params$parameters_props$set_seed)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, family = plot_props$font_family, size = plot_props$font_size))
  
  # now plot all elements which are == 0
  # for one call per Up reg./Down reg./... type
  for (variable in unique(plot_df$variable)) {
    plot_df_variable = plot_df[((plot_df$value == 0) & (plot_df$variable == variable)), ]
    scatter_plot_normal = scatter_plot_normal + # all poosible colours (colours$brain_region) are needed here!
      geom_shadowtext(data=plot_df_variable, aes(x=value, y=tissue_factor, label=shapes, family=font_family, color=variable), size=3, show.legend=FALSE, bg.r=0.1, alpha = 1, position = position_jitter(width = 0, height = 0.25, seed = snakemake@params$parameters_props$set_seed))
  }
  #scatter_plot_normal
  
  rects_flipped = data.frame(ystart = half_yticks[-length(half_yticks)], yend = half_yticks[-1], col = rev(tissues_sorted))
  
  # if log_toggle == TURE we rename the x-axis ticks to the raw counts
  if (log_toggle) {
    scatter_plot_flipped = scatter_plot + scale_x_continuous(position = "top", limits=c(0, max(log_breaks)), breaks=log_breaks, labels=log_labels) 
  } else {
    scatter_plot_flipped = scatter_plot + scale_x_continuous(position = "top", limits=c(0, maxi_male_female)) # Move y-axis to the right and set X-limits again!
  }
  
  scatter_plot_flipped = scatter_plot_flipped + 
    # sort x-axis revertible
    scale_y_discrete(breaks=rev(plot_df$tissue_factor), limits=rev(levels_with_dividers), expand=c(0, 0)) +
    # sort coloured background revertible
    geom_rect(data=rects_flipped, aes(x=NULL, y=NULL, xmin=-Inf, xmax=Inf, ymin=ystart, ymax=yend, fill=col), alpha = 0.4, show.legend=FALSE) +
    scale_fill_manual(values=tissue_colours_list) +
    geom_text(aes(label=shapes, family=font_family, color=variable), size=3, show.legend=FALSE, alpha = 0.75, position = position_jitter(width = 0, height = 0.25, seed = snakemake@params$parameters_props$set_seed)) 

  # now plot all elements which are == 0
  # for one call per Up reg./Down reg./... type
  for (variable in unique(plot_df$variable)) {
    plot_df_variable = plot_df[((plot_df$value == 0) & (plot_df$variable == variable)), ]
    scatter_plot_flipped = scatter_plot_flipped +
      geom_shadowtext(data=plot_df_variable, aes(x=value, y=tissue_factor, label=shapes, family=font_family, color=variable), size=3, show.legend=FALSE, bg.r=0.1, alpha = 1, position = position_jitter(width = 0, height = 0.25, seed = snakemake@params$parameters_props$set_seed))
  }
  
  #scatter_plot_flipped
  
  return(list(normal = scatter_plot_normal, flip = scatter_plot_flipped))
}

#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
diff_exp_props = snakemake@params$diff_exp_props
corr_props = snakemake@params$corr_props
hclust_props = snakemake@params$hclust_props
plot_layout = snakemake@params$plot_layout
props = snakemake@params$props
xticks_names = snakemake@params$xticks_names
colours = snakemake@params$colors
plot_props = snakemake@params$plots_props
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_candidate_comparison_scatter_plot_tissue_split_male_female.rds")


#------------------------------------ Script ----------------------------------- 
for (prop in props) {
  
  tissue = prop[1]
  time = prop[2]

  output_folder_scatter = sprintf("%s/%s_%s/results_%s_%s/figures/candidate_comparison_%s_split/scatter_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set_male, data_input$data_sub_set_female, tissue)
  dir.create(output_folder_scatter, recursive=TRUE)
  
  for (num_dereg_combinations in diff_exp_props$num_dereg_combinations_th) {
    for (num_dereg_tissues in diff_exp_props$num_dereg_tissues_th) {
      for (m in corr_props$method) {
        for (th in corr_props$corr_th) {
          for (num_corr_tissues in corr_props$num_corr_tissues_th) {
            for (m_hclust in hclust_props$method) {
              # only first part (for example cv from cv_zscore)
              m_hclust_split = data.frame(stri_split_fixed(m_hclust, "_"))[1,]
              for (top in hclust_props$top) {
                for (th_hclust in hclust_props$zscores_th) {
                  
                  num_tissue_spec_features_df = c()
                  num_features_df = c()
                  for (sex_key in c("male", "female")) {
                    print(sex_key)
                    
                    if (sex_key == "male") {
                      data_sub_set = data_input$data_sub_set_male
                    } else {
                      data_sub_set = data_input$data_sub_set_female
                    }
                    
                    # diff exp
                    print(data_sub_set)
                    file_to_load_up = sprintf("%s/%s_%s/results_%s/matrices/aggregation_diff_exp/diff_exp_sig_up_%s_at_least_in_%s_%s_comps_for_%s_%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_sub_set, data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue)
                    file_to_load_down = sprintf("%s/%s_%s/results_%s/matrices/aggregation_diff_exp/diff_exp_sig_down_%s_at_least_in_%s_%s_comps_for_%s_%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_sub_set, data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue)
                    
                    diff_exp_list_up = fread(file_to_load_up)
                    # remove columns containing only 0 
                    tmp = colSums(diff_exp_list_up != 0) > 0
                    diff_exp_list_up = diff_exp_list_up[, ..tmp]
                    diff_exp_list_down = fread(file_to_load_down)
                    tmp = colSums(diff_exp_list_down != 0) > 0
                    diff_exp_list_down = diff_exp_list_down[, ..tmp]
  
                    diff_exp_list_up_rownames = diff_exp_list_up$V1
                    diff_exp_list_up$V1 = c()
                    #tmp = (rowSums(diff_exp_list_up) == 1)
                    #diff_exp_list_up_tissue_spec = diff_exp_list_up[tmp,]
                    #diff_exp_list_up_rownames_tissue_spec = diff_exp_list_up_rownames[tmp]
                    #tmp = colSums(diff_exp_list_up_tissue_spec != 0) > 0
                    #diff_exp_list_up_tissue_spec = diff_exp_list_up_tissue_spec[, ..tmp]
                    
                    diff_exp_list_down_rownames = diff_exp_list_down$V1
                    diff_exp_list_down$V1 = c()
                    #tmp = (rowSums(diff_exp_list_down) == 1)
                    #diff_exp_list_down_tissue_spec = diff_exp_list_down[tmp,]
                    #diff_exp_list_down_rownames_tissue_spec = diff_exp_list_down_rownames[tmp]
                    #tmp = colSums(diff_exp_list_down_tissue_spec != 0) > 0
                    #diff_exp_list_down_tissue_spec = diff_exp_list_down_tissue_spec[, ..tmp]
                    
                    
                    # corr
                    file_to_load_pos = sprintf("%s/%s_%s/results_%s/matrices/aggregation_corr/corr_method=%s_sig_pos_for_%s_%s_%s_corr_th=%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_sub_set, m, num_corr_tissues, tissue, data_input$rna_class, th)
                    file_to_load_neg = sprintf("%s/%s_%s/results_%s/matrices/aggregation_corr/corr_method=%s_sig_neg_for_%s_%s_%s_corr_th=%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_sub_set, m, num_corr_tissues, tissue, data_input$rna_class, th)
                    
                    corr_list_pos = fread(file_to_load_pos)
                    tmp = colSums(corr_list_pos != 0) > 0
                    corr_list_pos = corr_list_pos[, ..tmp]
                    corr_list_neg = fread(file_to_load_neg)
                    tmp = colSums(corr_list_neg != 0) > 0
                    corr_list_neg = corr_list_neg[, ..tmp]
                    
                    corr_list_pos_rownames = corr_list_pos$V1
                    corr_list_pos$V1 = c()
                    #tmp = (rowSums(corr_list_pos) == 1)
                    #corr_list_pos_tissue_spec = corr_list_pos[tmp,]
                    #corr_list_pos_rownames_tissue_spec = corr_list_pos_rownames[tmp]
                    #tmp = colSums(corr_list_pos_tissue_spec != 0) > 0
                    #corr_list_pos_tissue_spec = corr_list_pos_tissue_spec[, ..tmp]
                    
                    corr_list_neg_rownames = corr_list_neg$V1
                    corr_list_neg$V1 = c()
                    #tmp = (rowSums(corr_list_neg) == 1)
                    #corr_list_neg_tissue_spec = corr_list_neg[tmp,]
                    #corr_list_neg_rownames_tissue_spec = corr_list_neg_rownames[tmp]
                    #tmp = colSums(corr_list_neg_tissue_spec != 0) > 0
                    #corr_list_neg_tissue_spec = corr_list_neg_tissue_spec[, ..tmp]
                    
                    
                    # hclust
                    file_to_load = sprintf("%s/%s_%s/results_%s/matrices/%s_clustering/%s_zscore_th=%s_clustered_by_cv_removed_zero_rows_plot_df.csv", results_folder, data_input$rna_class, data_input$detection_rate, sprintf("%s_%s", hclust_props$subset, sex_key), m_hclust_split, top, th_hclust)
                    print(file_to_load)

                    hclust_table = fread(file_to_load)
                    tmp = colSums(hclust_table != 0) > 0
                    hclust_table = hclust_table[, ..tmp]
                    
                    hclust_table_rownames = hclust_table$V1#hclust_table[[data_input$rna_class]]
                    hclust_table$V1 = c() #hclust_table[[data_input$rna_class]] = c()
                    #tmp = (rowSums(hclust_table) == 1)
                    #hclust_table_tissue_spec = hclust_table[tmp]
                    #hclust_table_rownames_tissue_spec = hclust_table_rownames[tmp]
                    #tmp = colSums(hclust_table_tissue_spec != 0) > 0
                    #hclust_table_tissue_spec = hclust_table_tissue_spec[, ..tmp]
                    
                    
                    # venn plot
                    # only if not empty
                    if (dim(diff_exp_list_up)[1] == 0) {
                      diff_exp_list_up = c()
                    }
                    if (dim(diff_exp_list_down)[1] == 0) {
                      diff_exp_list_down = c()
                    }
                    if (dim(corr_list_pos)[1] == 0) {
                      corr_list_pos = c()
                    } 
                    if (dim(corr_list_neg)[1] == 0) {
                      corr_list_neg = c()
                    }
                    if (dim(hclust_table)[1] == 0) {
                      hclust_table = c()
                    }
                    
                    #maximal_value = 0
                    #for (sub_list in list(diff_exp_list_up, diff_exp_list_down, corr_list_pos, corr_list_neg, hclust_table)) {
                    #  # only go into if clause if it is really a matrix
                    #  if (!is.null(nrow(sub_list)) && !is.null(ncol(sub_list))) {
                    #    maximal_value = max(maximal_value, colSums(sub_list))
                    #  }
                    #
                    
                    maximal_value = 0
                    #for (sub_list in list(diff_exp_list_up[rowSums(diff_exp_list_up[,..tmp]) > 0], corr_list_pos[rowSums(corr_list_pos[,..tmp]) > 0], diff_exp_list_down[rowSums(diff_exp_list_down[,..tmp]) > 0], corr_list_neg[rowSums(corr_list_neg[,..tmp]) > 0], hclust_table[rowSums(hclust_table[,..tmp]) > 0])) {
                    for (sub_list in list(diff_exp_list_up[rowSums(diff_exp_list_up) > 0], corr_list_pos[rowSums(corr_list_pos) > 0], diff_exp_list_down[rowSums(diff_exp_list_down) > 0], corr_list_neg[rowSums(corr_list_neg) > 0], hclust_table[rowSums(hclust_table) > 0])) {
                      # only go into if clause if it is really a matrix
                      if (!is.null(nrow(sub_list)) && !is.null(ncol(sub_list))) {
                        maximal_value = max(maximal_value, colSums(sub_list))
                      }
                    }
                    
                    # tissue list
                    tissue_list_up_pos = union(colnames(diff_exp_list_up), colnames(corr_list_pos))
                    tissue_list_down_neg = union(colnames(diff_exp_list_down), colnames(corr_list_neg))
                    tissue_list_complete = unique(c(tissue_list_up_pos, tissue_list_down_neg, colnames(hclust_table)))
                    
                    # tissue specific
                    p_both_tissue_spec_labels = list()
                    p_both_tissue_spec = list()
                    
                    p_both_labels = list()
                    p_both = list()
                    
                    i = 1
                    j = 1
                    
                    diff_exp_list_up_tissue_spec_features = c()
                    corr_list_pos_tissue_spec_features = c()
                    diff_exp_list_down_tissue_spec_features = c()
                    corr_list_neg_tissue_spec_features = c()
                    hclust_table_tissue_spec_features = c()
                    
                    diff_exp_list_up_features = c()
                    corr_list_pos_features = c()
                    diff_exp_list_down_features = c()
                    corr_list_neg_features = c()
                    hclust_table_features = c()
                    
                    diff_exp_list_tissue_spec_features = c()
                    corr_list_tissue_spec_features = c()
                    hclust_table_tissue_spec_features = c()
                    
                    diff_exp_list_features = c()
                    corr_list_features = c()
                    hclust_table_features = c()
                    
                    tissue_specific = c()
                    age_specific = c()
                    
                    tissue_cand = c()
                    age_cand = c()
                    
                    feature_candidates_tissue_spec = c()
                    feature_candidates_tissue_spec_without_hclust = c()
                    
                    feature_candidates = c()
                    feature_candidates_without_hclust = c()
                    
                    for (ti in tissue_list_complete) {
                      print(ti)
                      #tissue_list = hclust_tissueing[hclust_tissueing[[tissue]] == ti, ][[tissue]]
                      
                      tmp = (colnames(diff_exp_list_up) == ti)
                      tmp_neg = !tmp
                      diff_exp_list_up_tissue_spec_features[[ti]] = diff_exp_list_up_rownames[((rowSums(diff_exp_list_up[,..tmp]) > 0) & (rowSums(diff_exp_list_up[,..tmp_neg]) == 0))]
                      diff_exp_list_up_features[[ti]] = diff_exp_list_up_rownames[rowSums(diff_exp_list_up[,..tmp]) > 0]
                      tmp = (colnames(corr_list_pos) == ti)
                      tmp_neg = !tmp
                      corr_list_pos_tissue_spec_features[[ti]] = corr_list_pos_rownames[((rowSums(corr_list_pos[,..tmp]) > 0) & (rowSums(corr_list_pos[,..tmp_neg]) == 0))]
                      corr_list_pos_features[[ti]] = corr_list_pos_rownames[rowSums(corr_list_pos[,..tmp]) > 0]
                      
                      tmp = (colnames(diff_exp_list_down) == ti)
                      tmp_neg = !tmp
                      diff_exp_list_down_tissue_spec_features[[ti]] = diff_exp_list_down_rownames[((rowSums(diff_exp_list_down[,..tmp]) > 0) & (rowSums(diff_exp_list_down[,..tmp_neg]) == 0))]
                      diff_exp_list_down_features[[ti]] = diff_exp_list_down_rownames[rowSums(diff_exp_list_down[,..tmp]) > 0]
                      tmp = (colnames(corr_list_neg) == ti)
                      tmp_neg = !tmp
                      corr_list_neg_tissue_spec_features[[ti]] = corr_list_neg_rownames[((rowSums(corr_list_neg[,..tmp]) > 0) & (rowSums(corr_list_neg[,..tmp_neg]) == 0))]
                      corr_list_neg_features[[ti]] = corr_list_neg_rownames[rowSums(corr_list_neg[,..tmp]) > 0]
                      
                      tmp = (colnames(hclust_table) == ti)
                      tmp_neg = !tmp
                      hclust_table_tissue_spec_features[[ti]] = hclust_table_rownames[((rowSums(hclust_table[,..tmp]) > 0) & (rowSums(hclust_table[,..tmp_neg]) == 0))]
                      hclust_table_features[[ti]] = hclust_table_rownames[rowSums(hclust_table[,..tmp]) > 0]
                      
                      if (length(intersect(diff_exp_list_up_features[[ti]], corr_list_neg_features[[ti]])) > 0) {
                        print(sprintf("!Same %s in sig. up regulated and sig negative correlated in tissue %s for given parameters!", data_input$rna_class, ti))
                        print(paste(intersect(diff_exp_list_up_features[[ti]], corr_list_neg_features[[ti]]), collapse = ", "))
                        #opt = options(show.error.messages = FALSE)
                        #on.exit(options(opt))
                        #stop()
                      } else if (length(intersect(diff_exp_list_down_features[[ti]], corr_list_pos_features[[ti]])) > 0) {
                        print(sprintf("!Same %s in sig. down regulated and sig postive correlated in tissue %s for given parameters!", data_input$rna_class, ti))
                        print(paste(intersect(diff_exp_list_down_features[[ti]], corr_list_pos_features[[ti]]), collapse = ", "))
                        #opt = options(show.error.messages = FALSE)
                        #on.exit(options(opt))
                        #stop()
                      }
                      
                      diff_exp_list_tissue_spec_features[[ti]] = unique(c(diff_exp_list_up_tissue_spec_features[[ti]], diff_exp_list_down_tissue_spec_features[[ti]]))
                      corr_list_tissue_spec_features[[ti]] = unique(c(corr_list_pos_tissue_spec_features[[ti]], corr_list_neg_tissue_spec_features[[ti]]))
                      hclust_table_tissue_spec_features[[ti]] = unique(hclust_table_tissue_spec_features[[ti]])
                      
                      diff_exp_list_features[[ti]] = unique(c(diff_exp_list_up_features[[ti]], diff_exp_list_down_features[[ti]]))
                      corr_list_features[[ti]] = unique(c(corr_list_pos_features[[ti]], corr_list_neg_features[[ti]]))
                      hclust_table_features[[ti]] = unique(hclust_table_features[[ti]])
                      
                      # remove NULL
                      diff_exp_list_tissue_spec_features[[ti]] = replace_null_with_empty_vector(diff_exp_list_tissue_spec_features[[ti]])
                      corr_list_tissue_spec_features[[ti]] = replace_null_with_empty_vector(corr_list_tissue_spec_features[[ti]])
                      hclust_table_tissue_spec_features[[ti]] = replace_null_with_empty_vector(hclust_table_tissue_spec_features[[ti]])
                      tissue_specific[[ti]] = hclust_table_tissue_spec_features[[ti]]
                      
                      diff_exp_list_features[[ti]] = replace_null_with_empty_vector(diff_exp_list_features[[ti]])
                      corr_list_features[[ti]] = replace_null_with_empty_vector(corr_list_features[[ti]])
                      hclust_table_features[[ti]] = replace_null_with_empty_vector(hclust_table_features[[ti]])
                      tissue_cand[[ti]] = hclust_table_features[[ti]]
                      
                      # if all three tissue_spec are empty
                      if ((length(diff_exp_list_tissue_spec_features[[ti]]) == 0) && (length(corr_list_tissue_spec_features[[ti]]) == 0) && (length(hclust_table_tissue_spec_features[[ti]]) == 0)) {
                        print(sprintf("no tissue specific sig corr, dregulated or hclust noticeable %ss", data_input$rna_class))
                        #next
                      } else {
                        age_specific[[ti]] = unique(c(diff_exp_list_tissue_spec_features[[ti]], corr_list_tissue_spec_features[[ti]]))
                        
                        #file_name = sprintf("candidate_comp_tissue_specific_tissue=%s_sig_diff_exp_both_%s_at_least_in_%s_%s_comps_for_%s_%s_sig_corr_method=%s_both_for_%s_%s_%s_corr_th=%s_hclust_method=%s_th=%s_top=%s_%ss", ti, data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th, m_hclust, th_hclust, top, data_input$rna_class)
                        #file_path_venn_plot = sprintf("%s/%s", output_folder_venn, file_name)
                        #file_path_venn_tab = sprintf("%s/%s", output_folder_venn_tab, file_name)
          
                        # first, we run the plot function to fill without set name for the variable p_output
                        # then, we plot it again with show_set_name = TRUE to save the right plot
                        # and fill the variable p_output_with_set_name
                        #p_output = create_venn_plot(age_specific[[ti]], tissue_specific[[ti]], maximal_value, ti, colours, plot_props, file_path_venn_plot, file_path_venn_tab, show_set_name = TRUE)
                        #p_output_with_set_name = create_venn_plot(age_specific[[ti]], tissue_specific[[ti]], maximal_value, ti, colours, plot_props, file_path_venn_plot, file_path_venn_tab, show_set_name = TRUE)
          
                        #p_both_tissue_spec_labels[[i]] = p_output$p_labels + labs(title = sprintf("%s", ti))
                        #p_both_tissue_spec[[i]] = p_output$p + labs(title = sprintf("%s", ti))
                        
                        # write the tissue names on the left side (vertical)
                        #if (i == 1) {
                        #  # first plot with set names
                        #  p_both_labels[[i]] = tissue_name_side(p_output_with_set_name$p_labels, ti)
                        #  p_both[[i]] = tissue_name_side(p_output_with_set_name$p, ti)
                        #} else {
                        # the rest without
                        #p_both_labels[[i]] = tissue_name_side(p_output$p_labels, tissue, ti, xticks_names)
                        #p_both[[i]] = tissue_name_side(p_output$p, tissue, ti, xticks_names)
                        #}
                        # on top (landscape)
                        # never print set name
                        #p_both_landscape_labels[[i]] = tissue_name_top(p_output$p_labels, tissue, ti, xticks_names) + reduce_landscape_spacing()
                        #p_both_landscape[[i]] = tissue_name_top(p_output$p, tissue, ti, xticks_names) + reduce_landscape_spacing()
                        
                        # euler plot
                        #file_name_euler = sprintf("%s/sig_tissue_specific_comparison_diff_exp_%s_%s_at_least_in_%s_%s_comps_for_%s_%s_corr_method=%s_%s_for_%s_%s_%s_corr_th=%s", output_folder_euler, diff_exp_dir, data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, corr_dir, num_corr_tissues, tissue, data_input$rna_class, th)
                        #euler_plot(diff_exp_tissue_spec_features, corr_tissue_spec_features, plot_props, file_name_euler)
                        
                        feature_candidates_tissue_spec[[ti]] = unique(c(age_specific[[ti]], hclust_table_tissue_spec_features[[ti]]))
                        feature_candidates_tissue_spec_without_hclust[[ti]] = unique(c(age_specific[[ti]]))
                        
                        i = i + 1
                      }
                      
                      if ((length(diff_exp_list_features[[ti]]) == 0) && (length(corr_list_features[[ti]]) == 0) && (length(hclust_table_features[[ti]]) == 0)) {
                        print(sprintf("no sig corr, dregulated or hclust noticeable %ss", data_input$rna_class))
                        #next
                      } else {
                        age_cand[[ti]] = unique(c(diff_exp_list_features[[ti]], corr_list_features[[ti]]))
                        
                        #file_name = sprintf("candidate_comp_tissue=%s_sig_diff_exp_both_%s_at_least_in_%s_%s_comps_for_%s_%s_sig_corr_method=%s_both_for_%s_%s_%s_corr_th=%s_hclust_method=%s_th=%s_top=%s_%ss", ti, data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th, m_hclust, th_hclust, top, data_input$rna_class)
                        #file_path_venn_plot = sprintf("%s/%s", output_folder_venn, file_name)
                        #file_path_venn_tab = sprintf("%s/%s", output_folder_venn_tab, file_name)
                        
                        # first, we run the plot function to fill without set name for the variable p_output
                        # then, we plot it again with show_set_name = TRUE to save the right plot
                        # and fill the variable p_output_with_set_name
                        #p_output = create_venn_plot(age_cand[[ti]], tissue_cand[[ti]], maximal_value, ti, colours, plot_props, file_path_venn_plot, file_path_venn_tab, show_set_name = TRUE)
                        #p_output_with_set_name = create_venn_plot(age_cand[[ti]], tissue_cand[[ti]], maximal_value, ti, colours, plot_props, file_path_venn_plot, file_path_venn_tab, show_set_name = TRUE)
                        
                        #p_both_labels[[j]] = p_output$p_labels + labs(title = sprintf("%s", ti))
                        #p_both[[j]] = p_output$p + labs(title = sprintf("%s", ti))
                        
                        # write the tissue names on the left side (vertical)
                        #if (i == 1) {
                        #  # first plot with set names
                        #  p_both_labels[[i]] = tissue_name_side(p_output_with_set_name$p_labels, ti)
                        #  p_both[[i]] = tissue_name_side(p_output_with_set_name$p, ti)
                        #} else {
                        # the rest without
                        #p_both_labels[[i]] = tissue_name_side(p_output$p_labels, tissue, ti, xticks_names)
                        #p_both[[i]] = tissue_name_side(p_output$p, tissue, ti, xticks_names)
                        #}
                        # on top (landscape)
                        # never print set name
                        #p_both_landscape_labels[[i]] = tissue_name_top(p_output$p_labels, tissue, ti, xticks_names) + reduce_landscape_spacing()
                        #p_both_landscape[[i]] = tissue_name_top(p_output$p, tissue, ti, xticks_names) + reduce_landscape_spacing()
                        
                        # euler plot
                        #file_name_euler = sprintf("%s/sig_comparison_diff_exp_%s_%s_at_least_in_%s_%s_comps_for_%s_%s_corr_method=%s_%s_for_%s_%s_%s_corr_th=%s", output_folder_euler, diff_exp_dir, data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, corr_dir, num_corr_tissues, tissue, data_input$rna_class, th)
                        #euler_plot(diff_exp_features, corr_features, plot_props, file_name_euler)
                        
                        feature_candidates[[ti]] = unique(c(age_cand[[ti]], hclust_table_features[[ti]]))
                        feature_candidates_without_hclust[[ti]] = unique(c(age_cand[[ti]]))
                            
                        j = j + 1
                      }
                    }
                    
                    # merged plots
                    #file_name_venn_merged = sprintf("candidate_comp_tissue_specific_sig_diff_exp_both_%s_at_least_in_%s_%s_comps_for_%s_%s_sig_corr_method=%s_both_for_%s_%s_%s_corr_th=%s_hclust_method=%s_th=%s_top=%s_%ss", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th, m_hclust, th_hclust, top, data_input$rna_class)
                    #file_path_venn_merged = sprintf("%s/%s", output_folder_venn, file_name_venn_merged)
                    
                    # save plots
                    #all_p_tissue_spec_labels = arrangeGrob(grobs=p_both_tissue_spec_labels, ncol=plot_layout$ncol, nrow=plot_layout$nrow, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
                    #all_p_tissue_spec = arrangeGrob(grobs=p_both_tissue_spec, ncol=plot_layout$ncol, nrow=plot_layout$nrow, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
                    
                    #ggsave(sprintf("%s_labels_%sx%s_9x3.png", file_path_venn_merged, plot_layout$nrow, plot_layout$ncol), all_p_tissue_spec_labels, dpi=plot_props$dpi, width=plot_props$image_width, height=plot_props$image_height/2, units = plot_props$image_units)
                    #ggsave(sprintf("%s_labels_%sx%s_9x3.svg", file_path_venn_merged, plot_layout$nrow, plot_layout$ncol), all_p_tissue_spec_labels, width=plot_props$image_width, height=plot_props$image_height/2, units =plot_props$image_units, limitsize = FALSE, scale = 1)
                    
                    #ggsave(sprintf("%s_%sx%s_9x3.png", file_path_venn_merged, plot_layout$nrow, plot_layout$ncol), all_p_tissue_spec, dpi=plot_props$dpi, width=plot_props$image_width, height=plot_props$image_height/2, units = plot_props$image_units)
                    #ggsave(sprintf("%s_%sx%s_9x3.svg", file_path_venn_merged, plot_layout$nrow, plot_layout$ncol), all_p_tissue_spec, width=plot_props$image_width, height=plot_props$image_height/2, units =plot_props$image_units, limitsize = FALSE, scale = 1)
                    
                    # merged plots
                    #file_name_venn_merged = sprintf("candidate_comp_sig_diff_exp_both_%s_at_least_in_%s_%s_comps_for_%s_%s_sig_corr_method=%s_both_for_%s_%s_%s_corr_th=%s_hclust_method=%s_th=%s_top=%s_%ss", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th, m_hclust, th_hclust, top, data_input$rna_class)
                    #file_path_venn_merged = sprintf("%s/%s", output_folder_venn, file_name_venn_merged)
                    
                    # save plots
                    #all_p_labels = arrangeGrob(grobs=p_both_labels, ncol=plot_layout$ncol, nrow=plot_layout$nrow, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
                    #all_p = arrangeGrob(grobs=p_both, ncol=plot_layout$ncol, nrow=plot_layout$nrow, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
                    
                    #ggsave(sprintf("%s_labels_%sx%s_9x3.png", file_path_venn_merged, plot_layout$nrow, plot_layout$ncol), all_p_labels, dpi=plot_props$dpi, width=plot_props$image_width, height=plot_props$image_height/2, units = plot_props$image_units)
                    #ggsave(sprintf("%s_labels_%sx%s_9x3.svg", file_path_venn_merged, plot_layout$nrow, plot_layout$ncol), all_p_labels, width=plot_props$image_width, height=plot_props$image_height/2, units =plot_props$image_units, limitsize = FALSE, scale = 1)
                    
                    #ggsave(sprintf("%s_%sx%s_9x3.png", file_path_venn_merged, plot_layout$nrow, plot_layout$ncol), all_p, dpi=plot_props$dpi, width=plot_props$image_width, height=plot_props$image_height/2, units = plot_props$image_units)
                    #ggsave(sprintf("%s_%sx%s_9x3.svg", file_path_venn_merged, plot_layout$nrow, plot_layout$ncol), all_p, width=plot_props$image_width, height=plot_props$image_height/2, units =plot_props$image_units, limitsize = FALSE, scale = 1)
                    
                    #print("venn done")
                    
                    
                    # upset plot
                    #file_name_upset = sprintf("candidate_comp_tissue_specific_sig_diff_exp_both_%s_at_least_in_%s_%s_comps_for_%s_%s_sig_corr_method=%s_both_for_%s_%s_%s_corr_th=%s_hclust_method=%s_th=%s_top=%s_%ss", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th, m_hclust, th_hclust, top, data_input$rna_class)
                    #file_name_upset = sprintf("%s/%s", output_folder_upset, file_name_upset)
                    #create_upset_plot(feature_candidates_tissue_spec, tissue, xticks_names, colours, plot_props, file_name_upset)
                    
                    #file_name_upset = sprintf("candidate_comp_sig_diff_exp_both_%s_at_least_in_%s_%s_comps_for_%s_%s_sig_corr_method=%s_both_for_%s_%s_%s_corr_th=%s_hclust_method=%s_th=%s_top=%s_%ss", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th, m_hclust, th_hclust, top, data_input$rna_class)
                    #file_name_upset = sprintf("%s/%s", output_folder_upset, file_name_upset)
                    #create_upset_plot(feature_candidates, tissue, xticks_names, colours, plot_props, file_name_upset)
                  
                    #file_name_upset = sprintf("candidate_comp_tissue_specific_sig_diff_exp_both_%s_at_least_in_%s_%s_comps_for_%s_%s_sig_corr_method=%s_both_for_%s_%s_%s_corr_th=%s", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th, m_hclust, th_hclust, top, data_input$rna_class)
                    #file_name_upset = sprintf("%s/%s", output_folder_upset, file_name_upset)
                    #create_upset_plot(feature_candidates_tissue_spec_without_hclust, tissue, xticks_names, colours, plot_props, file_name_upset)
                    
                    #file_name_upset = sprintf("candidate_comp_sig_diff_exp_both_%s_at_least_in_%s_%s_comps_for_%s_%s_sig_corr_method=%s_both_for_%s_%s_%s_corr_th=%s", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th, m_hclust, th_hclust, top, data_input$rna_class)
                    #file_name_upset = sprintf("%s/%s", output_folder_upset, file_name_upset)
                    #create_upset_plot(feature_candidates_without_hclust, tissue, xticks_names, colours, plot_props, file_name_upset)
                    
                    #print("upset done")
                    
                    
                    # tissue specific
                    overlapping_features = c()
                    num_tissue_spec_features = c()
                    for (ti in tissue_list_complete) {
                        all_lists = list(diff_exp_list_up_tissue_spec_features[[ti]], diff_exp_list_down_tissue_spec_features[[ti]], corr_list_pos_tissue_spec_features[[ti]], corr_list_neg_tissue_spec_features[[ti]], hclust_table_tissue_spec_features[[ti]])
                        overlapping_features[[ti]] = Reduce(intersect, all_lists)
                        
                        # overlap_inlcuded
                        #num_tissue_spec_features[[ti]] = c(length(diff_exp_list_up_tissue_spec_features[[ti]]), length(diff_exp_list_down_tissue_spec_features[[ti]]), length(corr_list_pos_tissue_spec_features[[ti]]), length(corr_list_neg_tissue_spec_features[[ti]]), length(hclust_table_tissue_spec_features[[ti]]), length(overlapping_features[[ti]]))
                      
                        # overlap not included
                        num_tissue_spec_features[[ti]] = c(length(diff_exp_list_up_tissue_spec_features[[ti]]), length(diff_exp_list_down_tissue_spec_features[[ti]]), length(corr_list_pos_tissue_spec_features[[ti]]), length(corr_list_neg_tissue_spec_features[[ti]]), length(hclust_table_tissue_spec_features[[ti]]))
                    }
                    
                    num_tissue_spec_features_df[[sex_key]] = as.data.frame(do.call(rbind, num_tissue_spec_features))
                    
                    # overlap_inlcuded
                    #colnames(num_tissue_spec_features_df[[sex_key]]) = c("Up reg.", "Down reg.", "Pos. corr.", "Neg. corr.", "Region spec.", "Overlap")
                    #suffix = "_overlap"
                    # overlap not included
                    colnames(num_tissue_spec_features_df[[sex_key]]) = c("Up reg.", "Down reg.", "Pos. corr.", "Neg. corr.", "Region spec.")
                    suffix = ""
                    
                    # all candidates
                    overlapping_features = c()
                    num_features = c()
                    for (ti in tissue_list_complete) {
                      all_lists = list(diff_exp_list_up_features[[ti]], diff_exp_list_down_features[[ti]], corr_list_pos_features[[ti]], corr_list_neg_features[[ti]], hclust_table_features[[ti]])
                      overlapping_features[[ti]] = Reduce(intersect, all_lists)
                      
                      # overlap_inlcuded
                      #num_features[[ti]] = c(length(diff_exp_list_up_features[[ti]]), length(diff_exp_list_down_features[[ti]]), length(corr_list_pos_features[[ti]]), length(corr_list_neg_features[[ti]]), length(hclust_table_features[[ti]]), length(overlapping_features[[ti]]))
                      
                      # overlap not included
                      num_features[[ti]] = c(length(diff_exp_list_up_features[[ti]]), length(diff_exp_list_down_features[[ti]]), length(corr_list_pos_features[[ti]]), length(corr_list_neg_features[[ti]]), length(hclust_table_features[[ti]]))
                    }
                    
                    num_features_df[[sex_key]] = as.data.frame(do.call(rbind, num_features))
                    
                    # overlap_inlcuded
                    #colnames(num_features_df[[sex_key]]) = c("Up reg.", "Down reg.", "Pos. corr.", "Neg. corr.", "Region spec.", "Overlap")
                    #suffix = "_overlap"
                    # overlap not included
                    colnames(num_features_df[[sex_key]]) = c("Up reg.", "Down reg.", "Pos. corr.", "Neg. corr.", "Region spec.")
                    suffix = ""
                  }
                  
                  maxi_tissue_spec = c()
                  maxi_all = c()
                  for (sex_key in c("male", "female")) {
                    maxi_tissue_spec[[sex_key]] = max(max(num_tissue_spec_features_df[[sex_key]]))
                    maxi_all[[sex_key]] = max(max(num_features_df[[sex_key]])) 
                  }
                  maxi_male_female_tissue_spec = max(unname(unlist(maxi_tissue_spec)))
                  maxi_male_female_all = max(unname(unlist(maxi_all)))
                    
                  scatter_plot_tissue_spec = list()
                  scatter_plot_tissue_spec_log10 = list()
                  scatter_plot_all_cands = list()
                  scatter_plot_all_cands_log10 = list()
                  for (sex_key in c("male", "female")) {
                    print(sex_key)
                    # tissue specific
                    # raw number
                    x_axis_label = ""
                    log_toggle = FALSE
                    
                    p_tissue_spec = create_scatter_plot(num_tissue_spec_features_df[[sex_key]], tissue, x_axis_label, log_toggle, maxi_male_female_tissue_spec, colours, plot_props)
                    
                    if (sex_key == "male") {
                      scatter_plot_tissue_spec[[sex_key]] = p_tissue_spec$flip
                    } else {
                      scatter_plot_tissue_spec[[sex_key]] = p_tissue_spec$normal
                    }
                    
                    log_toggle = TRUE
                    num_tissue_spec_features_df_log10 = log10(num_tissue_spec_features_df[[sex_key]] + 1)
                    num_tissue_spec_features_df_log10[num_tissue_spec_features_df_log10 == -Inf] = 0
                    
                    p_tissue_spec_log10 = create_scatter_plot(num_tissue_spec_features_df_log10, tissue, x_axis_label, log_toggle, log10(maxi_male_female_tissue_spec), colours, plot_props)
                    
                    if (sex_key == "male") {
                      scatter_plot_tissue_spec_log10[[sex_key]] = p_tissue_spec_log10$flip
                    } else {
                      scatter_plot_tissue_spec_log10[[sex_key]] = p_tissue_spec_log10$normal
                    }
                    
                    # all candidates
                    # raw number
                    x_axis_label = ""
                    log_toggle = FALSE
                    p_all_cands = create_scatter_plot(num_features_df[[sex_key]], tissue, x_axis_label, log_toggle, maxi_male_female_all, colours, plot_props)
                    
                    if (sex_key == "male") {
                      scatter_plot_all_cands[[sex_key]] = p_all_cands$flip
                    } else {
                      scatter_plot_all_cands[[sex_key]] = p_all_cands$normal
                    }
                                      
                    # number in log10
                    #x_axis_label = sprintf("Num. of %ss  (log10)", data_input$rna_class)
                    
                    log_toggle = TRUE
                    num_features_df_log10 = log10(num_features_df[[sex_key]] + 1)
                    num_features_df_log10[num_features_df_log10 == -Inf] = 0
                    
                    p_all_cands_log10 = create_scatter_plot(num_features_df_log10, tissue, x_axis_label, log_toggle, log10(maxi_male_female_all), colours, plot_props)
                    
                    if (sex_key == "male") {
                      scatter_plot_all_cands_log10[[sex_key]] = p_all_cands_log10$flip
                    } else {
                      scatter_plot_all_cands_log10[[sex_key]] = p_all_cands_log10$normal
                    }
                  }

                  
                  # save all plots separately
                  for (sex_key in c("male", "female")) {
                    width = plot_props$image_width
                    height = plot_props$image_height
                    
                    # tissue specific
                    file_name = sprintf("candidate_comp_tissue_specific_sig_diff_exp_both_%s_at_least_in_%s_%s_comps_for_%s_%s_sig_corr_method=%s_both_for_%s_%s_%s_corr_th=%s_hclust_method=%s_th=%s_top=%s_%ss%s", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th, m_hclust, th_hclust, top, data_input$rna_class, suffix)
                    file_name_log10 = sprintf("candidate_comp_tissue_specific_sig_diff_exp_both_%s_at_least_in_%s_%s_comps_for_%s_%s_sig_corr_method=%s_both_for_%s_%s_%s_corr_th=%s_hclust_method=%s_th=%s_top=%s_%ss%s_log10", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th, m_hclust, th_hclust, top, data_input$rna_class, suffix)
                    
                    # raw
                    ggsave(sprintf("%s/%s_%s.png", output_folder_scatter, file_name, sex_key), scatter_plot_tissue_spec[[sex_key]], dpi=plot_props$dpi, width=width, height=height, units = plot_props$image_units)
                    ggsave(sprintf("%s/%s_%s.svg", output_folder_scatter, file_name, sex_key), scatter_plot_tissue_spec[[sex_key]], width=width, height=height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
                    # log10
                    ggsave(sprintf("%s/%s_%s.png", output_folder_scatter, file_name_log10, sex_key), scatter_plot_tissue_spec_log10[[sex_key]], dpi=plot_props$dpi, width=width, height=height, units = plot_props$image_units)
                    ggsave(sprintf("%s/%s_%s.svg", output_folder_scatter, file_name_log10, sex_key), scatter_plot_tissue_spec_log10[[sex_key]], width=width, height=height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
                    
                    
                    # all candidates
                    file_name = sprintf("candidate_comp_sig_diff_exp_both_%s_at_least_in_%s_%s_comps_for_%s_%s_sig_corr_method=%s_both_for_%s_%s_%s_corr_th=%s_hclust_method=%s_th=%s_top=%s_%ss%s", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th, m_hclust, th_hclust, top, data_input$rna_class, suffix)
                    file_name_log10 = sprintf("candidate_comp_sig_diff_exp_both_%s_at_least_in_%s_%s_comps_for_%s_%s_sig_corr_method=%s_both_for_%s_%s_%s_corr_th=%s_hclust_method=%s_th=%s_top=%s_%ss%s_log10", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th, m_hclust, th_hclust, top, data_input$rna_class, suffix)
                    
                    # raw
                    ggsave(sprintf("%s/%s_%s.png", output_folder_scatter, file_name, sex_key), scatter_plot_all_cands[[sex_key]], dpi=plot_props$dpi, width=width, height=height, units = plot_props$image_units)
                    ggsave(sprintf("%s/%s_%s.svg", output_folder_scatter, file_name, sex_key), scatter_plot_all_cands[[sex_key]], width=width, height=height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
                    # log10
                    ggsave(sprintf("%s/%s_%s.png", output_folder_scatter, file_name_log10, sex_key), scatter_plot_all_cands_log10[[sex_key]], dpi=plot_props$dpi, width=width, height=height, units = plot_props$image_units)
                    ggsave(sprintf("%s/%s_%s.svg", output_folder_scatter, file_name_log10, sex_key), scatter_plot_all_cands_log10[[sex_key]], width=width, height=height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
                    
                  }
                  
                  # cleanup
                  rm(scatter_plot_tissue_spec)
                  rm(scatter_plot_tissue_spec_log10)
                  
                  rm(scatter_plot_all_cands)
                  rm(scatter_plot_all_cands_log10)
                  
                  
                }
                # cleanup
                gc()
              }
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

