suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))

# save rdata
#saveRDS(snakemake, file = "snakemale_files/snakemake_external_protein_plot.rds")
#stop()

#snakemake = readRDS("snakemale_files/snakemake_external_protein_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("external_protein_plot")


# ---------------------- Load data ----------------------
raw_data_path = snakemake@params$raw_data_path
age_colour = snakemake@params$age_colour
intensity_col = snakemake@params$intensity_col
plot_parameters = snakemake@params$plot_parameters
xticks_names = snakemake@params$xticks_names
plots_props = snakemake@params$plots_props
colours = snakemake@params$colors
results_folder = snakemake@params$results_folder

output_folder = sprintf("%s/figures/box_plots", results_folder)
dir.create(output_folder, recursive=TRUE)


#------------------------------------ Script ----------------------------------- 
tbl = fread(raw_data_path, sep='\t', header=T)

for (i in 1:length(plot_parameters)) {
  plot_param = plot_parameters[[i]]
  
  print(sprintf("Gene: %s", plot_param$gene))
  print(sprintf("Tissue: %s", plot_param$tissue))
  
  split_prop = plot_param$split_prop
  plot_prop = plot_param$plot_prop
  plot_titles = plot_param$plot_titles
    
  tbl_ti_gene = tbl[((tbl$tissue == plot_param$tissue) & (tbl$symbol == plot_param$gene)),]
  
  p_box = list()
  fold_change = c()
  for (split_p in unique(tbl_ti_gene[[split_prop]])) {
    tbl_filtered = tbl_ti_gene[tbl_ti_gene[[split_prop]] == split_p,]
    
    tmp = colnames(tbl_filtered) %in% c(plot_prop, intensity_col)
    plot_df = tbl_filtered[,..tmp]
    
    colnames(plot_df) = c("value", "time")
    
    p_box[[split_p]] = ggplot(plot_df, aes(x = time, y = value, group = time, fill = as.factor(time))) + #plot_df%>%group_by(params$plot_group)%>%summarise(value=sum(value))%>%
      geom_boxplot(aes(fill = as.factor(time)), width = 1.75, color = "black", outlier.colour= NA, outlier.shape=16, outlier.size=1, notch=FALSE, fatten=TRUE) + 
      geom_boxplot(fill=NA, color="black", width = 1.75, fatten=FALSE, coef = 0, outlier.shape = NA, outlier.alpha = 0) +
      #geom_jitter(data = plot_df, aes(y = value, x = time), color = "black", shape = 19, size = 0.75, alpha = 0.75) +
      geom_boxplot(width = 1.75, outlier.colour="grey", outlier.shape=16, outlier.size=1, notch=FALSE, outlier.alpha = 1) +
      geom_point(data = plot_df, aes(y = value, x = time), color = "black", shape = 1, size = 2) +
      scale_fill_manual(values = age_colour) +
      ggtitle(split_p) +
      xlab("") +
      ylab("") +
      #ylab(sprintf("expr (%s, detection rate = %s%%\nper %s)", data_input$norm, expr_input$detection_rate, params$plot_group)) +
      theme_classic() +
      theme(
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position="none", 
        text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
        plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size, margin = margin(t = unit(0, "pt"), r = unit(5, "pt"), b = unit(2, "pt"), l = unit(5, "pt"))),
        legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
        #panel.ontop = TRUE,  # Ensure points are drawn on top
        #panel.border = element_rect(fill = NA) # Make border visible
      ) +
    coord_cartesian(clip = "off")  # Prevents clipping
    
    #hight = 3
    #width = 3
    #ggsave(sprintf("%s/protein_intensity_plot_per_%s_over_%s_%s_%s_%sx%s.png", output_folder, split_p, plot_prop, plot_param$gene, plot_param$tissue, height, width), p_box[[split_p]], dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
    #ggsave(sprintf("%s/protein_intensity_plot_per_%s_over_%s_%s_%s_%sx%s.svg", output_folder, split_p, plot_prop, plot_param$gene, plot_param$tissue, height, width), p_box[[split_p]], width=width, height=height, units = plots_props$image_units, limitsize = FALSE, scale = 1)
  
    # fold change calculation
    plot_df$value_log2 = log2(plot_df$value + 1)
    geom_median_list = c()
    for (tp in sort(unique(plot_df$time))) {
      geom_median_list = append(geom_median_list, rep(median(plot_df[plot_df$time == tp,]$value), dim(plot_df[plot_df$time == tp,])[1]))
    }
    plot_df$geom_median = geom_median_list
    fold_change[[split_p]] = 2^(plot_df[plot_df$time == sort(unique(plot_df$time))[2],][1]$geom_median - plot_df[plot_df$time == sort(unique(plot_df$time))[1],][1]$geom_median)
  }
  
  # sort the figures given the order from the config
  p_box = p_box[c("Male", "Female")]
  
  plot_nrows = 1
  plot_ncols = 2
  
  # determine which figures are not on the left side
  left_col = c()
  for (row in 1:plot_nrows - 1) {
    left_col = append(left_col, plot_ncols*row + 1)
  }
  # determine which figures are not on the left side
  right_col = c()
  for (row in 1:plot_nrows) {
    right_col = append(right_col, plot_ncols*row)
  }
  
  p_box_mod_margins = p_box
  for (img_i in 1:length(p_box_mod_margins)) {
    if ((!(img_i %in% left_col)) & (!(img_i %in% right_col))) {
      p_box_mod_margins[[img_i]] = p_box_mod_margins[[img_i]] + 
        theme(plot.margin = margin(-2, -1, -10, -13, "pt")
        )
    }
    else if (img_i %in% left_col) {
      p_box_mod_margins[[img_i]] = p_box_mod_margins[[img_i]] + 
        theme(plot.margin = margin(-2, -4, -10, -10, "pt")
        )
    } else if (img_i %in% right_col) {
      p_box_mod_margins[[img_i]] = p_box_mod_margins[[img_i]] + 
        theme(plot.margin = margin(-2, 4, -10, -18, "pt")
        )
    }
    if (!(img_i %in% left_col)) {
      p_box_mod_margins[[img_i]] = p_box_mod_margins[[img_i]] + 
        #scale_y_continuous(labels = function(breaks) {rep_along(breaks, "  ")}, limits = c(0, 1.15*maxi_expr_rna), expand = c(0,0)) +
        theme(axis.ticks.y = element_line(color = "transparent"),
              axis.text.y = element_text(color = "transparent")
        )
    }
  }
  
  y_axis_title = "Intensity"
  x_axis_title = xticks_names$categories[[sprintf("%s_capital", plot_prop)]]
  plot_title = plot_titles[1]
  y_axis_right_title = plot_titles[2]
  all_p_bar_manual = arrangeGrob(grobs=p_box_mod_margins, ncol=plot_ncols, nrow=plot_nrows, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), top = textGrob(plot_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90), right = textGrob(y_axis_right_title, gp=gpar(fontsize=plots_props$font_size_header, fontfamily=plots_props$font_family), rot = 270))
  
  width = 8
  height = 4
  ggsave(sprintf("%s/protein_intensity_plot_per_%s_over_%s_%s_%s_%sx%s.png", output_folder, split_prop, plot_prop, plot_param$gene, plot_param$tissue, height, width), all_p_bar_manual, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
  ggsave(sprintf("%s/protein_intensity_plot_per_%s_over_%s_%s_%s_%sx%s.svg", output_folder, split_prop, plot_prop, plot_param$gene, plot_param$tissue, height, width), all_p_bar_manual, width=width, height=height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
  
  print(fold_change)
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
