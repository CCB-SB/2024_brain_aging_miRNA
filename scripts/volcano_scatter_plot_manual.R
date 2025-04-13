suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))

# save rdata
#saveRDS(snakemake, file = "snakemake_volcano_scatter_plot_manual.rds")
#stop()

#snakemake = readRDS("snakemake_volcano_scatter_plot_manual.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("volcano_scatter_plot")


#---------------------------------- Functions ---------------------------------- 
plot_fix_layout = function(plot_mod_margins, max_value, top_margin_left_col, right_margin_left_col, bottom_margin_left_col, left_margin_left_col,
                           top_margin_middle_cols, right_margin_middle_cols, bottom_margin_middle_cols, left_margin_middle_cols,
                           top_margin_right_col, right_margin_right_col, bottom_margin_right_col, left_margin_right_col){
  for (img_i in 1:length(plot_mod_margins)) {
    if ((!(img_i %in% left_col)) & (!(img_i %in% right_col))) {
      plot_mod_margins[[img_i]] = plot_mod_margins[[img_i]] + 
        theme(plot.margin = margin(top_margin_middle_cols, right_margin_middle_cols, bottom_margin_middle_cols, left_margin_middle_cols, "pt")
        )
    }
    else if (img_i %in% left_col) {
      plot_mod_margins[[img_i]] = plot_mod_margins[[img_i]] + 
        theme(plot.margin = margin(top_margin_left_col, right_margin_left_col, bottom_margin_left_col, left_margin_left_col, "pt")
        )
    } else if (img_i %in% right_col) {
      plot_mod_margins[[img_i]] = plot_mod_margins[[img_i]] + 
        theme(plot.margin = margin(top_margin_right_col, right_margin_right_col, bottom_margin_right_col, left_margin_right_col, "pt")
        )
    }
    if (!(img_i %in% left_col)) {
      plot_mod_margins[[img_i]] = plot_mod_margins[[img_i]] + 
              ylim(0, max(unlist(max_value))) +
              theme(axis.ticks.y = element_line(color = "transparent"),
              axis.text.y = element_text(color = "transparent")
        )
    }
  }
  return(plot_mod_margins)
}


#---------------------------------- Load data ---------------------------------- 
data_input = snakemake@params$data_input
diff_log = snakemake@params$diff_log
comps = snakemake@params$comparisons
colours = snakemake@params$colors
xticks_names = snakemake@params$xticks_names
plots_props = snakemake@params$plots_props
manual_plots_parameters = snakemake@params$manual_plots_parameters
feature_path = manual_plots_parameters$heighlight_features
results_folder = snakemake@params$results_folder
log_fc_thres = log2(snakemake@params$parameters$updownregulated)
sig_lvl = snakemake@params$parameters$significance
effect_size_thres = snakemake@params$parameters$mineffectsize

tbl = fread(sprintf("%s/%s_%s/results_%s/matrices/diff_exp/diff_exp%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, diff_log), sep='\t') 

repel_max_time =  10 

feature_list_highlighting = unique(fread(feature_path, sep='\t', header = FALSE))

output_folder_merged = sprintf("%s/%s_%s/results_%s/figures/volcano/scatter_plots/merged", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_merged, recursive=TRUE)


# --------------------------------- Script -------------------------------------
s_string = sprintf("Significant (P<%s)", sig_lvl)
ns_string = sprintf("Not significant (P\u2265%s)", sig_lvl)
up_string = sprintf("Sig. upregulated (FC\u2265%s)", snakemake@params$parameters$updownregulated)
down_string = sprintf("Sig. downregulated (FC\u22641/%s)", snakemake@params$parameters$updownregulated)
e_string = sprintf("Considerable effect (d\u2265%s)", effect_size_thres)
ne_string = sprintf("Neglectable effect (d<%s)", effect_size_thres)
e_up_string = sprintf("Considerably upregulated (FC\u2265%s)", snakemake@params$parameters$updownregulated)
e_down_string = sprintf("Considerably downregulated (FC\u22641/%s)", snakemake@params$parameters$updownregulated)

colors = c()
# grey style
#colors[ns_string] = "lightgrey"
#colors[up_string] = "#31363B"
#colors[down_string] = "#31363B"
#colors[s_string] = "lightgrey"
#colors[e_string] = "lightgrey"
#colors[ne_string] = "lightgrey"
#colors[e_up_string] = "#31363B"
#colors[e_down_string] = "#31363B"

colors[ns_string] = "lightgrey"
colors[up_string] = colours$direction$up
colors[down_string] = colours$direction$down
colors[s_string] = "lightgrey" 
colors[e_string] = "lightgrey" 
colors[ne_string] = "lightgrey"
colors[e_up_string] = colours$direction$up
colors[e_down_string] = colours$direction$down

# set font-family for all plots
theme_set(theme_cowplot(font_family=snakemake@params$plots_props$font_family))

groups = list()
group_y_max = list()
group_x_max = list()
group_x_min = list()
group_plots = list()

ttest_raw_max_value = c()
ttest_adj_max_value = c()
wilcox_raw_max_value = c()
wilcox_adj_max_value = c()
cohend_abs_max_value = c()

p1_plots = list()
p2_plots = list()
p3_plots = list()
p4_plots = list()
p5_plots = list()
p6_plots = list()
p7_plots = list()
p8_plots = list()
p9_plots = list()
p10_plots = list()

for(i in 1:length(comps)){
  c = comps[[i]]
  c$g1 = c$g[1]
  c$g2 = c$g[2]
  
  group <- NULL
  if("volcano" %in% names(c)) {
    group <- c$volcano
    if(!(group %in% names(groups))) {
      groups[[group]] <- TRUE
      group_y_max[[group]] <- list("ttest_raw_pval" = -Inf,
                                    "ttest_adj_pval" = -Inf,
                                    "wilcox_raw_pval" = -Inf,
                                    "wilcox_adj_pval" = -Inf,
                                    "cohend_estimate" = -Inf)
      group_x_max[[group]] <- list("logfc" = -Inf)
      group_x_min[[group]] <- list("logfc" = Inf)
      group_plots[[group]] <- list()
    }
  }
  
  prop = names(comps)[i]
 
  paired_string = ifelse(!is.null(c$paired), sprintf("_paired_%s", c$paired), "")
  sub_string = ifelse(!is.null(c$subset), sprintf("_%s=%s", c$subset[1], c$subset[2]), "")
 
  plot_df = data.table(feature=tbl[[data_input$rna_class]],
                       ttest_raw_pval=tbl[[sprintf("ttest_rawp_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]],
                       ttest_adj_pval=tbl[[sprintf("ttest_adjp_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]],
                       wilcox_raw_pval=tbl[[sprintf("wilcox_rawp_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]],
                       wilcox_adj_pval=tbl[[sprintf("wilcox_adjp_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]],
                       fc=tbl[[sprintf("fc_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]],
                       cohend=tbl[[sprintf("cohend_estimate_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]])
  plot_df[,ttest_raw_pval_log10 := -log10(ttest_raw_pval)]
  plot_df[,ttest_adj_pval_log10 := -log10(ttest_adj_pval)]
  plot_df[,wilcox_raw_pval_log10 := -log10(wilcox_raw_pval)]
  plot_df[,wilcox_adj_pval_log10 := -log10(wilcox_adj_pval)]
  plot_df[,logfc := log2(fc)]
  plot_df[,cohend_abs := abs(cohend)]
  
  
  if(!is.null(group)) {
    group_y_max[[group]] <- list("ttest_raw_pval" = max(plot_df[,ttest_raw_pval_log10], group_y_max[[group]][["ttest_raw_pval"]], na.rm=T),
               "ttest_adj_pval" = max(plot_df[,ttest_adj_pval_log10], group_y_max[[group]][["ttest_adj_pval"]], na.rm=T),
               "wilcox_raw_pval" = max(plot_df[,wilcox_raw_pval_log10], group_y_max[[group]][["wilcox_raw_pval"]], na.rm=T),
               "wilcox_adj_pval" = max(plot_df[,wilcox_adj_pval_log10], group_y_max[[group]][["wilcox_adj_pval"]], na.rm=T),
               "cohend_abs" = max(plot_df[,cohend_abs], group_y_max[[group]][["cohend_abs"]], na.rm=T))
    group_x_max[[group]] <- list("logfc" = max(plot_df[,logfc], group_x_max[[group]][["logfc"]], na.rm=T))
    group_x_min[[group]] <- list("logfc" = min(plot_df[,logfc], group_x_min[[group]][["logfc"]], na.rm=T))
  }
  
  plot_df[ttest_raw_pval >= sig_lvl, ttest_raw_cat:=ns_string]
  plot_df[ttest_raw_pval < sig_lvl, ttest_raw_cat:=s_string]
  plot_df[logfc >= log_fc_thres & ttest_raw_pval < sig_lvl, ttest_raw_cat:=up_string]
  plot_df[logfc <= -log_fc_thres & ttest_raw_pval < sig_lvl, ttest_raw_cat:=down_string]
  
  plot_df[ttest_adj_pval >= sig_lvl, ttest_adj_cat:=ns_string]
  plot_df[ttest_adj_pval < sig_lvl, ttest_adj_cat:=s_string]
  plot_df[logfc >= log_fc_thres & ttest_adj_pval < sig_lvl, ttest_adj_cat:=up_string]
  plot_df[logfc <= -log_fc_thres & ttest_adj_pval < sig_lvl, ttest_adj_cat:=down_string]
  
  plot_df[wilcox_raw_pval >= sig_lvl, wilcox_raw_cat:=ns_string]
  plot_df[wilcox_raw_pval < sig_lvl, wilcox_raw_cat:=s_string]
  plot_df[logfc >= log_fc_thres & wilcox_raw_pval < sig_lvl, wilcox_raw_cat:=up_string]
  plot_df[logfc <= -log_fc_thres & wilcox_raw_pval < sig_lvl, wilcox_raw_cat:=down_string]
  
  plot_df[wilcox_adj_pval >= sig_lvl, wilcox_adj_cat:=ns_string]
  plot_df[wilcox_adj_pval < sig_lvl, wilcox_adj_cat:=s_string]
  plot_df[logfc >= log_fc_thres & wilcox_adj_pval < sig_lvl, wilcox_adj_cat:=up_string]
  plot_df[logfc <= -log_fc_thres & wilcox_adj_pval < sig_lvl, wilcox_adj_cat:=down_string]

  plot_df[cohend_abs >= effect_size_thres, cohend_abs_cat:=e_string]
  plot_df[cohend_abs < effect_size_thres, cohend_abs_cat:=ne_string]
  plot_df[logfc >= log_fc_thres & cohend_abs >= effect_size_thres, cohend_abs_cat:=e_up_string]
  plot_df[logfc <= -log_fc_thres & cohend_abs >= effect_size_thres, cohend_abs_cat:=e_down_string]
  
  colors_ttest_raw = colors[names(colors) %in% plot_df$ttest_raw_cat]
  colors_ttest_adj = colors[names(colors) %in% plot_df$ttest_adj_cat]
  colors_wilcox_raw = colors[names(colors) %in% plot_df$wilcox_raw_cat]
  colors_wilcox_adj = colors[names(colors) %in% plot_df$wilcox_adj_cat]
  colors_cohend_estimate = colors[names(colors) %in% plot_df$cohend_abs_cat]
  sub_plot_df = plot_df[!is.na(ttest_raw_pval_log10) & !is.na(logfc) & !is.na(cohend_abs)]

  sub_plot_df$ttest_raw_cat_show = sub_plot_df$feature
  sub_plot_df$ttest_adj_cat_show = sub_plot_df$feature
  sub_plot_df$wilcox_raw_cat_show = sub_plot_df$feature
  sub_plot_df$wilcox_adj_cat_show = sub_plot_df$feature
  sub_plot_df$cohend_abs_cat_show = sub_plot_df$feature
  mask_ttest_raw = ((sub_plot_df$ttest_raw_cat != up_string) & (sub_plot_df$ttest_raw_cat != down_string)) 
  mask_ttest_adj = ((sub_plot_df$ttest_adj_cat != up_string) & (sub_plot_df$ttest_adj_cat != down_string)) 
  mask_wilcox_raw = ((sub_plot_df$wilcox_raw_cat != up_string) & (sub_plot_df$wilcox_raw_cat != down_string)) 
  mask_wilcox_adj = ((sub_plot_df$wilcox_adj_cat != up_string) & (sub_plot_df$wilcox_adj_cat != down_string)) 
  mask_cohend_abs = ((sub_plot_df$cohend_abs_cat != e_up_string) & (sub_plot_df$cohend_abs_cat != e_down_string)) 
  sub_plot_df[mask_ttest_raw,]$ttest_raw_cat_show = ""
  sub_plot_df[mask_ttest_adj,]$ttest_adj_cat_show = ""
  sub_plot_df[mask_wilcox_raw,]$wilcox_raw_cat_show = ""
  sub_plot_df[mask_wilcox_adj,]$wilcox_adj_cat_show = ""
  sub_plot_df[mask_cohend_abs,]$cohend_abs_cat_show = ""
  
  sub_plot_df$ttest_raw_cat_show_manual = sub_plot_df$feature
  sub_plot_df$ttest_adj_cat_show_manual = sub_plot_df$feature
  sub_plot_df$wilcox_raw_cat_show_manual = sub_plot_df$feature
  sub_plot_df$wilcox_adj_cat_show_manual = sub_plot_df$feature
  sub_plot_df$cohend_abs_cat_show_manual = sub_plot_df$feature
  mask = !(sub_plot_df$feature %in% feature_list_highlighting$V1)
  sub_plot_df[mask,]$ttest_raw_cat_show_manual = ""
  sub_plot_df[mask,]$ttest_adj_cat_show_manual = ""
  sub_plot_df[mask,]$wilcox_raw_cat_show_manual = ""
  sub_plot_df[mask,]$wilcox_adj_cat_show_manual = ""
  sub_plot_df[mask,]$cohend_abs_cat_show_manual = ""
  
  x_axis_name_p1 = "Fold change (log2)"
  y_axis_name_p1 = "Raw p-value\nfrom t-test (-log10)"
  p1 = ggplot(sub_plot_df, aes(x=logfc, y=ttest_raw_pval_log10, color=ttest_raw_cat)) + 
    geom_point(size = 0.75) +
    geom_hline(yintercept = -log10(sig_lvl), color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_vline(xintercept = -log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_vline(xintercept = log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) +
    xlab(x_axis_name_p1) + ylab(y_axis_name_p1) +
    scale_color_manual(name="", values=colors_ttest_raw) + 
    theme(legend.position = "none", 
          plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent"),
          text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size, margin = margin(b=0.2, unit="cm")), 
          legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size, face = "bold"),
          plot.margin = margin(b=-0.4, unit = "cm")) + 
    guides(col=guide_legend(nrow=length(colors_ttest_raw)))

  x_axis_name_p2 = "Fold change (log2)"
  y_axis_name_p2 = "Adj. p-value\nfrom t-test (-log10)"
  p2 = ggplot(sub_plot_df, aes(x=logfc, y=ttest_adj_pval_log10, color=ttest_adj_cat)) + 
    geom_point(size = 0.75) +
    geom_hline(yintercept = -log10(sig_lvl), color="lightgrey", linetype="dashed", size = 0.25) +
    geom_vline(xintercept = -log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_vline(xintercept = log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) +
    xlab(x_axis_name_p2) + ylab(y_axis_name_p2) +
    scale_color_manual(name="", values=colors_ttest_adj) + 
    theme(legend.position = "none", 
          plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent"),
          text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size, margin = margin(b=0.2, unit="cm")), 
          legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size, face = "bold"),
          plot.margin = margin(b=-0.4, unit = "cm")) + 
    guides(col=guide_legend(nrow=length(colors_ttest_adj)))

  x_axis_name_p3 = "Fold change (log2)"
  y_axis_name_p3 = "Raw p-value\nfrom Wilcoxon test(-log10)"
  p3 = ggplot(sub_plot_df, aes(x=logfc, y=wilcox_raw_pval_log10, color=wilcox_raw_cat)) + 
    geom_point(size = 0.75) + 
    geom_hline(yintercept = -log10(sig_lvl), color="lightgrey", linetype="dashed", size = 0.25) +
    geom_vline(xintercept = -log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_vline(xintercept = log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) +
    xlab(x_axis_name_p3) + ylab(y_axis_name_p3) +
    scale_color_manual(name="", values=colors_wilcox_raw) + 
    theme(legend.position = "none", 
          plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent"),
          text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size, margin = margin(b=0.2, unit="cm")), 
          legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size, face = "bold"),
          plot.margin = margin(b=-0.4, unit = "cm")) + 
    guides(col=guide_legend(nrow=length(colors_wilcox_raw)))

  x_axis_name_p4 = "Fold change (log2)"
  y_axis_name_p4 = "Adj. p-value\nfrom Wilcoxon test (-log10)"
  p4 = ggplot(sub_plot_df, aes(x=logfc, y=wilcox_adj_pval_log10, color=wilcox_adj_cat)) + 
    geom_point(size = 0.75) + 
    geom_hline(yintercept = -log10(sig_lvl), color="lightgrey", linetype="dashed", size = 0.25) +
    geom_vline(xintercept = -log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_vline(xintercept = log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) +
    xlab(x_axis_name_p4) + ylab(y_axis_name_p4) + 
    scale_color_manual(name="", values=colors_wilcox_adj) + 
    theme(legend.position = "none", 
          plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent"),
          text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size, margin = margin(b=0.2, unit="cm")), 
          legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size, face = "bold"),
          plot.margin = margin(b=-0.4, unit = "cm")) + 
    guides(col=guide_legend(nrow=length(colors_wilcox_adj)))

  x_axis_name_p5 = "Fold change (log2)"
  y_axis_name_p5 = "Effect size (Cohen's d)"
  p5 = ggplot(sub_plot_df, aes(x=logfc, y=cohend_abs, color=cohend_abs_cat)) + 
    geom_point(size = 0.75) + 
    geom_hline(yintercept = abs(effect_size_thres), color="lightgrey", linetype="dashed", size = 0.25) +
    geom_vline(xintercept = -log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_vline(xintercept = log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) +
    xlab(x_axis_name_p5) + ylab(y_axis_name_p5) +
    scale_color_manual(name="", values=colors_cohend_estimate) + 
    theme(legend.position = "none", 
          plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent"),
          text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size, margin = margin(b=0.2, unit="cm")), 
          legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size, face = "bold"),
          plot.margin = margin(b=-0.4, unit = "cm")) + 
    guides(col=guide_legend(nrow=length(colors_cohend_estimate)))
  
  p6 = ggplot(sub_plot_df, aes(x=logfc, y=ttest_raw_pval_log10, color=ttest_raw_cat)) + 
    geom_point(size = 0.75) +
    geom_hline(yintercept = -log10(sig_lvl), color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_vline(xintercept = -log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_vline(xintercept = log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) +
    xlab(x_axis_name_p1) + ylab(y_axis_name_p1) +
    scale_color_manual(name="", values=colors_ttest_raw) + 
    theme(legend.position = "none", 
          plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent"),
          text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size, margin = margin(b=0.2, unit="cm")), 
          legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size, face = "bold"),
          plot.margin = margin(b=-0.4, unit = "cm")) + 
    guides(col=guide_legend(nrow=length(colors_ttest_raw)))
  
  p7 = ggplot(sub_plot_df, aes(x=logfc, y=ttest_adj_pval_log10, color=ttest_adj_cat)) + 
    geom_point(size = 0.75) +
    geom_hline(yintercept = -log10(sig_lvl), color="lightgrey", linetype="dashed", size = 0.25) +
    geom_vline(xintercept = -log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_vline(xintercept = log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) +
    xlab(x_axis_name_p2) + ylab(y_axis_name_p2) +
    scale_color_manual(name="", values=colors_ttest_adj) + 
    theme(legend.position = "none", 
          plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent"),
          text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size, margin = margin(b=0.2, unit="cm")), 
          legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size, face = "bold"),
          plot.margin = margin(b=-0.4, unit = "cm")) + 
    guides(col=guide_legend(nrow=length(colors_ttest_adj)))
  
  p8 = ggplot(sub_plot_df, aes(x=logfc, y=wilcox_raw_pval_log10, color=wilcox_raw_cat)) + 
    geom_point(size = 0.75) + 
    geom_hline(yintercept = -log10(sig_lvl), color="lightgrey", linetype="dashed", size = 0.25) +
    geom_vline(xintercept = -log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_vline(xintercept = log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) +
    xlab(x_axis_name_p3) + ylab(y_axis_name_p3) +
    scale_color_manual(name="", values=colors_wilcox_raw) + 
    theme(legend.position = "none", 
          plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent"),
          text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size, margin = margin(b=0.2, unit="cm")), 
          legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size, face = "bold"),
          plot.margin = margin(b=-0.4, unit = "cm")) + 
    guides(col=guide_legend(nrow=length(colors_wilcox_raw)))
  
  p9 = ggplot(sub_plot_df, aes(x=logfc, y=wilcox_adj_pval_log10, color=wilcox_adj_cat)) + 
    geom_point(size = 0.75) + 
    geom_hline(yintercept = -log10(sig_lvl), color="lightgrey", linetype="dashed", size = 0.25) +
    geom_vline(xintercept = -log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_vline(xintercept = log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) +
    xlab(x_axis_name_p4) + ylab(y_axis_name_p4) + 
    scale_color_manual(name="", values=colors_wilcox_adj) + 
    theme(legend.position = "none", 
          plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent"),
          text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size, margin = margin(b=0.2, unit="cm")), 
          legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size, face = "bold"),
          plot.margin = margin(b=-0.4, unit = "cm")) + 
    guides(col=guide_legend(nrow=length(colors_wilcox_adj)))
  
  p10 = ggplot(sub_plot_df, aes(x=logfc, y=cohend_abs, color=cohend_abs_cat)) + 
    geom_point(size = 0.75) + 
    geom_hline(yintercept = abs(effect_size_thres), color="lightgrey", linetype="dashed", size = 0.25) +
    geom_vline(xintercept = -log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_vline(xintercept = log_fc_thres, color="lightgrey", linetype="dashed", size = 0.25) +
    xlab(x_axis_name_p5) + ylab(y_axis_name_p5) +
    scale_color_manual(name="", values=colors_cohend_estimate) + 
    theme(legend.position = "none", 
          plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent"),
          text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size, margin = margin(b=0.2, unit="cm")), 
          legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size, face = "bold"),
          plot.margin = margin(b=-0.4, unit = "cm")) + 
    guides(col=guide_legend(nrow=length(colors_cohend_estimate)))
  
  ttest_raw_max_value[[c$subset[2]]] = max(sub_plot_df$ttest_raw_pval_log10)
  ttest_adj_max_value[[c$subset[2]]] = max(sub_plot_df$ttest_adj_pval_log10)
  wilcox_raw_max_value[[c$subset[2]]] = max(sub_plot_df$wilcox_raw_pval_log10)
  wilcox_adj_max_value[[c$subset[2]]] = max(sub_plot_df$wilcox_adj_pval_log10)
  cohend_abs_max_value[[c$subset[2]]] = max(sub_plot_df$cohend_abs)
  
  #if (prop == manual_plots_parameters$comparison[1]) {
  #  if ((c$g1 == manual_plots_parameters$comparison[2]) & (c$g2 == manual_plots_parameters$comparison[3])) {
  #    if (c$subset[1] == manual_plots_parameters$subset[1]) {
  #      if (c$subset[2] %in% manual_plots_parameters$subset[2:length(manual_plots_parameters$subset)]) {
  plot_title = unlist(unname(xticks_names[[c$subset[1]]][paste(c$subset[2])]))
  p1_plots[[c$subset[2]]] = p1 + xlab("") + ylab("") + ggtitle(plot_title)
  p2_plots[[c$subset[2]]] = p2 + xlab("") + ylab("") + ggtitle(plot_title)
  p3_plots[[c$subset[2]]] = p3 + xlab("") + ylab("") + ggtitle(plot_title)
  p4_plots[[c$subset[2]]] = p4 + xlab("") + ylab("") + ggtitle(plot_title)
  p5_plots[[c$subset[2]]] = p5 + xlab("") + ylab("") + ggtitle(plot_title)
  p6_plots[[c$subset[2]]] = p1 + xlab("") + ylab("") + ggtitle(plot_title) + geom_text_repel(aes(label = ttest_raw_cat_show_manual),
                                                                                                     max.time = repel_max_time, max.iter = 1e6,
                                                                                                     max.overlaps = Inf,
                                                                                                     box.padding = 0.1,
                                                                                                     point.patting = 0.1,
                                                                                                     force_pull = 2,
                                                                                                     force = 100,
                                                                                                     min.segment.length = 0, segment.size = 0.25,
                                                                                                     seed = snakemake@params$parameters_porps$set_seed,
                                                                                                     size = 6/11.04*3.88, family = plots_props$font_family, color = "black")
  p7_plots[[c$subset[2]]] = p2 + xlab("") + ylab("") + ggtitle(plot_title) + geom_text_repel(aes(label = ttest_adj_cat_show_manual),
                                                                                                     max.time = repel_max_time, max.iter = 1e6,
                                                                                                     max.overlaps = Inf,
                                                                                                     box.padding = 0.1,
                                                                                                     point.patting = 0.1,
                                                                                                     force_pull = 2,
                                                                                                     force = 100,
                                                                                                     min.segment.length = 0, segment.size = 0.25,
                                                                                                     seed = snakemake@params$parameters_porps$set_seed,
                                                                                                     size = 6/11.04*3.88, family = plots_props$font_family, color = "black")
  p8_plots[[c$subset[2]]] = p3 + xlab("") + ylab("") + ggtitle(plot_title) + geom_text_repel(aes(label = wilcox_raw_cat_show_manual),
                                                                                                     max.time = repel_max_time, max.iter = 1e6,
                                                                                                     max.overlaps = Inf,
                                                                                                     box.padding = 0.1,
                                                                                                     point.patting = 0.1,
                                                                                                     force_pull = 2,
                                                                                                     force = 100,
                                                                                                     min.segment.length = 0, segment.size = 0.25,
                                                                                                     seed = snakemake@params$parameters_porps$set_seed,
                                                                                                     size = 6/11.04*3.88, family = plots_props$font_family, color = "black")
  p9_plots[[c$subset[2]]] = p4 + xlab("") + ylab("") + ggtitle(plot_title) + geom_text_repel(aes(label = wilcox_adj_cat_show_manual),
                                                                                                     max.time = repel_max_time, max.iter = 1e6,
                                                                                                     max.overlaps = Inf, 
                                                                                                     box.padding = 0.1, 
                                                                                                     point.patting = 0.1, 
                                                                                                     force_pull = 2, 
                                                                                                     force = 100,
                                                                                                     min.segment.length = 0, segment.size = 0.25,
                                                                                                     seed = snakemake@params$parameters_porps$set_seed,
                                                                                                     size = 6/11.04*3.88, family = plots_props$font_family, color = "black")
  p10_plots[[c$subset[2]]] = p5 + xlab("") + ylab("") + ggtitle(plot_title) + geom_text_repel(aes(label = cohend_abs_cat_show_manual),
                                                                                                      max.time = repel_max_time, max.iter = 1e6,
                                                                                                      max.overlaps = Inf, 
                                                                                                      box.padding = 0.1, 
                                                                                                      point.patting = 0.1, 
                                                                                                      force_pull = 2, 
                                                                                                      force = 100,
                                                                                                      min.segment.length = 0, segment.size = 0.25,
                                                                                                      seed = snakemake@params$parameters_porps$set_seed,
                                                                                                      size = 6/11.04*3.88, family = plots_props$font_family, color = "black")
    #    }
    #  }
    #} else if ((c$g1 == manual_plots_parameters$comparison[3]) & (c$g2 == manual_plots_parameters$comparison[2])) {
    #  if (c$subset[1] == manual_plots_parameters$subset[1]) {
    #    if (c$subset[2] %in% manual_plots_parameters$subet[2:length(manual_plots_parameters$subset)]) {
    #      p1_plots[[c$subset[2]]] = p1 + xlab("") + ylab("") + ggtitle(paste(c$subset[2]))
    #      p2_plots[[c$subset[2]]] = p2 + xlab("") + ylab("") + ggtitle(paste(c$subset[2]))
    #      p3_plots[[c$subset[2]]] = p3 + xlab("") + ylab("") + ggtitle(paste(c$subset[2]))
    #      p4_plots[[c$subset[2]]] = p4 + xlab("") + ylab("") + ggtitle(paste(c$subset[2]))
    #      p5_plots[[c$subset[2]]] = p5 + xlab("") + ylab("") + ggtitle(paste(c$subset[2]))
    #      p6_plots[[c$subset[2]]] = p1 + xlab("") + ylab("") + ggtitle(paste(c$subset[2])) + geom_text_repel(aes(label = ttest_raw_cat_show_manual),
    #                                                                           max.time = repel_max_time, max.iter = 1e6,
    #                                                                           max.overlaps = Inf, 
    #                                                                           box.padding = 0.1, 
    #                                                                           point.patting = 0.1, 
    #                                                                           force_pull = 2, 
    #                                                                           force = 100,
    #                                                                           min.segment.length = 0, segment.size = 0.25,
    #                                                                           seed = snakemake@params$parameters_porps$set_seed,
    #                                                                           size = 6/11.04*3.88, family = plots_props$font_family, color = "black")
    #      p7_plots[[c$subset[2]]] = p2 + xlab("") + ylab("") + ggtitle(paste(c$subset[2])) + geom_text_repel(aes(label = ttest_adj_cat_show_manual),
    #                                                                           max.time = repel_max_time, max.iter = 1e6,
    #                                                                           max.overlaps = Inf, 
    #                                                                           box.padding = 0.1, 
    #                                                                           point.patting = 0.1, 
    #                                                                           force_pull = 2, 
    #                                                                           force = 100,
    #                                                                           min.segment.length = 0, segment.size = 0.25,
    #                                                                           seed = snakemake@params$parameters_porps$set_seed,
    #                                                                           size = 6/11.04*3.88, family = plots_props$font_family, color = "black")
    #      p8_plots[[c$subset[2]]] = p3 + xlab("") + ylab("") + ggtitle(paste(c$subset[2])) + geom_text_repel(aes(label = wilcox_raw_cat_manual),
    #                                                                           max.time = repel_max_time, max.iter = 1e6,
    #                                                                           max.overlaps = Inf, 
    #                                                                           box.padding = 0.1, 
    #                                                                           point.patting = 0.1, 
    #                                                                           force_pull = 2, 
    #                                                                           force = 100,
    #                                                                           min.segment.length = 0, segment.size = 0.25,
    #                                                                           seed = snakemake@params$parameters_porps$set_seed,
    #                                                                           size = 6/11.04*3.88, family = plots_props$font_family, color = "black")
    #      p9_plots[[c$subset[2]]] = p4 + xlab("") + ylab("") + ggtitle(paste(c$subset[2])) + geom_text_repel(aes(label = wilcox_adj_cat_show_manual),
    #                                                                           max.time = repel_max_time, max.iter = 1e6,
    #                                                                           max.overlaps = Inf, 
    #                                                                           box.padding = 0.1, 
    #                                                                           point.patting = 0.1, 
    #                                                                           force_pull = 2, 
    #                                                                           force = 100,
    #                                                                           min.segment.length = 0, segment.size = 0.25,
    #                                                                           seed = snakemake@params$parameters_porps$set_seed,
    #                                                                           size = 6/11.04*3.88, family = plots_props$font_family, color = "black")
    #      p10_plots[[c$subset[2]]] = p5 + xlab("") + ylab("") + ggtitle(paste(c$subset[2])) + geom_text_repel(aes(label = cohend_abs_cat_show_manual),
    #                                                                            max.time = repel_max_time, max.iter = 1e6,
    #                                                                            max.overlaps = Inf, 
    #                                                                            box.padding = 0.1, 
    #                                                                            point.patting = 0.1, 
    #                                                                            force_pull = 2, 
    #                                                                            force = 100,
    #                                                                            min.segment.length = 0, segment.size = 0.25,
    #                                                                            seed = snakemake@params$parameters_porps$set_seed,
    #                                                                            size = 6/11.04*3.88, family = plots_props$font_family, color = "black")
    #     }
    #   }
    #}
  #}
}

iteration_list_names = c()
for(i in 1:length(comps)){
  c = comps[[i]]
  
  p1_plots[[c$subset[2]]] = p1_plots[[c$subset[2]]] + ylim(0, max(unlist(ttest_raw_max_value)))
  p2_plots[[c$subset[2]]] = p2_plots[[c$subset[2]]] + ylim(0, max(unlist(ttest_adj_max_value)))
  p3_plots[[c$subset[2]]] = p3_plots[[c$subset[2]]] + ylim(0, max(unlist(wilcox_raw_max_value)))
  p4_plots[[c$subset[2]]] = p4_plots[[c$subset[2]]] + ylim(0, max(unlist(wilcox_adj_max_value)))
  p5_plots[[c$subset[2]]] = p5_plots[[c$subset[2]]] + ylim(0, max(unlist(cohend_abs_max_value)))
  p6_plots[[c$subset[2]]] = p6_plots[[c$subset[2]]] + ylim(0, max(unlist(ttest_raw_max_value)))
  p7_plots[[c$subset[2]]] = p7_plots[[c$subset[2]]] + ylim(0, max(unlist(ttest_adj_max_value)))
  p8_plots[[c$subset[2]]] = p8_plots[[c$subset[2]]] + ylim(0, max(unlist(wilcox_raw_max_value)))
  p9_plots[[c$subset[2]]] = p9_plots[[c$subset[2]]] + ylim(0, max(unlist(wilcox_adj_max_value)))
  p10_plots[[c$subset[2]]] = p10_plots[[c$subset[2]]] + ylim(0, max(unlist(cohend_abs_max_value)))
  
  iteration_list_names = append(iteration_list_names, c$subset[2])
}

# sort the plots like specified in the config
tmp = names(p1_plots) %in% iteration_list_names
p1_plots = p1_plots[tmp]
p2_plots = p2_plots[tmp]
p3_plots = p3_plots[tmp]
p4_plots = p4_plots[tmp]
p5_plots = p5_plots[tmp]
p6_plots = p6_plots[tmp]
p7_plots = p7_plots[tmp]
p8_plots = p8_plots[tmp]
p9_plots = p9_plots[tmp]
p10_plots = p10_plots[tmp]

plot_nrows = manual_plots_parameters$plot_nrows
plot_ncols = manual_plots_parameters$plot_ncols

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

p1_plots_mod_margins = p1_plots
p2_plots_mod_margins = p2_plots
p3_plots_mod_margins = p3_plots
p4_plots_mod_margins = p4_plots
p5_plots_mod_margins = p5_plots
p6_plots_mod_margins = p6_plots
p7_plots_mod_margins = p7_plots
p8_plots_mod_margins = p8_plots
p9_plots_mod_margins = p9_plots
p10_plots_mod_margins = p10_plots

#if ((plot_nrows == 1) & (plot_ncols ==2)) {
top_margin_left_col = 0 #
right_margin_left_col = -6 # -14
bottom_margin_left_col = -10 #!
left_margin_left_col = -8 #!

top_margin_middle_cols = top_margin_left_col
right_margin_middle_cols = -4
bottom_margin_middle_cols = bottom_margin_left_col
left_margin_middle_cols = -8

top_margin_right_col = top_margin_left_col
right_margin_right_col = 2 #!
bottom_margin_right_col = bottom_margin_left_col
left_margin_right_col = -16 
#}

p1_plots = plot_fix_layout(p1_plots_mod_margins, ttest_raw_max_value,
                top_margin_left_col, right_margin_left_col, bottom_margin_left_col, left_margin_left_col,
                top_margin_middle_cols, right_margin_middle_cols, bottom_margin_middle_cols, left_margin_middle_cols,
                top_margin_right_col, right_margin_right_col, bottom_margin_right_col, left_margin_right_col
                )
p2_plots = plot_fix_layout(p2_plots_mod_margins, ttest_adj_max_value,
                top_margin_left_col, right_margin_left_col, bottom_margin_left_col, left_margin_left_col,
                top_margin_middle_cols, right_margin_middle_cols, bottom_margin_middle_cols, left_margin_middle_cols,
                top_margin_right_col, right_margin_right_col, bottom_margin_right_col, left_margin_right_col
                )
p3_plots = plot_fix_layout(p3_plots_mod_margins, wilcox_raw_max_value,
                top_margin_left_col, right_margin_left_col, bottom_margin_left_col, left_margin_left_col,
                top_margin_middle_cols, right_margin_middle_cols, bottom_margin_middle_cols, left_margin_middle_cols,
                top_margin_right_col, right_margin_right_col, bottom_margin_right_col, left_margin_right_col
                )
p4_plots = plot_fix_layout(p4_plots_mod_margins, wilcox_adj_max_value,
                top_margin_left_col, right_margin_left_col, bottom_margin_left_col, left_margin_left_col,
                top_margin_middle_cols, right_margin_middle_cols, bottom_margin_middle_cols, left_margin_middle_cols,
                top_margin_right_col, right_margin_right_col, bottom_margin_right_col, left_margin_right_col
                )
p5_plots = plot_fix_layout(p5_plots_mod_margins, cohend_abs_max_value,
                top_margin_left_col, right_margin_left_col, bottom_margin_left_col, left_margin_left_col,
                top_margin_middle_cols, right_margin_middle_cols, bottom_margin_middle_cols, left_margin_middle_cols,
                top_margin_right_col, right_margin_right_col, bottom_margin_right_col, left_margin_right_col
                )
p6_plots = plot_fix_layout(p6_plots_mod_margins, ttest_raw_max_value,
                top_margin_left_col, right_margin_left_col, bottom_margin_left_col, left_margin_left_col,
                top_margin_middle_cols, right_margin_middle_cols, bottom_margin_middle_cols, left_margin_middle_cols,
                top_margin_right_col, right_margin_right_col, bottom_margin_right_col, left_margin_right_col
                )
p7_plots = plot_fix_layout(p7_plots_mod_margins, ttest_adj_max_value,
                top_margin_left_col, right_margin_left_col, bottom_margin_left_col, left_margin_left_col,
                top_margin_middle_cols, right_margin_middle_cols, bottom_margin_middle_cols, left_margin_middle_cols,
                top_margin_right_col, right_margin_right_col, bottom_margin_right_col, left_margin_right_col
                )
p8_plots = plot_fix_layout(p8_plots_mod_margins, wilcox_raw_max_value,
                top_margin_left_col, right_margin_left_col, bottom_margin_left_col, left_margin_left_col,
                top_margin_middle_cols, right_margin_middle_cols, bottom_margin_middle_cols, left_margin_middle_cols,
                top_margin_right_col, right_margin_right_col, bottom_margin_right_col, left_margin_right_col
                )
p9_plots = plot_fix_layout(p9_plots_mod_margins, wilcox_adj_max_value,
                top_margin_left_col, right_margin_left_col, bottom_margin_left_col, left_margin_left_col,
                top_margin_middle_cols, right_margin_middle_cols, bottom_margin_middle_cols, left_margin_middle_cols,
                top_margin_right_col, right_margin_right_col, bottom_margin_right_col, left_margin_right_col
                )
p10_plots = plot_fix_layout(p10_plots_mod_margins, cohend_abs_max_value,
                top_margin_left_col, right_margin_left_col, bottom_margin_left_col, left_margin_left_col,
                top_margin_middle_cols, right_margin_middle_cols, bottom_margin_middle_cols, left_margin_middle_cols,
                top_margin_right_col, right_margin_right_col, bottom_margin_right_col, left_margin_right_col
                )

all_p1_plots = arrangeGrob(grobs=p1_plots, ncol=manual_plots_parameters$plot_ncols, nrow=manual_plots_parameters$plot_nrows, bottom = textGrob(x_axis_name_p1, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_name_p1, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
all_p2_plots = arrangeGrob(grobs=p2_plots, ncol=manual_plots_parameters$plot_ncols, nrow=manual_plots_parameters$plot_nrows, bottom = textGrob(x_axis_name_p2, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_name_p2, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
all_p3_plots = arrangeGrob(grobs=p3_plots, ncol=manual_plots_parameters$plot_ncols, nrow=manual_plots_parameters$plot_nrows, bottom = textGrob(x_axis_name_p3, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_name_p3, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
all_p4_plots = arrangeGrob(grobs=p4_plots, ncol=manual_plots_parameters$plot_ncols, nrow=manual_plots_parameters$plot_nrows, bottom = textGrob(x_axis_name_p4, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_name_p4, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
all_p5_plots = arrangeGrob(grobs=p5_plots, ncol=manual_plots_parameters$plot_ncols, nrow=manual_plots_parameters$plot_nrows, bottom = textGrob(x_axis_name_p5, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_name_p5, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
all_p6_plots = arrangeGrob(grobs=p6_plots, ncol=manual_plots_parameters$plot_ncols, nrow=manual_plots_parameters$plot_nrows, bottom = textGrob(x_axis_name_p1, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_name_p1, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
all_p7_plots = arrangeGrob(grobs=p7_plots, ncol=manual_plots_parameters$plot_ncols, nrow=manual_plots_parameters$plot_nrows, bottom = textGrob(x_axis_name_p2, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_name_p2, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
all_p8_plots = arrangeGrob(grobs=p8_plots, ncol=manual_plots_parameters$plot_ncols, nrow=manual_plots_parameters$plot_nrows, bottom = textGrob(x_axis_name_p3, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_name_p3, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
all_p9_plots = arrangeGrob(grobs=p9_plots, ncol=manual_plots_parameters$plot_ncols, nrow=manual_plots_parameters$plot_nrows, bottom = textGrob(x_axis_name_p4, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_name_p4, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
all_p10_plots = arrangeGrob(grobs=p10_plots, ncol=manual_plots_parameters$plot_ncols, nrow=manual_plots_parameters$plot_nrows, bottom = textGrob(x_axis_name_p5, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_name_p5, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))

height = manual_plots_parameters$plot_height
width = manual_plots_parameters$plot_width
ggsave(sprintf("%s/%s__%s_%s=%s.ttest.raw_%sx%s_%sx%s.png", output_folder_merged, names(comps[1]), paste(comps[[1]]$g, collapse= "_vs_"), comps[[1]]$subset[1], paste(iteration_list_names, collapse= "_"), manual_plots_parameters$plot_nrows, manual_plots_parameters$plot_ncols, width, height), all_p1_plots, dpi=plots_props$dpi, width= width, height=height, units = plots_props$image_units)
ggsave(sprintf("%s/%s__%s_%s=%s.ttest.raw_%sx%s_%sx%s.svg", output_folder_merged, names(comps[1]), paste(comps[[1]]$g, collapse= "_vs_"), comps[[1]]$subset[1], paste(iteration_list_names, collapse= "_"), manual_plots_parameters$plot_nrows, manual_plots_parameters$plot_ncols, width, height), all_p1_plots, width=width, height=height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
ggsave(sprintf("%s/%s__%s_%s=%s.ttest.adj_%sx%s_%sx%s.png", output_folder_merged, names(comps[1]), paste(comps[[1]]$g, collapse= "_vs_"), comps[[1]]$subset[1], paste(iteration_list_names, collapse= "_"), manual_plots_parameters$plot_nrows, manual_plots_parameters$plot_ncols, width, height), all_p2_plots, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
ggsave(sprintf("%s/%s__%s_%s=%s.ttest.adj_%sx%s_%sx%s.svg", output_folder_merged, names(comps[1]), paste(comps[[1]]$g, collapse= "_vs_"), comps[[1]]$subset[1], paste(iteration_list_names, collapse= "_"), manual_plots_parameters$plot_nrows, manual_plots_parameters$plot_ncols, width, height), all_p2_plots, width= width, height=height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
ggsave(sprintf("%s/%s__%s_%s=%s.wilcox.raw_%sx%s_%sx%s.png", output_folder_merged, names(comps[1]), paste(comps[[1]]$g, collapse= "_vs_"), comps[[1]]$subset[1], paste(iteration_list_names, collapse= "_"), manual_plots_parameters$plot_nrows, manual_plots_parameters$plot_ncols, width, height), all_p3_plots, dpi=plots_props$dpi, width= width, height=height, units = plots_props$image_units)
ggsave(sprintf("%s/%s__%s_%s=%s.wilcox.raw_%sx%s_%sx%s.svg", output_folder_merged, names(comps[1]), paste(comps[[1]]$g, collapse= "_vs_"), comps[[1]]$subset[1], paste(iteration_list_names, collapse= "_"), manual_plots_parameters$plot_nrows, manual_plots_parameters$plot_ncols, width, height), all_p3_plots, width= width, height=height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
ggsave(sprintf("%s/%s__%s_%s=%s.wilcox.adj_%sx%s_%sx%s.png", output_folder_merged, names(comps[1]), paste(comps[[1]]$g, collapse= "_vs_"), comps[[1]]$subset[1], paste(iteration_list_names, collapse= "_"), manual_plots_parameters$plot_nrows, manual_plots_parameters$plot_ncols, width, height), all_p4_plots, dpi=plots_props$dpi, width= width, height=height, units = plots_props$image_units)
ggsave(sprintf("%s/%s__%s_%s=%s.wilcox.adj_%sx%s_%sx%s.svg", output_folder_merged, names(comps[1]), paste(comps[[1]]$g, collapse= "_vs_"), comps[[1]]$subset[1], paste(iteration_list_names, collapse= "_"), manual_plots_parameters$plot_nrows, manual_plots_parameters$plot_ncols, width, height), all_p4_plots, width= width, height=height, units=plots_props$image_units, limitsize = FALSE, scale = 1)
ggsave(sprintf("%s/%s__%s_%s=%s.cohend_estimate_%sx%s_%sx%s.png", output_folder_merged, names(comps[1]), paste(comps[[1]]$g, collapse= "_vs_"), comps[[1]]$subset[1], paste(iteration_list_names, collapse= "_"), manual_plots_parameters$plot_nrows, manual_plots_parameters$plot_ncols, width, height), all_p5_plots, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
ggsave(sprintf("%s/%s__%s_%s=%s.cohend_estimate_%sx%s_%sx%s.svg", output_folder_merged, names(comps[1]), paste(comps[[1]]$g, collapse= "_vs_"), comps[[1]]$subset[1], paste(iteration_list_names, collapse= "_"), manual_plots_parameters$plot_nrows, manual_plots_parameters$plot_ncols, width, height), all_p5_plots, width= width, height=height, units =plots_props$image_units, limitsize = FALSE, scale = 1)

ggsave(sprintf("%s/%s__%s_%s=%s.ttest.raw_%sx%s_labels_%sx%s.png", output_folder_merged, names(comps[1]), paste(comps[[1]]$g, collapse= "_vs_"), comps[[1]]$subset[1], paste(iteration_list_names, collapse= "_"), manual_plots_parameters$plot_nrows, manual_plots_parameters$plot_ncols, width, height), all_p6_plots, dpi=plots_props$dpi, width= width, height=height, units = plots_props$image_units)
ggsave(sprintf("%s/%s__%s_%s=%s.ttest.raw_%sx%s_labels_%sx%s.svg", output_folder_merged, names(comps[1]), paste(comps[[1]]$g, collapse= "_vs_"), comps[[1]]$subset[1], paste(iteration_list_names, collapse= "_"), manual_plots_parameters$plot_nrows, manual_plots_parameters$plot_ncols, width, height), all_p6_plots, width=width, height=height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
ggsave(sprintf("%s/%s__%s_%s=%s.ttest.adj_%sx%s_labels_%sx%s.png", output_folder_merged, names(comps[1]), paste(comps[[1]]$g, collapse= "_vs_"), comps[[1]]$subset[1], paste(iteration_list_names, collapse= "_"), manual_plots_parameters$plot_nrows, manual_plots_parameters$plot_ncols, width, height), all_p7_plots, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
ggsave(sprintf("%s/%s__%s_%s=%s.ttest.adj_%sx%s_labels_%sx%s.svg", output_folder_merged, names(comps[1]), paste(comps[[1]]$g, collapse= "_vs_"), comps[[1]]$subset[1], paste(iteration_list_names, collapse= "_"), manual_plots_parameters$plot_nrows, manual_plots_parameters$plot_ncols, width, height), all_p7_plots, width= width, height=height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
ggsave(sprintf("%s/%s__%s_%s=%s.wilcox.raw_%sx%s_labels_%sx%s.png", output_folder_merged, names(comps[1]), paste(comps[[1]]$g, collapse= "_vs_"), comps[[1]]$subset[1], paste(iteration_list_names, collapse= "_"), manual_plots_parameters$plot_nrows, manual_plots_parameters$plot_ncols, width, height), all_p8_plots, dpi=plots_props$dpi, width= width, height=height, units = plots_props$image_units)
ggsave(sprintf("%s/%s__%s_%s=%s.wilcox.raw_%sx%s_labels_%sx%s.svg", output_folder_merged, names(comps[1]), paste(comps[[1]]$g, collapse= "_vs_"), comps[[1]]$subset[1], paste(iteration_list_names, collapse= "_"), manual_plots_parameters$plot_nrows, manual_plots_parameters$plot_ncols, width, height), all_p8_plots, width= width, height=height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
ggsave(sprintf("%s/%s__%s_%s=%s.wilcox.adj_%sx%s_labels_%sx%s.png", output_folder_merged, names(comps[1]), paste(comps[[1]]$g, collapse= "_vs_"), comps[[1]]$subset[1], paste(iteration_list_names, collapse= "_"), manual_plots_parameters$plot_nrows, manual_plots_parameters$plot_ncols, width, height), all_p9_plots, dpi=plots_props$dpi, width= width, height=height, units = plots_props$image_units)
ggsave(sprintf("%s/%s__%s_%s=%s.wilcox.adj_%sx%s_labels_%sx%s.svg", output_folder_merged, names(comps[1]), paste(comps[[1]]$g, collapse= "_vs_"), comps[[1]]$subset[1], paste(iteration_list_names, collapse= "_"), manual_plots_parameters$plot_nrows, manual_plots_parameters$plot_ncols, width, height), all_p9_plots, width= width, height=height, units=plots_props$image_units, limitsize = FALSE, scale = 1)
ggsave(sprintf("%s/%s__%s_%s=%s.cohend_estimate_%sx%s_labels_%sx%s.png", output_folder_merged, names(comps[1]), paste(comps[[1]]$g, collapse= "_vs_"), comps[[1]]$subset[1], paste(iteration_list_names, collapse= "_"), manual_plots_parameters$plot_nrows, manual_plots_parameters$plot_ncols, width, height), all_p10_plots, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
ggsave(sprintf("%s/%s__%s_%s=%s.cohend_estimate_%sx%s_labels_%sx%s.svg", output_folder_merged, names(comps[1]), paste(comps[[1]]$g, collapse= "_vs_"), comps[[1]]$subset[1], paste(iteration_list_names, collapse= "_"), manual_plots_parameters$plot_nrows, manual_plots_parameters$plot_ncols, width, height), all_p10_plots, width= width, height=height, units =plots_props$image_units, limitsize = FALSE, scale = 1)

#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])

