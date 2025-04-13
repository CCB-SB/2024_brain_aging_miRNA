suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))

#snakemake = readRDS("/local/23-04_brain_aging_sncrna_data_analysis/annika/brain_aging/snakemake_fc_versus_expr_scatter_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("fc_versus_expr_scatter_plot")


# ------------------------------- Functions -----------------------------------
scatter_plot = function(plot_df_all, thresholds, property, shape, legend_title, data_input, diff_log, plots_props, output_folder) {
  test_value = "logfc"
  #for (test_value in c("ttest_adj_pval_log10", "logfc", "auc")) {
  xlab_name = sprintf("Expression value (%s%s)", gsub("_norm", "", data_input$norm), sprintf(", %s", gsub("_", "", diff_log)))
  #if (test_value == "ttest_adj_pval_log10") {
  #  ylab_name = "Adj. p-value from ttest (-log10)"
  #} else if (test_value == "logfc") {
  ylab_name = "Fold change (log2)"
  #} else {
  #  ylab_name = "Area under curve"
  #}
    
  if (shape == TRUE) {
    p = ggplot(plot_df_all, aes(x=expr, y=.data[[test_value]], color=group, shape=time, text="ID")) + #label=RNA,
      geom_point(size = 2, alpha = 0.75) + 
      #geom_text(hjust=-0.1) +
      scale_colour_manual(values = unlist(colors[[property]])) +
      scale_shape() +
      #geom_label_repel(aes(label=RNA), show.legend = FALSE) +
      xlab(xlab_name) +
      ylab(ylab_name) + 
      theme_classic() +
      theme(legend.position="bottom", 
            text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
            axis.text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
            axis.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
            plot.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size_header),
            legend.text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
            legend.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size, face = "bold"),
            legend.margin=margin(t = -8), legend.key.size = unit(0.3, "cm"), legend.spacing.y = unit(0.4, "cm")) +
      guides(col="none", shape=guide_legend(nrow=2, title=legend_title))
  } else {
    if (length(unique(annot[[property]])) == 2) {
      alpha_value = 1
    } else {
      alpha_value = 0.75
    }
    p = ggplot(plot_df_all, aes(x=expr, y=.data[[test_value]], color=group, text="ID")) + #label=RNA
      geom_point(size = 2, alpha = alpha_value) + 
      #geom_text(hjust=-0.1) +
      # scale_colour_manual(values = unlist(colors[[property]])) +
      scale_colour_manual(values = unlist(colors[[property]])[names(unlist(colors[[property]])) %in% plot_df_all$group]) +
      #geom_label_repel(aes(label=RNA), show.legend = FALSE) +
      xlab(xlab_name) +
      ylab(ylab_name) + 
      theme_classic() +
      theme(legend.position="none", 
            text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
            axis.text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
            axis.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
            plot.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size_header),
            legend.text = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size),
            legend.title = element_text(family = snakemake@params$plots_props$font_family, size = snakemake@params$plots_props$font_size, face = "bold"),
            legend.margin=margin(t = -8), legend.key.size = unit(0.3, "cm"), legend.spacing.y = unit(0.4, "cm"))
  }
  p
  #if (test_value == "ttest_adj_pval_log10") {
  #th = -log10(thresholds$adj_p_value)
  #p = p + geom_hline(yintercept=th, linetype="dashed", color="lightgrey", size=0.5)
  #} else if (test_value == "logfc") {
  th_1 = log2(thresholds$fc)
  th_2 = log2(1/thresholds$fc)
  p = p + 
    geom_hline(yintercept=th_1, linetype="dashed", color="lightgrey", size=0.5) + 
    geom_hline(yintercept=th_2, linetype="dashed", color="lightgrey", size=0.5) +
    geom_vline(xintercept=log10(15), linetype="dashed", color="lightgrey", size=0.5) 
  #}

  ggsave(sprintf("%s/%s_vs_expr_coloured=%s.png", output_folder, test_value, property), plot = p, dpi = snakemake@params$plots_props$dpi, width = snakemake@params$plots_props$image_width, height = snakemake@params$plots_props$image_height, units = snakemake@params$plots_props$image_units)
  ggsave(sprintf("%s/%s_vs_expr_coloured=%s.svg", output_folder, test_value, property), plot = p, width = snakemake@params$plots_props$image_width, height = snakemake@params$plots_props$image_height, units = snakemake@params$plots_props$image_units, limitsize = FALSE, scale = 1)
  #}
}

# ------------------------------- Input data -----------------------------------
data_input = snakemake@params$data_input
diff_log = snakemake@params$diff_exp_log
thresholds = snakemake@params$thresholds
data_sub_set = data_input$data_sub_set
property = snakemake@params$prop
xticks_names = snakemake@params$xticks_names
colors = snakemake@params$colors
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "/local/23-04_brain_aging_sncrna_data_analysis/annika/brain_aging/snakemake_fc_versus_expr_scatter_plot.rds")

if (data_sub_set == "CA1_aging") {
  comps = snakemake@params$comparisons_CA1
} else if (data_sub_set == "CA1_aging_without_age=26m_28m"){
  comps = snakemake@params$comparisons_CA1_without_age_26m_28m
} else{ 
  stop(sprintf("cannot get comparisons for data subset %s", data_sub_set))
}

annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
diff_exp = fread(sprintf("%s/%s_%s/results_%s/matrices/diff_exp/diff_exp%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, diff_log), sep='\t') 

output_folder = sprintf("%s/%s_%s/results_%s/figures/fc_vs_expr/scatter_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder, recursive=TRUE)

# ------------------------- plot fc versus expr median -------------------------
prop = names(comps)[1]
if (prop != property) {
  plot_df_all = c()
  for(i in 1:length(comps)) {
    print(i)
    # which comp do we investigate and which groups are considered
    prop = names(comps)[i]
    c = comps[[i]]
    print(c)
    groups = c$g
    
    c$g1 = c$g[1]
    c$g2 = c$g[2]

    print(prop)
    
      shape = TRUE
      for (ti in unique(annot[[property]])) {

        
        paired_string = ifelse(!is.null(c$paired), sprintf("_paired_%s", c$paired), "")
        sub_string = sprintf("_%s=%s", property, ti)
        
        plot_df = data.table(feature=diff_exp[[data_input$rna_class]],
                            expr=diff_exp[[sprintf("median_%s__%s_%s=%s", prop, c$g1, property, ti)]],
                            ttest_adj_pval=diff_exp[[sprintf("ttest_adjp_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]],
                            fc=diff_exp[[sprintf("fc_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]],
                            auc=diff_exp[[sprintf("auc_value_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]],
                            group=rep(sprintf("%s", ti), length(diff_exp[[data_input$rna_class]])),
                            time=rep(sprintf("%s", c$g1), length(diff_exp[[data_input$rna_class]]))
                            )
        # plot_df[,expr := 2^expr_log2]
        plot_df[,expr := expr]
        plot_df[,ttest_adj_pval_log10 := -log10(ttest_adj_pval)]
        plot_df[,logfc := log2(fc)]
        
        plot_df_all = rbind(plot_df_all, plot_df)
      }
  }

  legend_title = xticks_names$categories[[sprintf("%s_capital", prop)]]
  scatter_plot(plot_df_all, thresholds, property, shape, legend_title, data_input, diff_log, plots_props, output_folder) 

} else {
  plot_df_all = c()
  for(i in 1:length(comps)) {
    shape = FALSE
    for (ti in groups[1:(length(unique(annot[[property]]))-1)]) {
      # which comp do we investigate and which groups are considered
      prop = names(comps)[i]
      c = comps[[i]]
      
      c$g1 = c$g[1]
      c$g2 = c$g[2]
      
      paired_string = ifelse(!is.null(c$paired), sprintf("_paired_%s", c$paired), "")

      plot_df = data.table(feature=diff_exp[[data_input$rna_class]],
                            expr=diff_exp[[sprintf("median_%s__%s", prop, c$g1)]],
                            ttest_adj_pval=diff_exp[[sprintf("ttest_adjp_%s__%s_vs_%s%s", prop, c$g1, c$g2, paired_string)]],
                            fc=diff_exp[[sprintf("fc_%s__%s_vs_%s%s", prop, c$g1, c$g2, paired_string)]],
                            auc=diff_exp[[sprintf("auc_value_%s__%s_vs_%s%s", prop, c$g1, c$g2, paired_string)]],
                            group=rep(sprintf("%s", ti), length(diff_exp[[data_input$rna_class]]))
      )
      # plot_df[,expr := 2^expr_log2]
      plot_df[,expr := expr]
      plot_df[,ttest_adj_pval_log10 := -log10(ttest_adj_pval)]
      plot_df[,logfc := log2(fc)]
      
      plot_df_all = rbind(plot_df_all, plot_df)
    }
  }

  scatter_plot(plot_df_all, thresholds, property, shape, "", data_input, diff_log, plots_props, output_folder) 
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
