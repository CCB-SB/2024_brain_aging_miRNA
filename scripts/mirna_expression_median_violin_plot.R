suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(foreach))  # for parallelization
suppressPackageStartupMessages(library(doParallel))  # for parallelization

# save rdata
#saveRDS(snakemake, file = "snakemake_mirna_expression_median_violin_plot.rds")
#stop()

#snakemake = readRDS("snakemake_mirna_expression_median_violin_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("mirna_expression_median_violin_plot")

#---------------------------------- Functions ---------------------------------- 
ploting_violin = function(df, annot, diff_exp_feature, rna, feature_config, plot_group, comp_group, params, order_of_x_axis, data_input, x_axis_title, y_axis_title_big, y_axis_title_small, xticks_names, plots_props, output_folder_violin_single, output_folder_box_single){
  p_violin = list()
  p_box = list()
  p_violin_manual = list()
  p_box_manual = list()
  j = 1
 
  if (plot_group == "brain_region") {
    iteration_list = c("all", unique(annot[[plot_group]]))
  } else {
    iteration_list = unique(annot[[plot_group]])
  }
  
  for (tissue in iteration_list) {
    #print(tissue)
    if (tissue != "all") {
      id_group = annot[annot[[plot_group]] == tissue,][[data_input$identifier_column]] 
      tmp = colnames(df) %in% id_group
      expr_feature_tissue = df[,..tmp]
      plot_title = sprintf("%s", unlist(xticks_names[[plot_group]][[tissue]]))
    } else {
      expr_feature_tissue = df
      plot_title = sprintf("%s", tissue)
    }
    
    plot_df = melt(expr_feature_tissue)
    # sort plot_df according to annot
    # only select rows from annot which are in plot_df
    target_order = annot[annot[[data_input$identifier_column]] %in% plot_df$variable][[data_input$identifier_column]]
    plot_df = plot_df[match(target_order, plot_df$variable), ]
    plot_df[[plot_group]] = annot[annot[[data_input$identifier_column]] %in% plot_df$variable][[plot_group]]

    comp_category_list = c()
    for (id in plot_df$variable) {
      comp_category_list = append(comp_category_list, annot[annot[[data_input$identifier_column]] == id,][[comp_group]])
    }
    plot_df$comp = comp_category_list
    
    if (unique(unique(plot_df$comp) %in% c("young", "old")) == TRUE) {
      plot_df$comp = factor(plot_df$comp, levels = c("young","old"),ordered = TRUE)
    }
    
    #x_axis_title = xticks_names$categories[[sprintf("%s_capital", comp)]]
    #if (grepl(" ", xticks_names$categories[[data_input$detection_group]]) == TRUE) {
    #  detection_group = strsplit(xticks_names$categories[[data_input$detection_group]], " ")[[1]][1]
    #}
    #y_axis_title_big = sprintf("Expr (%s, detection rate = %s%% per %s)", gsub("_norm", "" ,data_input$norm), gsub("p", "", data_input$detection_rate), detection_group)
    #y_axis_title_small = sprintf("Expr (%s,\ndetection rate = %s%%\nper %s)", gsub("_norm", "" ,data_input$norm), gsub("p", "", data_input$detection_rate), detection_group)
    
    if (length(order_of_x_axis) != 0) {
      time_median = c()
      time_sd = c()
      significance = c()
      for (time in order_of_x_axis) {
        #time_median = append(time_median, median(plot_df[plot_df$comp == time,]$value))
        #time_sd = append(time_sd, sd(plot_df[plot_df$comp == time,]$value))
        time_median = ""
        time_sd = ""
        
        if (time != order_of_x_axis[1]) {
          if (tissue == "all") {
            ttest_adj_pval = diff_exp_feature[[sprintf("ttest_adjp_%s__%s_vs_%s", comp_group, time, order_of_x_axis[1])]]  #ttest_adjp_yy__12_vs_3
            fc = log2(diff_exp_feature[[sprintf("fc_%s__%s_vs_%s", comp_group, time, order_of_x_axis[1])]])
          } else {
            if (plot_group == "brain_region") {
              ttest_adj_pval = diff_exp_feature[[sprintf("ttest_adjp_%s__%s_vs_%s_%s=%s", comp_group, time, order_of_x_axis[1], plot_group, tissue)]]  #ttest_adjp_yy__12_vs_3_brain_region_fixed=cc
              fc = log2(diff_exp_feature[[sprintf("fc_%s__%s_vs_%s_%s=%s", comp_group, time, order_of_x_axis[1], plot_group, tissue)]])  #ttest_adjp_yy__12_vs_3_brain_region_fixed=cc
            } else {
              ttest_adj_pval = diff_exp_feature[[sprintf("ttest_adjp_%s__%s_vs_%s", comp_group, time, order_of_x_axis[1])]]
              fc = log2(diff_exp_feature[[sprintf("fc_%s__%s_vs_%s", comp_group, time, order_of_x_axis[1])]]) 
            }
          }
          #print(sprintf("ttest_adjp_%s__%s_vs_%s", comp_group, time, order_of_x_axis[1]))
          
          if (is.na(fc)) {
            significance = append(significance, "")
          } else if (fc >= log2(params$fc_th)) {
            if (is.na(ttest_adj_pval)) { 
              significance = append(significance, "Up")
            } else if (ttest_adj_pval < 0.001) { 
              significance = append(significance, "***\nUp")
            } else if (ttest_adj_pval < 0.01) {
              significance = append(significance, "**\nUp")
            } else if (ttest_adj_pval < 0.05) {
              significance = append(significance, "*\nUp")
            } else {
              significance = append(significance, "Up")
            }
          } else if (fc <= -log2(params$fc_th)) {
            if (is.na(ttest_adj_pval)) { 
              significance = append(significance, "Down")
            } else if (ttest_adj_pval < 0.001) { 
              significance = append(significance, "***\nDown")
            } else if (ttest_adj_pval < 0.01) {
              significance = append(significance, "**\nDown")
            } else if (ttest_adj_pval < 0.05) {
              significance = append(significance, "*\nDown")
            } else {
              significance = append(significance, "Down")
            }
          } else {
            significance = append(significance, "")
          }
        } else {
          significance = append(significance, "")
        }
      }
      time = order_of_x_axis
    } else {
      time_median = ""
      time_sd = ""
      significance = ""
      
      time = sort(unique(annot[[comp_group]]))
    }

    plot_df_violin = data.frame(comp = time, time_median = time_median, time_sd = time_sd, significance = significance)
    
    # add the significance * to the timepoint only for which the value is the highest
    plot_df = merge(plot_df, plot_df_violin, by="comp")
    #sig_samples = plot_df[plot_df$significance != "",]$variable
    #max_value = max(plot_df[plot_df$variable %in% sig_samples]$value)
    #plot_df[plot_df$value != max_value,]$significance = ""
    for (t in plot_df$comp) {
      sig_samples = plot_df[(plot_df$comp == t) & (plot_df$significance != ""),]$variable
      max_value = max(plot_df[(plot_df$comp == t) & (plot_df$variable %in% sig_samples),]$value)
      plot_df[(plot_df$comp == t) & (plot_df$value != max_value),]$significance = ""
    }
    
    # replace all not max with "" for each params$comp_group
    #for (t in plot_df[[comp_group]]) {
    #  max_value = max(plot_df[(plot_df$variable %in% sig_samples) & (plot_df[[comp_group]] == t)]$value)
    #  plot_df[(plot_df$value != max_value) & (plot_df[[comp_group]] == t),]$significance = ""
    #}
    
    # corr line
    # if (tissue != "all") {
    #   mask = (corr[[data_input$rna_class]] == rna)
    #   corr_value = corr[mask,][[tissue]] 
    #   p = p_value[mask,][[tissue]]
    #   # %.3f makes 3 decimal places
    #   #corr_text = sprintf("corr: %.3f\np: %.3f", corr_value, p)
    #   if (p < params$sig_niveau) {
    #     corr_text = sprintf("sig. corr.: %.3f", corr_value) 
    #   } else {
    #     corr_text = sprintf("corr.: %.3f", corr_value) 
    #   }
    # } else {
    #   corr_text = ""
    # }
    
    tmp = unname(unlist(xticks_names[[comp_group]])[unique(annot[[comp_group]])])
    if (sum(nchar(tmp)) > 20) {
      x_axis_ticks_rot = TRUE
    } else {
      x_axis_ticks_rot = FALSE
    }
    
    if (is.numeric(plot_df$comp)) {
      tmp = max(plot_df$comp) + 2
    } else {
      tmp = length(unique(plot_df$comp)) + 0.15
    }
    
    font_size_correction = 2.845
    
    # sort x axis
    if (length(order_of_x_axis) != 0) {
      plot_df$comp = factor(plot_df$comp, levels=order_of_x_axis)
    }
    
    # create violin plot
    p_violin[[j]] = ggplot(plot_df, aes(x = comp, y = value, colour = comp, fill = comp)) + 
      geom_violin(trim = TRUE) +
      scale_fill_manual(values = alpha(unlist(colours[[comp_group]]), 0.3)) +
      scale_colour_manual(values = alpha(unlist(colours[[comp_group]]), 1)) +
      geom_jitter(data = plot_df, aes(y = value, x = comp), color = "darkgray", size = 0.75) +
      stat_summary(fun=function(x) 2^(median(log2(x))), geom="crossbar", size=.2, aes(color = comp)) + 
      ggtitle(plot_title) +
      xlab("") +
      ylab("") +
      ##ylab(sprintf("expr (%s, detection rate = %s%%\nper %s)", data_input$norm, expr_input$detection_rate, params$plot_group)) +
      ##scale_x_discrete(labels = plot_x_ticks, limits = names(plot_x_ticks)) +
      scale_x_discrete(labels = unlist(xticks_names[[comp_group]])) +
      scale_y_continuous(limits = c(0, 1.1*max(plot_df$value)), expand = c(0,0)) +
      geom_text(aes(label = significance), vjust=-0.3, size=3, color="black") +
      theme_classic() +
      theme(##plot.margin = unit(c(0.1,-1.9,0.25,-1.75), plots_props$image_units), #top, right, bottom, left
        plot.margin = unit(c(0,0.1,0,0.1), plots_props$image_units), #top, right, bottom, left
        ##axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9), 
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position="none", 
        text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
        plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size, margin = margin(b=0, unit="cm")),
        legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size)
      ) +
    coord_cartesian(clip = "off")

    if (x_axis_ticks_rot) {
      p_violin[[j]] =  p_violin[[j]] + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    }
    
    if (!is.null(manual_plots_parameters$features)) {
      if ((rna == manual_plots_parameters$features) & (tissue %in% manual_plots_parameters$tissues)) {
        #print(tissue)
        #print(sprintf("fc_log2: %s", fc))
        #print(sprintf("fc_th_log2: %s", log2(params$fc_th)))
        xpos = 1.1*tmp
        ypos = 1.1*maxi_expr_rna
        p_violin_manual[[tissue]] = p_violin[[j]] + 
          #geom_text(x = xpos, y = ypos-0.5, label = corr_text, lineheight = 0.75, hjust = 1, vjust = 1, family = plots_props$font_family, size = 6/font_size_correction, color = "black") +
          #ylim(0, ypos) + 
          scale_y_continuous(limits = c(0, ypos), expand = c(0,0)) +
          ylab("") +
          xlab("")
      }
    }

    #xpos = 1.05*tmp
    #ypos = 1.1*max(plot_df$value)
    #p_violin[[j]] = p_violin[[j]] + geom_text(x = xpos, y = ypos, label = corr_text, lineheight = 0.75, hjust = 1, vjust = 1, family = plots_props$font_family, size = 6/font_size_correction, color = "black")
    
    if (rna %in% feature_config) {
      p_violin_group = p_violin[[j]] + 
        ylab(y_axis_title_small) + 
        xlab(x_axis_title) 
      
      ggsave(sprintf("%s/%s.png", output_folder_violin_single, tissue), p_violin_group, dpi=plots_props$dpi, width=plots_props$image_width / 2, height=(plots_props$image_height*1.5)/2, units = plots_props$image_units)
      ggsave(sprintf("%s/%s.svg", output_folder_violin_single, tissue), p_violin_group, width=plots_props$image_width / 2, height=(plots_props$image_height*1.5)/2, units = plots_props$image_units, limitsize = FALSE, scale = 1)
    } 
    
    # create box plot
    p_box[[j]] = ggplot(plot_df, aes(x = comp, y = value, group = comp, fill = comp)) + 
      geom_boxplot(width = 0.9, outlier.colour="grey", outlier.shape=16, outlier.size=1, notch=FALSE) +
      scale_fill_manual(values = unlist(colours[[comp_group]])) +
      ggtitle(plot_title) +
      xlab("") +
      ylab("") +
      #ylab(sprintf("expr (%s, detection rate = %s%%\nper %s)", data_input$norm, expr_input$detection_rate, params$plot_group)) +
      scale_x_discrete(labels = unlist(xticks_names[[comp_group]])) +
      scale_y_continuous(limits = c(0, 1.1*max(plot_df$value)), expand = c(0,0)) +
      geom_text(aes(label = significance), vjust=-0.3, size=3) +
      theme_classic() +
      theme(#plot.margin = unit(c(0.1,-1.9,0.25,-1.75), plots_props$image_units), #top, right, bottom, left
        plot.margin = unit(c(0,0.1,0,0.1), plots_props$image_units), #top, right, bottom, left
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9), 
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position="none", 
        text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
        #plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
        plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size, margin = margin(b=0, unit="cm")), # Using pt for margins      
        #legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size)
      ) + 
      coord_cartesian(clip = "off")
    
    if (x_axis_ticks_rot) {
      p_box[[j]] = p_box[[j]] + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    }
    
    if (!is.null(manual_plots_parameters$features)) {
      if ((rna == manual_plots_parameters$features) & (tissue %in% manual_plots_parameters$tissues)) {
        xpos = 1.1*tmp
        ypos = 1.1*maxi_expr_rna
        p_box_manual[[tissue]] = p_box[[j]] + 
          #geom_text(x = xpos, y = ypos, label = corr_text, lineheight = 0.75, hjust = 1, vjust = 1, family = plots_props$font_family, size = 6/font_size_correction) +
          scale_y_continuous(limits = c(0, ypos), expand = c(0,0)) +
          ylab("") +
          xlab("")
      }
    }
    
    #xpos = 1.05*tmp
    #ypos = 1.1*max(plot_df$value) 
    #p_box[[j]] = p_box[[j]] + geom_text(x = xpos, y = ypos, label = corr_text, lineheight = 0.75, hjust = 1, vjust = 1, family = plots_props$font_family, size = 6/font_size_correction)
    
    if (!is.null(manual_plots_parameters$features)) {
      if ((rna %in% feature_config) | (rna == manual_plots_parameters$features)){
        p_box_group = p_box[[j]] + 
          ylab(y_axis_title_small) + 
          xlab(x_axis_title) 
  
        ggsave(sprintf("%s/%s.png", output_folder_box_single, tissue), p_box_group, dpi=plots_props$dpi, width=plots_props$image_width / 2, height=(plots_props$image_height*1.5)/2, units = plots_props$image_units)
        ggsave(sprintf("%s/%s.svg", output_folder_box_single, tissue), p_box_group, width=plots_props$image_width / 2, height=(plots_props$image_height*1.5)/2, units =plots_props$image_units, limitsize = FALSE, scale = 1)
      }
    }
    
    if (tissue != "all") {
      j = j + 1
    }
  }
  
  return(list("violin"=p_violin, "violin_manual"=p_violin_manual, "box"=p_box, "box_manual"=p_box_manual))
}

rep_along = function(x, values) {
  rep(values, length.out = length(x))
}


#---------------------------------- Load data ---------------------------------- 
data_input = snakemake@params$data_input
input_files = snakemake@params$files
feature_config = snakemake@params$features
comparisons = snakemake@params$comparisons
params = snakemake@params$parameters_props 
manual_plots_parameters = snakemake@params$manual_plots_parameters
colours = snakemake@params$colors
xticks_names = snakemake@params$xticks_names
plots_props = snakemake@params$plots_props
results_folder = snakemake@params$results_folder

expr = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
parts = strsplit(data_input$data_sub_set, "_")[[1]]
suffix = paste(parts[-1], collapse = "_")
raw_counts = fread(sprintf("%s_%s/%s_%s_expression_raw_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, suffix), sep='\t', header=T)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 

diff_exp = fread(sprintf("%s/%s_%s/results_%s/matrices/diff_exp/diff_exp%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, params$diff_log), sep='\t') 


#------------------------------------ Script ----------------------------------- 
feature_list = c()
#for (input_file in input_files) {
#  feature_list = c(feature_list, fread(sprintf(input_file, data_input$data_sub_set))[[data_input$rna_class]])
#}
feature_list = unique(c(feature_list, feature_config))

#for (m in params$corr_methods) {
#  print(m)
for(i in 1:length(comparisons)) {
  # which comparisons do we investigate and which groups are considered
  comparison = comparisons[[i]]
  #print(comparison$comp)
  
  order_of_x_axis = comparison$order_of_x_axis
  
  plot_group = comparison$comp[1]
  comp_group = comparison$comp[2]
  
  plot_nrows = comparison$plot_layout[1]
  plot_ncols = comparison$plot_layout[2]
  
  plot_height = comparison$plot_size[1] 
  plot_width = comparison$plot_size[2]
  
  # check if files exist
  #print(file.exists(sprintf("%s/%s_%s/results_%s/matrices/aggregation_corr/corr_method=%s_%ss_with_%s_per_%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m, data_input$rna_class, comp_group, plot_group)))
  #print(file.exists(sprintf("%s/%s_%s/results_%s/matrices/aggregation_corr/corr_method=%s_adj_pvalues_%ss_with_%s_per_%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m, data_input$rna_class, comp_group, plot_group)))
  
  # old
  #corr = fread(sprintf("%s/%s_%s/results_%s/matrices/aggregation_corr/corr_method=%s_%ss_with_%s_per_%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m, data_input$rna_class, comp_group, params$plot_group), sep='\t') 
  #p_value = fread(sprintf("%s/%s_%s/results_%s/matrices/aggregation_corr/corr_method=%s_adj_pvalues_%ss_with_%s_per_%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m, data_input$rna_class, comp_group, params$plot_group), sep='\t') 
  
  #corr = fread(sprintf("%s/%s_%s/results_%s/matrices/aggregation_corr/corr_method=%s_%ss_with_%s_per_%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m, data_input$rna_class, comp_group, plot_group), sep='\t', header = TRUE, check.names = FALSE) 
  #p_value = fread(sprintf("%s/%s_%s/results_%s/matrices/aggregation_corr/corr_method=%s_adj_pvalues_%ss_with_%s_per_%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m, data_input$rna_class, comp_group, plot_group), sep='\t', header = TRUE, check.names = FALSE) 
  
  #colnames(corr) = gsub("\\.", " ", colnames(corr))
  #colnames(p_value) = gsub("\\.", " ", colnames(p_value))
  
  # create output folders
  output_folder_violin = sprintf("%s/%s_%s/results_%s/figures/expr/violin_plots/median_%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, plot_group)
  dir.create(output_folder_violin, recursive=TRUE)
  output_folder_box = sprintf("%s/%s_%s/results_%s/figures/expr/box_plots/median_%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, plot_group)
  dir.create(output_folder_box, recursive=TRUE)
  
  n_cores = 64
  ## print(sprintf("running in %s threads", n_cores))
  ## FOR PARALLEL EXECUTION
  #cl <- makeCluster(n_cores, type="FORK")
  #registerDoParallel(cl)
  ##
  
  #print(length(feature_list))
  for (rna in feature_list) {
  ## FOR PARALLEL EXECUTION
  #foreach (rna=feature_list) %dopar% {

    print(rna)
    number_of_plots_possible = plot_nrows * plot_ncols
    
    if ((rna %in% feature_config) | (rna %in% manual_plots_parameters$features)) {
      output_folder_violin_single = sprintf("%s/%s_%s/results_%s/figures/expr/violin_plots/median_%s_single/%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, plot_group, rna)
      dir.create(output_folder_violin_single, recursive=TRUE)
      output_folder_box_single = sprintf("%s/%s_%s/results_%s/figures/expr/box_plots/median_%s_single/%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, plot_group, rna)
      dir.create(output_folder_box_single, recursive=TRUE)
    } else {
      output_folder_violin_single = ""
      output_folder_box_single = ""
    }
    
    expr_feature = expr[expr[[data_input$rna_class]] == rna,]
    diff_exp_feature = diff_exp[diff_exp[[data_input$rna_class]] == rna,]
    
    if (rna %in% manual_plots_parameters$features) {
      tmp = (expr_feature != rna)
      expr_feature_tmp = expr_feature[,..tmp]
      #samples_for_plot = annot[annot[[params$plot_group]] %in% manual_plots_parameters$tissues][[data_input$identifier_column]]
      samples_for_plot = annot[annot[[plot_group]] %in% manual_plots_parameters$tissues][[data_input$identifier_column]]
      maxi_expr_rna = max(expr_feature_tmp[,..samples_for_plot])
    }
    
    # axis titles
    x_axis_title = xticks_names$categories[[sprintf("%s_capital", comp_group)]]
    if (grepl(" ", xticks_names$categories[[data_input$detection_group]]) == TRUE) {
      detection_group = strsplit(xticks_names$categories[[data_input$detection_group]], " ")[[1]][1]
    } else {
      detection_group = data_input$detection_group
    }
    y_axis_title_small = sprintf("Expr (%s,\ndetection rate = %s%%\nper %s)", gsub("_norm", "" ,data_input$norm), gsub("p", "", data_input$detection_rate), detection_group)
    if (plot_height <= 9) {
      y_axis_title_big = y_axis_title_small
    } else {
      y_axis_title_big = sprintf("Expr (%s, detection rate = %s%% per %s)", gsub("_norm", "" ,data_input$norm), gsub("p", "", data_input$detection_rate), detection_group)
    }
    
    p = ploting_violin(expr_feature, annot, diff_exp_feature, rna, feature_config, plot_group, comp_group, params, order_of_x_axis, data_input, x_axis_title, y_axis_title_big, y_axis_title_small, xticks_names, plots_props, output_folder_violin_single, output_folder_box_single)
    
    # creating concatinated plot
    all_p_violin = arrangeGrob(grobs=p$violin, ncol=plot_ncols, nrow=plot_nrows, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title_big, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
    all_p_box = arrangeGrob(grobs=p$box, ncol=plot_ncols, nrow=plot_nrows, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title_big, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
    
    # save violin plot
    #ggsave(sprintf("%s/%s.png", output_folder_violin, rna), all_p_violin, dpi=plots_props$dpi, width=plot_width, height=plot_height, units = plots_props$image_units)
    #ggsave(sprintf("%s/%s.svg", output_folder_violin, rna), all_p_violin, width=plot_width, height=plot_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
    
    # save box plot
    #ggsave(sprintf("%s/%s.png", output_folder_box, rna), all_p_box, dpi=plots_props$dpi, width=plot_width, height=plot_height, units = plots_props$image_units)
    #ggsave(sprintf("%s/%s.svg", output_folder_box, rna), all_p_box, width=plot_width, height=plot_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
    
    # manual part 
    if (!is.null(manual_plots_parameters$features)) {
      if ((rna == manual_plots_parameters$features) & (comparison$comp[1] == "brain_region")){
        # sort the figures given the order from the config
        p_violin_manual = p$violin_manual[manual_plots_parameters$tissues]
        p_box_manual = p$box_manual[manual_plots_parameters$tissues]
        
        plot_nrows_manual = manual_plots_parameters$plot_layout[1]
        plot_ncols_manual = manual_plots_parameters$plot_layout[2]
        
        plot_height_manual = manual_plots_parameters$plot_size[1]
        plot_width_manual = manual_plots_parameters$plot_size[2]
        
        # determine which figures are not on the left side
        left_col = c()
        for (row in 1:plot_nrows_manual - 1) {
          left_col = append(left_col, plot_ncols_manual*row + 1)
        }
        # determine which figures are not on the left side
        right_col = c()
        for (row in 1:plot_nrows_manual) {
          right_col = append(right_col, plot_ncols_manual*row)
        }
        
        #f (data_input$data_sub_set == "CA2_Experiment=diet restriction_quantity_control") {
        top_margin_left_col = 7 #!
        right_margin_left_col = -4 # -9
        bottom_margin_left_col = -10 #!
        left_margin_left_col = -5 #!
        
        top_margin_middle_cols = top_margin_left_col
        right_margin_middle_cols = -1
        bottom_margin_middle_cols = bottom_margin_left_col
        left_margin_middle_cols = -8
        
        top_margin_right_col = top_margin_left_col
        right_margin_right_col = 2 #!
        bottom_margin_right_col = bottom_margin_left_col
        left_margin_right_col = -11 
        #}
        
        p_violin_manual_mod_margins = p_violin_manual
        p_box_manual_mod_margins = p_box_manual
        for (img_i in 1:length(p_violin_manual_mod_margins)) {
          if ((!(img_i %in% left_col)) & (!(img_i %in% right_col))) {
            p_violin_manual_mod_margins[[img_i]] = p_violin_manual_mod_margins[[img_i]] + 
              theme(plot.margin = margin(top_margin_middle_cols, right_margin_middle_cols, bottom_margin_middle_cols, left_margin_middle_cols, "pt")
              )
            p_box_manual_mod_margins[[img_i]] = p_box_manual_mod_margins[[img_i]] + 
              theme(plot.margin = margin(top_margin_middle_cols, right_margin_middle_cols, bottom_margin_middle_cols, left_margin_middle_cols, "pt")
              )
          }
          else if (img_i %in% left_col) {
            p_violin_manual_mod_margins[[img_i]] = p_violin_manual_mod_margins[[img_i]] + 
              theme(plot.margin = margin(top_margin_left_col, right_margin_left_col, bottom_margin_left_col, left_margin_left_col, "pt")
              )
            p_box_manual_mod_margins[[img_i]] = p_box_manual_mod_margins[[img_i]] + 
              theme(plot.margin = margin(top_margin_left_col, right_margin_left_col, bottom_margin_left_col, left_margin_left_col, "pt")
              )
          } else if (img_i %in% right_col) {
            p_violin_manual_mod_margins[[img_i]] = p_violin_manual_mod_margins[[img_i]] + 
              theme(plot.margin = margin(top_margin_right_col, right_margin_right_col, bottom_margin_right_col, left_margin_right_col, "pt")
              )
            p_box_manual_mod_margins[[img_i]] = p_box_manual_mod_margins[[img_i]] + 
              theme(plot.margin = margin(top_margin_right_col, right_margin_right_col, bottom_margin_right_col, left_margin_right_col, "pt")
              )
          }
          if (!(img_i %in% left_col)) {
            p_violin_manual_mod_margins[[img_i]] = p_violin_manual_mod_margins[[img_i]] + 
              scale_y_continuous(labels = function(breaks) {rep_along(breaks, "  ")}, limits = c(0, 1.15*maxi_expr_rna), expand = c(0,0)) +
              theme(axis.ticks.y = element_line(color = "transparent"),
                    axis.text.y = element_text(color = "transparent")
              )
            p_box_manual_mod_margins[[img_i]] = p_box_manual_mod_margins[[img_i]] + 
              #scale_y_continuous(labels = function(breaks) {rep_along(breaks, "  ")}, limits = c(0, 1.15*maxi_expr_rna), expand = c(0,0)) +
              theme(axis.ticks.y = element_line(color = "transparent"),
                    axis.text.y = element_text(color = "transparent")
              )
          }
        }
        
        all_p_violin_manual = arrangeGrob(grobs=p_violin_manual_mod_margins, ncol=plot_ncols_manual, nrow=plot_nrows_manual, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title_big, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90), right = textGrob(rna, gp=gpar(fontsize=plots_props$font_size_header, fontfamily=plots_props$font_family), rot = 270))
        # = arrangeGrob(grobs=p_box_manual_mod_margins, ncol=plot_ncols_manual, nrow=plot_nrows_manual, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title_big, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90), right = textGrob(rna, gp=gpar(fontsize=plots_props$font_size_header, fontfamily=plots_props$font_family), rot = 270))
        
        ggsave(sprintf("%s/%s.png", output_folder_violin_single, paste(manual_plots_parameters$tissues, collapse = "_")), all_p_violin_manual, dpi = plots_props$dpi, width = plot_width_manual, height = plot_height_manual, units = plots_props$image_units) #width=manual_plots_parameters$plot_ncols*(plots_props$image_width/2), height=(plots_props$image_height*1.5)/2
        ggsave(sprintf("%s/%s.svg", output_folder_violin_single, paste(manual_plots_parameters$tissues, collapse = "_")), all_p_violin_manual, width = plot_width_manual, height = plot_height_manual, units = plots_props$image_units, limitsize = FALSE, scale = 1) #width=manual_plots_parameters$plot_ncols*(plots_props$image_width/2), height=(plots_props$image_height*1.5)/2
        
        #ggsave(sprintf("%s/%s.png", output_folder_box_single, paste(manual_plots_parameters$tissues, collapse = "_")), all_p_box_manual, dpi = plots_props$dpi, width = plot_width_manual, height = plot_height_manual, units = plots_props$image_units) #width=manual_plots_parameters$plot_ncols*(plots_props$image_width/2), height=(plots_props$image_height*1.5)/2
        #ggsave(sprintf("%s/%s.svg", output_folder_box_single, paste(manual_plots_parameters$tissues, collapse = "_")), all_p_box_manual, width = plot_width_manual, height = plot_height_manual, units = plots_props$image_units, limitsize = FALSE, scale = 1) #width=manual_plots_parameters$plot_ncols*(plots_props$image_width/2), height=(plots_props$image_height*1.5)/2
        
        #rm(all_p_violin_manual)
        #rm(all_p_box_manual)
        
        # if (manual_plots_parameters$plot_nrows == 1){
        #   p_violin_manual_adj = list()
        #   p_box_manual_adj = list()
        #   for (ti in manual_plots_parameters$tissues) {
        #     #if (ti != manual_plots_parameters$tissues[1]){
        #     if (ti == manual_plots_parameters$tissues[length(manual_plots_parameters$tissues)]) {
        #       print(ti)
        #       p_violin_manual_adj[[ti]] = p_bar_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,0.25,-0.5,-0.1), "cm"))
        #       p_box_manual_adj[[ti]] = p_box_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,0.25,-0.5,-0.1), "cm"))
        #       #} else if (ti == manual_plots_parameters$tissues[length(manual_plots_parameters$tissues)]) {
        #     } else if (ti != manual_plots_parameters$tissues[1]){
        #       p_violin_manual_adj[[ti]] = p_bar_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,-0.1,-0.5,-0.1), "cm"))
        #       p_box_manual_adj[[ti]] = p_box_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,-0.1,-0.5,-0.1), "cm"))
        #     } else {
        #       p_violin_manual_adj[[ti]] = p_bar_manual[[ti]] + theme(plot.margin = unit(c(0,-0.1,-0.5,-0.35), "cm"))
        #       p_box_manual_adj[[ti]] = p_box_manual[[ti]] + theme(plot.margin = unit(c(0,-0.1,-0.5,-0.35), "cm"))
        #     } 
        #     #p_violin_manual_adj[[ti]] = p_violin_manual_adj[[ti]] + theme(plot.margin = unit(c(0,-0.15,-0.5,-0.15), "cm"))
        #     #p_box_manual_adj[[ti]] = p_box_manual_adj[[ti]] + theme(plot.margin = unit(c(0,-0.15,-0.5,-0.15), "cm"))
        #   }
        #   
        #   all_p_violin_manual = arrangeGrob(grobs=p_violin_manual_adj, ncol=manual_plots_parameters$plot_ncols, nrow=manual_plots_parameters$plot_nrows, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title_small, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
        #   all_p_box_manual = arrangeGrob(grobs=p_box_manual_adj, ncol=manual_plots_parameters$plot_ncols, nrow=manual_plots_parameters$plot_nrows, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title_small, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
        #   
        #   width = 2 * plots_props$image_width
        #   height = (plots_props$image_height*1.5)/2
        #   ggsave(sprintf("%s/%s_%sx%s.png", output_folder_bar_single_rna, paste(manual_plots_parameters$tissues, collapse = "_"), width, height), all_p_violin_manual, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
        #   ggsave(sprintf("%s/%s_%sx%s.svg", output_folder_bar_single_rna, paste(manual_plots_parameters$tissues, collapse = "_"), width, height), all_p_violin_manual, width=width, height=height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
        #   
        #   ggsave(sprintf("%s/%s_%sx%s.png", output_folder_box_single_rna, paste(manual_plots_parameters$tissues, collapse = "_"), width, height), all_p_box_manual, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
        #   ggsave(sprintf("%s/%s_%sx%s.svg", output_folder_box_single_rna, paste(manual_plots_parameters$tissues, collapse = "_"), width, height), all_p_box_manual, width=width, height=height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
        # }
      }
    }
    
    # cleanup an call garbage cleaner
    rm(all_p_violin)
    rm(all_p_box)
    rm(all_p_violin_manual)
    rm(all_p_box_manual)
    rm(p)
    
    
    
    # raw counts
    # create output folders
    output_folder_violin_raw = sprintf("%s/%s_%s/results_%s/figures/expr/violin_plots/median_%s_raw", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, plot_group)
    dir.create(output_folder_violin_raw, recursive=TRUE)
    output_folder_box_raw = sprintf("%s/%s_%s/results_%s/figures/expr/box_plots/median_%s_raw", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, plot_group)
    dir.create(output_folder_box_raw, recursive=TRUE)
    output_folder_violin_raw_single = sprintf("%s/%s_%s/results_%s/figures/expr/violin_plots/median_%s_raw_single/%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, plot_group, rna)
    dir.create(output_folder_violin_raw_single, recursive=TRUE)
    output_folder_box_raw_single = sprintf("%s/%s_%s/results_%s/figures/expr/box_plots/median_%s_raw_single/%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, plot_group, rna)
    dir.create(output_folder_box_raw_single, recursive=TRUE)
    
    raw_counts_feature = raw_counts[raw_counts[[data_input$rna_class]] == rna,]
    
    # axis titles
    x_axis_title = xticks_names$categories[[sprintf("%s_capital", comp_group)]]
    y_axis_title_big = "Raw counts"
    y_axis_title_small = "Raw counts"
    
    p_raw = ploting_violin(raw_counts_feature, annot, diff_exp_feature, rna, feature_config, plot_group, comp_group, params, order_of_x_axis, data_input, x_axis_title, y_axis_title_big, y_axis_title_small, xticks_names, plots_props, output_folder_violin_raw_single, output_folder_box_raw_single)
      
    # creating concatinated plot
    all_p_violin_raw = arrangeGrob(grobs=p_raw$violin, ncol=plot_ncols, nrow=plot_nrows, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title_big, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
    all_p_box_raw = arrangeGrob(grobs=p_raw$box, ncol=plot_ncols, nrow=plot_nrows, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title_big, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
    
    # save violin plot
    ggsave(sprintf("%s/%s.png", output_folder_violin_raw, rna), all_p_violin_raw, dpi=plots_props$dpi, width=plot_width, height=plot_height, units = plots_props$image_units)
    ggsave(sprintf("%s/%s.svg", output_folder_violin_raw, rna), all_p_violin_raw, width=plot_width, height=plot_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
    
    # save box plot
    ggsave(sprintf("%s/%s.png", output_folder_box_raw, rna), all_p_box_raw, dpi=plots_props$dpi, width=plot_width, height=plot_height, units = plots_props$image_units)
    ggsave(sprintf("%s/%s.svg", output_folder_box_raw, rna), all_p_box_raw, width=plot_width, height=plot_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
    
    # manual part 
    if (!is.null(manual_plots_parameters$features)) {
      if (rna == manual_plots_parameters$features) {
        # sort the figures given the order from the config
        p_violin_raw_manual = p_raw$violin_manual[manual_plots_parameters$tissues]
        p_box_raw_manual = p_raw$box_manual[manual_plots_parameters$tissues]
        
        plot_nrows_manual = manual_plots_parameters$plot_layout[1]
        plot_ncols_manual = manual_plots_parameters$plot_layout[2]
        
        plot_height_manual = manual_plots_parameters$plot_size[1]
        plot_width_manual = manual_plots_parameters$plot_size[2]
        
        all_p_violin_raw_manual = arrangeGrob(grobs=p_violin_raw_manual, ncol=plot_ncols_manual, nrow=plot_nrows_manual, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title_small, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
        all_p_box_raw_manual = arrangeGrob(grobs=p_box_raw_manual, ncol=plot_ncols_manual, nrow=plot_nrows_manual, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title_small, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
        
        ggsave(sprintf("%s/%s.png", output_folder_violin_raw_single, paste(manual_plots_parameters$tissues, collapse = "_")), all_p_violin_raw_manual, dpi = plots_props$dpi, width = plot_width_manual, height = plot_height_manual, units = plots_props$image_units) #width=manual_plots_parameters$plot_ncols*(plots_props$image_width/2), height=(plots_props$image_height*1.5)/2
        ggsave(sprintf("%s/%s.svg", output_folder_violin_raw_single, paste(manual_plots_parameters$tissues, collapse = "_")), all_p_violin_raw_manual, width = plot_width_manual, height = plot_height_manual, units = plots_props$image_units, limitsize = FALSE, scale = 1) #width=manual_plots_parameters$plot_ncols*(plots_props$image_width/2), height=(plots_props$image_height*1.5)/2
        
        ggsave(sprintf("%s/%s.png", output_folder_box_raw_single, paste(manual_plots_parameters$tissues, collapse = "_")), all_p_box_raw_manual, dpi = plots_props$dpi, width = plot_width_manual, height = plot_height_manual, units = plots_props$image_units) #width=manual_plots_parameters$plot_ncols*(plots_props$image_width/2), height=(plots_props$image_height*1.5)/2
        ggsave(sprintf("%s/%s.svg", output_folder_box_raw_single, paste(manual_plots_parameters$tissues, collapse = "_")), all_p_box_raw_manual, width = plot_width_manual, height = plot_height_manual, units = plots_props$image_units, limitsize = FALSE, scale = 1) #width=manual_plots_parameters$plot_ncols*(plots_props$image_width/2), height=(plots_props$image_height*1.5)/2
        
        rm(all_p_violin_raw_manual)
        rm(all_p_box_raw_manual)
        
        # if (manual_plots_parameters$plot_nrows == 1){
        #   p_violin_manual_adj = list()
        #   p_box_manual_adj = list()
        #   for (ti in manual_plots_parameters$tissues) {
        #     #if (ti != manual_plots_parameters$tissues[1]){
        #     if (ti == manual_plots_parameters$tissues[length(manual_plots_parameters$tissues)]) {
        #       print(ti)
        #       p_violin_manual_adj[[ti]] = p_bar_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,0.25,-0.5,-0.1), "cm"))
        #       p_box_manual_adj[[ti]] = p_box_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,0.25,-0.5,-0.1), "cm"))
        #       #} else if (ti == manual_plots_parameters$tissues[length(manual_plots_parameters$tissues)]) {
        #     } else if (ti != manual_plots_parameters$tissues[1]){
        #       p_violin_manual_adj[[ti]] = p_bar_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,-0.1,-0.5,-0.1), "cm"))
        #       p_box_manual_adj[[ti]] = p_box_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,-0.1,-0.5,-0.1), "cm"))
        #     } else {
        #       p_violin_manual_adj[[ti]] = p_bar_manual[[ti]] + theme(plot.margin = unit(c(0,-0.1,-0.5,-0.35), "cm"))
        #       p_box_manual_adj[[ti]] = p_box_manual[[ti]] + theme(plot.margin = unit(c(0,-0.1,-0.5,-0.35), "cm"))
        #     } 
        #     #p_violin_manual_adj[[ti]] = p_violin_manual_adj[[ti]] + theme(plot.margin = unit(c(0,-0.15,-0.5,-0.15), "cm"))
        #     #p_box_manual_adj[[ti]] = p_box_manual_adj[[ti]] + theme(plot.margin = unit(c(0,-0.15,-0.5,-0.15), "cm"))
        #   }
        #   
        #   all_p_violin_manual = arrangeGrob(grobs=p_violin_manual_adj, ncol=manual_plots_parameters$plot_ncols, nrow=manual_plots_parameters$plot_nrows, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title_small, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
        #   all_p_box_manual = arrangeGrob(grobs=p_box_manual_adj, ncol=manual_plots_parameters$plot_ncols, nrow=manual_plots_parameters$plot_nrows, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title_small, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
        #   
        #   width = 2 * plots_props$image_width
        #   height = (plots_props$image_height*1.5)/2
        #   ggsave(sprintf("%s/%s_%sx%s.png", output_folder_bar_single_rna, paste(manual_plots_parameters$tissues, collapse = "_"), width, height), all_p_violin_manual, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
        #   ggsave(sprintf("%s/%s_%sx%s.svg", output_folder_bar_single_rna, paste(manual_plots_parameters$tissues, collapse = "_"), width, height), all_p_violin_manual, width=width, height=height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
        #   
        #   ggsave(sprintf("%s/%s_%sx%s.png", output_folder_box_single_rna, paste(manual_plots_parameters$tissues, collapse = "_"), width, height), all_p_box_manual, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
        #   ggsave(sprintf("%s/%s_%sx%s.svg", output_folder_box_single_rna, paste(manual_plots_parameters$tissues, collapse = "_"), width, height), all_p_box_manual, width=width, height=height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
        # }
      }
    }
    
    # cleanup an call garbage cleaner
    rm(all_p_violin_raw)
    rm(all_p_box_raw)
    rm(all_p_violin_raw_manual)
    rm(all_p_box_raw_manual)
    rm(p_raw)
    
    gc()
  }
  ## FOR PARALLEL EXECUTION
  #stopCluster(cl)
  ##
}
#}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
