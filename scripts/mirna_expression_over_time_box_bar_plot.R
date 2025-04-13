suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(foreach))  # for parallelization
suppressPackageStartupMessages(library(doParallel))  # for parallelization

#snakemake = readRDS("snakemake_mirna_expression_over_time_box_bar_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("mirna_expression_over_time_box_bar_plot")


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
input_files = snakemake@params$files
feature_config = snakemake@params$features
manual_plots_parameter_list = snakemake@params$manual_plots_parameters
params = snakemake@params$parameters 
corr_methods = snakemake@params$corr_methods
plots_props = snakemake@params$plots_props
xticks_names = snakemake@params$xticks_names
colors = snakemake@params$colors
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_mirna_expression_over_time_box_bar_plot.rds")
#stop()

expr = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
diff_exp = fread(sprintf("%s/%s_%s/results_%s/matrices/diff_exp/diff_exp%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, params$diff_log), sep='\t') 

#features = fread(snakemake@input$feature_list, sep='\t')[[data_input$rna_class]]  # , colClasses=c(ID="character")

poly_degs = c(2,3,0)


#------------------------------------ Script ----------------------------------- 
manual_plots_parameter_mirna_list = unique(unname(sapply(manual_plots_parameter_list, function(x) x$feature)))

feature_list = c()
for (input_file in input_files) {
  feature_list = c(feature_list, fread(input_file)[[data_input$feature_column]])
}
if (is.null(feature_config)) {
  if (manual_plots_parameter_mirna_list == "") {
    if (is.null(feature_list)) {
      print("No features selected. Stop!")
      break
    } else {
    feature_list = unique(c(feature_list))
    }
  } else {
    if (is.null(feature_list)) {
      feature_list = unique(c(manual_plots_parameter_mirna_list))
    } else {
      feature_list = unique(c(feature_list, manual_plots_parameter_mirna_list))
    }
  }
} else {
  if (manual_plots_parameter_mirna_list == "") {
    if (is.null(feature_list)) {
      feature_list = unique(c(feature_config))
    } else {
      feature_list = unique(c(feature_list, feature_config))
    }
  }
  if (is.null(feature_list)) {
    feature_list = unique(c(feature_config, manual_plots_parameter_mirna_list))
  } else {
    feature_list = unique(c(feature_list, feature_config, manual_plots_parameter_mirna_list))
  }
}

for (i in 1:length(manual_plots_parameter_list)) {
  manual_plots_parameters = manual_plots_parameter_list[[i]]
  
  #feature_list = unique(c(manual_plots_parameters$features, feature_config))
  #print(feature_list)
  print(sprintf("features to be processed: %s", length(feature_list)))
  
  for (m in corr_methods) {
    for (poly_deg in poly_degs) {
      # load correlation values
      corr = fread(sprintf("%s/%s_%s/results_%s/matrices/aggregation_corr/corr_method=%s_%ss_with_%s_per_%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m, data_input$rna_class, params$comb_group, params$plot_group), sep='\t') 
      p_value = fread(sprintf("%s/%s_%s/results_%s/matrices/aggregation_corr/corr_method=%s_adj_pvalues_%ss_with_%s_per_%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m, data_input$rna_class, params$comb_group, params$plot_group), sep='\t') 
      
      # create output folders
      output_folder_bar = sprintf("%s/%s_%s/results_%s/figures/expr/bar_plots/median_time_corr_method=%s_poly_deg=%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m, poly_deg)
      dir.create(output_folder_bar, recursive=TRUE)
      output_folder_box = sprintf("%s/%s_%s/results_%s/figures/expr/box_plots/time_corr_method=%s_poly_deg=%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m, poly_deg)
      dir.create(output_folder_box, recursive=TRUE)
      output_folder_bar_single = sprintf("%s/%s_%s/results_%s/figures/expr/bar_plots/median_time_single_corr_method=%s_poly_deg=%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m, poly_deg)
      dir.create(output_folder_bar_single, recursive=TRUE)
      output_folder_box_single = sprintf("%s/%s_%s/results_%s/figures/expr/box_plots/time_single_corr_method=%s_poly_deg=%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m, poly_deg)
      dir.create(output_folder_box_single, recursive=TRUE)
       
      # Set up parallel computing
      n_cores = 64
      # print(sprintf("running in %s threads", n_cores))
      #cl <- makeCluster(n_cores, type="FORK")
      #registerDoParallel(cl)
      
      #print(length(feature_list))
      for (rna in feature_list) {
      #foreach (rna=feature_list) %dopar% {
        print(rna)
        expr_feature = expr[expr[[data_input$feature_column]] == rna,]
        diff_exp_feature = diff_exp[diff_exp[[data_input$rna_class]] == rna,]
        
        #number_of_plots_possible = params$plot_nrows * params$plot_ncols
        
        if (rna %in% feature_list) {
          output_folder_bar_single_rna = sprintf("%s/%s", output_folder_bar_single, rna)
          dir.create(output_folder_bar_single_rna, recursive=TRUE)
          output_folder_box_single_rna = sprintf("%s/%s", output_folder_box_single, rna)
          dir.create(output_folder_box_single_rna, recursive=TRUE)
        }
        if (rna == manual_plots_parameters$features) {
          output_folder_bar_single_rna = sprintf("%s/%s", output_folder_bar_single, rna)
          dir.create(output_folder_bar_single_rna, recursive=TRUE)
          output_folder_box_single_rna = sprintf("%s/%s", output_folder_box_single, rna)
          dir.create(output_folder_box_single_rna, recursive=TRUE)
          
          tmp = (expr_feature != rna)
          expr_feature_tmp = expr_feature[,..tmp]
          samples_for_plot = annot[annot[[params$plot_group]] %in% manual_plots_parameters$tissues][[data_input$identifier_column]]
          maxi_expr_rna = max(expr_feature_tmp[,..samples_for_plot])
        }
        
        tmp = unique(names(colors[[params$plot_group]])) %in% unique(annot$brain_region)
        tissue_list = unique(names(colors[[params$plot_group]]))[tmp]
  
        p_bar = list()
        p_box = list()
        p_bar_manual = list()
        p_box_manual = list()
        i = 1
        for (tissue in tissue_list){
          id_group = annot[annot[[params$plot_group]] == tissue,][[data_input$identifier_column]] 
          tmp = colnames(expr_feature) %in% id_group
          expr_feature_tissue = expr_feature[,..tmp]
          
          plot_df = melt(expr_feature_tissue)
          plot_df[[params$plot_group]] = str_split_fixed(plot_df$variable, "_", 3)[,2]
          
          time_list = c()
          for (id in plot_df$variable) {
            time_list = append(time_list, annot[annot[[data_input$identifier_column]] == id,][[params$comb_group]])
          }
          plot_df[[params$comb_group]] = as.numeric(time_list)
      
          time_median = c()
          time_sd = c()
          significance = c()
          for (time in unique(annot[[params$comb_group]])) {
            time_median = append(time_median, median(plot_df[plot_df[[params$comb_group]] == time,]$value))
            time_sd = append(time_sd, sd(plot_df[plot_df[[params$comb_group]] == time,]$value))
            
            if (time != min(as.numeric(unique(annot[[params$comb_group]])))) {
              ttest_adj_pval = diff_exp_feature[[sprintf("ttest_adjp_%s__%s_vs_%s_%s=%s", params$comb_group, time, sort(as.numeric(unique(annot[[params$comb_group]])))[1], params$plot_group, tissue)]]  #ttest_adjp_age__12_vs_3_brain_region_fixed=cc
              fc = log2(diff_exp_feature[[sprintf("fc_%s__%s_vs_%s_%s=%s", params$comb_group, time, sort(as.numeric(unique(annot[[params$comb_group]])))[1], params$plot_group, tissue)]])  #ttest_adjp_yy__12_vs_3_brain_region_fixed=cc
              if (is.na(fc)) {
                significance = append(significance, "")
              } else if (abs(fc) >= log2(params$fc_th)) {
                if (is.na(ttest_adj_pval)) { 
                  #significance = append(significance, "+/-")
                  significance = append(significance, "")
                } else if (ttest_adj_pval < 0.001) { 
                  #significance = append(significance, "***\n+/-")
                  significance = append(significance, "***")
                } else if (ttest_adj_pval < 0.01) {
                  #significance = append(significance, "**\n+/-")
                  significance = append(significance, "**")
                } else if (ttest_adj_pval < 0.05) {
                  #significance = append(significance, "*\n+/-")
                  significance = append(significance, "*")
                } else {
                  #significance = append(significance, "+/-")
                  significance = append(significance, "")
                }
              } else {
                significance = append(significance, "")
              }
            } else {
              significance = append(significance, "")
            }
          }
          
          plot_df_bar = data.frame(time = unique(annot[[params$comb_group]]), time_median = time_median, time_sd = time_sd, significance = significance)
          significance = data.frame(time = unique(annot[[params$comb_group]]), significance = significance)
            
          # corr line
          mask = (corr[[data_input$rna_class]] == rna)
          corr_value = corr[mask,][[tissue]] 
          p = p_value[mask,][[tissue]]
          # %.3f makes 3 decimal places
          #corr_text = sprintf("corr: %.3f\np: %.3f", corr_value, p)
          if (p < params$sig_niveau) {
            corr_text = sprintf("sig. corr.: %.3f", corr_value) 
          } else {
            corr_text = sprintf("corr.: %.3f", corr_value) 
          }
          
          x_axis_title = xticks_names$categories[[sprintf("%s_capital", params$comb_group)]]
          y_axis_title_big = sprintf("Expr (%s, detection rate = %s%% per %s)", gsub("_norm", "" ,data_input$norm), gsub("p", "", data_input$detection_rate), xticks_names$categories[[params$plot_group]])
          y_axis_title_small = sprintf("Expr (%s,\ndetection rate = %s%%\nper %s)", gsub("_norm", "" ,data_input$norm), gsub("p", "", data_input$detection_rate), xticks_names$categories[[params$plot_group]])
          plot_title = sprintf("%s", unlist(xticks_names[[params$plot_group]][[tissue]]))
          
          font_size_correction = 2.845
          
          x_end = max(unique(annot[[params$comb_group]]))
          x_end = 1.05*x_end
          
          # create bar plot
          p_bar[[i]] = ggplot(plot_df_bar, aes(x = time, y = time_median, fill = as.factor(time))) + 
            geom_bar(stat="identity",  width=1.75) +
            geom_linerange(aes(ymin=time_median, ymax=time_median+time_sd), position=position_dodge(.9)) +
            geom_text(aes(label = significance), vjust=-0.75, size=5) +
            scale_fill_manual(values = unlist(colors[[params$comb_group]])) +
            geom_jitter(data = plot_df, aes_string(y = "value", x = params$comb_group), color = "darkgray", size = 0.25)
          if(poly_deg > 0) {
            p_bar[[i]] = p_bar[[i]] +
              # fitting polynom
              geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, poly_deg), se = FALSE, fullrange=TRUE, color = "black", size = 0.75) 
          } else {
            p_bar[[i]] = p_bar[[i]] +
              geom_smooth(aes(group = 1), method = "loess", se = FALSE, color = "black", size = 0.75) 
            # display correlation value
            #geom_abline(slope = slope, intercept = intercept, color = "black", size = 0.75) 
          }
          p_bar[[i]] = p_bar[[i]] +
            ggtitle(plot_title) +
            xlab("") +
            ylab("") +
            xlim(0, max(unique(annot[[params$comb_group]]))) +
            #ylab(sprintf("expr (%s, detection rate = %s%%\nper %s)", data_input$norm, expr_input$detection_rate, params$plot_group)) +
            scale_x_continuous(limits = c(0, x_end), breaks=c(3, 10, 20, 30),) +
            scale_y_continuous(limits = c(0, 1.3*max(plot_df$value)), expand = c(0,0)) +
            theme_classic() +
            theme(#plot.margin = unit(c(0.1,-1.9,0.25,-1.75), plots_props$image_units), #top, right, bottom, left
                  plot.margin = unit(c(0,0.1,0,0.1), plots_props$image_units), #top, right, bottom, left
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
          
          xpos <- x_end
          ypos <- 1.25*max(plot_df$value) 
  
          if ((rna == manual_plots_parameters$features) & (tissue %in% manual_plots_parameters$tissues)) {
            ypos <- 1.25*maxi_expr_rna
            p_bar_manual[[tissue]] = p_bar[[i]] + ylim(0, ypos)
            p_bar_manual[[tissue]] = p_bar_manual[[tissue]] + geom_text(x = xpos, y = ypos, label = corr_text, hjust = 1, vjust = 1, family = plots_props$font_family, size = 6/font_size_correction)
          }
                  
          p_bar[[i]] = p_bar[[i]] + geom_text(x = xpos, y = ypos, label = corr_text, hjust = 1, vjust = 1, family = plots_props$font_family, size = 6/font_size_correction)
          
          if (rna %in% feature_list) {
            p_bar_rna = p_bar[[i]] + xlab(x_axis_title)
            p_bar_rna = p_bar_rna + ylab(y_axis_title_small)
            #ggsave(sprintf("%s/%s.png", output_folder_bar_single_rna, tissue), p_bar_rna, dpi=plots_props$dpi, width=plots_props$image_width / 2, height=(plots_props$image_height*1.5)/2, units = plots_props$image_units)
            #ggsave(sprintf("%s/%s.svg", output_folder_bar_single_rna, tissue), p_bar_rna, width=plots_props$image_width / 2, height=(plots_props$image_height*1.5)/2, units = plots_props$image_units, limitsize = FALSE, scale = 1)
          }
          
          
          
          # create box plot
          # this line just adds a column with a fixed name "time" to plot_df with the same values as plot_df[[cparams$omb_group]]
          plot_df[["time"]] = plot_df[[params$comb_group]]
          
          # add the significance * to the timepoint only for which the value is the highest
          plot_df = merge(plot_df, significance, by="time")
          for (t in plot_df$time) {
            sig_samples = plot_df[(plot_df$time == t) & (plot_df$significance != ""),]$variable
            max_value = max(plot_df[(plot_df$time == t) & (plot_df$variable %in% sig_samples),]$value)
            plot_df[(plot_df$time == t) & (plot_df$value != max_value),]$significance = ""
          }
          
          # p_box[[i]] = ggplot(plot_df, aes_string(x = params$comb_group, y = "value", group = params$comb_group, fill = as.factor(params$comb_group))) + #plot_df%>%group_by(params$plot_group)%>%summarise(value=sum(value))%>%
          p_box[[i]] = ggplot(plot_df, aes(x = time, y = value, group = time, fill = as.factor(time))) + #plot_df%>%group_by(params$plot_group)%>%summarise(value=sum(value))%>%
            geom_boxplot(aes(fill = as.factor(time)), width = 1.75, color = "black", outlier.colour= NA, outlier.shape=16, outlier.size=1, notch=FALSE, fatten=TRUE) + 
            geom_boxplot(fill=NA, color="black", width = 1.75, fatten=FALSE, coef = 0, outlier.shape = NA, outlier.alpha = 0) +
            #geom_jitter(data = plot_df, aes(y = value, x = time), color = "darkgray", shape = 1, size = 0.75, alpha = 0.75) +
            geom_boxplot(width = 1.75, outlier.colour="grey", outlier.shape=16, outlier.size=1, notch=FALSE, outlier.alpha = 1) +
            scale_fill_manual(values = unlist(colors[[params$comb_group]])) 
          if (poly_deg > 0) {
            p_box[[i]] = p_box[[i]] +
            # fitting polynom
            geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, poly_deg), se = FALSE, fullrange=TRUE, color = "black", size = 0.75) 
          } else {
            p_box[[i]] = p_box[[i]] +
              geom_smooth(aes(group = 1), method = "loess", se = FALSE, color = "black", size = 0.75) 
          }
          p_box[[i]] = p_box[[i]] +
            geom_text(aes(label = significance), color="black", vjust=-0.35, size=5) +
            ggtitle(plot_title) +
            xlab("") +
            ylab("") +
            #ylab(sprintf("expr (%s, detection rate = %s%%\nper %s)", data_input$norm, expr_input$detection_rate, params$plot_group)) +
            scale_x_continuous(limits = c(0, x_end), breaks=c(3, 10, 20, 30)) +
            scale_y_continuous(limits = c(0, 1.3*max(plot_df$value)), expand = c(0,0)) +
            theme_classic() +
            theme(#plot.margin = unit(c(0.1,-1.9,0.25,-1.75), plots_props$image_units), #top, right, bottom, left
                  plot.margin = unit(c(0,0.1,0,0.1), plots_props$image_units), #top, right, bottom, left
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
          #p_box[[i]]
          xpos <- x_end
  
          if ((rna == manual_plots_parameters$features) & (tissue %in% manual_plots_parameters$tissues)) {
            ypos <- 1.25*maxi_expr_rna
            p_box_manual[[tissue]] = p_box[[i]] + ylim(0, ypos)
            p_box_manual[[tissue]] = p_box_manual[[tissue]] + geom_text(x = xpos, y = ypos, label = corr_text, hjust = 1, vjust = 1, family = plots_props$font_family, size = 6/font_size_correction)
          }
          
          ypos <- 1.25*max(plot_df$value) 
          
          p_box[[i]] = p_box[[i]] + ylim(0, ypos) +
                                    geom_text(x = xpos, y = ypos, label = corr_text, hjust = 1, vjust = 1, family = plots_props$font_family, size = 6/font_size_correction)
          
          if (rna %in% feature_list) {
            p_box_rna = p_box[[i]] + xlab(x_axis_title)
            p_box_rna = p_box_rna + ylab(y_axis_title_small)
            #ggsave(sprintf("%s/%s.png", output_folder_box_single_rna, tissue), p_box_rna, dpi=plots_props$dpi, width=plots_props$image_width / 2, height=(plots_props$image_height*1.5)/2, units = plots_props$image_units)
            #ggsave(sprintf("%s/%s.svg", output_folder_box_single_rna, tissue), p_box_rna, width=plots_props$image_width / 2, height=(plots_props$image_height*1.5)/2, units = plots_props$image_units, limitsize = FALSE, scale = 1)
          }
          
          i = i + 1
        }
        # creating concatinated plot
        all_p_bar = arrangeGrob(grobs=p_bar, ncol=params$plot_ncols, nrow=params$plot_nrows, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title_big, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90), right = textGrob(rna, gp=gpar(fontsize=plots_props$font_size_header, fontfamily=plots_props$font_family), rot = 270))
        all_p_box = arrangeGrob(grobs=p_box, ncol=params$plot_ncols, nrow=params$plot_nrows, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title_big, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90), right = textGrob(rna, gp=gpar(fontsize=plots_props$font_size_header, fontfamily=plots_props$font_family), rot = 270))
      
        # save bar plot
        ggsave(sprintf("%s/%s.png", output_folder_bar, rna), all_p_bar, dpi=plots_props$dpi, width=2*plots_props$image_width, height=2*plots_props$image_height, units = plots_props$image_units)
        ggsave(sprintf("%s/%s.svg", output_folder_bar, rna), all_p_bar, width=2*plots_props$image_width, height=2*plots_props$image_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
              
        # save box plot
        ggsave(sprintf("%s/%s.png", output_folder_box, rna), all_p_box, dpi=plots_props$dpi, width=2*plots_props$image_width, height=2*plots_props$image_height, units = plots_props$image_units)
        ggsave(sprintf("%s/%s.svg", output_folder_box, rna), all_p_box, width=2*plots_props$image_width, height=2*plots_props$image_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
      
        if (rna == manual_plots_parameters$features) {
          # sort the figures given the order from the config
          p_bar_manual = p_bar_manual[manual_plots_parameters$tissues]
          p_box_manual = p_box_manual[manual_plots_parameters$tissues]
          
          plot_nrows = manual_plots_parameters$plot_layout[1]
          plot_ncols = manual_plots_parameters$plot_layout[2]
          
          if (manual_plots_parameters$plot_size[2] > 6) {
            y_axis_title = y_axis_title_big
          } else {
            y_axis_title = y_axis_title_small
          }
          
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
          
          p_bar_manual_mod_margins = p_bar_manual
          p_box_manual_mod_margins = p_box_manual
          for (img_i in 1:length(p_box_manual_mod_margins)) {
            if ((!(img_i %in% left_col)) & (!(img_i %in% right_col))) {
              p_box_manual_mod_margins[[img_i]] = p_box_manual_mod_margins[[img_i]] + 
                theme(plot.margin = margin(10, -1, -7, -13, "pt")
                )
              p_bar_manual_mod_margins[[img_i]] = p_bar_manual_mod_margins[[img_i]] + 
                theme(plot.margin = margin(10, -1, -7, -13, "pt")
                )
            }
            else if (img_i %in% left_col) {
              p_box_manual_mod_margins[[img_i]] = p_box_manual_mod_margins[[img_i]] + 
                theme(plot.margin = margin(10, -4, -7, -10, "pt")
                )
              p_bar_manual_mod_margins[[img_i]] = p_bar_manual_mod_margins[[img_i]] + 
                theme(plot.margin = margin(10, -4, -7, -10, "pt")
                )
            } else if (img_i %in% right_col) {
              p_box_manual_mod_margins[[img_i]] = p_box_manual_mod_margins[[img_i]] + 
                theme(plot.margin = margin(10, 4, -7, -18, "pt")
                )
              p_bar_manual_mod_margins[[img_i]] = p_bar_manual_mod_margins[[img_i]] + 
                theme(plot.margin = margin(10, 4, -7, -18, "pt")
                )
            }
            if (!(img_i %in% left_col)) {
              p_box_manual_mod_margins[[img_i]] = p_box_manual_mod_margins[[img_i]] + 
                #scale_y_continuous(labels = function(breaks) {rep_along(breaks, "  ")}, limits = c(0, 1.15*maxi_expr_rna), expand = c(0,0)) +
                theme(axis.ticks.y = element_line(color = "transparent"),
                      axis.text.y = element_text(color = "transparent")
                )
              p_bar_manual_mod_margins[[img_i]] = p_bar_manual_mod_margins[[img_i]] + 
                #scale_y_continuous(labels = function(breaks) {rep_along(breaks, "  ")}, limits = c(0, 1.15*maxi_expr_rna), expand = c(0,0)) +
                theme(axis.ticks.y = element_line(color = "transparent"),
                      axis.text.y = element_text(color = "transparent")
                )
            }
          }
          
          # col_distribution = c()
          # for (col_i in 1:plot_ncols) {
          #   col_tissues = c()
          #   for (row_i in 0:(plot_nrows-1)) {
          #     col_tissues = append(col_tissues, manual_plots_parameters$tissues[col_i + row_i * plot_ncols])
          #   }
          #   col_distribution[[sprintf("col_%s", col_i)]] = col_tissues
          # }
          # 
          # p_bar_manual_adj = list()
          # p_box_manual_adj = list()
          # for (ti in manual_plots_parameters$tissues) {
          #   if (plot_ncols == 3) {
          #     #if (ti != manual_plots_parameters$tissues[1]){
          #     if (ti %in% col_distribution$col_1) {
          #       p_bar_manual_adj[[ti]] = p_bar_manual[[ti]] + theme(plot.margin = unit(c(0,-0.25,0,-0.4), "cm"))
          #       p_box_manual_adj[[ti]] = p_box_manual[[ti]] + theme(plot.margin = unit(c(0,-0.25,0,-0.4), "cm"))
          #       #} else if (ti == manual_plots_parameters$tissues[length(manual_plots_parameters$tissues)]) {
          #     } else if (ti %in% col_distribution$col_2){
          #       p_bar_manual_adj[[ti]] = p_bar_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,0,0,0), "cm"))
          #       p_box_manual_adj[[ti]] = p_box_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,0,0,0), "cm"))
          #     } else if (ti %in% col_distribution$col_3){
          #       p_bar_manual_adj[[ti]] = p_bar_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,0.15,0,-0.15), "cm"))
          #       p_box_manual_adj[[ti]] = p_box_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,0.15,0,-0.15), "cm"))
          #     } 
          #     #p_bar_manual_adj[[ti]] = p_bar_manual_adj[[ti]] + theme(plot.margin = unit(c(0,-0.15,-0.5,-0.15), "cm"))
          #     #p_box_manual_adj[[ti]] = p_box_manual_adj[[ti]] + theme(plot.margin = unit(c(0,-0.15,-0.5,-0.15), "cm"))
          #   } else if (plot_ncols == 4) {
          #     if (ti %in% col_distribution$col_1) {
          #       p_bar_manual_adj[[ti]] = p_bar_manual[[ti]] + theme(plot.margin = unit(c(0,-0.25,0,-0.4), "cm"))
          #       p_box_manual_adj[[ti]] = p_box_manual[[ti]] + theme(plot.margin = unit(c(0,-0.25,0,-0.4), "cm"))
          #       #} else if (ti == manual_plots_parameters$tissues[length(manual_plots_parameters$tissues)]) {
          #     } else if (ti %in% col_distribution$col_2){
          #       p_bar_manual_adj[[ti]] = p_bar_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,0,0,0), "cm"))
          #       p_box_manual_adj[[ti]] = p_box_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,0,0,0), "cm"))
          #     } else if (ti %in% col_distribution$col_3){
          #       p_bar_manual_adj[[ti]] = p_bar_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,0,0,0), "cm"))
          #       p_box_manual_adj[[ti]] = p_box_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,0,0,0), "cm"))
          #     }else if (ti %in% col_distribution$col_4){
          #       p_bar_manual_adj[[ti]] = p_bar_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,0.15,0,-0.15), "cm"))
          #       p_box_manual_adj[[ti]] = p_box_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,0.15,0,-0.15), "cm"))
          #     } 
          #   } else if (plot_ncols == 2) {
          #     if (ti %in% col_distribution$col_1) {
          #       p_bar_manual_adj[[ti]] = p_bar_manual[[ti]] + theme(plot.margin = unit(c(0,-0.25,0,-0.4), "cm"))
          #       p_box_manual_adj[[ti]] = p_box_manual[[ti]] + theme(plot.margin = unit(c(0,-0.25,0,-0.4), "cm"))
          #       #} else if (ti == manual_plots_parameters$tissues[length(manual_plots_parameters$tissues)]) {
          #     } else if (ti %in% col_distribution$col_2){
          #       p_bar_manual_adj[[ti]] = p_bar_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,0.5,0,0), "cm"))
          #       p_box_manual_adj[[ti]] = p_box_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,0.5,0,0), "cm"))
          #     } 
          #   }    
          # }
          
          all_p_bar_manual = arrangeGrob(grobs=p_bar_manual_mod_margins, ncol=plot_ncols, nrow=plot_nrows, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90), right = textGrob(rna, gp=gpar(fontsize=plots_props$font_size_header, fontfamily=plots_props$font_family), rot = 270))
          all_p_box_manual = arrangeGrob(grobs=p_box_manual_mod_margins, ncol=plot_ncols, nrow=plot_nrows, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90), right = textGrob(rna, gp=gpar(fontsize=plots_props$font_size_header, fontfamily=plots_props$font_family), rot = 270))
          
          plot_width = manual_plots_parameters$plot_size[1]
          plot_height = manual_plots_parameters$plot_size[2]
          
          ggsave(sprintf("%s/%s_%s_%sx%s.png", output_folder_bar_single_rna, rna, paste(manual_plots_parameters$tissues, collapse = "_"), plot_width, plot_height), all_p_bar_manual, dpi=plots_props$dpi, width=plot_width, height=plot_height, units = plots_props$image_units)
          ggsave(sprintf("%s/%s_%s_%sx%s.svg", output_folder_bar_single_rna, rna, paste(manual_plots_parameters$tissues, collapse = "_"), plot_width, plot_height), all_p_bar_manual, width=plot_width, height=plot_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
          
          ggsave(sprintf("%s/%s_%s_%sx%s.png", output_folder_box_single_rna, rna, paste(manual_plots_parameters$tissues, collapse = "_"), plot_width, plot_height), all_p_box_manual, dpi=plots_props$dpi, width=plot_width, height=plot_height, units = plots_props$image_units)
          ggsave(sprintf("%s/%s_%s_%sx%s.svg", output_folder_box_single_rna, rna, paste(manual_plots_parameters$tissues, collapse = "_"), plot_width, plot_height), all_p_box_manual, width=plot_width, height=plot_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
          
          rm(all_p_bar_manual)
          rm(p_bar_manual_mod_margins)
          rm(all_p_box_manual)
          rm(p_box_manual_mod_margins)
          
          if (plot_nrows == 1){
            p_bar_manual_adj = list()
            p_box_manual_adj = list()
            for (ti in manual_plots_parameters$tissues) {
              #if (ti != manual_plots_parameters$tissues[1]){
              if (ti == manual_plots_parameters$tissues[length(manual_plots_parameters$tissues)]) {
                print(ti)
                p_bar_manual_adj[[ti]] = p_bar_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,0.25,-0.5,-0.1), "cm"))
                p_box_manual_adj[[ti]] = p_box_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,0.25,-0.5,-0.1), "cm"))
              #} else if (ti == manual_plots_parameters$tissues[length(manual_plots_parameters$tissues)]) {
              } else if (ti != manual_plots_parameters$tissues[1]){
                p_bar_manual_adj[[ti]] = p_bar_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,-0.1,-0.5,-0.1), "cm"))
                p_box_manual_adj[[ti]] = p_box_manual[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,-0.1,-0.5,-0.1), "cm"))
              } else {
                p_bar_manual_adj[[ti]] = p_bar_manual[[ti]] + theme(plot.margin = unit(c(0,-0.1,-0.5,-1), "cm"))
                p_box_manual_adj[[ti]] = p_box_manual[[ti]] + theme(plot.margin = unit(c(0,-0.1,-0.5,-1), "cm"))
              } 
              #p_bar_manual_adj[[ti]] = p_bar_manual_adj[[ti]] + theme(plot.margin = unit(c(0,-0.15,-0.5,-0.15), "cm"))
              #p_box_manual_adj[[ti]] = p_box_manual_adj[[ti]] + theme(plot.margin = unit(c(0,-0.15,-0.5,-0.15), "cm"))
            }
            
            all_p_bar_manual = arrangeGrob(grobs=p_bar_manual_adj, ncol=plot_ncols, nrow=plot_nrows, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
            all_p_box_manual = arrangeGrob(grobs=p_box_manual_adj, ncol=plot_ncols, nrow=plot_nrows, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
            
            # plot_width = 2 * plots_props$image_width
            # plot_height = (plots_props$image_height*1.5)/2
            # ggsave(sprintf("%s/expr_%s_%s_%sx%s_jitter_adj.png", output_folder_bar_single_rna, rna, paste(manual_plots_parameters$tissues, collapse = "_"), plot_width, plot_height), all_p_bar_manual, dpi=plots_props$dpi, width=plot_width, height=plot_height, units = plots_props$image_units)
            # ggsave(sprintf("%s/expr_%s_%s_%sx%s_jitter_adj.svg", output_folder_bar_single_rna, rna, paste(manual_plots_parameters$tissues, collapse = "_"), plot_width, plot_height), all_p_bar_manual, width=plot_width, height=plot_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
            # 
            # ggsave(sprintf("%s/expr_%s_%s_%sx%s_jitter_adj.png", output_folder_box_single_rna, rna, paste(manual_plots_parameters$tissues, collapse = "_"), plot_width, plot_height), all_p_box_manual, dpi=plots_props$dpi, width=plot_width, height=plot_height, units = plots_props$image_units)
            # ggsave(sprintf("%s/expr_%s_%s_%sx%s_jitter_adj.svg", output_folder_box_single_rna, rna, paste(manual_plots_parameters$tissues, collapse = "_"), plot_width, plot_height), all_p_box_manual, width=plot_width, height=plot_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
            # 
            # plot_width = manual_plots_parameters$plot_size[1]
            # plot_height = manual_plots_parameters$plot_size[2]
            # ggsave(sprintf("%s/expr_%s_%s_%sx%s_jitter_adj.png", output_folder_bar_single_rna, rna, paste(manual_plots_parameters$tissues, collapse = "_"), plot_width, plot_height), all_p_bar_manual, dpi=plots_props$dpi, width=plot_width, height=plot_height, units = plots_props$image_units)
            # ggsave(sprintf("%s/expr_%s_%s_%sx%s_jitter_adj.svg", output_folder_bar_single_rna, rna, paste(manual_plots_parameters$tissues, collapse = "_"), plot_width, plot_height), all_p_bar_manual, width=plot_width, height=plot_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
            # 
            # ggsave(sprintf("%s/expr_%s_%s_%sx%s_jitter_adj.png", output_folder_box_single_rna, rna, paste(manual_plots_parameters$tissues, collapse = "_"), plot_width, plot_height), all_p_box_manual, dpi=plots_props$dpi, width=plot_width, height=plot_height, units = plots_props$image_units)
            # ggsave(sprintf("%s/expr_%s_%s_%sx%s_jitter_adj.svg", output_folder_box_single_rna, rna, paste(manual_plots_parameters$tissues, collapse = "_"), plot_width, plot_height), all_p_box_manual, width=plot_width, height=plot_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
           }
          rm(all_p_bar_manual)
          rm(all_p_box_manual)
        }
        
        # cleanup an call garbage cleaner
        rm(p_bar)
        rm(p_box)
        
        rm(p_bar_manual)
        rm(p_box_manual)
        
        gc()
      }
      
      #stopCluster(cl)
    }
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
