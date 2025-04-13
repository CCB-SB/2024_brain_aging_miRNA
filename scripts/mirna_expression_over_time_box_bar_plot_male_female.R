suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))

#snakemake = readRDS("snakemake_mirna_expression_over_time_box_bar_plot_male_female.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("mirna_expression_over_time_box_bar_plot_male_female")


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
input_files = snakemake@params$files
manual_plots_parameter_list = snakemake@params$manual_plots_parameters
params = snakemake@params$parameters 
plots_props = snakemake@params$plots_props
xticks_names = snakemake@params$xticks_names
colors = snakemake@params$colors
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_mirna_expression_over_time_box_bar_plot_male_female.rds")
#stop()


#------------------------------------ Script ----------------------------------- 
for (i in 1:length(manual_plots_parameter_list)) {
  manual_plots_parameters = manual_plots_parameter_list[[i]]
  plot_nrows_male = manual_plots_parameters$plot_layout_male[1]
  plot_ncols_male = manual_plots_parameters$plot_layout_male[2]
  plot_nrows_female = manual_plots_parameters$plot_layout_female[1]
  plot_ncols_female = manual_plots_parameters$plot_layout_female[2]
  
  plot_height = manual_plots_parameters$plot_size[1]
  plot_width = manual_plots_parameters$plot_size[2]
  
  for (m in params$corr_methods) {
    for (poly_deg in params$polynom_degree) {
  
      p_box = list()
      for (sex_key in c("male", "female")) {
        
        if (sex_key == "male") {
          data_sub_set = data_input$data_sub_set_male
          rna = manual_plots_parameters$feature_male
          tissues_sex = manual_plots_parameters$tissues_male
        } else {
          data_sub_set = data_input$data_sub_set_female
          rna = manual_plots_parameters$feature_female
          tissues_sex = manual_plots_parameters$tissues_female
        }
        
        # load data files
        expr = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
        annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
        # check if diff_exp file exists
        diff_exp_filename = sprintf("%s/%s_%s/results_%s/matrices/diff_exp/diff_exp%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_sub_set, params$diff_log)
        diff_exp = fread(diff_exp_filename, sep='\t') 
        
        
        # load correlation values
        corr = fread(sprintf("%s/%s_%s/results_%s/matrices/aggregation_corr/corr_method=%s_%ss_with_%s_per_%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_sub_set, m, data_input$rna_class, params$comp_group, params$plot_group), sep='\t') 
        p_value = fread(sprintf("%s/%s_%s/results_%s/matrices/aggregation_corr/corr_method=%s_adj_pvalues_%ss_with_%s_per_%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_sub_set, m, data_input$rna_class, params$comp_group, params$plot_group), sep='\t') 
        
        # create output folders
        output_folder_box_single = sprintf("%s/%s_%s/results_%s_%s/figures/expr/box_plots/time_single_corr_method=%s_poly_deg=%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set_male, data_input$data_sub_set_female, m, poly_deg)
        #print(output_folder_box_single)
        dir.create(output_folder_box_single, recursive=TRUE)
  
        expr_feature = expr[expr[[data_input$feature_column]] == rna,]
        diff_exp_feature = diff_exp[diff_exp[[data_input$rna_class]] == rna,]
          
        tmp = (expr_feature != rna)
        expr_feature_tmp = expr_feature[,..tmp]
        samples_for_plot = annot[annot[[params$plot_group]] %in% tissues_sex][[data_input$identifier_column]]
        maxi_expr_rna = max(expr_feature_tmp[,..samples_for_plot])
        
        for (tissue in tissues_sex) {
          id_group = annot[annot[[params$plot_group]] == tissue,][[data_input$identifier_column]] 
          tmp = colnames(expr_feature) %in% id_group
          expr_feature_tissue = expr_feature[,..tmp]
          
          plot_df = melt(expr_feature_tissue)
          plot_df[[params$plot_group]] = annot[annot[[data_input$identifier_column]] %in% plot_df$variable][[params$plot_group]]
          
          time_list = c()
          for (id in plot_df$variable) {
            time_list = append(time_list, annot[annot[[data_input$identifier_column]] == id,][[params$comp_group]])
          }
          plot_df[[params$comp_group]] = as.numeric(time_list)
      
          time_median = c()
          time_sd = c()
          significance = c()
          for (time in unique(annot[[params$comp_group]])) {
            time_median = append(time_median, median(plot_df[plot_df[[params$comp_group]] == time,]$value))
            time_sd = append(time_sd, sd(plot_df[plot_df[[params$comp_group]] == time,]$value))
            
            if (time != min(as.numeric(unique(annot[[params$comp_group]])))) {
              ttest_adj_pval = diff_exp_feature[[sprintf("ttest_adjp_%s__%s_vs_%s_%s=%s", params$comp_group, time, sort(as.numeric(unique(annot[[params$comp_group]])))[1], params$plot_group, tissue)]]  #ttest_adjp_age__12_vs_3_brain_region_fixed=cc
              fc = log2(diff_exp_feature[[sprintf("fc_%s__%s_vs_%s_%s=%s", params$comp_group, time, sort(as.numeric(unique(annot[[params$comp_group]])))[1], params$plot_group, tissue)]])  #ttest_adjp_yy__12_vs_3_brain_region_fixed=cc
              if (is.na(fc)) {
                significance = append(significance, "")
              } else if (abs(fc) >= log2(params$fc_th)) {
                if (is.na(ttest_adj_pval)) { 
                  #significance = append(significance, "+/-")
                  significance = append(significance, "")
                } else if (ttest_adj_pval < 0.001) { 
                  significance = append(significance, "***")
                } else if (ttest_adj_pval < 0.01) {
                  significance = append(significance, "**")
                } else if (ttest_adj_pval < 0.05) {
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
          
          significance = data.frame(time = unique(annot[[params$comp_group]]), significance = significance)
            
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
          
          x_axis_title = xticks_names$categories[[sprintf("%s_capital", params$comp_group)]]
          y_axis_title_big = sprintf("Expr (%s, detection rate = %s%% per %s)", gsub("_norm", "" ,data_input$norm), gsub("p", "", data_input$detection_rate), xticks_names$categories[[params$plot_group]])
          y_axis_title_small = sprintf("Expr (%s,\ndetection rate = %s%%\nper %s)", gsub("_norm", "" ,data_input$norm), gsub("p", "", data_input$detection_rate), xticks_names$categories[[params$plot_group]])
          plot_title = sprintf("%s", unlist(xticks_names[[params$plot_group]][[tissue]]))
          #plot_title = sprintf("%s: %s", sex_key, unlist(xticks_names[[params$plot_group]][[tissue]]))
          
          font_size_correction = 2.845
          
          x_end = max(unique(annot[[params$comp_group]]))
          x_end = 1.05*x_end
          
          # plot box plot
          # this line just adds a column with a fixed name "time" to plot_df with the same values as plot_df[[cparams$omb_group]]
          plot_df[["time"]] = plot_df[[params$comp_group]]
          
          # add the significance * to the timepoint only for which the value is the highest
          plot_df = merge(plot_df, significance, by="time")
          sig_samples = plot_df[plot_df$significance != "",]$variable
          
          #max_value = max(plot_df[plot_df$variable %in% sig_samples]$value)
          #plot_df[plot_df$value != max_value,]$significance = ""
          
          # replace all not max with "" for each params$comp_group
          for (t in plot_df[[params$comp_group]]) {
            max_value = max(plot_df[(plot_df$variable %in% sig_samples) & (plot_df[[params$comp_group]] == t)]$value)
            plot_df[(plot_df$value != max_value) & (plot_df[[params$comp_group]] == t),]$significance = ""
          }
          
          xpos = x_end+0.1
          ypos = 1.25*maxi_expr_rna
          
          # p_box[[sprintf("%s_%s", sex_key, tissue)]] = ggplot(plot_df, aes_string(x = params$comp_group, y = "value", group = params$comp_group, fill = as.factor(params$comp_group))) + #plot_df%>%group_by(params$plot_group)%>%summarise(value=sum(value))%>%
          p_box[[sprintf("%s_%s", sex_key, tissue)]] = ggplot(plot_df, aes(x = time, y = value, group = time, fill = as.factor(time))) + #plot_df%>%group_by(params$plot_group)%>%summarise(value=sum(value))%>%
            geom_boxplot(aes(fill = as.factor(time)), width = 1.75, color = "black", outlier.colour= NA, outlier.shape=16, outlier.size=1, notch=FALSE, fatten=TRUE) + 
            geom_boxplot(fill=NA, color="black", width = 1.75, fatten=FALSE, coef = 0, outlier.shape = NA, outlier.alpha = 0) +
            #geom_jitter(data = plot_df, aes(y = value, x = time), color = "darkgray", size = 0.75, shape = 1, alpha = 0.75) + 
            geom_boxplot(width = 1.75, outlier.colour="grey", outlier.shape=16, outlier.size=1, notch=FALSE, outlier.alpha = 1) +
            scale_fill_manual(values = unlist(colors[[params$comp_group]])) 
          if (poly_deg > 0) {
            p_box[[sprintf("%s_%s", sex_key, tissue)]] = p_box[[sprintf("%s_%s", sex_key, tissue)]] +
            # fitting polynom
            geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, poly_deg), se = FALSE, fullrange=TRUE, color = "black", size = 0.75) 
          } else {
            p_box[[sprintf("%s_%s", sex_key, tissue)]] = p_box[[sprintf("%s_%s", sex_key, tissue)]] +
              geom_smooth(aes(group = 1), method = "loess", se = FALSE, color = "black", size = 0.75) 
          }
          p_box[[sprintf("%s_%s", sex_key, tissue)]] = p_box[[sprintf("%s_%s", sex_key, tissue)]] +
            geom_text(aes(label = significance), color="black", vjust=-0.35, size=3) +
            ggtitle(plot_title) +
            xlab("") +
            ylab("") +
            #ylab(sprintf("expr (%s, detection rate = %s%%\nper %s)", data_input$norm, expr_input$detection_rate, params$plot_group)) +
            scale_x_continuous(limits = c(0, x_end), breaks=c(3, 10, 20, 30)) +
            scale_y_continuous(limits = c(0, 1.3*max(plot_df$value)), expand = c(0,0)) +
            ylim(0, ypos) +
            geom_text(x = xpos, y = ypos, label = corr_text, hjust = 1, vjust = 1, family = plots_props$font_family, size = 6/font_size_correction) +
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
        }
      }
  
      # sort the figures given the order from the config
      p_box_male = p_box[c(sprintf("male_%s", manual_plots_parameters$tissues_male))]
      p_box_female = p_box[c(sprintf("female_%s", manual_plots_parameters$tissues_female))]
      if(data_input$rna_class == "miRNA") {
        if (plot_nrows_male == 1) {
          left_text = y_axis_title_small
          # -18
          left_col_margin_left = -2
          left_col_margin_right = -16
          right_col_margin_left = -28
          right_col_margin_right = 10
          middle_col_margin_left = -15
          middle_col_margin_right = -3
        } else {
          left_text = y_axis_title_big
          # -18
          left_col_margin_left = -4
          left_col_margin_right = -14
          right_col_margin_left = -28
          right_col_margin_right = 10
          middle_col_margin_left = -16
          middle_col_margin_right = -2
        }
      } else {
        if (plot_nrows_male == 1) {
          left_text = y_axis_title_small
          # -16
          left_col_margin_left = -6
          left_col_margin_right = -10
          right_col_margin_left = -24
          right_col_margin_right = 8
          middle_col_margin_left = -15
          middle_col_margin_right = -1
        } else {
          left_text = y_axis_title_big
          # -16
          left_col_margin_left = -8
          left_col_margin_right = -8
          right_col_margin_left = -26
          right_col_margin_right = 10
          middle_col_margin_left = -16
          middle_col_margin_right = 0
        }
      }
      
      # determine which figures are not on the left side
      left_col = c()
      for (row in 1:plot_nrows_male - 1) {
        left_col = append(left_col, plot_ncols_male*row + 1)
      }
      # determine which figures are not on the left side
      right_col = c()
      for (row in 1:plot_nrows_male) {
        right_col = append(right_col, plot_ncols_male*row)
      }
      
      p_box_male_mod_margins = p_box_male
      for (img_i in 1:length(p_box_male_mod_margins)) {
        if ((!(img_i %in% left_col)) & (!(img_i %in% right_col))) {
          p_box_male_mod_margins[[img_i]] = p_box_male_mod_margins[[img_i]] + 
            theme(plot.margin = margin(10, middle_col_margin_right, -7, middle_col_margin_left, "pt")
            )
        }
        else if (img_i %in% left_col) {
          p_box_male_mod_margins[[img_i]] = p_box_male_mod_margins[[img_i]] + 
            theme(plot.margin = margin(10, left_col_margin_right, -7, left_col_margin_left, "pt")
            )
        } else if (img_i %in% right_col) {
          p_box_male_mod_margins[[img_i]] = p_box_male_mod_margins[[img_i]] + 
            theme(plot.margin = margin(10, right_col_margin_right, -7, right_col_margin_left, "pt")
            )
        }
        if (!(img_i %in% left_col)) {
          p_box_male_mod_margins[[img_i]] = p_box_male_mod_margins[[img_i]] + 
            #scale_y_continuous(labels = function(breaks) {rep_along(breaks, "  ")}, limits = c(0, 1.15*maxi_expr_rna), expand = c(0,0)) +
            theme(axis.ticks.y = element_line(color = "transparent"),
                  axis.text.y = element_text(color = "transparent")
            )
        }
      }
      
      # determine which figures are not on the left side
      left_col = c()
      for (row in 1:plot_nrows_female - 1) {
        left_col = append(left_col, plot_ncols_female*row + 1)
      }
      # determine which figures are not on the left side
      right_col = c()
      for (row in 1:plot_nrows_female) {
        right_col = append(right_col, plot_ncols_female*row)
      }
      
      p_box_female_mod_margins = p_box_female
      for (img_i in 1:length(p_box_female_mod_margins)) {
        if ((!(img_i %in% left_col)) & (!(img_i %in% right_col))) {
          p_box_female_mod_margins[[img_i]] = p_box_female_mod_margins[[img_i]] + 
            theme(plot.margin = margin(10, middle_col_margin_right, -7, middle_col_margin_left, "pt")
            )
        }
        else if (img_i %in% left_col) {
          p_box_female_mod_margins[[img_i]] = p_box_female_mod_margins[[img_i]] + 
            theme(plot.margin = margin(10, left_col_margin_right, -7, left_col_margin_left, "pt")
            )
        } else if (img_i %in% right_col) {
          p_box_female_mod_margins[[img_i]] = p_box_female_mod_margins[[img_i]] + 
            theme(plot.margin = margin(10, right_col_margin_right, -7, right_col_margin_left, "pt")
            )
        }
        if (!(img_i %in% left_col)) {
          p_box_female_mod_margins[[img_i]] = p_box_female_mod_margins[[img_i]] + 
            #scale_y_continuous(labels = function(breaks) {rep_along(breaks, "  ")}, limits = c(0, 1.15*maxi_expr_rna), expand = c(0,0)) +
            theme(axis.ticks.y = element_line(color = "transparent"),
                  axis.text.y = element_text(color = "transparent")
            )
        }
      }
        
      # creating concatinated plot
      all_p_box_sex = list()
      all_p_box_sex[["male"]] = arrangeGrob(grobs=p_box_male_mod_margins, ncol=plot_ncols_male, nrow=plot_nrows_male, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)))
      all_p_box_sex[["female"]] = arrangeGrob(grobs=p_box_female_mod_margins, ncol=plot_ncols_female, nrow=plot_nrows_female, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)))
      
      if (manual_plots_parameters$feature_male == manual_plots_parameters$feature_female) {
        feature_name = manual_plots_parameters$feature_male
      } else {
        feature_name = ""
      }
  
      all_p_box = arrangeGrob(grobs=all_p_box_sex, ncol=2, nrow=1, left = textGrob(left_text, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90), right = textGrob(feature_name, gp=gpar(fontsize=plots_props$font_size_header, fontfamily=plots_props$font_family), rot = 270))
      
      ggsave(sprintf("%s/%s_%s_%s_%s.png", output_folder_box_single, manual_plots_parameters$feature_male, paste(manual_plots_parameters$tissues_male, collapse = "_"), manual_plots_parameters$feature_female, paste(manual_plots_parameters$tissues_female, collapse = "_")), all_p_box, dpi=plots_props$dpi, width=plot_width, height=plot_height, units = plots_props$image_units)
      ggsave(sprintf("%s/%s_%s_%s_%s.svg", output_folder_box_single, manual_plots_parameters$feature_male, paste(manual_plots_parameters$tissues_male, collapse = "_"), manual_plots_parameters$feature_female, paste(manual_plots_parameters$tissues_female, collapse = "_")), all_p_box, width=plot_width, height=plot_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
  
      
      # if (plot_nrows_male == 1){
      #   p_box_adj = list()
      #   for (ti in names(p_box)) {
      #     #if (ti != manual_plots_parameters$tissues[1]){
      #     if (ti == sprintf("female_%s", manual_plots_parameters$tissues_female[length(manual_plots_parameters$tissues_female)])) {
      #       print(ti)
      #       p_box_adj[[ti]] = p_box[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,0.25,-0.3,-0.1), "cm"))
      #     #} else if (ti == manual_plots_parameters$tissues[length(manual_plots_parameters$tissues)]) {
      #     } else if (ti != manual_plots_parameters$tissues_female[1]){
      #       p_box_adj[[ti]] = p_box[[ti]] + theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,-0.1,-0.3,-0.1), "cm"))
      #     } else {
      #       p_box_adj[[ti]] = p_box[[ti]] + theme(plot.margin = unit(c(0,-0.1,-0.3,-0.25), "cm"))
      #     }
      #     #p_bar_manual_adj[[ti]] = p_bar_manual_adj[[ti]] + theme(plot.margin = unit(c(0,-0.15,-0.5,-0.15), "cm"))
      #     #p_box_manual_adj[[ti]] = p_box_manual_adj[[ti]] + theme(plot.margin = unit(c(0,-0.15,-0.5,-0.15), "cm"))
      #   }
      # 
      #   all_p_box_adj = arrangeGrob(grobs=p_box_adj, ncol=(plot_ncols_male+plot_ncols_female), nrow=plot_nrows_male, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title_small, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
      # 
      #   ggsave(sprintf("%s/expr_%s_%s_%s_%s_%s_jitter_adj.png", output_folder_box_single, rna, manual_plots_parameters$feature_male, paste(manual_plots_parameters$tissues_male, collapse = "_"), manual_plots_parameters$feature_female, paste(manual_plots_parameters$tissues_female, collapse = "_")), all_p_box_adj, dpi=plots_props$dpi, width=plot_width, height=plot_height, units = plots_props$image_units)
      #   ggsave(sprintf("%s/expr_%s_%s_%s_%s_%s_jitter_adj.svg", output_folder_box_single, rna, manual_plots_parameters$feature_male, paste(manual_plots_parameters$tissues_male, collapse = "_"), manual_plots_parameters$feature_female, paste(manual_plots_parameters$tissues_female, collapse = "_")), all_p_box_adj, width=plot_width, height=plot_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
      #  }
      rm(p_box)
      rm(p_box_male)
      rm(p_box_female)
      rm(p_box_male_mod_margins)
      rm(p_box_female_mod_margins)
      rm(p_box_adj)
    }
    
    # cleanup an call garbage cleaner
    rm(all_p_box)
    rm(all_p_box_adj)
    rm(all_p_box_sex)
    
    #gc()
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
