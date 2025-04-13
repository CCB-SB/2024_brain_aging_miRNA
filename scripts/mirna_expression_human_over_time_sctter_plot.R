suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(foreach))  # for parallelization
suppressPackageStartupMessages(library(doParallel))  # for parallelization

#snakemake = readRDS("snakemake_mirna_expression_human_over_time_sctter_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("mirna_expression_human_over_time_sctter_plot")


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
miRNA_human_path = snakemake@params$miRNA_human_path
feature_config = snakemake@params$features
manual_plots_parameters = snakemake@params$manual_plots_parameters
params = snakemake@params$parameters 
corr_methods = snakemake@params$corr_methods
plots_props = snakemake@params$plots_props
xticks_names = snakemake@params$xticks_names
colors = snakemake@params$colors
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_mirna_expression_human_over_time_sctter_plot.rds")
#stop()

# load annot
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
expr = fread(sprintf("%s/data_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)

poly_degs = c(2,3,0)


#------------------------------------ Script ----------------------------------- 
feature_list = unique(append(feature_config, manual_plots_parameters$features))

#feature_list = unique(c(manual_plots_parameters$features, feature_config))
#print(feature_list)
print(sprintf("features to be processed: %s", length(feature_list)))

for (m in corr_methods) {
  for (poly_deg in poly_degs) {
    # load correlation values
    corr = fread(sprintf("%s/%s_%s/results_%s/matrices/corr_%s_%s_with_%s/%s/correlation_with_%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, data_input$rna_class, data_input$data_sub_set, params$comb_group, m, params$comb_group), sep='\t') 
    p_value = fread(sprintf("%s/%s_%s/results_%s/matrices/corr_%s_%s_with_%s/%s/padj_with_%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, data_input$rna_class, data_input$data_sub_set, params$comb_group, m, params$comb_group), sep='\t') 
    
    # create output folders
    #output_folder_scatter = sprintf("%s/%s_%s/results_%s/figures/expr/scatter_plots/corr_method=%s_poly_deg=%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m, poly_deg)
    output_folder_scatter = sprintf("%s/%s_%s/results_%s/figures/expr/scatter_plots/corr_method=%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m)
    dir.create(output_folder_scatter, recursive=TRUE)
     
    # Set up parallel computing
    #n_cores = 64
    # print(sprintf("running in %s threads", n_cores))
    #cl <- makeCluster(n_cores, type="FORK")
    #registerDoParallel(cl)
    
    p_scatter = list()
    p_scatter_manual = list()
    feature_list_new = c()
    #print(length(feature_list))
    for (rna in feature_list) {
    #foreach (rna=feature_list) %dopar% {
      print(rna)
      expr_feature = expr[expr[[data_input$feature_column]] == rna,]
      
      if (nrow(expr_feature) == 0) {
        print(sprintf("%s is not included in the filtered dataset. Skip this feature", rna))
        next
      } else {
        feature_list_new = append(feature_list_new, rna)
      }
      
      plot_df = melt(expr_feature)

      time_list = c()
      sex_status_list = c()
      for (id in plot_df$variable) {
        time_list = append(time_list, annot[annot[[data_input$identifier_column]] == id,][[params$comb_group]])
        sex_status_list = append(sex_status_list, ifelse(annot[annot[[data_input$identifier_column]] == id,][[params$colour_group]] == 1, "Male", "Female"))
      }
      plot_df$plot_group = as.numeric(time_list)
      plot_df$colour_group = sex_status_list
    
      # corr line
      mask = (corr$RNA == rna)
      corr_value_male = corr[mask,]$`disease_status=0_samples=male`
      corr_value_female = corr[mask,]$`disease_status=0_samples=female`
      p_male = p_value[mask,]$`disease_status=0_samples=male`
      p_female = p_value[mask,]$`disease_status=0_samples=female`
      # %.3f makes 3 decimal places
      #corr_text = sprintf("corr: %.3f\np: %.3f", corr_value, p)
      if (p_male < params$sig_niveau) {
        corr_text_male = sprintf("sig. corr.: %.3f", corr_value_male)
      } else {
        corr_text_male = sprintf("corr.: %.3f", corr_value_male)
      }
      if (p_female < params$sig_niveau) {
        corr_text_female = sprintf("sig. corr.: %.3f", corr_value_female)
      } else {
        corr_text_female = sprintf("corr.: %.3f", corr_value_female)
      }
      corr_text = sprintf("Male: %s\nFemale: %s", corr_text_male, corr_text_female)
        
      x_axis_title = xticks_names$categories[[sprintf("%s_capital", params$comb_group)]]
      y_axis_title_big = sprintf("Expr (%s, detection rate = %s%% per age)", gsub("_norm", "" ,data_input$norm), gsub("p", "", data_input$detection_rate))
      y_axis_title_small = sprintf("Expr (%s,\ndetection rate = %s%%\nper age)", gsub("_norm", "" ,data_input$norm), gsub("p", "", data_input$detection_rate))
      plot_title = rna
      font_size_correction = 2.845
        
      x_end = max(unique(annot[[params$comb_group]]))
      x_end = 1.05*x_end

      # create box plot
      # this line just adds a column with a fixed name "time" to plot_df with the same values as plot_df[[cparams$omb_group]]
      #plot_df[["time"]] = plot_df[[params$comb_group]]
      
      #colour_disease_status = c("#CC6928", "#459BAB") 
      #names(colour_disease_status) = c("1", "0")
      plot_df$colour_group_fac = as.factor(plot_df$colour_group)
      
      # p_box[[i]] = ggplot(plot_df, aes_string(x = params$comb_group, y = "value", group = params$comb_group, fill = as.factor(params$comb_group))) + #plot_df%>%group_by(params$plot_group)%>%summarise(value=sum(value))%>%
      p_scatter[[rna]] = ggplot(plot_df, aes(x = plot_group, y = value, color = colour_group_fac )) + #plot_df%>%group_by(params$plot_group)%>%summarise(value=sum(value))%>%
        geom_point(size = 0.5) +
        geom_smooth(method = "loess", se = FALSE, size = 0.8, alpha = 1) + # Add a smoothed trend line
        scale_color_manual(values = unlist(colors$sex))
      # if (poly_deg > 0) {
      #   p_scatter[[rna]] = p_scatter[[rna]] +
      #   # fitting polynom
      #   geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, poly_deg), se = FALSE, fullrange=TRUE, color = "black", size = 0.75) 
      # } else {
      #   p_scatter[[rna]] = p_scatter[[rna]] +
      #     geom_smooth(aes(group = 1), method = "loess", se = FALSE, color = "black", size = 0.75) 
      # }
      p_scatter[[rna]] = p_scatter[[rna]] +
        ggtitle(plot_title) +
        xlab("") +
        ylab("") +
        #ylab(sprintf("expr (%s, detection rate = %s%%\nper %s)", data_input$norm, expr_input$detection_rate, params$plot_group)) +
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
              #plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header),
              plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size, hjust = 2), #  margin = margin(t = unit(10, "pt"), r = unit(5, "pt"), b = unit(5, "pt"), l = unit(-40, "pt"))   
              legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
              legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size)
        )
      #p_scatter[[rna]]
      maxi_expr_rna = max(plot_df$value)
      xpos = max(plot_df$plot_group)
      #if (rna %in% manual_plots_parameters$features) {
        #ypos = 1.25*maxi_expr_rna
        #p_scatter_manual[[rna]] = p_scatter[[rna]] + ylim(0, ypos)
        #p_scatter_manual[[rna]] = p_scatter_manual[[rna]] + geom_text(x = xpos, y = ypos, label = corr_text, hjust = 1, vjust = 1, family = plots_props$font_family, size = 6/font_size_correction)
      #}
      
      ypos = 1.5*max(plot_df$value) 
      
      p_scatter[[rna]] = p_scatter[[rna]] + ylim(0, ypos) +
                                geom_text(x = xpos, y = ypos, label = corr_text, hjust = 1, vjust = 1, family = plots_props$font_family, size = 4/font_size_correction, color = "black") +
                                coord_cartesian(clip = "off") 
      
      if (rna %in% feature_list) {
        p_scatter_rna = p_scatter[[rna]] + xlab(x_axis_title)
        p_scatter_rna = p_scatter_rna + ylab(y_axis_title_small)
        width = 3.5 #plots_props$image_width / 2
        height = 3.5 #(plots_props$image_height*1.5)/2
        ggsave(sprintf("%s/%s.png", output_folder_scatter, rna), p_scatter_rna, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
        ggsave(sprintf("%s/%s.svg", output_folder_scatter, rna), p_scatter_rna, width=width, height=height, units = plots_props$image_units, limitsize = FALSE, scale = 1)
      }
    }
    
    # plot_nrows = 1
    # plot_ncols = length(p_scatter)
    #   
    # # determine which figures are not on the left side
    # left_col = c()
    # for (row in 1:plot_nrows - 1) {
    #   left_col = append(left_col, plot_ncols*row + 1)
    # }
    # # determine which figures are not on the left side
    # right_col = c()
    # for (row in 1:plot_nrows) {
    #   right_col = append(right_col, plot_ncols*row)
    # }
    # 
    # p_scatter_mod_margins = p_scatter
    # for (img_i in 1:length(p_scatter_mod_margins)) {
    #   if ((!(img_i %in% left_col)) & (!(img_i %in% right_col))) {
    #     p_scatter_mod_margins[[img_i]] = p_scatter_mod_margins[[img_i]] +
    #       theme(plot.margin = margin(10, -1, -7, -13, "pt")
    #       )
    #   }
    #   else if (img_i %in% left_col) {
    #     p_scatter_mod_margins[[img_i]] = p_scatter_mod_margins[[img_i]] +
    #       theme(plot.margin = margin(10, -4, -7, -10, "pt")
    #       )
    #   } else if (img_i %in% right_col) {
    #     p_scatter_mod_margins[[img_i]] = p_scatter_mod_margins[[img_i]] +
    #       theme(plot.margin = margin(10, 4, -7, -18, "pt")
    #       )
    #   }
    #   if (!(img_i %in% left_col)) {
    #     p_scatter_mod_margins[[img_i]] = p_scatter_mod_margins[[img_i]] +
    #       #scale_y_continuous(labels = function(breaks) {rep_along(breaks, "  ")}, limits = c(0, 1.15*maxi_expr_rna), expand = c(0,0)) +
    #       theme(axis.ticks.y = element_line(color = "transparent"),
    #             axis.text.y = element_text(color = "transparent")
    #       )
    #   }
    # }
    #   
    # all_p_scatter_manual = arrangeGrob(grobs=p_scatter_mod_margins, ncol=plot_ncols, nrow=plot_nrows, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90), right = textGrob(rna, gp=gpar(fontsize=plots_props$font_size_header, fontfamily=plots_props$font_family), rot = 270))
    # 
    # plot_width = 9
    # plot_height = 3.5
    #ggsave(sprintf("%s/%s_%sx%s.png", output_folder_scatter, paste(feature_list_new, collapse = "_"), plot_width, plot_height), all_p_scatter_manual, dpi=plots_props$dpi, width=plot_width, height=plot_height, units = plots_props$image_units)
    #ggsave(sprintf("%s/%s_%sx%s.png", output_folder_scatter, paste(feature_list_new, collapse = "_"), plot_width, plot_height), all_p_scatter_manual, width=plot_width, height=plot_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)


    # cleanup an call garbage cleaner
    rm(p_scatter)
    rm(p_scatter_manual)
    
    # rm(p_bar_manual)
    # rm(p_box_manual)
    # 
    # gc()
  }
  
  #stopCluster(cl)
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
