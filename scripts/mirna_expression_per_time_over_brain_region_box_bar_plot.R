suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(foreach))  # for parallelization
suppressPackageStartupMessages(library(doParallel))  # for parallelization

#snakemake = readRDS("snakemake_mirna_expression_per_time_over_brain_region_box_bar_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("mirna_expression_per_time_over_brain_region_box_bar_plot")


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
input_files = snakemake@params$files
feature_config = snakemake@params$features
params = snakemake@params$parameters 
plots_props = snakemake@params$plots_props
xticks_names = snakemake@params$xticks_names
colors = snakemake@params$colors
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_mirna_expression_per_time_over_brain_region_box_bar_plot.rds")

expr = fread(sprintf("%s/data_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
diff_exp = fread(sprintf("%s/%s_%s/results_%s/matrices/diff_exp/diff_exp%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, params$diff_log), sep='\t') 

output_folder_bar = sprintf("%s/%s_%s/results_%s/figures/expr/bar_plots/median_%s_per_%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, params$comb_group, params$plot_group)
dir.create(output_folder_bar, recursive=TRUE)
output_folder_box = sprintf("%s/%s_%s/results_%s/figures/expr/box_plots/%s_per_%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, params$comb_group, params$plot_group)
dir.create(output_folder_box, recursive=TRUE)


#------------------------------------ Script ----------------------------------- 
feature_list = c()
for (input_file in input_files) {
  feature_list = c(feature_list, fread(input_file)[[data_input$rna_class]])
}
feature_list = unique(c(feature_list, feature_config))

# Set up parallel computing
n_cores = 64
# print(sprintf("running in %s threads", n_cores))
cl <- makeCluster(n_cores, type="FORK")
registerDoParallel(cl)

#print(length(feature_list))
# for (rna in feature_list) {
foreach (rna=feature_list) %dopar% {
  #print(rna)
  expr_feature = expr[expr[[data_input$rna_class]] == rna,]
  diff_exp_feature = diff_exp[diff_exp[[data_input$rna_class]] == rna,]
  
  number_of_plots_possible = params$plot_nrows * params$plot_ncols
  
  plot_df_time = c()
  plot_df_bar = c()
  for (time in sort(unique(annot[[params$plot_group]]))){
    #print(time)
    id_group = annot[annot[[params$plot_group]] == time,][[data_input$identifier_column]] 
    tmp = colnames(expr_feature) %in% id_group
    expr_feature_time = expr_feature[,..tmp]
    
    plot_df = melt(expr_feature_time)
    plot_df[[params$plot_group]] = str_split_fixed(plot_df$variable, "_", 3)[,2]
    
    tissue_list = c()
    for (id in plot_df$variable) {
      tissue_list = append(tissue_list, annot[annot[[data_input$identifier_column]]== id,][[params$comb_group]])
    }
    plot_df[[params$comb_group]] = tissue_list
    plot_df_time[[paste(time)]] = plot_df
    # this line just adds a column with a fixed name "comb_group_name" to plot_df_time with the same values as plot_df_time[[paste(time)]][[params$comb_group]]
    plot_df_time[[paste(time)]][["comb_group_name"]] = plot_df_time[[paste(time)]][[params$comb_group]]
    
    tissue_median = c()
    tissue_sd = c()
    for (tissue in unique(annot[[params$comb_group]])) {
      tissue_median = append(tissue_median, median(plot_df[plot_df[[params$comb_group]] == tissue,]$value))
      tissue_sd = append(tissue_sd, sd(plot_df[plot_df[[params$comb_group]] == tissue,]$value))
      tissue_sd[is.na(tissue_sd)] = 0
    }

    plot_df_bar[[paste(time)]] = data.frame(tissue = unique(annot[[params$comb_group]]), tissue_median = tissue_median, tissue_sd = tissue_sd)
  }
  
  x_axis_title = xticks_names$categories[[sprintf("%s_capital", params$comb_group)]]
  #y_axis_title_big = sprintf("expr (%s, detection rate = %s%% per %s)", gsub("_norm", "", data_input$norm), gsub("p", "", data_input$detection_rate), xticks_names$categories[[params$comb_group]])
  y_axis_title_small = sprintf("Expr (%s,\ndetection rate = %s%% per %s)", gsub("_norm", "", data_input$norm), gsub("p", "", data_input$detection_rate), xticks_names$categories[[params$comb_group]])
  
  p_bar = list()
  p_box = list()
  i = 1
  for (time in sort(unique(annot[[params$plot_group]]))){
    # create bar plot
    p_bar[[i]] = ggplot(plot_df_bar[[paste(time)]], aes(x = tissue, y = tissue_median, fill = as.factor(tissue))) + 
      geom_bar(stat="identity",  width=0.9) +
      geom_linerange(aes(ymin=tissue_median, ymax=tissue_median+tissue_sd), position=position_dodge(.9)) +
      #geom_text(aes(label = significance), vjust=-2, size=5) +
      scale_fill_manual(values = unlist(colors[[params$comb_group]])) +
      geom_jitter(data = plot_df_time[[paste(time)]], aes_string(y = "value", x = params$comb_group), color = "darkgray", size = 0.25) +
      #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 3), se = FALSE, fullrange=TRUE, color = "black", size = 0.75) +
      ggtitle(xticks_names[[params$plot_group]][[paste(time)]]) +
      xlab("") +
      ylab("") +
      #ylab(sprintf("expr (%s, detection rate = %s%%\nper %s)", data_input$norm, expr_input$detection_rate, params$plot_group)) +
      #scale_x_continuous(limits = c(0, 30), breaks=c(3, 10, 20, 30),) +
      scale_y_continuous(limits = c(0, 1.2*max(do.call("rbind", plot_df_time)$value)), expand = c(0,0)) +
      theme_classic() +
      theme(#plot.margin = unit(c(0.1,-1.9,0.25,-1.75), plots_props$image_units), #top, right, bottom, left
            plot.margin = unit(c(0.1,0.1,0,0.1), plots_props$image_units), #top, right, bottom, left
            #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9), 
            plot.background = element_rect(fill='transparent', color=NA),
            legend.position="none", 
            text = element_text(family = plots_props$font_family, size = plots_props$font_size),
            axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
            axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
            axis.text.x = element_text(angle=90, vjust=.5, hjust=1, size = 6),
            plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header),
            legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
            legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size)
      )

    # create box plot
    # p_box[[i]] = ggplot(plot_df_time, aes_string(x = params$comb_group, y = "value", group = params$comb_group, fill = as.factor(params$comb_group))) + #plot_df[[paste(time)]]%>%group_by(params$comb_group)%>%summarise(value=sum(value))%>%
    p_box[[i]] = ggplot(plot_df_time[[paste(time)]], aes(x = comb_group_name, y = value, group = comb_group_name, fill = as.factor(comb_group_name))) + #plot_df[[paste(time)]]%>%group_by(params$comb_group)%>%summarise(value=sum(value))%>%
      geom_boxplot(width = 0.9, outlier.colour="grey", outlier.shape=16, outlier.size=1, notch=FALSE) +
      scale_fill_manual(values = unlist(colors[[params$comb_group]])) +
      #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 3), se = FALSE, fullrange=TRUE, color = "black", size = 0.75) +
      # geom_text(aes(label = significance), vjust=-2, size=5) +. # <- was ist significance hier? wie kommt die Spalte in plot_df_time?
      ggtitle(xticks_names[[params$plot_group]][[paste(time)]]) +
      xlab("") +
      ylab("") +
      #ylab(sprintf("expr (%s, detection rate = %s%%\nper %s)", data_input$norm, expr_input$detection_rate, params$plot_group)) +
      #scale_x_continuous(limits = c(0, 30), breaks=c(3, 10, 20, 30)) +
      scale_y_continuous(limits = c(0, 1.2*max(do.call("rbind", plot_df_time)$value)), expand = c(0,0)) +
      theme_classic() +
      theme(#plot.margin = unit(c(0.1,-1.9,0.25,-1.75), plots_props$image_units), #top, right, bottom, left
            plot.margin = unit(c(0.1,0.1,0,0.1), plots_props$image_units), #top, right, bottom, left
            plot.background = element_rect(fill='transparent', color=NA),
            legend.position="none", 
            text = element_text(family = plots_props$font_family, size = plots_props$font_size),
            axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
            axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6), 
            plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header),
            legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
            legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size)
      )
    i = i + 1
  }
  # creating concatinated plot
  all_p_bar = arrangeGrob(grobs=p_bar, ncol=params$plot_ncols, nrow=params$plot_nrows, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title_small, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
  all_p_box = arrangeGrob(grobs=p_box, ncol=params$plot_ncols, nrow=params$plot_nrows, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title_small, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))

  # save bar plot
  ggsave(sprintf("%s/%s.png", output_folder_bar, rna), all_p_bar, dpi=plots_props$dpi, width=2*plots_props$image_width, height=plots_props$image_height, units = plots_props$image_units)
  ggsave(sprintf("%s/%s.svg", output_folder_bar, rna), all_p_bar, width=2*plots_props$image_width, height=plots_props$image_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
  
  # save box plot
  ggsave(sprintf("%s/%s.png", output_folder_box, rna), all_p_box, dpi=plots_props$dpi, width=2*plots_props$image_width, height=plots_props$image_height, units = plots_props$image_units)
  ggsave(sprintf("%s/%s.svg", output_folder_box, rna), all_p_box, width=2*plots_props$image_width, height=plots_props$image_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)

  # cleanup an call garbage cleaner
  rm(p_bar)
  rm(p_box)
  
  rm(all_p_bar)
  rm(all_p_box)
  
  gc()
}

stopCluster(cl)


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
