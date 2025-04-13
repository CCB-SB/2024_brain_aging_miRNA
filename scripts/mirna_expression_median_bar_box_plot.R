suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(foreach))  # for parallelization
suppressPackageStartupMessages(library(doParallel))  # for parallelization


#snakemake = readRDS("snakemake_mirna_expression_median_bar_box_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("mirna_expression_median_bar_box_plot")


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
#saveRDS(snakemake, file = "snakemake_mirna_expression_median_bar_box_plot.rds")

expr = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 

output_folder_bar = sprintf("%s/%s_%s/results_%s/figures/expr/bar_plots/median_%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, params$comb_group)
dir.create(output_folder_bar, recursive=TRUE)
output_folder_box = sprintf("%s/%s_%s/results_%s/figures/expr/box_plots/%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, params$comb_group)
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


# for (rna in feature_list) {
foreach (rna=feature_list) %dopar% {
  expr_feature = expr[expr[[data_input$rna_class]] == rna,]
  
  plot_df = melt(expr_feature)
  plot_df[[params$comb_group]] = str_split_fixed(plot_df$variable, "_", 3)[,2]
  
  tissue_median = c()
  tissue_sd = c()
  for (tissue in unique(annot[[params$comb_group]])) {
    tissue_median = append(tissue_median, median(plot_df[plot_df[[params$comb_group]] == tissue,]$value))
    tissue_sd = append(tissue_sd, sd(plot_df[plot_df[[params$comb_group]] == tissue,]$value))
  }
  
  plot_df_bar = data.frame(tissue = unique(annot[[params$comb_group]]), tissue_median = tissue_median, tissue_sd = tissue_sd)
  
  # get the right order of tissues
  tissue_order = plot_df_bar[order(-plot_df_bar$tissue_median),]$tissue
  tissue_order_names = xticks_names[[params$comb_group]][tissue_order]
  
  # create bar plot
  p_bar = ggplot(plot_df_bar, aes(x = reorder(tissue, - tissue_median), y = tissue_median, fill = tissue)) + #plot_df_bar%>%group_by(params$comb_group)%>%summarise(value=sum(value))%>%
    geom_bar(stat="identity",  width=0.75) +
    scale_fill_manual(values = unlist(colors[[params$comb_group]])) +
    #geom_point(data = plot_df, aes_string(y = "value", x = comb_group), color = "darkgray") +
    geom_linerange(aes(ymin=tissue_median, ymax=tissue_median+tissue_sd), position=position_dodge(.9)) +
    scale_x_discrete(limits=tissue_order, labels=tissue_order_names) + # <- we need this to fix the order
    ggtitle(sprintf("%s", rna)) +
    xlab(sprintf("Median over all samples in a %s", xticks_names$categories[[params$comb_group]])) +
    ylab(sprintf("Expr (%s,\ndetection rate = %s%%\nper %s)", gsub("_norm", "", data_input$rna_class), gsub("p", "", data_input$detection_rate), xticks_names$categories[[params$comb_group]])) +
    #scale_y_continuous(limits = c(0, 1.1*(max(plot_df_bar$tissue_median) + plot_df_bar[plot_df_bar$time_median == which.max(plot_df_bar$tissue_median),]$tissue_sd)), expand = c(0,0)) +
    theme_classic() +
    theme(aspect.ratio=plots_props$image_height / plots_props$image_width, 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9), 
          legend.position="none", 
          text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
          plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header),
          legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size)
    )

  # save bar plot
  ggsave(sprintf("%s/%s.png", output_folder_bar, rna), p_bar, dpi=plots_props$dpi, width=plots_props$image_width, height=plots_props$image_height, units = plots_props$image_units)
  ggsave(sprintf("%s/%s.svg", output_folder_bar, rna), p_bar, width=plots_props$image_width, height=plots_props$image_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
  
  
  # create box plot
  # p_box = ggplot(plot_df, aes(x = reorder(comb_group, - value), y = value, fill = comb_group)) + #plot_df%>%group_by(params$comb_group)%>%summarise(value=sum(value))%>%
  p_box = ggplot(plot_df, aes_string(x = params$comb_group, y = "value", fill = params$comb_group)) + #plot_df%>%group_by(params$comb_group)%>%summarise(value=sum(value))%>%
    geom_boxplot(outlier.colour="darkgray", outlier.shape=16, outlier.size=1, notch=FALSE) +
    scale_fill_manual(values = unlist(colors[[params$comb_group]])) +
    scale_x_discrete(limits=tissue_order, labels=tissue_order_names) + # <- we need this to fix the order
    ggtitle(sprintf("%s", rna)) +
    xlab(sprintf("%s", xticks_names$categories[[params$comb_group]])) +
    ylab(sprintf("Expr (%s,\ndetection rate = %s%%\nper %s)", gsub("_norm", "", data_input$rna_class), gsub("p", "", data_input$detection_rate), xticks_names$categories[[params$comb_group]])) +
    theme_classic() +
    theme(aspect.ratio=plots_props$image_height / plots_props$image_width, 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9), 
          legend.position="none", 
          text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
          plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header),
          legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size)
    )
  
  # save box plot
  ggsave(sprintf("%s/%s.png", output_folder_box, rna), p_box, dpi=plots_props$dpi, width=plots_props$image_width, height=plots_props$image_height, units = plots_props$image_units)
  ggsave(sprintf("%s/%s.svg", output_folder_box, rna), p_box, width=plots_props$image_width, height=plots_props$image_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
  
  # cleanup an call garbage cleaner
  rm(p_bar)
  rm(p_box)
  
  gc()
}

# stop parallel stuff
stopCluster(cl)


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
