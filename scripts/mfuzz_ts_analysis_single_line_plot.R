options(bitmapType='cairo')
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
#suppressPackageStartupMessages(library(ComplexHeatmap))
#suppressPackageStartupMessages(library(viridisLite))
#suppressPackageStartupMessages(library(proxy))
suppressPackageStartupMessages(library(Mfuzz))
suppressPackageStartupMessages(library(data.table))
#suppressPackageStartupMessages(library(clValid))
suppressPackageStartupMessages(library(pbapply))
#suppressPackageStartupMessages(library(fcvalid))
suppressPackageStartupMessages(library(gridBase))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(openxlsx))

#snakemake = readRDS("snakemake_mfuzz_ts_analysis_single_line_plot.rds")

data.table::setDTthreads(1)
set.seed(snakemake@params$parameters_props$set_seed)

print("mfuzz_ts_analysis_single_line_plot")


#---------------------------------- Functions ----------------------------------
cluster_line_plots = function(median_table_df, cluster_result, acore_dt, cl, membership_plot_th, data_input, output_folder_path_fig, num_of_clusters, plot_props){
  
  font_size_scaling_factor = 0.12725
  
  rownames(median_table_df) = median_table_df$feature
  median_table_rownames = median_table_df$feature
  median_table_df$feature = c()
  
  eset = ExpressionSet(assayData=as.matrix(median_table_df))
  eset = standardise(eset)
  
  # if the variance for a feature_tissue row is 0 (for example all values = 1)
  # then the wir is NaN after standardization
  # we replace the NaN with 0
  # first, we have to get it out the ExpressionSet
  expr_zscore = exprs(eset)
  # replace values
  #na_mask = is.na(expr_zscore)
  #expr_zscore[na_mask]= 0
  expr_zscore = na.omit(expr_zscore)
  # then into the ExpressionSet again
  eset = ExpressionSet(assayData=as.matrix(expr_zscore))
  
  #pdf plots
  num_cluster = c()
  for (cluster in sort(unique(acore_dt$CLUSTER))){
    num_cluster[cluster] = length(acore_dt[acore_dt$CLUSTER == cluster]$CLUSTER)
  }

  if (membership_plot_th == 0) {
    #cairo_pdf(file.path(output_folder_path_fig, sprintf("cluster_%s_expression.pdf", cl)), height=29.7, width=21, onefile=TRUE)  
    #par(family = plot_props$font_family, mgp=c(1.5,0.5,0), mar = c(2.5,3,1,0.5))
    #mfuzz.plot2(eset, cl=cluster_result, mfrow=c(1,1), colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
    #dev.off()
    
    png(file.path(output_folder_path_fig, sprintf("cluster_%s_expression.png", cl)), width = 2/3 * plot_props$image_width, height = plot_props$image_height, units = plot_props$image_units, res = plot_props$dpi)
    par(family = plot_props$font_family, mgp=c(1.5,0.5,0), mar = c(2.5,3,1,0.5))
    mfuzz.plot2(eset, cl=cluster_result, single=cl, mfrow=c(1,1), colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
    dev.off()
    
    svglite(file.path(output_folder_path_fig, sprintf("cluster_%s_expression.svg", cl)), width = 2/3 * plot_props$image_width / 2.54, height = plot_props$image_height / 2.54)
    par(family = plot_props$font_family, mgp=c(1.5,0.5,0), mar = c(2.5,3,1,0.5))
    mfuzz.plot2(eset, cl=cluster_result, single=cl, mfrow=c(1,1), colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
    dev.off()
    
    # without_ylab
    png(file.path(output_folder_path_fig, sprintf("cluster_%s_expression_without_ylab.png", cl)), width = 2/3 * plot_props$image_width, height = plot_props$image_height, units = plot_props$image_units, res = plot_props$dpi)
    par(family = plot_props$font_family, mgp=c(1.5,0.5,0), mar = c(2.5,3,1,0.5))
    mfuzz.plot2(eset, cl=cluster_result, single=cl, mfrow=c(1,1), ylab="", colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
    dev.off()
    
    svglite(file.path(output_folder_path_fig, sprintf("cluster_%s_expression_without_ylab.svg", cl)), width = 2/3 * plot_props$image_width / 2.54, height = plot_props$image_height / 2.54)
    par(family = plot_props$font_family, mgp=c(1.5,0.5,0), mar = c(2.5,3,1,0.5))
    mfuzz.plot2(eset, cl=cluster_result, single=cl, mfrow=c(1,1),  ylab="", colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
    dev.off()
  } else {
  #cairo_pdf(file.path(output_folder_path_fig, sprintf("cluster_%s_expression.pdf", cl)), height=29.7, width=21, onefile=TRUE)  
  #par(family = plot_props$font_family, mgp=c(1.5,0.5,0), mar = c(2.5,3,1,0.5))
  #mfuzz.plot2(eset, cl=cluster_result, mfrow=c(1,1), colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, min.mem=membership_plot_th, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
  #mfuzz.plot2(eset, cl=cluster_results[[i]], colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, min.mem=0.5)
  #dev.off()
  
  width = 2/3 * plot_props$image_width
  height = plot_props$image_height
  png(file.path(output_folder_path_fig, sprintf("cluster_%s_expression_%sx%s.png", cl, width, height)), width = width, height = height, units = plot_props$image_units, res = plot_props$dpi)
  par(family = plot_props$font_family, mgp=c(1.5,0.5,0), mar = c(2.5,3,1,0.5))
  mfuzz.plot2(eset, cl=cluster_result, single=cl, mfrow=c(1,1), colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, min.mem=membership_plot_th, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
  dev.off()
  
  svglite(file.path(output_folder_path_fig, sprintf("cluster_%s_expression_%sx%s.svg", cl, width, height)), width = width / 2.54, height = height / 2.54)
  par(family = plot_props$font_family, mgp=c(1.5,0.5,0), mar = c(2.5,3,1,0.5))
  mfuzz.plot2(eset, cl=cluster_result, single=cl, mfrow=c(1,1), colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, min.mem=membership_plot_th, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
  dev.off()
  
  width = 2/3 * plot_props$image_width
  height = 2/3 * plot_props$image_height
  png(file.path(output_folder_path_fig, sprintf("cluster_%s_expression_%sx%s.png", cl, width, height)), width = width, height = height, units = plot_props$image_units, res = plot_props$dpi)
  par(family = plot_props$font_family, mgp=c(1.5,0.5,0), mar = c(2.5,3,1,0.5))
  mfuzz.plot2(eset, cl=cluster_result, single=cl, mfrow=c(1,1), colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, min.mem=membership_plot_th, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
  dev.off()
  
  svglite(file.path(output_folder_path_fig, sprintf("cluster_%s_expression_%sx%s.svg", cl, width, height)), width = width / 2.54, height = height / 2.54)
  par(family = plot_props$font_family, mgp=c(1.5,0.5,0), mar = c(2.5,3,1,0.5))
  mfuzz.plot2(eset, cl=cluster_result, single=cl, mfrow=c(1,1), colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, min.mem=membership_plot_th, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
  dev.off()
  
  # without_ylab
  width = 2/3 * plot_props$image_width
  height = plot_props$image_height
  png(file.path(output_folder_path_fig, sprintf("cluster_%s_expression_%sx%s_without_ylab.png", cl, width, height)), width = width, height = height, units = plot_props$image_units, res = plot_props$dpi)
  par(family = plot_props$font_family, mgp=c(1.5,0.5,0), mar = c(2.5,3,1,0.5))
  mfuzz.plot2(eset, cl=cluster_result, single=cl, mfrow=c(1,1), ylab="", colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, min.mem=membership_plot_th, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
  dev.off()
  
  svglite(file.path(output_folder_path_fig, sprintf("cluster_%s_expression_%sx%s_without_ylab.svg", cl, width, height)), width = width / 2.54, height = height / 2.54)
  par(family = plot_props$font_family, mgp=c(1.5,0.5,0), mar = c(2.5,3,1,0.5))
  mfuzz.plot2(eset, cl=cluster_result, single=cl, mfrow=c(1,1), ylab="", colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, min.mem=membership_plot_th, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
  dev.off()
  
  width = 2/3 * plot_props$image_width
  height = 2/3 * plot_props$image_height
  png(file.path(output_folder_path_fig, sprintf("cluster_%s_expression_%sx%s_without_ylab.png", cl, width, height)), width = width, height = height, units = plot_props$image_units, res = plot_props$dpi)
  par(family = plot_props$font_family, mgp=c(1.5,0.5,0), mar = c(2.5,3,1,0.5))
  mfuzz.plot2(eset, cl=cluster_result, single=cl, mfrow=c(1,1), ylab="", colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, min.mem=membership_plot_th, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
  dev.off()
  
  svglite(file.path(output_folder_path_fig, sprintf("cluster_%s_expression_%sx%s_without_ylab.svg", cl, width, height)), width = width / 2.54, height = height / 2.54)
  par(family = plot_props$font_family, mgp=c(1.5,0.5,0), mar = c(2.5,3,1,0.5))
  mfuzz.plot2(eset, cl=cluster_result, single=cl, mfrow=c(1,1), ylab="", colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, min.mem=membership_plot_th, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
  dev.off()
  }
}


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
cluster_props = snakemake@params$cluster_props
plot_category = snakemake@params$plot_category
diff_exp_log = snakemake@params$diff_log
thresholds = snakemake@params$thresholds
interesting_cluster = snakemake@params$interesting_cluster
colours = snakemake@params$colors
plot_props = snakemake@params$plots_props
results_folder = snakemake@params$results_folder
xticks_names = snakemake@params$xticks_names

# save rdata
#saveRDS(snakemake, file = "snakemake_mfuzz_ts_analysis_single_line_plot.rds")

expr = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
diff_exp = fread(sprintf("%s/%s_%s/results_%s/matrices/diff_exp/diff_exp%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, diff_exp_log), sep='\t') 

expr_rownames = expr[[data_input$rna_class]]
expr[[data_input$rna_class]] = c()

input_folder_path_tab = sprintf("%s/%s_%s/results_%s/matrices/mfuzz_clustering/k=%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, cluster_props$num_of_clusters)
cluster_result = readRDS(file.path(snakemake@input$cluster_results_folder_path, sprintf("cluster_result_k=%s.rds", cluster_props$num_of_clusters)))
acore_dt = fread(file.path(snakemake@input$cluster_results_folder_path, sprintf("cluster_members_k=%s.csv", cluster_props$num_of_clusters)))

output_folder_path_fig = sprintf("%s/%s_%s/results_%s/figures/mfuzz_clustering/k=%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, cluster_props$num_of_clusters)

#------------------------------------ Script ----------------------------------- 
median_table_df = data.frame()
for (tissue in unique(annot[[plot_category[1]]])){

  #print(tissue)
  
  median_table = data.table(feature=sprintf("%s_%s", diff_exp[[data_input$rna_class]], tissue))
  for (tp in paste(sort(unique(annot[[plot_category[2]]])))) {
    median_table[[tp]] = 2^(as.numeric(diff_exp[[sprintf("median_%s__%s_%s=%s", plot_category[2], tp, plot_category[1], tissue)]]))
  }
  
  median_table_df = rbind(median_table_df, data.frame(median_table, check.names = FALSE))
}

output_folder_path_fig_all = sprintf("%s/all/line_plots/single", output_folder_path_fig)
dir.create(output_folder_path_fig_all, recursive=TRUE)

for (cl in interesting_cluster) {
  cluster_line_plots(median_table_df, cluster_result, acore_dt, cl, 0, data_input, output_folder_path_fig_all, cluster_props$num_of_clusters, plot_props)
}

for (membership_plot_th in thresholds$membership_plot_ths) {
  output_folder_path_fig_th = sprintf("%s/%s/line_plots/single", output_folder_path_fig, membership_plot_th)
  dir.create(output_folder_path_fig_th, recursive=TRUE)
  
  for (cl in interesting_cluster) {
    cluster_line_plots(median_table_df, cluster_result, acore_dt, cl, membership_plot_th, data_input, output_folder_path_fig_th, cluster_props$num_of_clusters, plot_props)
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
