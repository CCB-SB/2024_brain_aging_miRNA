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

source("scripts/mfuzz_ggplot.R")

#snakemake = readRDS("snakemake_mfuzz_ts_analysis_overview_line_plot.rds")

data.table::setDTthreads(1)
set.seed(snakemake@params$parameters_props$set_seed)

print("mfuzz_ts_analysis_overview_line_plot")


#---------------------------------- Functions ----------------------------------
cluster_line_plots_overview = function(median_table_df, cluster_result, acore_dt, membership_plot_ths, data_input, output_folder_path_fig, output_folder_path_tab, num_of_clusters, plot_props){
  
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
  
  num_cluster = c()
  for (cluster in sort(unique(acore_dt$CLUSTER))){
    num_cluster[cluster] = length(acore_dt[acore_dt$CLUSTER == cluster]$CLUSTER)
  }
  ovl = overlap(cluster_result)
  
  # pdf plot containing all clusters and a membership overview
  cairo_pdf(file.path(sprintf(output_folder_path_fig), "cluster_expression_and_overlap.pdf"), height=29.7, width=21, onefile=TRUE)
  par(family = plot_props$font_family)
  mfuzz.plot2(eset, cl=cluster_result, mfrow=c(9,7), colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
  #mfuzz.plot2(eset, cl=cluster_results[[i]], colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE)
  overlap.plot(cluster_result, over=ovl, thres=0.1)
  boxplot(MEM.SHIP~CLUSTER, acore_dt)
  for (membership_plot_th in membership_plot_ths) {
    barplot(V1~CLUSTER, acore_dt[, mean(MEM.SHIP >= membership_plot_th)*100, by="CLUSTER"], ylab = sprintf("Genes with membership >= %s (%%)", membership_plot_th))
  }
  b = barplot(num_cluster, xlab = "CLUSTER", ylab = "Genes with cluster membership", width = 2, names = sort(unique(acore_dt$CLUSTER))) #ylim = c(0, max(num_cluster + 10)
  text(x=b, y= num_cluster+5, labels=as.character(num_cluster))
  
  dev.off()
  
  #par(family = plot_props$font_family)
  #num_cluster = c()
  #for (cluster in sort(unique(acore_dt$CLUSTER))){
  #  num_cluster[cluster] = length(acore_dt[acore_dt$CLUSTER == cluster]$CLUSTER)
  #}
  #ovl = overlap(cluster_result)
  
  # only cluster plots
  cairo_pdf(file.path(output_folder_path_fig, "cluster_expression.pdf"), height=29.7, width=21, onefile=TRUE)
  par(family = plot_props$font_family)
  mfuzz.plot2(eset, cl=cluster_result, mfrow=c(9,7), colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
  #mfuzz.plot2(eset, cl=cluster_results[[i]], colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE)
  boxplot(MEM.SHIP~CLUSTER, acore_dt)
  dev.off()
  
  png(file.path(output_folder_path_fig, "cluster_expression.png"), width = 2 * plot_props$image_width, height = ceiling(num_of_clusters * 1.5), units = plot_props$image_units, res = plot_props$dpi)
  par(family = plot_props$font_family)
  mfuzz.plot2(eset, cl=cluster_result, mfrow=c(16,4), colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
  boxplot(MEM.SHIP~CLUSTER, acore_dt)
  dev.off()
  svglite(file.path(output_folder_path_fig, "cluster_expression.svg"), width = 2 * plot_props$image_width / 2.54, height = ceiling(num_of_clusters * 1.5) / 2.54)
  par(family = plot_props$font_family)
  mfuzz.plot2(eset, cl=cluster_result, mfrow=c(16,4), colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
  boxplot(MEM.SHIP~CLUSTER, acore_dt)
  dev.off()
  
  svglite(file.path(output_folder_path_fig, "cluster_expression_adj.svg"), width = 2 * plot_props$image_width / 2.54, height = ceiling(num_of_clusters * 1.5) / 2.54)
  par(family = plot_props$font_family)
  mfuzz.plot2(eset, cl=cluster_result, mfrow=c(10,6), colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
  #boxplot(MEM.SHIP~CLUSTER, acore_dt)
  dev.off()
}

cluster_line_plots = function(median_table_df, cluster_result, acore_dt, membership_plot_th, data_input, output_folder_path_fig, output_folder_path_tab, num_of_clusters, plot_props){
  
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
  
  ovl = overlap(cluster_result)
  
  # pdf with all clusters showing only the ones with a membership higher than a given threshold
  cairo_pdf(file.path(output_folder_path_fig, sprintf("cluster_expression_and_overlap.pdf", membership_plot_th)), height=29.7, width=21, onefile=TRUE)  
  par(family = plot_props$font_family)
  mfuzz.plot2(eset, cl=cluster_result, mfrow=c(9,7), colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, min.mem=membership_plot_th, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
  barplot(V1~CLUSTER, acore_dt[, mean(MEM.SHIP >= membership_plot_th)*100, by="CLUSTER"], ylab = sprintf("Genes with membership >= %s (%%)", membership_plot_th))
  #mfuzz.plot2(eset, cl=cluster_results[[i]], colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, min.mem=0.5)
  overlap.plot(cluster_result, over=ovl, thres=0.1)
  dev.off()
  
  cairo_pdf(file.path(output_folder_path_fig, sprintf("cluster_expression.pdf", membership_plot_th)),
            height=29.7, width=21, onefile=TRUE)  
  par(family = plot_props$font_family)
  mfuzz.plot2(eset, cl=cluster_result, mfrow=c(9,7), colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, min.mem=membership_plot_th, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
  barplot(V1~CLUSTER, acore_dt[, mean(MEM.SHIP >= membership_plot_th)*100, by="CLUSTER"], ylab = sprintf("Genes with membership >= %s (%%)", membership_plot_th))
  #mfuzz.plot2(eset, cl=cluster_results[[i]], colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, min.mem=0.5)
  dev.off()
  
  png(file.path(output_folder_path_fig, sprintf("cluster_expression.png", membership_plot_th)), width = 2 * plot_props$image_width, height = ceiling(num_of_clusters * 1.5), units = plot_props$image_units, res = plot_props$dpi)
  par(family = plot_props$font_family)
  mfuzz.plot2(eset, cl=cluster_result, mfrow=c(16,4), colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, min.mem=membership_plot_th, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
  barplot(V1~CLUSTER, acore_dt[, mean(MEM.SHIP >= membership_plot_th)*100, by="CLUSTER"], ylab = sprintf("Genes with membership >= %s (%%)", membership_plot_th))
  dev.off()
  svglite(file.path(output_folder_path_fig, sprintf("cluster_expression.svg", membership_plot_th)), width = 2 * plot_props$image_width / 2.54, height = ceiling(num_of_clusters * 1.5) / 2.54)
  par(family = plot_props$font_family)
  mfuzz.plot2(eset, cl=cluster_result, mfrow=c(16,4), colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, min.mem=membership_plot_th, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
  barplot(V1~CLUSTER, acore_dt[, mean(MEM.SHIP >= membership_plot_th)*100, by="CLUSTER"], ylab = sprintf("Genes with membership >= %s (%%)", membership_plot_th))
  dev.off()
  
  svglite(file.path(output_folder_path_fig, sprintf("cluster_expression_adj.svg", membership_plot_th)), width = 2 * plot_props$image_width / 2.54, height = ceiling(num_of_clusters * 1.5) / 2.54)
  par(family = plot_props$font_family)
  mfuzz.plot2(eset, cl=cluster_result, mfrow=c(16,4), colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, min.mem=membership_plot_th, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
  #barplot(V1~CLUSTER, acore_dt[, mean(MEM.SHIP >= membership_plot_th)*100, by="CLUSTER"], ylab = sprintf("Genes with membership >= %s (%%)", membership_plot_th))
  dev.off()
}

cluster_overview = function(median_table_df, acore_dt, feature, feature_list, annot, plot_category, num_of_clusters, data_input, output_folder_path_tab, file_name){
  
  acore_dt[[data_input$rna_class]] = str_split_fixed(acore_dt$NAME, "_", 2)[,1]
  acore_dt[[plot_category[1]]] = str_split_fixed(acore_dt$NAME, "_", 2)[,2]
  
  cluster_info_df = c()
  num_of_features = c()
  for (cluster in 1:num_of_clusters) {
    acore_dt_cluster = acore_dt[acore_dt$CLUSTER == cluster,]
    
    cluster_info = c()
    for (feat in unique(feature_list)) {
      mask = (acore_dt_cluster[[feature]] == feat)
      acore_dt_cluster_feat = acore_dt_cluster[mask,]
      cluster_info = rbind(cluster_info, c(dim(acore_dt_cluster_feat)[1], dim(acore_dt_cluster_feat)[1] / dim(acore_dt_cluster)[1]))
    }
    cluster_info_df = rbind(cluster_info_df, data.frame(feature = unique(feature_list), absolute_features = cluster_info[,1], relative_features = cluster_info[,2], CLUSTER = cluster))
  
    num_of_features = append(num_of_features, sum(cluster_info[,1]))
  }
  
  # save cluster overviews
  write.table(cluster_info_df, sprintf("%s/%s.csv", output_folder_path_tab, file_name), sep = "\t", row.names = FALSE, col.names = TRUE)
  write.xlsx(cluster_info_df, sprintf("%s/%s.xlsx", output_folder_path_tab, file_name), colNames = TRUE, rowNames = FALSE, append = FALSE)
  
  num_of_features_df = data.frame(CLUSTER=1:num_of_clusters, num_features=num_of_features)
  
  write.table(num_of_features_df, sprintf("%s/cluster_overview.csv", output_folder_path_tab), sep = "\t", row.names = FALSE, col.names = TRUE)
  write.xlsx(num_of_features_df, sprintf("%s/cluster_overview.xlsx", output_folder_path_tab), colNames = TRUE, rowNames = FALSE, append = FALSE)
  return(cluster_info_df)
}

cluster_line_plots_with_tissue_specific_info = function(median_table_df, cluster_result, acore_dt, membership_plot_th, specific_tissues, data_input, output_folder_path_fig, output_folder_path_tab, num_of_clusters, plot_props){
  
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
  
  label_list = c()
  for (cl in 1:num_of_clusters) {
    mask  =  (specific_tissues$CLUSTER == cl)
    specific_tissues_cl = specific_tissues[mask,]
    if (dim(specific_tissues_cl)[1] == 1) {
      label_list = append(label_list, sprintf("%s: %s%%", specific_tissues_cl$feature, specific_tissues_cl$relative_features))
    } else if (dim(specific_tissues_cl)[1] > 1) {
      label_elem = ""
      for (row in 1:dim(specific_tissues_cl)[1]) {
        specific_tissues_cl_row = specific_tissues_cl[row,]
        label_elem = sprintf("%s%s: %s%%", label_elem, specific_tissues_cl_row$feature, specific_tissues_cl_row$relative_features)
        if (row != dim(specific_tissues_cl)[1]) {
          label_elem = sprintf("%s\n", label_elem)
        }
      }
      label_list = append(label_list, label_elem)
    } else {
      label_list = append(label_list, "")
    }
  }
  
  #pdf plots
  cairo_pdf(file.path(output_folder_path_fig, "cluster_expression_tissue_specific_info.pdf"),
            height=29.7, width=21, onefile=TRUE)
  par(family = plot_props$font_family)
  num_cluster = c()
  for (cluster in sort(unique(acore_dt$CLUSTER))){
    num_cluster[cluster] = length(acore_dt[acore_dt$CLUSTER == cluster]$CLUSTER)
  }
  #for (i in 1:num_of_clusters) {
  mfuzz.plot2(eset, cl=cluster_result, mfrow=c(9,7), colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor) #main = label_list
  #  text(x = c(i, 2), y = c(3, 2), labels = label_list[i], col = "red")
  #}
  dev.off()
  
  cairo_pdf(file.path(output_folder_path_fig, sprintf("cluster_expression_tissue_specific_info.pdf", membership_plot_th)),
            height=29.7, width=21, onefile=TRUE)  
  par(family = plot_props$font_family)
  mfuzz.plot2(eset, cl=cluster_result, mfrow=c(9,7), colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, min.mem=membership_plot_th, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
  #mfuzz.plot2(eset, cl=cluster_results[[i]], colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, min.mem=0.5)
  dev.off()
  
  
  par(family = plot_props$font_family)
  png(file.path(output_folder_path_fig, "cluster_expression_tissue_specific_info.png"), width = 2 * plot_props$image_width, height = ceiling(num_of_clusters * 1.5), units = plot_props$image_units, res = plot_props$dpi)
  mfuzz.plot2(eset, cl=cluster_result, mfrow=c(16,4), colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
  dev.off()
  svglite(file.path(output_folder_path_fig, "cluster_expression_tissue_specific_info.svg"), width = 2 * plot_props$image_width / 2.54, height = ceiling(num_of_clusters * 1.5) / 2.54)
  mfuzz.plot2(eset, cl=cluster_result, mfrow=c(16,4), colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
  dev.off()
  png(file.path(output_folder_path_fig, sprintf("cluster_expression_tissue_specific_info.png", membership_plot_th)), width = 2 * plot_props$image_width, height = ceiling(num_of_clusters * 1.5), units = plot_props$image_units, res = plot_props$dpi)
  mfuzz.plot2(eset, cl=cluster_result, mfrow=c(16,4), colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, min.mem=membership_plot_th, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
  dev.off()
  svglite(file.path(output_folder_path_fig, sprintf("cluster_expression_tissue_specific_info.svg", membership_plot_th)), width = 2 * plot_props$image_width / 2.54, height = ceiling(num_of_clusters * 1.5) / 2.54)
  mfuzz.plot2(eset, cl=cluster_result, mfrow=c(16,4), colo="fancy", time.labels=colnames(median_table_df), x11=FALSE, centre=FALSE, min.mem=membership_plot_th, cex.main = plot_props$font_size_header * font_size_scaling_factor, cex.lab = plot_props$font_size * font_size_scaling_factor, cex.axis = plot_props$font_size * font_size_scaling_factor)
  dev.off()
}


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
cluster_props = snakemake@params$cluster_props
plot_category = snakemake@params$plot_category
diff_exp_log = snakemake@params$diff_log
thresholds = snakemake@params$thresholds
colours = snakemake@params$colors
plot_props = snakemake@params$plots_props
results_folder = snakemake@params$results_folder
xticks_names = snakemake@params$xticks_names

# save rdata
#saveRDS(snakemake, file = "snakemake_mfuzz_ts_analysis_overview_line_plot.rds")

expr = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
diff_exp = fread(sprintf("%s/%s_%s/results_%s/matrices/diff_exp/diff_exp%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, diff_exp_log), sep='\t') 

feature_list = fread(snakemake@input$feature_list, sep='\t')  # , colClasses=c(ID="character")

expr_rownames = expr[[data_input$rna_class]]
expr[[data_input$rna_class]] = c()

input_folder_path_tab = sprintf("%s/%s_%s/results_%s/matrices/mfuzz_clustering/k=%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, cluster_props$num_of_clusters)
cluster_result = readRDS(file.path(snakemake@input$cluster_results_folder_path, sprintf("cluster_result_k=%s.rds", cluster_props$num_of_clusters)))
acore_dt = fread(file.path(snakemake@input$cluster_results_folder_path, sprintf("cluster_members_k=%s.csv", cluster_props$num_of_clusters)))

output_folder_path_fig = sprintf("%s/%s_%s/results_%s/figures/mfuzz_clustering/k=%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, cluster_props$num_of_clusters)
output_folder_path_tab = input_folder_path_tab

  
#------------------------------------ Script ----------------------------------- 
# get medians from the diff exp table
median_table_df = data.frame()
for (tissue in unique(annot[[plot_category[1]]])){
  median_table = data.table(feature=sprintf("%s_%s", diff_exp[[data_input$rna_class]], tissue))
  for (tp in paste(sort(unique(annot[[plot_category[2]]])))) {
    median_table[[tp]] = 2^(as.numeric(diff_exp[[sprintf("median_%s__%s_%s=%s", plot_category[2], tp, plot_category[1], tissue)]]))
  }
  median_table_df = rbind(median_table_df, data.frame(median_table, check.names = FALSE))
}

output_folder_path_fig_all = sprintf("%s/all/line_plots", output_folder_path_fig)
dir.create(output_folder_path_fig_all, recursive=TRUE)
output_folder_path_tab_all = sprintf("%s/all", output_folder_path_tab)
dir.create(output_folder_path_tab_all, recursive=TRUE)

# create cluster overview for all clustered elments and give an overview of the membership
cluster_line_plots_overview(median_table_df, cluster_result, acore_dt, thresholds$membership_plot_ths, data_input, output_folder_path_fig_all, output_folder_path_tab_all, cluster_props$num_of_clusters, plot_props)

# tissue containment in each cluster
file_name = sprintf("cluster_overview_%s", plot_category[1])
cluster_info_tissue_df = cluster_overview(median_table_df, acore_dt, plot_category[1], annot[[plot_category[1]]], annot, plot_category, cluster_props$num_of_clusters, data_input, output_folder_path_tab_all, file_name)

# feature occurances in each cluster
median_table_rownames = median_table_df$feature
all_features = str_split_fixed(median_table_rownames, "_", 2)[,1]
file_name = sprintf("cluster_overview_%s_occ", data_input$rna_class)
cluster_info_feature_df = cluster_overview(median_table_df, acore_dt, paste(data_input$rna_class), all_features, annot, plot_category, cluster_props$num_of_clusters, data_input, output_folder_path_tab_all, file_name)

# NOT WORKING because the plotting packages is not working in that way
# cluster overview plots showing the aggregated tissue and rna occ in the left and right corners resp.
#cluster_line_plots_with_tissue_specific_info(median_table_df, cluster_result, acore_dt, th, specific_tissues, data_input, output_folder_path_fig_all, output_folder_path_tab_all, cluster_props$num_of_clusters, plot_props)

# write tables for the given rnas if their occ fulfills the given tissue threshold
for (rna in feature_list[[data_input$rna_class]]) {
  print(rna)
  for(th in thresholds$rna_occ_ths) {
    mask = ((cluster_info_feature_df$feature == rna) & (cluster_info_feature_df$absolute_features > th))
    cluster_info_features_df_rna = cluster_info_feature_df[mask,]
    
    specific_tissues = c()
    for (cluster in 1:cluster_props$num_of_clusters) {
      mask = ((cluster_info_tissue_df$CLUSTER == cluster) & (cluster_info_tissue_df$relative_features >= th))
      cluster_info_tissue_df_cluster = cluster_info_tissue_df[mask,]
      specific_tissues = rbind(specific_tissues, cluster_info_tissue_df_cluster)
    }
    
    cluster_info_features_df_rna_spec_ti = merge(cluster_info_features_df_rna, specific_tissues, by = "CLUSTER")
    
    write.table(cluster_info_features_df_rna_spec_ti, sprintf("%s/%s_occurance_greater_%s_in_cluster.csv", output_folder_path_tab_all, rna, th), sep = "\t", row.names = FALSE, col.names = TRUE)
    write.xlsx(cluster_info_features_df_rna_spec_ti, sprintf("%s/%s_occurance_greater_%s_in_cluster.xlsx", output_folder_path_tab_all, rna, th), colNames = TRUE, rowNames = FALSE, append = FALSE)
  }
}

for (membership_plot_th in thresholds$membership_plot_ths) {
  # create output folders
  output_folder_path_fig_th = sprintf("%s/%s/line_plots", output_folder_path_fig, membership_plot_th)
  dir.create(output_folder_path_fig_th, recursive=TRUE)
  output_folder_path_tab_th = sprintf("%s/%s", output_folder_path_tab, membership_plot_th)
  dir.create(output_folder_path_tab_th, recursive=TRUE)
  
  # create cluster plots only for the elements fulfilling the membership threshold
  cluster_line_plots(median_table_df, cluster_result, acore_dt, membership_plot_th, data_input, output_folder_path_fig_th, output_folder_path_tab_th, cluster_props$num_of_clusters, plot_props)

  # filter for membership threshold
  acore_dt_filtered = acore_dt[acore_dt$MEM.SHIP >= membership_plot_th]
  
  # tissue containment in each cluster
  file_name = sprintf("cluster_overview_%s", plot_category[1])
  cluster_info_tissue_df = cluster_overview(median_table_df, acore_dt_filtered, plot_category[1], annot[[plot_category[1]]], annot, plot_category, cluster_props$num_of_clusters, data_input, output_folder_path_tab_th, file_name)
  
  # feature occurances in each cluster
  file_name = sprintf("cluster_overview_%s_occ", data_input$rna_class)
  cluster_info_feature_df = cluster_overview(median_table_df, acore_dt_filtered, paste(data_input$rna_class), all_features, annot, plot_category, cluster_props$num_of_clusters, data_input, output_folder_path_tab_th, file_name)

  # NOT WORKING because the plotting packages is not working in that way
  # cluster overview plots showing the aggregated tissue and rna occ in the left and right corners resp.
  #cluster_line_plots_with_tissue_specific_info(median_table_df, cluster_result, acore_dt, membership_plot_th, specific_tissues, data_input, output_folder_path_fig_th, output_folder_path_tab_th, cluster_props$num_of_clusters, plot_props)
  
  # write tables for the given rnas if their occ fulfills the given tissue threshold
  for (rna in feature_list[[data_input$rna_class]]) {
    #print(rna)
    for(th in thresholds$rna_occ_ths) {
      mask = ((cluster_info_feature_df$feature == rna) & (cluster_info_feature_df$absolute_features > th))
      cluster_info_features_df_rna = cluster_info_feature_df[mask,]
      
      cluster_info_features_df_rna_spec_ti = merge(cluster_info_features_df_rna, specific_tissues, by = "CLUSTER")
      
      write.table(cluster_info_features_df_rna_spec_ti, sprintf("%s/%s_occurance_greater_%s_in_cluster.csv", output_folder_path_tab_th, rna, th), sep = "\t", row.names = FALSE, col.names = TRUE)
      write.xlsx(cluster_info_features_df_rna_spec_ti, sprintf("%s/%s_occurance_greater_%s_in_cluster.xlsx", output_folder_path_tab_th, rna, th), colNames = TRUE, rowNames = FALSE, append = FALSE)
    }
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
