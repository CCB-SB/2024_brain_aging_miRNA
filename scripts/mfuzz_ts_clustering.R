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

#snakemake = readRDS("snakemake_mfuzz_ts_clustering.rds")

data.table::setDTthreads(1)
set.seed(snakemake@params$parameters_props$set_seed)

print("mfuzz_ts_clustering")


#---------------------------------- Functions ----------------------------------
find_elbow = function(scores){
  # as suggested here: https://www.datasciencecentral.com/profiles/blogs/how-to-automatically-determine-the-number-of-clusters-in-your-dat
  d1 = diff(scores)
  d2 = c(0, diff(d1))
  strength = d2-d1
  #plot(strength)
  return(which.max(strength))
}

# mieaa_go_gmt = strsplit(readLines("data/mieaa__GO_Annotations_indirect_mature.gmt"), split='\t')
# mieaa_go_list = lapply(1:length(mieaa_go_gmt), function(i){
#   mieaa_go_gmt[[i]][3:length(mieaa_go_gmt[[i]])]
# })
# names(mieaa_go_list) = unlist(lapply(1:length(mieaa_go_gmt), function(i) mieaa_go_gmt[[i]][1]))
# # drop categories with > 100 entries
# mieaa_go_list = mieaa_go_list[unlist(lapply(mieaa_go_list, length)) <= 100]
# mieaa_go_annot = annotationListToMatrix(mieaa_go_list, genenames = rownames(input_expr_mtx))
# bhis = lapply(cluster_range, function(k){
#   mf_cls = mfuzz(eset, c=k, m=m1)
#   clValid::BHI(mf_cls$cluster, mieaa_go_annot)
# })


compute_cluster_overview = function(median_table_df, data_input, output_folder_path_fig, output_folder_path_tab, cluster_range, threads, plots_props){
  
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
  
  m1 = mestimate(eset)
  
  avg_min_centroid_dist = unlist(pblapply(cluster_range, function(c) {
    set.seed(snakemake@params$parameters_props$set_seed)
    Dmin(eset, m1, crange=c, visu=FALSE, repeats = 5)
  }, cl=threads))
  
  cluster_results = pblapply(cluster_range, function(k){
    set.seed(snakemake@params$parameters_props$set_seed)
    mf_cls = mfuzz(eset, c=k, m=m1)
  }, cl=threads)
  
  other_metrics = data.table(do.call(rbind, lapply(cluster_results, function(cl_res){
    suppressWarnings(e1071::fclustIndex(cl_res, exprs(eset), index=c("xie.beni"))) #, "fukuyama.sugeno", "partition.coefficient", "partition.entropy")))
  })))
  other_metrics[, avg_min_centroid_dist:=avg_min_centroid_dist]
  other_metrics[, k:=cluster_range]
  #other_metrics[, penalized_avg_min_centroid_dist:=avg_min_centroid_dist]
  
  fwrite(other_metrics, file.path(output_folder_path_tab, "cluster_validity_measures.csv"))
  other_metrics_melt = melt(other_metrics, id.vars = c("k"))
  p = ggplot(other_metrics_melt, aes(x=k, y=value)) + facet_wrap(~variable, scales="free_y") + geom_point() + geom_line() + xlab("k") + ylab("Value") + theme_cowplot()
  save_plot(file.path(output_folder_path_fig, "cluster_validity_measures.pdf"), p, base_width = 210, base_height = 160, unit="mm")
  
  print(other_metrics)
  print(sprintf("highest gradient: %s", find_elbow(other_metrics$avg_min_centroid_dist)))
}

compute_clusters_and_plots = function(median_table_df, data_input, output_folder_path_fig, output_folder_path_tab, num_clusters, threads, plots_props){
  
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
  
  m1 = mestimate(eset)
  
  avg_min_centroid_dist = unlist(pblapply(num_clusters, function(c) {
    set.seed(snakemake@params$parameters_props$set_seed)
    Dmin(eset, m1, crange=c, visu=FALSE, repeats = 5)
  }, cl=threads))
  
  cluster_results = pblapply(num_clusters, function(k){
    set.seed(snakemake@params$parameters_props$set_seed)
    mf_cls = mfuzz(eset, c=k, m=m1)
  }, cl=threads)
  
  other_metrics = data.table(do.call(rbind, lapply(cluster_results, function(cl_res){
    suppressWarnings(e1071::fclustIndex(cl_res, exprs(eset), index=c("xie.beni"))) #, "fukuyama.sugeno", "partition.coefficient", "partition.entropy")))
  })))
  other_metrics[, avg_min_centroid_dist:=avg_min_centroid_dist]
  other_metrics[, k:=num_clusters]
  #other_metrics[, penalized_avg_min_centroid_dist:=avg_min_centroid_dist]
  
  for(i in 1:length(num_clusters)){
    members = acore(eset, cluster_results[[i]], min.acore=0)
    acore_dt = rbindlist(lapply(seq_along(members), function(i){ data.table(CLUSTER=i, members[[i]])}))
    k = num_clusters[i]
    fwrite(acore_dt, file=file.path(output_folder_path_tab, sprintf("cluster_members_k=%s.csv", k)))
    saveRDS(cluster_results[[i]], file = file.path(output_folder_path_tab, sprintf("cluster_result_k=%s.rds", k)))
  }
  
  cairo_pdf(file.path(output_folder_path_fig, "overlap_plots.pdf"),
            height=max(num_clusters) * 1.5, width=9)
  par(mfrow=c(ceiling(max(num_clusters)/4),4))
  for(i in 1:length(num_clusters)){
    ovl = overlap(cluster_results[[i]])
    overlap.plot(cluster_results[[i]], over=ovl, thres=0.1)
  }
  dev.off()
}


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
cluster_props = snakemake@params$cluster_props
plot_category = snakemake@params$plot_category
diff_exp_log = snakemake@params$diff_log
colours = snakemake@params$colors
plot_props = snakemake@params$plots_props
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_mfuzz_ts_clustering.rds")

expr = fread(sprintf("%s/data_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
diff_exp = fread(sprintf("%s/%s_%s/results_%s/matrices/diff_exp/diff_exp%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, diff_exp_log), sep='\t') 

expr_rownames = expr[[data_input$rna_class]]
expr[[data_input$rna_class]] = c()

output_folder_path_fig_overview = sprintf("%s/%s_%s/results_%s/figures/mfuzz_clustering/k=%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, cluster_props$max_num_of_clusters)
dir.create(output_folder_path_fig_overview, recursive=TRUE)
output_folder_path_tab_overview = sprintf("%s/%s_%s/results_%s/matrices/mfuzz_clustering/k=%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, cluster_props$max_num_of_clusters)
dir.create(output_folder_path_tab_overview, recursive=TRUE)

output_folder_path_fig = sprintf("%s/%s_%s/results_%s/figures/mfuzz_clustering/k=%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, cluster_props$num_of_clusters)
dir.create(output_folder_path_fig, recursive=TRUE)
output_folder_path_tab = sprintf("%s/%s_%s/results_%s/matrices/mfuzz_clustering/k=%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, cluster_props$num_of_clusters)
dir.create(output_folder_path_tab, recursive=TRUE)

threads = 20


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

if (FALSE) {
  cluster_range = seq(4,cluster_props$max_num_of_clusters,1)

  compute_cluster_overview(median_table_df, data_input, output_folder_path_fig_overview, output_folder_path_tab_overview, cluster_range, threads, plot_props)
}

num_clusters = cluster_props$num_of_clusters
compute_clusters_and_plots(median_table_df, data_input, output_folder_path_fig, output_folder_path_tab, num_clusters, threads, plot_props)


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
