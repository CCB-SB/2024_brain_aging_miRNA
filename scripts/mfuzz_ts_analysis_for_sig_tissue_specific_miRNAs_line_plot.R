options(bitmapType='cairo')
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
##suppressPackageStartupMessages(library(ComplexHeatmap))
##suppressPackageStartupMessages(library(viridisLite))
##suppressPackageStartupMessages(library(proxy))
suppressPackageStartupMessages(library(Mfuzz)) #v2.62.0
suppressPackageStartupMessages(library(data.table))
##suppressPackageStartupMessages(library(clValid))
suppressPackageStartupMessages(library(pbapply))
##suppressPackageStartupMessages(library(fcvalid))
suppressPackageStartupMessages(library(gridBase)) #v0.4_7
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(ggrepel))

#snakemake = readRDS("snakemake_mfuzz_ts_analysis_for_sig_tissue_specific_miRNAs_line_plot.rds")

data.table::setDTthreads(1)
set.seed(snakemake@params$parameters_props$set_seed)

print("mfuzz_ts_analysis_for_sig_tissue_specific_miRNAs_line_plot")


#---------------------------------- Functions ----------------------------------
datatable2mtx = function(dt){
  res = as.matrix(dt[, 2:ncol(dt)])
  rownames(res) = dt[[1]]
  return(res)
}

find_elbow = function(scores){
  # as suggested here: https://www.datasciencecentral.com/profiles/blogs/how-to-automatically-determine-the-number-of-clusters-in-your-dat
  d1 = diff(scores)
  d2 = c(0, diff(d1))
  strength = d2-d1
  plot(strength)
  which.max(strength)
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

reashape_features = function(data) {

  diff_exp_clean = data$diff_exp[data$diff_exp != ""]
  corr_clean = data$corr[data$corr != ""]
  
  merged_list = unique(c(diff_exp_clean, corr_clean))
  result = data.frame("feature" = merged_list, "method"="tmp")

  for (rna in merged_list) {
    if ((rna %in% diff_exp_clean == TRUE) & (rna %in% corr_clean == TRUE)) {
      result[result[[data_input$rna_class]] == rna,]$method = "both"
    } else if ((rna %in% diff_exp_clean == TRUE) & (rna %in% corr_clean == FALSE)) {
      # If only one value is empty, add the non-empty value to the result list
      result[result[[data_input$rna_class]] == rna,]$method = "diff_exp"
    } else {
      result[result[[data_input$rna_class]] == rna,]$method = "corr"
    }
  }
  
  return(result)
}

compute_clusters_and_plots = function(median_table, feature_distribution_reshaped, data_input, output_folder_path_fig, output_folder_path_tab, cluster_range, threads, plots_props){
  
  colnames(feature_distribution_reshaped) = c("Var1", "method")
  
  rownames(median_table) = median_table[[data_input$rna_class]]
  median_table_rownames = median_table[[data_input$rna_class]]
  median_table[[data_input$rna_class]] = c()
  
  eset = ExpressionSet(assayData=as.matrix(median_table))
  eset = standardise(eset)
  
  # if the variance for a feature_tissue row is 0 (for example all values = 1)
  # then the wir is NaN after standardization
  # we replace the NaN with 0
  # first, we have to get it out the ExpressionSet
  expr_zscore = exprs(eset)
  # replace values
  #na_mask = is.na(expr_zscore)
  #expr_zscore[na_mask]= 0
  if (dim(na.omit(expr_zscore))[1] < dim(expr_zscore)[1]) {
    print(expr_zscore)
  }
  expr_zscore = na.omit(expr_zscore)
  # then into the ExpressionSet again
  eset = ExpressionSet(assayData=as.matrix(expr_zscore))
  
  m1 = mestimate(eset)
  
  avg_min_centroid_dist = unlist(pblapply(cluster_range, function(c) {
    set.seed(42)
    Dmin(eset, m1, crange=c, visu=FALSE, repeats = 5)
  }, cl=threads))
  
  cluster_results = pblapply(cluster_range, function(k){
    set.seed(42)
    mf_cls = mfuzz(eset, c=k, m=m1)
  }, cl=threads)
  
  #other_metrics = data.table(do.call(rbind, lapply(cluster_results, function(cl_res){
  #  suppressWarnings(e1071::fclustIndex(cl_res, exprs(eset), index=c("xie.beni"))) #"xie.beni", "fukuyama.sugeno", "partition.coefficient", "partition.entropy"
  #})))
  #other_metrics[, avg_min_centroid_dist:=avg_min_centroid_dist]
  #other_metrics[, k:=cluster_range]
  ##other_metrics[, penalized_avg_min_centroid_dist:=avg_min_centroid_dist]
  
  # Plotting
  #------------------------------------
  # one plot for all clusters
  # clusters in different colours 
  # properties in different line types
  #------------------------------------
  for (k in 1:length(cluster_range)) {
    cluster_results_k = cluster_results[k] 
    
    # set colours for clusters
    plot_colors = rainbow(length(unique(cluster_results_k)))
    
    # Combine data for ggplot
    median_table_zscored_melted = melt(eset@assayData[["exprs"]], id.vars = "row_names")
    
    members = acore(eset, cluster_results[[k]], min.acore=0)
    acore_dt = rbindlist(lapply(seq_along(members), function(i){ data.table(CLUSTER=i, members[[i]])}))
    colnames(acore_dt) = c("cluster", "Var1", "cluster_membership_props") 
    
    #data_plot = merge(median_table_zscored_melted, acore_dt, by="Var1")
    data_plot = merge(merge(median_table_zscored_melted, acore_dt, by = "Var1", all = TRUE), feature_distribution_reshaped, by = "Var1", all = TRUE)
    
    data_plot$cluster = factor(data_plot$cluster)
    
    p = ggplot(data_plot, aes(x = Var2, y = value, group = Var1)) + 
      geom_line(aes(color=cluster, linetype = method)) +
      #geom_text_repel(aes(label = Var1), data = data_plot[data_plot$Var2 == sort(unique(data_plot$Var2))[1],], color = "black", size = 3) +
      scale_linetype_manual(values=c("diff_exp" = "solid", "both" = "dashed", "corr" = "dotted")) + 
      scale_x_continuous(breaks = sort(unique(data_plot$Var2))) +#, labels = 1:10) 
      xlab("") + 
      ylab(sprintf("Standardised expression\nvalues (%s)", data_input$norm)) +
      theme_classic() +
      theme(legend.position="bottom", 
            legend.margin = margin(t = -15, r = 0, b = 0, l = 0),
            text = element_text(family = plots_props$font_family, size = plots_props$font_size),
            axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
            axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
            plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header),
            legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
            legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size)
      ) +
      guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
    p
    
    # save bar plot
    ggsave(sprintf("%s/line_plot_k=%s.png", output_folder_path_fig, cluster_range[k]), p, dpi=plots_props$dpi, width=plots_props$image_width, height=plots_props$image_height, units = plots_props$image_units)
    ggsave(sprintf("%s/line_plot_k=%s.svg", output_folder_path_fig, cluster_range[k]), p, width=plots_props$image_width, height=plots_props$image_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
  }
  
  ### 
  for(i in 1:length(cluster_range)){
    members = acore(eset, cluster_results[[i]], min.acore=0)
    acore_dt = rbindlist(lapply(seq_along(members), function(i){ data.table(CLUSTER=i, members[[i]])}))
    k = cluster_range[i]
    fwrite(acore_dt, file=file.path(output_folder_path_tab, sprintf("cluster_members_k=%s.csv", k)))
  }
}


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
cluster_props = snakemake@params$cluster_props
plot_category = snakemake@params$plot_category
diff_exp_props = snakemake@params$diff_exp_props
corr_props = snakemake@params$corr_props
colours = snakemake@params$colors
plot_props = snakemake@params$plots_props
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_mfuzz_ts_analysis_for_sig_tissue_specific_miRNAs_line_plot.rds")

expr = fread(sprintf("%s/data_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
diff_exp = fread(sprintf("%s/%s_%s/results_%s/matrices/diff_exp/diff_exp%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, diff_exp_props$diff_log), sep='\t') 

expr_rownames = expr[[data_input$rna_class]]
expr[[data_input$rna_class]] = c()

folder_path_tab = sprintf("%s/%s_%s/results_%s/matrices/sig_tissue_specific_comp", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)

output_folder_path_fig = sprintf("%s/%s_%s/results_%s/figures/sig_tissue_specific_comp/line_plots/clustering", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_path_fig, recursive=TRUE)
output_folder_path_tab = sprintf("%s/line_plot/clustering", folder_path_tab)
dir.create(output_folder_path_tab, recursive=TRUE)


#------------------------------------ Script ----------------------------------- 
direction = list(c("both", "both"), c("up", "pos"), c("down", "neg"))

for (tissue in unique(annot[[plot_category[1]]])){
  for (num_dereg_combinations in diff_exp_props$num_dereg_combinations_th) {
    for (num_dereg_tissues in diff_exp_props$num_dereg_tissues_th) {
      for (m in corr_props$method) {
        for (th in corr_props$corr_th) {
          for (num_corr_tissues in corr_props$num_corr_tissues_th) {
            for (dir in direction) {
              file_name = sprintf("sig_%s_specific_comparison_diff_exp_%s_%s_at_least_in_%s_%s_comps_for_%s_%s_corr_method=%s_%s_for_%s_%s_%s_corr_th=%s", tissue, dir[1], data_input$rna_class, num_dereg_combinations, plot_category[2], num_dereg_tissues, plot_category[1], m, dir[2], num_corr_tissues, plot_category[1], data_input$rna_class, th)
              matching_files = list.files(sprintf("%s/venn_plots", folder_path_tab), pattern = sprintf("%s.csv", file_name))
              
              print(file_name)
              
              if (length(matching_files) > 0) {
                feature_distribution = fread(sprintf("%s/venn_plots/%s.csv", folder_path_tab, file_name), sep='\t')
                feature_list = unique(c(feature_distribution$diff_exp, feature_distribution$corr)[nzchar(c(feature_distribution$diff_exp, feature_distribution$corr))])
                
                # filter expr with feature_list
                diff_exp_filtered = diff_exp[diff_exp[[data_input$rna_class]] %in% feature_list,]
                
                median_table = data.table(feature=diff_exp_filtered[[data_input$rna_class]])
                for (tp in paste(sort(unique(annot[[plot_category[2]]])))) {
                  median_table[[tp]] = 2^(as.numeric(diff_exp_filtered[[sprintf("median_%s__%s_%s=%s", plot_category[2], tp, plot_category[1], tissue)]]))
                }
                
                if (dim(median_table)[1] < 4) {
                  print(sprintf("Not enough features for %s: %s", plot_category[1], tissue))
                  next
                } else if ((unique(is.na(feature_distribution$diff_exp)) == TRUE) || (unique(is.na(feature_distribution$corr)) == TRUE)) {
                  if (unique(is.na(feature_distribution$diff_exp)) == TRUE) {
                    diff_exp_or_corr = "diff_exp"
                  } else {
                    diff_exp_or_corr = "corr"
                  }
                  print(sprintf("Not enough features for %s: %s and method: %s", plot_category[1], tissue, diff_exp_or_corr))
                  next
                }
                
                median_table_df = data.frame(median_table, check.names = FALSE)
                
                # rearrange features
                feature_distribution_reshaped = reashape_features(feature_distribution)
                
                #max_clusters = 12
                #cluster_range = seq(4,max_clusters,1)
                cluster_range = cluster_props$range_num_of_clusters
                threads = 20
                
                final_folder_name = sprintf("sig_%s_specific_comparison_diff_exp_%s_%s_at_least_in_%s_%s_comps_for_%s_%s_corr_method=%s_%s_for_%s_%s_%s_corr_th=%s", tissue, dir[1], data_input$rna_class, num_dereg_combinations, plot_category[2], num_dereg_tissues, plot_category[1], m, dir[2], num_corr_tissues, plot_category[1], data_input$rna_class, th)
                output_folder_path_fig_tissue = sprintf("%s/%s", output_folder_path_fig, final_folder_name)
                dir.create(output_folder_path_fig_tissue, recursive=TRUE)
                output_folder_path_tab_tissue = sprintf("%s/%s", output_folder_path_tab, final_folder_name)
                dir.create(output_folder_path_tab_tissue, recursive=TRUE)
                compute_clusters_and_plots(median_table_df, feature_distribution_reshaped, data_input, output_folder_path_fig_tissue, output_folder_path_tab_tissue, cluster_range, threads, plot_props)
                
              } else {
                print(sprintf("No features for %s: %s", plot_category[1], tissue))
                next
              }
            }
          }
        }
      }
    }
  }  
}
  

#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])

