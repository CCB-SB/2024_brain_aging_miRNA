suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))

#snakemake = readRDS("snakemake_correlation_mirna_target_genes_cas_scatter_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("correlation_mirna_target_genes_cas_scatter_plot")


# ---------------------- Load data ----------------------
data_input = snakemake@params$data_input
feature = snakemake@params$feature
split_prop = snakemake@params$split_prop
corr_params = snakemake@params$corr_params
manual_plots_parameters = snakemake@params$manual_plots_parameters
xticks_names = snakemake@params$xticks_names
colours = snakemake@params$colors
plots_props = snakemake@params$plots_props
results_folder = snakemake@params$results_folder
save_folder = snakemake@params$save_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_correlation_mirna_target_genes_cas_scatter_plot.rds")
#stop()

# load annot
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character"))
# load mrna expr
tbl = read.csv(sprintf("%s_%s/mRNA_normalized_counts.csv", data_input$raw_data_folder, data_input$data_sub_set), sep='\t', header=TRUE, check.names = FALSE)
names_expressed_mrna = rownames(tbl)
rownames(tbl) = c()
names_mrna = toupper(names_expressed_mrna)
expr_mrna = as.data.table(tbl)
print(sprintf("Shape of gene expression matrix: %s", paste(dim(expr_mrna), collapse = "x")))

# load mirna data
tbl = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
expr_mirna = tbl[, sapply(tbl, is.numeric), with=F]
names_mirna = tbl[[data_input$rna_class]]

output_folder = sprintf("%s/%s_%s/results_%s/figures/%s_%s/scatter_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, feature, save_folder)
dir.create(output_folder, recursive=TRUE)

brain_region_th = 4


#------------------------------------ Script ----------------------------------- 
# filter samples
# get intersection of column names
common_samples = intersect(colnames(expr_mrna), colnames(expr_mirna))
if (length(setdiff(colnames(expr_mirna), colnames(expr_mrna))) > 0) {
  print(sprintf("There are %s (%s) samples missing in the expression of the genes compared to the miRNAs", length(setdiff(colnames(expr_mirna), colnames(expr_mrna))), paste(setdiff(colnames(expr_mirna), colnames(expr_mrna)), collapse = ", ")))
}

# select columns from intersection
expr_mrna = expr_mrna[, ..common_samples]
expr_mirna = expr_mirna[, ..common_samples]
print(sprintf("Shape of gene expression matrix: %s", dim(expr_mrna)))
print(sprintf("Shape of miRNA expression matrix: %s", dim(expr_mirna)))

# filter annot also to only keep the samples which are in both expression matrices
print(sprintf("Number of samples in annot before filtering: %s", dim(annot)[1]))
tmp = (annot[[data_input$identifier_column]] %in% common_samples)
annot_intersection = annot[tmp,]
print(sprintf("Shape of annot after filtering: %s", dim(annot_intersection)[1]))

for (target_genes_source in c("target_genes_all", "target_genes_without_weak")) {
  for (method in corr_params$method) {
    for (corr_th in corr_params$corr_th) {
      if (save_folder == "target_genes_corr") {
        dir.create(sprintf("%s/%s/%s", output_folder, target_genes_source, method), recursive=TRUE)
        
        input_file_path = sprintf("%s/%s_%s/results_%s/matrices/%s_%s/upset_plots/%s/%s/upset_%s_corr_m=%s_corr_th=%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, feature, save_folder, target_genes_source, method, feature, method, corr_th)
        brain_region_gene_distribution = fread(input_file_path, sep='\t')
        brain_region_gene_distribution_rownames = brain_region_gene_distribution$V1
        brain_region_gene_distribution$V1 = c()
        
        tmp = (colSums(brain_region_gene_distribution) >= brain_region_th)
        gene_list = colnames(brain_region_gene_distribution[,..tmp])
        
        if (!length(gene_list)) {
          next
        }
        
        # if there are some genes not included in the expression file
        if (length(setdiff(gene_list, names_mrna)) > 0) {
          print(sprintf("There are %s (%s) genes missing in the expression which a contained in the inserted target gene file", length(setdiff(gene_list, names_mrna)), paste(setdiff(gene_list, names_mrna), collapse = ", ")))
        }
        #print(sprintf("Shape of gene expression matrix: %s", dim(expr_mrna_genes)))
        
        gene_list_filtered = intersect(gene_list, names_mrna)
      } else {
        dir.create(sprintf("%s/%s/%s", output_folder, target_genes_source, method), recursive=TRUE)
        
        input_file_path = sprintf("%s/%s_%s/results_%s/matrices/%s_%s/target_genes_martin/%s/sig_neg_correlation_corr_the=0.3.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, feature, save_folder, method)
        brain_region_gene_distribution = fread(input_file_path, sep='\t')
        
        gene_list_filtered = intersect(manual_plots_parameters$genes, names_mrna)
      }
      
      plot_df = c()
      for (gene in gene_list_filtered) {
        print(gene)
        tissue = split_prop[1]
        time = split_prop[2]
        
        # filter features
        # filter gene expression by genes from the genes list
        tmp = (names_mrna == gene)
        expr_mrna_genes = expr_mrna[tmp,]
        names_mrna_genes = names_mrna[tmp]
        
        # filter miRNA expression by given feature from the config
        expr_mirna_feature = expr_mirna[names_mirna == feature,]
        
        if (save_folder == "target_genes_corr") {
          # for the plot we only want to show the tissues for which the gene and mirna is sig. corr.
          tissue_list = brain_region_gene_distribution_rownames[as.logical(brain_region_gene_distribution[[gene]])]
        } else {
          # use all tissues
          #tissue_list = unique(annot_intersection[[tissue]])
          tissue_list = brain_region_gene_distribution[brain_region_gene_distribution$mRNA == gene,][[tissue]]
        }
        group_id = annot_intersection[annot_intersection[[tissue]] %in% tissue_list][[data_input$identifier_column]]
        expr_mrna_genes_filtered = expr_mrna_genes[, ..group_id]
        expr_mirna_feature_filtered = expr_mirna_feature[, ..group_id]
        
        expr_mrna_genes_median = c()
        expr_mirna_feature_median = c()
        for (ti in tissue_list) {
          for (t in unique(annot_intersection[[time]])) {
            group_ID = annot_intersection[((annot_intersection[[tissue]] == ti) & (annot_intersection[[time]] == t)),][[data_input$identifier_column]]
            expr_mrna_genes_tissue_time = expr_mrna_genes_filtered[,..group_ID]
            expr_mirna_feature_tissue_time = expr_mirna_feature_filtered[,..group_ID]
              
            expr_mrna_genes_median[[sprintf("%s_%s", ti, t)]] = rowMedians(as.matrix(expr_mrna_genes_tissue_time), useNames = FALSE)
            expr_mirna_feature_median[[sprintf("%s_%s", ti, t)]] = rowMedians(as.matrix(expr_mirna_feature_tissue_time), useNames = FALSE)
          }
        }
        
        # concat the list elements to a matrix
        expr_mrna_genes_median_merged = data.frame(expr_mrna_genes_median)
        expr_mrna_genes_median_merged$names = names_mrna_genes
        expr_mrna_genes_median_merged_melt = melt(expr_mrna_genes_median_merged)
        #expr_mrna_genes_median_merged_melt$type = 0
        colnames(expr_mrna_genes_median_merged_melt) = gsub("value", "value_gene", colnames(expr_mrna_genes_median_merged_melt))
        
        expr_mirna_feature_median_merged = data.frame(expr_mirna_feature_median)
        #expr_mirna_feature_median_merged$names = feature
        expr_mirna_feature_median_melt = melt(expr_mirna_feature_median_merged)
        #expr_mirna_feature_median_melt$type = 1
        colnames(expr_mirna_feature_median_melt) = gsub("value", "value_mirna", colnames(expr_mirna_feature_median_melt))
        
        plot_df[[gene]] = merge(expr_mrna_genes_median_merged_melt, expr_mirna_feature_median_melt, by = "variable")
      }
      
      maxi_x = c()
      mini_x = c()
      maxi_y = c()
      mini_y = c()
      for (gene in manual_plots_parameters$genes) {
        maxi_x[[gene]] = max(plot_df[[gene]]$value_mirna)
        mini_x[[gene]] = min(plot_df[[gene]]$value_mirna)
        
        maxi_y[[gene]] = max(plot_df[[gene]]$value_gene)
        mini_y[[gene]] = min(plot_df[[gene]]$value_gene)
      }
      mirna_maxi = max(unlist(maxi_x))
      mirna_mini = min(unlist(mini_x))
      
      gene_maxi = max(unlist(maxi_y))
      gene_mini = min(unlist(mini_y))
      
      p_scatter_manual = list()
      for (gene in gene_list_filtered) {
        #plot_df$show_names = plot_df$names
        #plot_df$show_names = ifelse(plot_df$value_gene < hline_value, "", plot_df$show_names)
        
        #plot_df$shape_type = ifelse(plot_df$names %in% gene_list$Gene[gene_list$Origin == "target_gene"], "diamond", "circle")
              
        x_lab_title_small = sprintf("Expression\nvalues of %s\n(median per brain region and age)\n(%s, detection rate = %s%%\nper %s)", feature, gsub("_norm", "" ,data_input$norm), gsub("p", "", data_input$detection_rate), xticks_names$categories[[tissue]])
        x_lab_title_big = sprintf("Expression values of %s\n(median per brain region and age)\n(%s, detection rate = %s%% per %s)", feature, gsub("_norm", "" ,data_input$norm), gsub("p", "", data_input$detection_rate), xticks_names$categories[[tissue]])
        y_lab_title_big = "Gene expression values (Hahn et al.) for sig. neg. corr. target genes"
        y_lab_title_small = "Gene expression values (Hahn et al.)\nfor sig. neg. corr. target genes"
        #plot_title_big = sprintf("For %s", feature)
        plot_title = gene
        #plot_title_small = sprintf("For %s\nand target gene %s", feature, gene)
        font_correction = 6/17.07000017
        
        tmp = str_split_fixed(plot_df[[gene]]$variable, "_", 2)
        plot_df[[gene]]$var1 = tmp[,1]
        plot_df[[gene]]$var2 = tmp[,2]
        
        # coloured by tissue
        colours_plot = unlist(colours[[tissue]])
        # scatter plot
        scatter_plot_tissue = ggplot(plot_df[[gene]], aes(x = value_mirna, y = value_gene, color = var1)) + 
          geom_point(size = 1) +
          scale_color_manual(values=colours_plot) +
          scale_shape_manual(values=c("circle" = 19, "diamond" = 18)) +
          #geom_text(aes(label = show_names), vjust = -1, hjust = 0.5, color = "black", size=6*font_correction) +
          #geom_hline(yintercept=hline_value, linetype="dashed", color = "#727272", size=0.3) +
          ggtitle(plot_title) +
          xlab(x_lab_title_small) +
          ylab(y_lab_title_big) +
          theme_classic() +
          theme(legend.position="none", 
                text = element_text(family = plots_props$font_family, size = plots_props$font_size),
                axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
                axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
                plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header),
                legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
                legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
                #axis.text.y = element_text(angle=90, hjust=0.5, vjust=0.5),
                #plot.margin = unit(c(-0.25,0.25,-0.25,-1.25), "cm"),
                plot.background = element_rect(fill = "transparent", color = NA),   # Remove plot background
          )
        #scatter_plot_tissue

        file_name = sprintf("scatter_plot_%s_vs_%s_corr_m=%s_corr_th=%s_coloured_by_%s", gene, feature, method, corr_th, tissue)
        #height = 10
        #width = 4
        #ggsave(sprintf("%s/%s/%s/%s_%sx%s.png", output_folder, target_genes_source, method, file_name, width, height), scatter_plot_tissue, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
        #ggsave(sprintf("%s/%s/%s/%s_%sx%s.svg", output_folder, target_genes_source, method, file_name, width, height), scatter_plot_tissue, width=width, height=height, units = plots_props$image_units, limitsize = FALSE, scale = 1)
        
        scatter_plot_tissue = scatter_plot_tissue +
          xlab(x_lab_title_big) + 
          ylab(y_lab_title_small) 
        height = plots_props$image_height
        width = 6
        ggsave(sprintf("%s/%s/%s/%s_%sx%s.png", output_folder, target_genes_source, method, file_name, width, height), scatter_plot_tissue, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
        ggsave(sprintf("%s/%s/%s/%s_%sx%s.svg", output_folder, target_genes_source, method, file_name, width, height), scatter_plot_tissue, width=width, height=height, units = plots_props$image_units, limitsize = FALSE, scale = 1)
        
        if ((gene %in% manual_plots_parameters$genes) & (target_genes_source == manual_plots_parameters$target_genes_source) & (method == manual_plots_parameters$method) & (corr_th == manual_plots_parameters$corr_th)) {
          p_scatter_manual[[gene]] = scatter_plot_tissue + xlim(mirna_mini, mirna_maxi) + xlab("") + ylab("") # ylim(gene_mini, gene_maxi)
        }
  
        # coloured by time
        #colours_plot = unlist(colours[[time]])
        #scatter_plot_time = ggplot(plot_df, aes(x = value_gene, y = value_mirna, color = var2)) + 
        #  geom_point(size = 1) +
        #  scale_color_manual(values=colours_plot) +
        #  scale_shape_manual(values=c("circle" = 19, "diamond" = 18)) +
        #  #geom_text(aes(label = show_names), vjust = -1, hjust = 0.5, color = "black", size=6*font_correction) +
        #  #geom_hline(yintercept=hline_value, linetype="dashed", color = "#727272", size=0.3) +
        #  ggtitle(plot_title) +
        #  xlab(x_lab_title_small) +
        #  ylab(y_lab_title_big) +
        #  theme_classic() +
        #  theme(legend.position="none", 
        #        text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        #        axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
        #        axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        #        plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header),
        #        legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        #        legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
        #        axis.text.y = element_text(angle=90, hjust=0.5, vjust=0.5),
        #        #plot.margin = unit(c(-0.25,0.25,-0.25,-1.25), "cm")
        #  )
        ##scatter_plot_time
        #  
        #file_name = sprintf("scatter_plot_%s_vs_%s_corr_m=%s_corr_th=%s_coloured_by_%s", gene, feature, method, corr_th, time)
        ##height = 10
        ##width = 4
        ##ggsave(sprintf("%s/%s/%s/%s_%sx%s.png", output_folder, target_genes_source, method, file_name, width, height), scatter_plot_time, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
        ##ggsave(sprintf("%s/%s/%s/%s_%sx%s.svg", output_folder, target_genes_source, method, file_name, width, height), scatter_plot_time, width=width, height=height, units = plots_props$image_units, limitsize = FALSE, scale = 1)
        #
        #scatter_plot_time = scatter_plot_time +
        #  xlab(x_lab_title_big) + 
        #  ylab(y_lab_title_small) 
        #height = plots_props$image_height
        #width = 6
        #ggsave(sprintf("%s/%s/%s/%s_%sx%s.png", output_folder, target_genes_source, method, file_name, width, height), scatter_plot_time, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
        #ggsave(sprintf("%s/%s/%s/%s_%sx%s.svg", output_folder, target_genes_source, method, file_name, width, height), scatter_plot_time, width=width, height=height, units = plots_props$image_units, limitsize = FALSE, scale = 1)
      }
      
      if ((target_genes_source == manual_plots_parameters$target_genes_source) & (method == manual_plots_parameters$method) & (corr_th == manual_plots_parameters$corr_th)) {
        
        # sort the figures given the order from the config
        p_scatter_manual = p_scatter_manual[manual_plots_parameters$genes]
  
        plot_nrows = manual_plots_parameters$plot_layout[1]
        plot_ncols = manual_plots_parameters$plot_layout[2]
        
        y_axis_title = y_lab_title_small
        
        #col_distribution = c()
        #for (col_i in 1:plot_ncols) {
        #  col_tissues = c()
        #  for (row_i in 0:(plot_nrows-1)) {
        #    col_tissues = append(col_tissues, manual_plots_parameters$genes[col_i + row_i * plot_ncols])
        #  }
        #  col_distribution[[sprintf("col_%s", col_i)]] = col_tissues
        #}
        #  
        #p_scatter_manual_adj = list()
        #for (gene in manual_plots_parameters$genes) {
        #  if (gene %in% col_distribution$col_1) {
        #    p_scatter_manual_adj[[gene]] = p_scatter_manual[[gene]] + theme(plot.margin = unit(c(0.1,0.1,0,-0.4 ), "cm"))
        #  } else if (gene %in% col_distribution$col_2){ 
        #    p_scatter_manual_adj[[gene]] = p_scatter_manual[[gene]] + theme(plot.margin = unit(c(0.1,0.1,0,-0.1), "cm")) #theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,0,0,0), "cm"))
        #  } else if (gene %in% col_distribution$col_3){
        #    p_scatter_manual_adj[[gene]] = p_scatter_manual[[gene]] + theme(plot.margin = unit(c(0.1,0.15,0,-0.15), "cm")) #theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(0,0.15,0,-0.15), "cm"))
        #  } 
        #}
        
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
        
        p_scatter_manual_adj_mod_margins = p_scatter_manual
        for (img_i in 1:length(p_scatter_manual_adj_mod_margins)) {
          if ((!(img_i %in% left_col)) & (!(img_i %in% right_col))) {
            p_scatter_manual_adj_mod_margins[[img_i]] = p_scatter_manual_adj_mod_margins[[img_i]] + 
              theme(plot.margin = margin(5, -1, -7, -13, "pt")
              )
          }
          else if (img_i %in% left_col) {
            p_scatter_manual_adj_mod_margins[[img_i]] = p_scatter_manual_adj_mod_margins[[img_i]] + 
              theme(plot.margin = margin(5, -4, -7, -10, "pt")
              )
          } else if (img_i %in% right_col) {
            p_scatter_manual_adj_mod_margins[[img_i]] = p_scatter_manual_adj_mod_margins[[img_i]] + 
              theme(plot.margin = margin(5, 4, -7, -18, "pt")
              )
          }
          if (!(img_i %in% left_col)) {
            p_scatter_manual_adj_mod_margins[[img_i]] = p_scatter_manual_adj_mod_margins[[img_i]] + 
              #scale_y_continuous(labels = function(breaks) {rep_along(breaks, "  ")}, limits = c(0, 1.15*maxi_expr_rna), expand = c(0,0)) +
              theme(axis.ticks.y = element_line(color = "transparent"),
                    axis.text.y = element_text(color = "transparent")
              )
          }
        }
        
        # Function to extract y-limits of each plot
        get_ylims = function(plot) {
          ggplot_build(plot)$layout$panel_params[[1]]$y.range
        }
        
        # Split the plots into row-wise groups
        # Plot matrix must be rectangular, therefore fill with NA if
        # names(p_scatter_manual_adj_mod_margins) < plot_nrows * 2 (columns)
        plot_names = names(p_scatter_manual_adj_mod_margins)
        plot_names = c(plot_names, rep(NA, plot_nrows * 2 - length(plot_names)))
        plot_matrix = matrix(plot_names, nrow = plot_nrows, byrow = TRUE, ncol=2)
        
        # Compute row-wise y-max safely
        row_ymax = sapply(1:plot_nrows, function(i) {
          # Extract plots from the row
          tmp = names(p_scatter_manual_adj_mod_margins) %in% plot_matrix[i,]
          plots_in_row = p_scatter_manual_adj_mod_margins[tmp]
          # Remove NULLs
          plots_in_row = plots_in_row[!sapply(plots_in_row, is.null)]
          
          # Get the upper y-limit for each plot and remove any NA values
          ymax_values = unlist(lapply(plots_in_row, function(p) get_ylims(p)[2]))
          # Remove NA values
          ymax_values = ymax_values[!is.na(ymax_values)]
          
          if (length(ymax_values) > 0) {
            # Return the max if values exist
            return(max(ymax_values))
          } else {
            # Return NA if no valid values found
            return(NA)
          }
        })

        # Adjust ylim for each plot row-wise
        for (i in 1:plot_nrows) {
          for (j in 1:plot_ncols) {
            if (!is.null(plot_matrix[[i,j]])) {
              #print(sprintf("i:%s j:%s title:%s", i, j, p_scatter_manual_adj_mod_margins[[plot_matrix[i,j]]]$labels$title))
              p_scatter_manual_adj_mod_margins[[plot_matrix[i,j]]] = p_scatter_manual_adj_mod_margins[[plot_matrix[i,j]]] + ylim(gene_mini, row_ymax[i])
            }
          }
        }

        all_p_scatter_manual = arrangeGrob(grobs=p_scatter_manual_adj_mod_margins, ncol=plot_ncols, nrow=plot_nrows, bottom = textGrob(x_lab_title_big, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90))
  
        plot_width = manual_plots_parameters$plot_size[1]
        plot_height = manual_plots_parameters$plot_size[2]
        
        ggsave(sprintf("%s/%s/%s/%s_%s_%sx%s.png", output_folder, target_genes_source, method, file_name, paste(manual_plots_parameters$genes, collapse = "_"), plot_width, plot_height), all_p_scatter_manual, dpi=plots_props$dpi, width=plot_width, height=plot_height, units = plots_props$image_units)
        ggsave(sprintf("%s/%s/%s/%s_%s_%sx%s.svg", output_folder, target_genes_source, method, file_name, paste(manual_plots_parameters$genes, collapse = "_"), plot_width, plot_height), all_p_scatter_manual, width=plot_width, height=plot_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
      }
      
      rm(p_scatter_manual)
      rm(p_scatter_manual_adj)
      rm(p_scatter_manual_adj_mod_margins)
      
      gc()
    }
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
