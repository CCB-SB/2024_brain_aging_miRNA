suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(matrixStats))


#snakemake = readRDS("snakemake_target_genes_scatter_plot_CA2.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("target_genes_scatter_plot_CA2")


# ---------------------- Load data ----------------------
data_input = snakemake@params$data_input
feature = snakemake@params$feature
gene_list_folder_path = snakemake@params$gene_list_folder_path
split_prop = snakemake@params$split_prop
hline_value = snakemake@params$hline_th
xticks_names = snakemake@params$xticks_names
colours = snakemake@params$colors
plots_props = snakemake@params$plots_props
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_target_genes_scatter_plot_CA2.rds")

# load annot
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
 
# load mrna expr
# load CA2 data from GEo downloaded file
expr_mrna = fread("/local/22-02_brain_aging_mrna_data_aengel/GSE227689_BulkSeqRejuvenation_Counttable.txt", sep='\t', colClasses=c(ID="character")) 

# set rownames
rownames(expr_mrna) = expr_mrna$V1
expr_mrna$V1 = c()

# filter only the interesting ones from above
names_mrna = rownames(expr_mrna)
expr_mrna = as.data.table(expr_mrna)

# rename the colnames to the naming in our data
colnames(expr_mrna) = gsub("BulkSeq", "CA2", colnames(expr_mrna))

# load mirna data
tbl = fread(sprintf("%s/data_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
expr_mirna = tbl[, sapply(tbl, is.numeric), with=F]
names_mirna = tbl[[data_input$rna_class]]

# match columns
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

# load gene list
gene_list = fread(gene_list_folder_path, sep='\t', header=T)
print(sprintf("Length of target gene list: %s", length(gene_list$Gene)))

output_folder = sprintf("%s/%s_%s/results_%s/figures/target_genes_CAS/scatter_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder, recursive=TRUE)


#------------------------------------ Script ----------------------------------- 
# filter gene expression by genes from the genes list
tmp = (names_mrna %in% gene_list$Gene)
expr_mrna_genes = expr_mrna[tmp,]
names_mrna_genes = names_mrna[tmp]

# if there are some genes not included in the expression file
if (length(setdiff(gene_list$Gene, names_mrna)) > 0) {
  print(sprintf("There are %s (%s) genes missing in the expression which a contained in the inserted target gene file", length(setdiff(gene_list$Gene, names_mrna)), paste(setdiff(gene_list$Gene, names_mrna), collapse = ", ")))
}
print(sprintf("Shape of gene expression matrix: %s", dim(expr_mrna_genes)))

# filter miRNA expression by given feature from the config
expr_mirna_feature = expr_mirna[names_mirna == feature,]

expr_mrna_genes_median = c()
expr_mirna_feature_median = c()
for (ti in unique(annot_intersection[[split_prop]])) {
  group_ID = annot_intersection[annot_intersection[[split_prop]] == ti,][[data_input$identifier_column]]
  expr_mrna_genes_tissue = expr_mrna_genes[,..group_ID]
  expr_mirna_feature_tissue = expr_mirna_feature[,..group_ID]
  
  expr_mrna_genes_median[[ti]] = rowMedians(as.matrix(expr_mrna_genes_tissue), useNames = FALSE)
  expr_mirna_feature_median[[ti]] = rowMedians(as.matrix(expr_mirna_feature_tissue), useNames = FALSE)
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

plot_df = merge(expr_mrna_genes_median_merged_melt, expr_mirna_feature_median_melt, by = "variable")

plot_df$show_names = plot_df$names
plot_df$show_names = ifelse(plot_df$value_gene < hline_value, "", plot_df$show_names)

plot_df$shape_type = ifelse(plot_df$names %in% gene_list$Gene[gene_list$Origin == "target_gene"], "diamond", "circle")

colours_plot = unlist(colours[[split_prop]])

x_lab_title_thin = sprintf("MiRNA expression\nvalues per brain region\n(%s, detection rate = %s%%\nper %s)", gsub("_norm", "" ,data_input$norm), gsub("p", "", data_input$detection_rate), xticks_names$categories[[split_prop]])
x_lab_title_wide = sprintf("MiRNA expression values per brain region\n(%s, detection rate = %s%% per %s)", gsub("_norm", "" ,data_input$norm), gsub("p", "", data_input$detection_rate), xticks_names$categories[[split_prop]])
y_lab_title_tall = "Gene counts values (Hahn et al.) for CAS"
y_lab_title_small = "Gene counts values\n(Hahn et al.) for CAS"
plot_title = sprintf("For %s", feature)
  
font_correction = 6/17.07000017

# scatter plot
scatter_plot = ggplot(plot_df, aes(x = value_mirna, y = value_gene, color = variable, shape=shape_type)) + 
  geom_point(size = 0.5) +
  scale_color_manual(values=colours_plot) +
  scale_shape_manual(values=c("circle" = 19, "diamond" = 18)) +
  geom_text(aes(label = show_names), vjust = -1, hjust = 0.5, color = "black", size=6*font_correction) +
  geom_hline(yintercept=hline_value, linetype="dashed", color = "#727272", size=0.3) +
  ggtitle(plot_title) +
  xlab(x_lab_title_thin) +
  ylab(y_lab_title_tall) +
  theme_classic() +
  theme(legend.position="none", 
        text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
        axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header),
        legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
        axis.text.y = element_text(angle=90, hjust=0.5, vjust=0.5),
        #plot.margin = unit(c(-0.25,0.25,-0.25,-1.25), "cm")
  )
#scatter_plot

file_name = sprintf("%s", feature)
height = 10
width = 4
ggsave(sprintf("%s/%s_%sx%s.png", output_folder, file_name, width, height), scatter_plot, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
ggsave(sprintf("%s/%s_%sx%s.svg", output_folder, file_name, width, height), scatter_plot, width=width, height=height, units = plots_props$image_units, limitsize = FALSE, scale = 1)

scatter_plot = scatter_plot +
  xlab(x_lab_title_wide) + 
  ylab(y_lab_title_small) 
height = plots_props$image_height
width = 6
ggsave(sprintf("%s/%s_%sx%s.png", output_folder, file_name, width, height), scatter_plot, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
ggsave(sprintf("%s/%s_%sx%s.svg", output_folder, file_name, width, height), scatter_plot, width=width, height=height, units = plots_props$image_units, limitsize = FALSE, scale = 1)


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
