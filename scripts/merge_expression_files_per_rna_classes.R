suppressPackageStartupMessages(library(data.table))

#snakemake = readRDS("snakemake_merge_expression_files_per_rna_classes.rds")

print("merge_expression_files_per_rna_classes")


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input

# save rdata
#saveRDS(snakemake, file = "snakemake_merge_expression_files_per_rna_classes.rds")

file_name = sprintf("annotation_%s_%s", data_input$data_sub_set, data_input$feature_filtering)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 


#------------------------------------ Script ----------------------------------- 
tbl = c()
for (rna_cl in data_input$rna_classes) {
    tbl[[rna_cl]] = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, rna_cl, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
}

# check if all tables contain the same samples
# remove the first column from all matrices
tbl_without_feature_names = lapply(tbl, function(mat) mat[, -1, drop = FALSE])
# extract all colnames
colnames_list = lapply(tbl_without_feature_names, colnames)
# compare the colnames
same_colnames = all(sapply(colnames_list, function(x) setequal(x, colnames_list[[1]])))

if (length(unique(same_colnames)) != 1) {
  print("There are some files containing different samples")
  break
} else {
  if (unique(same_colnames) == FALSE) {
    print("There are some files containing different samples")
    break
  }
}

# replace "miRNA" to "RNA" in each matrix
tbl = lapply(tbl, function(mat) {
  colnames(mat) <- ifelse(colnames(mat) == "miRNA", "RNA", colnames(mat))
  return(mat)
})

# extract the reference column order (first matrix)
ref_colnames = colnames(tbl[[1]])

# reorder all matrices to match the reference column order
tbl = lapply(tbl, function(mat) mat[, ..ref_colnames])

# combine matrices by row
merged_tbl = do.call(rbind, tbl)


# save filtered table
fwrite(merged_tbl, sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, "ncRNA", data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t')


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
