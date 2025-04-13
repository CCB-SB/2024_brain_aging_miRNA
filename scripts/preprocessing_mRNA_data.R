suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(stringr))

data_input = snakemake@params$data_input

# meta data from https://github.com/OliInTheValley/SpatioTemporal_Analysis/blob/main/input_data/BulkSeq_Aging/BulkSeq_Aging_metadata.zip please unzip
annot_mrna = read.table(sprintf("%s/BulkSeq_Aging_metadata.txt", data_input$annotation_data_folder), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
annot_mrna$sampleID =  gsub("BulkSeq", "CA1", annot_mrna$sampleID)
tmp = str_split_fixed(annot_mrna$sampleID, "_", 3)
annot_mrna$tissue = sprintf("%s", tmp[,2])

# GEO data
expr_mrna_geo = read.table(sprintf("%s_%s/mR_expression_raw_aging.txt", data_input$raw_data_folder, data_input$data_sub_set), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
names_expressed_mrna_geo = expr_mrna_geo$X
expr_mrna_geo$X = c()
colnames(expr_mrna_geo) = gsub("BulkSeq", "CA1", colnames(expr_mrna_geo))
rownames(expr_mrna_geo) = names_expressed_mrna_geo

# Normalize with DeSeq2
dds_geo = DESeqDataSetFromMatrix(countData=expr_mrna_geo, colData=annot_mrna, design = ~age + tissue) 
design(dds_geo) = ~age + tissue
dds_geo = DESeq(dds_geo, fitType = 'local')
normalized_counts = counts(dds_geo, normalized = TRUE)

# load mrna expr
# this creates a variable with name "counttable_CA1_list" in the workspace
#load("../../data_external/data_CA1_CA2_mRNA/dds_CA1Publication.bin")
# this creates a variable with name "dds_CA1_list" in the workspace
#load("../../data_external/data_CA1_CA2_mRNA/dds_CA1_full_Deseqonly.bin")
# get names of "expressed" mRNA from counttable_CA1_list$All$expressed
#names_expressed_mrna = rownames(counttable_CA1_list$All$expressed)
# get dds (deseq data stuff) object from list
#dds = dds_CA1_list$All
# we need to do this Umweg because we can only execute the normalize on the Deseq object
#expr_mrna = counts(dds, normalized=TRUE)
# filter only the interesting ones from above
#expr_mrna = expr_mrna[names_expressed_mrna, ]

#expr_mrna["Xkr4", "CA1_th_59"]
#normalized_counts["Xkr4", "CA1_th_59"]

#expr_mrna_sorted = expr_mrna[,colnames(normalized_counts)]
#tmp = rownames(normalized_counts) %in% rownames(expr_mrna_sorted)
#normalized_counts_sorted = normalized_counts[tmp,]
#norm(expr_mrna_sorted - normalized_counts_sorted, type = "I")
#norm(expr_mrna_sorted - normalized_counts_sorted, type = "F")

# save 
file_path = sprintf("%s_%s/mRNA_normalized_counts.csv", data_input$raw_data_folder, data_input$data_sub_set)
write.table(normalized_counts, file_path, sep='\t', row.names = T)


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])

