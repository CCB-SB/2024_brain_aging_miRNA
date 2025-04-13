library(data.table)
library(stringr)

#snakemake = readRDS("snakemake_filter_samples_number_alignment.rds")


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
data_input$feature_filtering = "filtered"
norm_list = snakemake@params$norm_list
number_of_reads = snakemake@params$number_of_reads
id_col = snakemake@params$id_col
fastq_filename_col = snakemake@params$fastq_filename_col

# read mapping
mapping_info = fread(snakemake@input$mapping_info, sep = "\t")

# save rdata
#saveRDS(snakemake, file = "snakemake_filter_samples_number_alignment.rds")

if (data_input$data_sub_set != "TMS_brain") {
  data_set = data_input$data_sub_set #strsplit(data_input$data_sub_set, "_")[[1]][2]
} else {
  data_set = data_input$data_sub_set
}

# load annot
file_name = sprintf("annotation_%s", data_set)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set,data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
  
save_path = sprintf("%s_%s", data_input$preprocesing_data_folder, data_set)
dir.create(save_path, recursive=TRUE)
save_path_quant = sprintf("%s_%s", data_input$preprocesing_data_folder, data_set)
dir.create(save_path_quant, recursive=TRUE)

# create output path
#tmp = strsplit(data_input$preprocesing_data_folder, "/")[[1]]
#output_folder_prefix = sprintf("%s/%s/%s", tmp[1], tmp[2], tmp[3])
#output_folder_suffix = ""
#if (length(tmp) >= 5) {
#  output_folder_suffix = tmp[5]
#  if (length(tmp) >= 6) {
#    for (i in 6:length(tmp)) {
#      print(tmp[i])
#      output_folder_suffix = sprintf("%s/%s", output_folder_suffix, tmp[i])
#    }
#  }
#}

#output_folder_idefix = sprintf("samples_filtering_alignment_against_mmu10_over_%s_reads", number_of_reads)

#save_path = sprintf("%s/%s/%s_%s", output_folder_prefix, output_folder_idefix, output_folder_suffix, data_set)
#dir.create(save_path, recursive=TRUE)
#save_path_quant = sprintf("%s/%s/%s", output_folder_prefix, output_folder_idefix, output_folder_suffix)
#dir.create(save_path_quant, recursive=TRUE)


#------------------------------------ Script ----------------------------------- 

if (data_set == "TMS_brain") {
  # extract tissue
  mapping_info$tissue = str_split_fixed(mapping_info$Sample, "_", 2)[,1]
  # select CA1 samples
  mapping_info_filtered = mapping_info[mapping_info$tissue == "Brain",]
} else {#if ((data_set == "microglia")) {
  tmp = c(fastq_filename_col, id_col)
  mapping_info = merge(mapping_info, annot[,..tmp], by=fastq_filename_col, suffixes = c("","_drop"))
  # rename columns such that Sample is qual to id_col
  mapping_info[[fastq_filename_col]] = mapping_info$Sample
  mapping_info$Sample = mapping_info[[id_col]]

  mapping_info_filtered = mapping_info
} 
#else {
## extract cohort
#tmp = str_split_fixed(mapping_info$fastq_names, "_", 5)
#mapping_info$cohort = tmp[,2]
#mapping_info$Sample = sprintf("%s_%s_%s", tmp[,2], tmp[,3], tmp[,4])
## find first part (cohort) from data_set
#tmp = strsplit(data_set, "_")[[1]][1]
## select CA1 samples
#mapping_info_filtered = mapping_info[mapping_info$cohort == tmp,]
#}

# order for column reads_aligned
ordering_decreasing = order(as.numeric(mapping_info_filtered$reads_aligned))
mapping_info_filtered_sort = mapping_info_filtered[ordering_decreasing,]


# take only samples with more than 2e6 aligned reads
mapping_info_filtered_sort$alignment_over_th = rep(FALSE, dim(mapping_info_filtered_sort)[1])
mapping_info_filtered_sort[as.numeric(mapping_info_filtered_sort$reads_aligned) >= as.numeric(number_of_reads), alignment_over_th:=TRUE]

#print(mapping_info_filtered_sort)

print(sprintf("samples to be removed: %s", mapping_info_filtered_sort[mapping_info_filtered_sort$alignment_over_th == FALSE,]$Sample))
keep = mapping_info_filtered_sort[mapping_info_filtered_sort$alignment_over_th == TRUE,]$Sample

# select the samples to be kept
if (data_set == "TMS_brain") {
  tmp = annot$ID %in% keep
} else {
  tmp = annot[[id_col]]%in% keep
}

# annot
annot_filtered = annot[tmp,]
#fwrite(annot_filtered, sprintf("%s/%s_%s.csv", save_path, file_name, data_input$feature_filtering), sep = "\t")

# expression matrices
input_folder = sprintf("%s_%s", data_input$preprocesing_data_folder, data_set)
files = list.files(input_folder)

for (rna_class in data_input$rna_class_list) {
  file_prefix = sprintf("%s", rna_class)

  file_prefix_list = c()
  for (file in files) {
    if (startsWith(file, file_prefix)) {
      file_prefix_list = c(file_prefix_list, file)
    }
  }
  
  for (file in file_prefix_list) {
    # load expr
    tbl = fread(sprintf("%s/%s", input_folder, file), sep='\t', header=T)
    # filter expr
    samples_keep = c(colnames(tbl)[1], keep)
    tmp = (colnames(tbl) %in% samples_keep)
    tbl_filtered = tbl[,..tmp]
    # save filtered expr
    file_filtered = gsub("_complete", "", file)
    parts = strsplit(file_filtered, "_")[[1]]
    file_filtered = sprintf("%s_%s_%s", parts[1], "filtered", paste(parts[-1], collapse = "_"))
    file_filtered_csv = gsub("tsv", "csv", file_filtered)
    fwrite(tbl_filtered, sprintf("%s/%s", input_folder, file_filtered_csv), sep = "\t")
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])

