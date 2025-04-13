suppressPackageStartupMessages(library(data.table))

# save rdata
#saveRDS(snakemake, file = "snakemake_create_binned_human_age.rds")
#stop()

#snakemake = readRDS("snakemake_create_binned_human_age.rds")

print("create_binned_human_age")


#---------------------------------- Functions ----------------------------------
split_into_three_intervals = function(values, number_of_bins) {
  min_val = min(values)
  max_val = max(values)
  
  # Calculate interval size
  interval_size = (max_val - min_val) / number_of_bins
  
  # Define interval boundaries
  threshold1 = min_val + interval_size
  threshold2 = min_val + 2 * interval_size
  
  return(list(
    min_val = floor(min_val),
    threshold1 = round(threshold1),
    threshold2 = round(threshold2),
    max_val = ceiling(max_val) + 1
  ))
}


# ---------------------- Load data ----------------------
data_input = snakemake@params$data_input
property_to_be_binned = snakemake@params$property_to_be_binned
number_of_bins = snakemake@params$number_of_bins
miRNA_human_path = snakemake@params$miRNA_human_path

input_file_path = sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering)
annot = fread(input_file_path, sep='\t') 


#------------------------------------ Script ----------------------------------- 
intervall_ends = unlist(unname(split_into_three_intervals(floor(annot[[property_to_be_binned]]), number_of_bins)))

age_bin_groups = floor(annot[[property_to_be_binned]])

for (i in 1:(length(intervall_ends) -1)) {
  # age_bin_groups = ifelse(((age_bin_groups >= intervall_ends[i]) & (age_bin_groups <= intervall_ends[i+1])), floor(intervall_ends[i] / 10) * 10, age_bin_groups)
  age_bin_groups = ifelse(((age_bin_groups >= intervall_ends[i]) & (age_bin_groups < intervall_ends[i+1])), intervall_ends[i], age_bin_groups)
}
annot[[sprintf("age_bin_%s_groups", number_of_bins)]] = age_bin_groups

write.table(annot,input_file_path, sep='\t', row.names = F)


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
