suppressPackageStartupMessages(library(data.table))

#snakemake = readRDS("snakemake_filter_expression_files_per_group.rds")


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
data_input$feature_filtering = "filtered"
norm_list = snakemake@params$norm_list
min_detection_list = snakemake@params$min_detection_list
detection_rate_list = snakemake@params$detection_rate_list 
filter_variable = snakemake@params$filter_variable

# save rdata
#saveRDS(snakemake, file = "snakemake_filter_expression_files_per_group.rds")

annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set,data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 


#------------------------------------ Script ----------------------------------- 
for (min_detection in min_detection_list) {
  i = 1
  for (rna_class in data_input$rna_class_list) {
    for (norm_key in c("rpmm")) {
      print(norm_key)
      input_folder = sprintf("%s_%s", data_input$preprocesing_data_folder, data_input$data_sub_set)
      detect = fread(sprintf("%s/%s_filtered_detection_matrix_min_detect=%s.csv", input_folder, rna_class, min_detection), sep='\t', header=T)
      files = list.files(input_folder)
      
      file_prefix = sprintf("%s_filtered_quantification_%s", rna_class, norm_key)
      file_prefix_list = c()
      for (file in files) {
        if (startsWith(file, file_prefix)) {
          # exclude file if ".detection_rate_" is in file name
          if (!grepl(".detection_rate_", file, fixed = TRUE)) {
            file_prefix_list = c(file_prefix_list, file)
          }
        }
      }
    }

    for (file in file_prefix_list) {
      print(file)
      # load data table
      tbl = fread(sprintf("%s/%s", input_folder, file), sep='\t', header=T)
      
      # calculate features to keep
      
      if(colnames(detect)[1] == "V1"){
        setnames(detect, "V1", rna_class)
      }
      
      for (detection_rate in detection_rate_list) {
        to_keep = rep(F, nrow(detect))
        for(g in unique(annot[,get(filter_variable)])){
          if(g != ""){
            sub = annot[get(filter_variable) == g, ID]
            to_keep = to_keep | rowMeans(detect[,..sub]) >= detection_rate
          }
        }
        
        # filter table
        #print(tbl[,1])
        #print(detect[,1])
        stopifnot(assertthat::are_equal(tbl[,1], detect[,1]))
        
        # remove .csv from end
        filename_to_save = gsub('.csv$', '', file)

        # save filtered table
        fwrite(tbl[to_keep], sprintf("%s/%s.detection_rate_%sp_per_%s.csv", input_folder, filename_to_save, (detection_rate * 100), filter_variable), sep='\t')
      }
    }
    i = i + 1
  }
}

#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
