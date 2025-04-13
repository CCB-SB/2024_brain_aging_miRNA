suppressPackageStartupMessages(library(data.table))


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
log_list = snakemake@params$log_list


#------------------------------------ Script ----------------------------------- 
for (rna_class in data_input$rna_class_list) {
  for (norm in data_input$norm) {
    for (log_base in log_list) {
      # load expr in every loop
      tbl = fread(sprintf("%s_%s/%s_complete_quantification_%s.csv", data_input$preprocesing_data_folder, data_input$data_sub_set, rna_class, norm), sep='\t', header=T)
      
      # log-transformation 
      tbl[, (2:ncol(tbl)):=(log( tbl[,2:ncol(tbl)] + 1 , base = log_base))]
    
      # save log-transformed expr
      fwrite(tbl, sprintf("%s_%s/%s_complete_quantification_%s_log%s.csv", data_input$preprocesing_data_folder, data_input$data_sub_set, rna_class, norm, log_base), sep='\t')
    }
  }
}

#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
