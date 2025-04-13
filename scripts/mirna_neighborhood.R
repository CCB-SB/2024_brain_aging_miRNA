suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))

# save rdata
#saveRDS(snakemake, file = "snakemake_mirna_neighborhood.rds")
#stop()

#snakemake = readRDS("snakemake_mirna_neighborhood.rds")

set.seed(snakemake@params$parameters_porps$set_seed)

print("mirna_neighborhood")


#---------------------------------- Functions ---------------------------------- 
# Function to parse the file and generate the required sets
process_miRNA_data = function(file_path) {
  # Read the file and filter for miRNA lines
  raw_data = read.table(file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE, comment.char = "#")
  colnames(raw_data) = c("chromosome", "Source", "Type", "Start", "End", "Score", "strand", "Phase", "Attributes")
  
  # Filter only mature miRNA entries
  miRNA_data = raw_data %>% filter(Type == "miRNA")
  
  # Extract Name from the Attributes column
  miRNA_data$name = sapply(strsplit(miRNA_data$Attributes, ";"), function(x) {
    gsub("Name=", "", x[grep("Name=", x)])
  })
  
  # Calculate the mean coordinate for each miRNA
  miRNA_data$mean_coordinate = rowMeans(miRNA_data[, c("Start", "End")])
  
  # Create a matrix with necessary fields
  miRNA_matrix = miRNA_data %>%
    select(name, chromosome, mean_coordinate, strand) %>%
    as.matrix()
  
  # Function to find miRNAs within a +/- 10kb range
  find_nearby_miRNAs = function(target, miRNA_data, strand_condition) {
    target_chr = target["chromosome"]
    target_coord = as.numeric(target["mean_coordinate"])
    target_strand = target["strand"]
    
    miRNA_data %>%
      filter(
        chromosome == target_chr,
        abs(as.numeric(mean_coordinate) - target_coord) <= 10000,
        eval(parse(text = strand_condition))
      ) %>%
      pull(name)
  }
  
  # Generate the three sets for each miRNA
  results = lapply(1:nrow(miRNA_matrix), function(i) {
    target = miRNA_matrix[i, ]
    list(
      name = target["name"],
      same_strand = unique(find_nearby_miRNAs(target, miRNA_data, "strand == target_strand")),
      opposite_strand = unique(find_nearby_miRNAs(target, miRNA_data, "strand != target_strand")),
      both_strands = unique(find_nearby_miRNAs(target, miRNA_data, "TRUE"))
    )
  })
  
  return(list(Matrix = miRNA_matrix, Results = results))
}


make_mirna_chromosome_df = function(processed_matrix) {
  miRNA_matrix = as.data.frame(processed_matrix, stringsAsFactors = FALSE)
  colnames(miRNA_matrix)[1] = "miRNA"
  
  # Ensure mean_coordinate is numeric for plotting
  miRNA_matrix$mean_coordinate = as.numeric(miRNA_matrix$mean_coordinate)
  
  # Sort chromosomes in natural order (1 to largest, then X, Y)
  miRNA_matrix$chromosome = factor(
    miRNA_matrix$chromosome,
    levels = c(
      paste0("chr", 1:22),  # Adjust to the largest chromosome number in your data
      "chrX",
      "chrY"
    )
  )
  
  # Remove unused levels
  miRNA_matrix$chromosome = droplevels(miRNA_matrix$chromosome)
  
  # Calculate cumulative offsets for genomic coordinates
  chromosome_offsets = miRNA_matrix %>%
    group_by(chromosome) %>%
    summarize(chromosome_start = min(mean_coordinate), chromosome_end = max(mean_coordinate)) %>%
    mutate(cumulative_offset = cumsum(lag(chromosome_end, default = 0)))
  
  # Add cumulative coordinates to miRNA_matrix
  miRNA_matrix = miRNA_matrix %>%
    left_join(chromosome_offsets, by = "chromosome") %>%
    mutate(cumulative_coordinate = mean_coordinate + cumulative_offset)
  
  # Calculate midpoints for chromosome labels
  chromosome_positions = chromosome_offsets %>%
    mutate(chromosome_mid = cumulative_offset + (chromosome_end - chromosome_start) / 2)
  
  return(list(miRNA_matrix=miRNA_matrix, chromosome_positions=chromosome_positions))
}

get_sig_corr_mirnas = function(miRNA_matrix, corr_th=0, sig_th=0.05) {
  miRNA_matrix_sig = miRNA_matrix[miRNA_matrix$adj_p_value < sig_th,]
  miRNA_matrix_corr_sig = miRNA_matrix_sig[abs(miRNA_matrix_sig$corr) > corr_th,]
  return(miRNA_matrix_corr_sig$miRNA)
}


count_mirnas_on_strands = function(miRNA_sets, mirna_list_sig_corr) {
  # find the sig. corr. mirnas in the same and opposite strand for each sig. corr. mirna from mmu.gff3
  miRNA_sets_sig_corr = list()
  for (one_set in miRNA_sets) {
    if (one_set$name %in% mirna_list_sig_corr) {
      # get the ones from the same strand an look if the are sig. correlated
      same_strand = intersect(one_set$same_strand, mirna_list_sig_corr)
      opposite_strand = intersect(one_set$opposite_strand, mirna_list_sig_corr)
      miRNA_sets_sig_corr = append(miRNA_sets_sig_corr, list(list(name=unname(one_set$name), same_strand=same_strand, opposite_strand=opposite_strand)))
    }
  }
  # Group by miRNA name
  miRNA_sets_sig_corr_merged = lapply(unique(sapply(miRNA_sets_sig_corr, function(x) x$name)), function(miRNA_name) {
    # Filter all sub-lists with the same name
    sublists <- Filter(function(x) x$name == miRNA_name, miRNA_sets_sig_corr)
    # Merge same_strand and opposite_strand using unique union
    list(
      name = miRNA_name,
      same_strand = unique(unlist(lapply(sublists, function(x) x$same_strand))),
      opposite_strand = unique(unlist(lapply(sublists, function(x) x$opposite_strand)))
    )
  })
  # now count unique mirnas in same and opposite strand
  # output dataframe format: cols = c(miRNA, same, opposite)
  heatmap_df = data.frame(miRNA=character(), n_same=numeric(), n_opposite=numeric(), stringsAsFactors = FALSE)
  # get the correlated and significant mirnas
  for (one_set in miRNA_sets_sig_corr_merged) {
    # if (one_set$name == "mmu-miR-669a-5p") {
    #   print(one_set)
    # }
    # get the ones from the same strand an look if the are sig. correlated
    n_same = length(unique(one_set$same_strand)) - 1 # the original one is always included in "same_strand"
    n_opposite = length(unique(one_set$opposite_strand))
    heatmap_df = rbind(heatmap_df, c(one_set$name, n_same, n_opposite))
  }
  colnames(heatmap_df) = c("miRNA", "n_same", "n_opposite")
  heatmap_df$n_same = as.numeric(heatmap_df$n_same)
  heatmap_df$n_opposite = as.numeric(heatmap_df$n_opposite)

  # sum up rows with the same miRNA
  heatmap_df_sum = heatmap_df %>%
    group_by(miRNA) %>%
    summarise(
      n_same = sum(n_same, na.rm = TRUE), 
      n_opposite = sum(n_opposite, na.rm = TRUE),
      .groups = "drop"
    )
  return(heatmap_df_sum)
}


manhattan_plot = function(miRNA_matrix, chromosome_positions, corr_th, sig_th, plot_title="All", m, time, chr_colors) {
  # Filter only corr and sig rows
  miRNA_matrix_sig = miRNA_matrix[miRNA_matrix$miRNA %in% get_sig_corr_mirnas(miRNA_matrix, corr_th=0, sig_th=sig_th),]  # corr_th=0 to get all sig miRNAs
  miRNA_matrix_corr_sig = miRNA_matrix[miRNA_matrix$miRNA %in% get_sig_corr_mirnas(miRNA_matrix, corr_th=corr_th, sig_th=sig_th),]
  # Generate the Manhattan plot
  p = ggplot(miRNA_matrix, aes(x = cumulative_coordinate, y = corr, color = colour_code)) +
    geom_point(size = 0.75) + #, alpha = 0.7
    #geom_point(data = miRNA_matrix_sig, mapping = aes(x = cumulative_coordinate, y = corr),  size = 0.5, fill=NA, stroke = 0.5, color="black", shape = 21, ) + #alpha = 0.7
    #scale_color_manual(values = rep(chr_colors, length.out = length(levels(miRNA_matrix$chromosome)))) +
    scale_color_manual(values = chr_colors) +
    scale_x_continuous(limits = c(min(miRNA_matrix$cumulative_coordinate), max(miRNA_matrix$cumulative_coordinate)), expand = expansion(mult = 0.01, add = 0), breaks = chromosome_positions$chromosome_mid, labels = chromosome_positions$chromosome) + 
    #geom_text(data = chromosome_positions, aes(x = chromosome_mid, y = 1.2, label = chromosome), inherit.aes = FALSE, angle = 90, vjust = 0.5, size = 5 / 14 * 2 / 3 *plots_props$font_size) +
    geom_vline(data = chromosome_positions, aes(xintercept = cumulative_offset + chromosome_end), linetype = "dotted", color = "black", linewidth=0.25) +
    geom_hline(yintercept = corr_th, linetype = "dashed", color = "black", linewidth=0.25) +
    geom_hline(yintercept = -corr_th, linetype = "dashed", color = "black", linewidth=0.25) +
    labs(
      title = plot_title,
      x = "", #"Cumulative Genomic Coordinate",
      y = "", #sprintf("%s correlation (%s)", m, time), #"Correlation Value",
      color = "Chromosome"
    ) +
    ylim(-1, 1) +
    theme_classic() +
    theme(#plot.margin = unit(c(0.1,-1.9,0.25,-1.75), plots_props$image_units), #top, right, bottom, left
      plot.margin = unit(c(0,0.1,0,0.1), plots_props$image_units), #top, right, bottom, left
      #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9), 
      plot.background = element_rect(fill='transparent', color=NA),
      text = element_text(family = plots_props$font_family, size = plots_props$font_size),
      axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
      axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
      plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size, margin = margin(t = unit(10, "pt"), r = unit(5, "pt"), b = unit(2, "pt"), l = unit(5, "pt"))), # Using pt for margins      legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
      legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(family = plots_props$font_family, size = plots_props$font_size, angle = 90, vjust = 0.5, hjust=1),
      #axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )
  #p
  
  return(p)
}


heatmap_count_plot = function(heatmap_df, max_same=NA, max_opposite=NA) {
  # count the occurrences for each (n_same, n_opposite)
  heatmap_counts = as.data.frame(table(heatmap_df$n_same, heatmap_df$n_opposite))
  colnames(heatmap_counts) = c("n_same", "n_opposite", "num_miRNAs")
  # Convert numeric columns back from factors (if necessary)
  heatmap_counts$n_same = as.numeric(as.character(heatmap_counts$n_same))
  heatmap_counts$n_opposite = as.numeric(as.character(heatmap_counts$n_opposite))
  # View result
  #print(heatmap_counts)
  
  # Define the desired size of the matrix
  if (is.na(max_same)) {
    max_same = max(c(heatmap_counts$n_same, 1))
  }
  if (is.na(max_opposite)) {
    max_opposite = max(c(heatmap_counts$n_opposite, 1))
  }
  
  # Create all possible (n_same, n_opposite) combinations
  all_combinations = expand.grid(n_same = 0:max_same, n_opposite = 0:max_opposite)
  
  # Merge with existing data, filling missing values with 0
  full_heatmap_counts = all_combinations %>%
    left_join(heatmap_counts, by = c("n_same", "n_opposite")) %>%
    replace_na(list(num_miRNAs = 0))
  
  # Convert to matrix format
  heatmap_matrix = full_heatmap_counts %>%
    pivot_wider(names_from = n_opposite, values_from = num_miRNAs, values_fill = 0) %>%
    column_to_rownames(var = "n_same") %>%
    as.matrix()
  
  # Flip matrix row order so (0,0) is bottom-left
  # transpose to get opposite as y axis
  heatmap_matrix = t(heatmap_matrix)
  heatmap_matrix = heatmap_matrix[nrow(heatmap_matrix):1, ]
  
  # Create row annotations (labels on the left side)
  row_annotation = rowAnnotation(
    n_opposite = anno_text(rownames(heatmap_matrix), gp = gpar(fontsize = 9))
  )
  
  # Column annotation with 90-degree rotated and centered labels
  # col labels only for even
  even_cols = seq(2, ncol(heatmap_matrix), by = 2)  # Get column indices that are even
  col_labels = colnames(heatmap_matrix)
  col_labels[even_cols] = ""
  col_annotation = HeatmapAnnotation(
    n_same = anno_text(col_labels, rot = 0, just = "center", gp = gpar(fontsize = 9, vjust = 1), location = unit(0, "mm"))
  )
  
  # Define function to add text annotations in heatmap cells
  cell_fun = function(j, i, x, y, width, height, fill) {
    # black border
    grid.rect(x, y, width, height, gp = gpar(fill=fill, col = "grey", lwd = 1))  # Black cell border
    value = heatmap_matrix[i, j]
    if (value > 100) {
      text_col = "white"
    } else {
      text_col = "black"
    }
    if (value != 0) {
      grid.text(value, x, y, gp = gpar(fontsize = 6, col = text_col))
    }
  }
  
  # Plot the heatmap
  p_heat = ComplexHeatmap::Heatmap(
    heatmap_matrix,
    col = colorRamp2(c(0, 2, max(heatmap_matrix)), c("#FFFFFF","#DFDFB9", "#525200")),  # Gradient color
    cluster_rows = FALSE,  # No clustering
    cluster_columns = FALSE,  # No clustering
    show_row_names = FALSE,  # Hide row names since annotation is added
    show_column_names = FALSE,  # Show column labels
    left_annotation = row_annotation,  # Place row numbers on the left
    bottom_annotation = col_annotation,
    cell_fun = cell_fun,  # Add text labels for nonzero values
    show_heatmap_legend = FALSE,  # REMOVE COLOR BAR COMPLETELY
    column_title = "Number of sig. corr. miRNAs on the same strand",  # X-axis label
    row_title = "Number on\nopposite strand    ",   # Y-axis label
    column_title_side = "bottom",
    column_title_gp = gpar(fontsize = 9, just="left"),
    row_title_gp = gpar(fontsize = 9, just="left")
  )
  return(p_heat)
}


#---------------------------------- Load data ---------------------------------- 
data_input = snakemake@params$data_input
method = snakemake@params$method
props = snakemake@params$props
sig_th = snakemake@params$thresholds$adj_p_value
corr_ths = snakemake@params$thresholds$corr_th
file_path_mmu = "../../data_external/mmu.gff3"
xticks_names = snakemake@params$xticks_names
plots_props = snakemake@params$plots_props
colours = snakemake@params$colors
plot_layout = snakemake@params$plot_layout
plot_size = snakemake@params$plot_size
results_folder = snakemake@params$results_folder

chr_colors = c("chr1" = "#C6C6C6",
               "chr2" = "#A7A7A7",
               "chr3" = "#C6C6C6",
               "chr4" = "#A7A7A7",
               "chr5" = "#C6C6C6",
               "chr6" = "#A7A7A7",
               "chr7" = "#C6C6C6",
               "chr8" = "#A7A7A7",
               "chr9" = "#C6C6C6",
               "chr10" = "#A7A7A7",
               "chr11" = "#C6C6C6",
               "chr12" = "#A7A7A7",
               "chr13" = "#C6C6C6",
               "chr14" = "#A7A7A7",
               "chr15" = "#C6C6C6",
               "chr16" = "#A7A7A7",
               "chr17" = "#C6C6C6",
               "chr18" = "#A7A7A7",
               "chr19" = "#C6C6C6",
               "chrX" = "#A7A7A7",
               "pos" = colours$direction$pos,
               "neg" = colours$direction$neg,
               "sig_pos" = colours$direction$sig_pos,
               "sig_neg" = colours$direction$sig_neg
) #c("#459BAB", "#CC6928") #c("#E9BD5B", "#91C2AB")

# input folder
input_folder_tab = sprintf("%s/%s_%s/results_%s/matrices/corr_plots/age", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
input_folder_tab_per_tissue = sprintf("%s/%s_%s/results_%s/matrices/aggregation_corr", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)

# create output folder
output_folder_fig_manhattan = sprintf("%s/%s_%s/results_%s/figures/mirna_neighborhood/manhattan_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_fig_manhattan, recursive=TRUE)
output_folder_fig_heatmap = sprintf("%s/%s_%s/results_%s/figures/mirna_neighborhood/heatmap_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_fig_heatmap, recursive=TRUE)
output_folder_fig_line = sprintf("%s/%s_%s/results_%s/figures/mirna_neighborhood/line_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_fig_line, recursive=TRUE)

output_folder_tab = sprintf("%s/%s_%s/results_%s/matrices/mirna_neighborhood/heatmap_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_tab, recursive=TRUE)

# Specify the file path and run the function
result = process_miRNA_data(file_path_mmu)
# Access the sets for each miRNA
miRNA_sets = result$Results
# Example output for the first miRNA
#print(miRNA_sets[[1]])

# Convert miRNA_matrix to a data frame
processed_result = make_mirna_chromosome_df(result$Matrix)
miRNA_matrix = processed_result$miRNA_matrix
chromosome_positions = processed_result$chromosome_positions


#------------------------------------ Script ----------------------------------- 
for (m in method) {
  for (prop in props) {
    tissue = prop[1]
    time = prop[2]
    # load correlation and p values
    # load data
    #corr_values = fread(sprintf("%s/corr_method=%s_%ss_with_%s.csv", input_folder_tab, m, data_input$rna_class, time), sep='\t')[, -1, with = FALSE]
    #p_values = fread(sprintf("%s/corr_method=%s_adj_pvalues_%ss_with_%s.csv", input_folder_tab, m, data_input$rna_class, time), sep='\t')[, -1, with = FALSE]
    #corr_value_df = merge(corr_values, p_values, "miRNA")
    
    for (corr_th in corr_ths) {
      # # merge to keep only interesting mirnas
      # miRNA_matrix_with_corr = merge(miRNA_matrix, corr_value_df, by = "miRNA")
      # 
      # # create colour coding for different chromosoms and for pos., neg., sig. pos. and sig. neg. miRNAs
      # miRNA_matrix_with_corr$colour_code = as.character(miRNA_matrix_with_corr$chromosome)
      # #miRNA_matrix_with_corr$colour_code = ifelse(miRNA_matrix_with_corr$corr >= corr_th, "pos", miRNA_matrix_with_corr$colour_code)
      # #miRNA_matrix_with_corr$colour_code = ifelse(miRNA_matrix_with_corr$corr <= -corr_th, "neg", miRNA_matrix_with_corr$colour_code)
      # miRNA_matrix_with_corr$colour_code = ifelse(((miRNA_matrix_with_corr$corr >= corr_th) & (miRNA_matrix_with_corr$adj_p_value < sig_th)), "sig_pos", miRNA_matrix_with_corr$colour_code)
      # miRNA_matrix_with_corr$colour_code = ifelse(((miRNA_matrix_with_corr$corr <= -corr_th) & (miRNA_matrix_with_corr$adj_p_value < sig_th)), "sig_neg", miRNA_matrix_with_corr$colour_code)
      # 
      # # Generate the Manhattan plot
      # p = manhattan_plot(miRNA_matrix_with_corr, chromosome_positions, corr_th = corr_th, sig_th = sig_th, plot_title = "All", m, time, chr_colors)
      # file_path = sprintf("%s/manhattan_%s_corr_method=%s_corr_th=%s_%ss_with_%s", output_folder_fig_manhattan, "all", m, corr_th, data_input$rna_class, time)
      # ggsave(sprintf("%s.png", file_path), p)
      # 
      # # prepare the data for the heatmap
      # # get the correlated and significant mirnas
      # mirna_list_sig_corr = get_sig_corr_mirnas(miRNA_matrix_with_corr, corr_th = corr_th, sig_th = sig_th)
      # if (length(mirna_list_sig_corr) != 0) {
      #   # output dataframe format: cols = c(miRNA, same, opposite)
      #   heatmap_df = count_mirnas_on_strands(miRNA_sets, mirna_list_sig_corr)
      #   
      #   # Convert to matrix format for ComplexHeatmap
      #   p_heat = heatmap_count_plot(heatmap_df)
      #   file_path = sprintf("%s/count_heatmap_%s_corr_method=%s_corr_th=%s_%ss_with_%s", output_folder_fig_heatmap, "all", m, corr_th, data_input$rna_class, time)
      #   ggsave(sprintf("%s.png", file_path), grid.grabExpr(draw(p_heat)))
      #   
      # } else {
      #   print(sprintf("No sig. corr. miRNAs for tissue %s", "all"))
      # }
      
      
      ############
      # tissue specific
      # load correlation and p values
      corr_values = fread(sprintf("%s/corr_method=%s_%ss_with_%s_per_%s.csv", input_folder_tab_per_tissue, m, data_input$rna_class, time, tissue), sep = "\t",)[, -1, with = FALSE]
      p_values = fread(sprintf("%s/corr_method=%s_adj_pvalues_%ss_with_%s_per_%s.csv", input_folder_tab_per_tissue, m, data_input$rna_class, time, tissue), sep = "\t",)[, -1, with = FALSE]
      
      n_sig_corr_mirnas = data.frame()
      corr_value_dfs = list()
      heatmap_dfs = list()
      manhattan_plots = list()
      count_heatmap_plots = list()
      for (ti in colnames(corr_values)[-1]) {
        corr_values_ti = corr_values[, ..ti]
        colnames(corr_values_ti) = c("corr")
        corr_values_ti$miRNA = corr_values$miRNA
        p_values_ti = p_values[, ..ti]
        colnames(p_values_ti) = c("adj_p_value")
        p_values_ti$miRNA = p_values$miRNA
        corr_value_df_ti = merge(corr_values_ti, p_values_ti, "miRNA")
        corr_value_dfs[[ti]] = corr_value_df_ti
        
        # merge to keep only interesting mirnas
        miRNA_matrix_with_corr_ti = merge(miRNA_matrix, corr_value_df_ti, by = "miRNA")
        
        miRNA_matrix_with_corr_ti$colour_code = as.character(miRNA_matrix_with_corr_ti$chromosome)
        #miRNA_matrix_with_corr_ti$colour_code = ifelse(miRNA_matrix_with_corr_ti$corr >= corr_th, "pos", miRNA_matrix_with_corr_ti$colour_code)
        #miRNA_matrix_with_corr_ti$colour_code = ifelse(miRNA_matrix_with_corr_ti$corr <= -corr_th, "neg", miRNA_matrix_with_corr_ti$colour_code)
        miRNA_matrix_with_corr_ti$colour_code = ifelse(((miRNA_matrix_with_corr_ti$corr >= corr_th) & (miRNA_matrix_with_corr_ti$adj_p_value < sig_th)), "sig_pos", miRNA_matrix_with_corr_ti$colour_code)
        miRNA_matrix_with_corr_ti$colour_code = ifelse(((miRNA_matrix_with_corr_ti$corr <= -corr_th) & (miRNA_matrix_with_corr_ti$adj_p_value < sig_th)), "sig_neg", miRNA_matrix_with_corr_ti$colour_code)
        
        # Generate the Manhattan plot
        manhattan_plots[[ti]] = manhattan_plot(miRNA_matrix_with_corr_ti, chromosome_positions, corr_th = corr_th, sig_th = sig_th, plot_title=ti, m, time, chr_colors)
        file_path = sprintf("%s/manhattan_%s_corr_method=%s_corr_th=%s_%ss_with_%s", output_folder_fig_manhattan, ti, m, corr_th, data_input$rna_class, time, tissue)
        ggsave(sprintf("%s.png", file_path), manhattan_plots[[ti]])
        
        # prepare the data for the heatmap
        # get the correlated and significant mirnas
        mirna_list_sig_corr_ti = get_sig_corr_mirnas(miRNA_matrix_with_corr_ti, corr_th = corr_th, sig_th = sig_th)
        if (length(mirna_list_sig_corr_ti) != 0) {
          # output dataframe format: cols = c(miRNA, same, opposite)
          heatmap_df = count_mirnas_on_strands(miRNA_sets, mirna_list_sig_corr_ti)
          
          # Convert to matrix format for ComplexHeatmap
          count_heatmap_plots[[ti]] = heatmap_count_plot(heatmap_df)
          file_path = sprintf("%s/count_heatmap_%s_corr_method=%s_corr_th=%s_%ss_with_%s", output_folder_fig_heatmap, ti, m, corr_th, data_input$rna_class, time, tissue)
          ggsave(sprintf("%s.png", file_path), grid.grabExpr(draw(count_heatmap_plots[[ti]])))
          
          heatmap_dfs[[ti]] = heatmap_df
        } else {
          print(sprintf("No sig. corr. miRNAs for tissue %s", ti))
        }
      }
      
      # sort the figures given the order from the config
      manhattan_plots_sorted = manhattan_plots[names(colours[[tissue]])]
      count_heatmap_plots = count_heatmap_plots[names(colours[[tissue]])]
      
      plot_nrows = plot_layout$plot_nrows
      plot_ncols = plot_layout$plot_ncols
      
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
      
      manhattan_plots_sorted_mod_margins = manhattan_plots_sorted
      count_heatmap_plots_mod_margins = count_heatmap_plots
      for (img_i in 1:length(manhattan_plots_sorted_mod_margins)) {
        if ((!(img_i %in% left_col)) & (!(img_i %in% right_col))) {
          manhattan_plots_sorted_mod_margins[[img_i]] = manhattan_plots_sorted_mod_margins[[img_i]] + 
            theme(plot.margin = margin(-8, -1, -10, -13, "pt")
            )
          #count_heatmap_plots_mod_margins[[img_i]] = count_heatmap_plots_mod_margins[[img_i]] + 
            #theme(plot.margin = margin(10, -1, -7, -13, "pt")
            #)
        }
        else if (img_i %in% left_col) {
          manhattan_plots_sorted_mod_margins[[img_i]] = manhattan_plots_sorted_mod_margins[[img_i]] + 
            theme(plot.margin = margin(-8, -4, -10, -10, "pt")
            )
          #count_heatmap_plots_mod_margins[[img_i]] = count_heatmap_plots_mod_margins[[img_i]] + 
          #  theme(plot.margin = margin(10, -4, -7, -10, "pt")
          #  )
        } else if (img_i %in% right_col) {
          manhattan_plots_sorted_mod_margins[[img_i]] = manhattan_plots_sorted_mod_margins[[img_i]] + 
            theme(plot.margin = margin(-8, 4, -10, -18, "pt")
            )
          #count_heatmap_plots_mod_margins[[img_i]] = count_heatmap_plots_mod_margins[[img_i]] + 
          #  theme(plot.margin = margin(10, 4, -7, -18, "pt")
          #  )
        }
        if (!(img_i %in% left_col)) {
          manhattan_plots_sorted_mod_margins[[img_i]] = manhattan_plots_sorted_mod_margins[[img_i]] + 
            #scale_y_continuous(labels = function(breaks) {rep_along(breaks, "  ")}, limits = c(0, 1.15*maxi_expr_rna), expand = c(0,0)) +
            theme(axis.ticks.y = element_line(color = "transparent"),
                  axis.text.y = element_text(color = "transparent")
            )
          #count_heatmap_plots_mod_margins[[img_i]] = count_heatmap_plots_mod_margins[[img_i]] + 
            #scale_y_continuous(labels = function(breaks) {rep_along(breaks, "  ")}, limits = c(0, 1.15*maxi_expr_rna), expand = c(0,0)) +
            #theme(axis.ticks.y = element_line(color = "transparent"),
            #      axis.text.y = element_text(color = "transparent")
            #)
        }
      }
      
      if (m == "spearman") {
        y_axis_title = sprintf("Spearman correlation (%s)", time)
      } else if (m == "pearson") {
        y_axis_title = sprintf("Pearson correlation (%s)", time)
      } else {
        y_axis_title = sprintf("%s correlation (%s)", m, time)
      }
      
      x_axis_title = "Cumulative Genomic Coordinate"
      all_manhattan_plots_sorted_mod_margins = arrangeGrob(grobs=manhattan_plots_sorted_mod_margins, ncol=plot_ncols, nrow=plot_nrows, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90)) #, right = textGrob(rna, gp=gpar(fontsize=plots_props$font_size_header, fontfamily=plots_props$font_family), rot = 270)
      #all_count_heatmap_plots_mod_margins = arrangeGrob(grobs=count_heatmap_plots_mod_margins, ncol=plot_ncols, nrow=plot_nrows, bottom = textGrob(x_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family)), left = textGrob(y_axis_title, gp=gpar(fontsize=plots_props$font_size, fontfamily=plots_props$font_family), rot = 90)) #, right = textGrob(rna, gp=gpar(fontsize=plots_props$font_size_header, fontfamily=plots_props$font_family), rot = 270)
    
      plot_width = plot_size$image_width
      plot_height = plot_size$image_height
      
      file_name = sprintf("%s/manhattan_corr_method=%s_corr_th=%s_%ss_with_%s", output_folder_fig_manhattan, m, corr_th, data_input$rna_class, time)
      ggsave(sprintf("%s_%sx%s.png", file_name, plot_width, plot_height), all_manhattan_plots_sorted_mod_margins, dpi=plots_props$dpi, width=plot_width, height=plot_height, units = plots_props$image_units)
      ggsave(sprintf("%s_%sx%s.svg", file_name, plot_width, plot_height), all_manhattan_plots_sorted_mod_margins, width=plot_width, height=plot_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
      
      ##########################
      # make combined heatmap plot
      combined_heatmap_dfs = do.call(rbind, heatmap_dfs)
      combined_heatmap = heatmap_count_plot(combined_heatmap_dfs)
      
      height = 4.5
      width = 12
      ggsave(sprintf("%s/count_heatmap_corr_method=%s_corr_th=%s_%ss_with_%s_%sx%s.png", output_folder_fig_heatmap, m, corr_th, data_input$rna_class, time, width, height), grid.grabExpr(draw(combined_heatmap)), width = width, height = height, dpi = plots_props$dpi, units = plots_props$image_units)
      ggsave(sprintf("%s/count_heatmap_corr_method=%s_corr_th=%s_%ss_with_%s_%sx%s.svg", output_folder_fig_heatmap, m, corr_th, data_input$rna_class, time, width, height), grid.grabExpr(draw(combined_heatmap)), width = width, height = height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
      
      #file_name = sprintf("%s/count_heatmap_corr_method=%s_corr_th=%s_%ss_with_%s", output_folder_fig_heatmap, m, corr_th, data_input$rna_class, time, tissue)
      #ggsave(sprintf("%s_%sx%s.png", file_name, plot_width, plot_height), all_count_heatmap_plots_mod_margins, dpi=plots_props$dpi, width=plot_width, height=plot_height, units = plots_props$image_units)
      #ggsave(sprintf("%s_%sx%s.svg", file_name, plot_width, plot_height), all_count_heatmap_plots_mod_margins, width=plot_width, height=plot_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
      
      rm(all_manhattan_plots_sorted_mod_margins)
      rm(all_count_heatmap_plots_mod_margins)
      rm(manhattan_plots_sorted_mod_margins)
      rm(count_heatmap_plots_mod_margins)
      rm(manhattan_plots_sorted)
      rm(count_heatmap_plots)
      
      gc()
      
      
      # Combine all list elements into one dataframe
      heatmap_dfs_combined = do.call(rbind, lapply(names(heatmap_dfs), function(name) {
        df = heatmap_dfs[[name]]  # Extract dataframe
        df$tissue <- name  # Add a new column with the list element name
        df  # Return the modified dataframe
      }))
      # remove rows where both are 0
      # heatmap_dfs_combined = heatmap_dfs_combined[!(heatmap_dfs_combined$n_same == 0 & heatmap_dfs_combined$n_opposite == 0), ]
      
      fwrite(heatmap_dfs_combined, sprintf("%s/heatmap_dfs_combined.csv", output_folder_tab), sep = "\t", row.names = FALSE,)
      
      # Identify duplicated miRNAs
      #duplicated_miRNAs = heatmap_dfs_combined[heatmap_dfs_combined$miRNA %in% heatmap_dfs_combined$miRNA[duplicated(heatmap_dfs_combined$miRNA)], ]
      # View result
      #print(duplicated_miRNAs)
      
      
      ##########################
      # make line plot for each mirna
      # Summarize data: count miRNAs per (n_same, tissue) combination
      heatmap_summary_same = heatmap_dfs_combined %>%
        group_by(tissue, n_same) %>%
        summarise(num_miRNAs = n(), .groups = "drop")
      # Determine the full range of n_same values
      n_same_max = max(heatmap_summary_same$n_same, na.rm = TRUE)  # Find global max n_same
      
      # Summarize data: count miRNAs per (n_same, tissue) combination
      heatmap_summary_opposite = heatmap_dfs_combined %>%
        group_by(tissue, n_opposite) %>%
        summarise(num_miRNAs = n(), .groups = "drop")
      # Determine the full range of n_opposite values
      n_opposite_max = max(heatmap_summary_opposite$n_opposite, na.rm = TRUE)  # Find global max n_opposite
      
      n_max = max(c(n_same_max, n_opposite_max))
      
      # Create all possible (tissue, n_same) combinations
      all_n_same = expand.grid(n_same = 0:n_max, tissue = unique(heatmap_dfs_combined$tissue))
      all_n_opposite = expand.grid(n_opposite = 0:n_max, tissue = unique(heatmap_dfs_combined$tissue))
      
      # Merge with actual data, filling missing num_miRNAs with 0
      # same
      heatmap_summary_expanded_same = all_n_same %>%
        left_join(heatmap_summary_same, by = c("n_same", "tissue")) %>%
        replace_na(list(num_miRNAs = 0))
      # opposite
      heatmap_summary_expanded_opposite = all_n_opposite %>%
        left_join(heatmap_summary_opposite, by = c("n_opposite", "tissue")) %>%
        replace_na(list(num_miRNAs = 0))
      
      # divide the heatmap_summary_expanded by the respective sig. corr. mirnas per tissue
      heatmap_summary_expanded_same$percent_miRNAs = 0
      heatmap_summary_expanded_opposite$percent_miRNAs = 0
      for (ti in n_sig_corr_mirnas$tissue) {
        percent_same = heatmap_summary_expanded_same[heatmap_summary_expanded_same$tissue == ti,]$num_miRNAs / sum(heatmap_summary_expanded_same[heatmap_summary_expanded_same$tissue == ti,]$num_miRNAs)
        percent_opposite = heatmap_summary_expanded_opposite[heatmap_summary_expanded_opposite$tissue == ti,]$num_miRNAs / sum(heatmap_summary_expanded_opposite[heatmap_summary_expanded_opposite$tissue == ti,]$num_miRNAs)
        
        heatmap_summary_expanded_same[heatmap_summary_expanded_same$tissue == ti,]$percent_miRNAs = 100 * percent_same
        heatmap_summary_expanded_opposite[heatmap_summary_expanded_opposite$tissue == ti,]$percent_miRNAs = 100 * percent_opposite
      }
      
      # Create the smoothed line plot
      p_line_same = ggplot(heatmap_summary_expanded_same, aes(x = n_same, y = num_miRNAs, color = tissue, group = tissue)) +
        geom_point(size = 2) +  # Add points for raw data
        geom_line(linewidth = 1.2) +  # lines
        labs(x = "Number of sig. corr. miRNAs on the same strand", y = "Number of sig. corr. miRNAs", title = "") +
        theme_minimal() +
        theme(legend.title = element_blank())  # Remove legend title
      
      
      # now make the line plot in percent
      p_line_percent_same = ggplot(heatmap_summary_expanded_same, aes(x = n_same, y = percent_miRNAs, color = tissue, group = tissue)) +
        geom_point(size = 2) +  # Add points for raw data
        geom_line(linewidth = 1.2) +  # lines
        labs(x = "Number of sig. corr. miRNAs on the same strand", y = "Proportion of sig. corr. miRNAs (%)", title = "") +
        theme_minimal() +
        theme(legend.title = element_blank())  # Remove legend title
      
      
      # combined plot
      # x-axis = number of sig. corr. mirnas
      # y-axis pos same strand
      # y-axis neg opposite strand
      mask_same = heatmap_summary_expanded_same$num_miRNAs != 0
      mask_opposite = heatmap_summary_expanded_opposite$num_miRNAs != 0
      
      p_line_both = ggplot() +
        geom_point(data=heatmap_summary_expanded_same[mask_same, ], mapping=aes(x = num_miRNAs, y = n_same, color = tissue, group = tissue), size = 2) +
        geom_line(data=heatmap_summary_expanded_same[mask_same, ], mapping=aes(x = num_miRNAs, y = n_same, color = tissue, group = tissue), linewidth = 1.2) +  # lines
        geom_point(data=heatmap_summary_expanded_opposite[mask_opposite, ], mapping=aes(x = num_miRNAs, y = -n_opposite, color = tissue, group = tissue), size = 2) +
        geom_line(data=heatmap_summary_expanded_opposite[mask_opposite, ], mapping=aes(x = num_miRNAs, y = -n_opposite, color = tissue, group = tissue), linewidth = 1.2) +  # lines
        labs(x = "Number of sig. corr. miRNAs", y = "Number of sig. corr. miRNAs", title = "") +
        theme_minimal() +
        theme(legend.title = element_blank())  # Remove legend title
    }
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
