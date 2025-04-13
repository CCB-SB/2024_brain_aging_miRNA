suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(svglite))


#snakemake = readRDS("snakemake_alignment_stats_bar_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("alignment_stats_bar_plot")


# ---------------------- Load data ----------------------
data_input = snakemake@params$data_input
xticks_names = snakemake@params$xticks_names
plots_props = snakemake@params$plots_props
results_folder = snakemake@params$results_folder
grouping = snakemake@params$grouping[1]
subgrouping = snakemake@params$grouping[2]
number_of_reads = snakemake@params$number_of_reads

# save rdata
#saveRDS(snakemake, file = sprintf("snakemake_alignment_stats_bar_plot.rds"))

colors = c("aligned"="#788A38", "not aligned"="grey")

annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
mapping_info_ori = fread(snakemake@input$mapping_info, sep='\t')

output_folder = sprintf("%s/%s_%s/results_%s/figures/mapping_statistics/bar_plots/alignment", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder, recursive=TRUE)


# ----------------------- Script -----------------------
if (unique(grepl('BA10_CA1_cp_26_1', mapping_info_ori$fastq_name)) == FALSE) {
  mapping_info = mapping_info_ori
} else {
  mapping_info = mapping_info_ori[mapping_info_ori$fastq_name != "BA10_CA1_cp_26_1", ]
}
mapping_info_1 = mapping_info[mapping_info$Mismatches == 1, ]

# keep only samples which are in annot
tmp = read.table(text=mapping_info_1$fastq_name, sep = "_", fill = TRUE, as.is = TRUE)
if (startsWith(data_input$data_sub_set, "CA1") | startsWith(data_input$data_sub_set, "CA2")) {
  mapping_info_1[[data_input$identifier_column]] = sprintf("%s_%s_%s", tmp[,2], tmp[,3], tmp[,4])
} else if (startsWith(data_input$data_sub_set, "microglia")) {
  tmp = (colnames(annot) %in% c("fastq_name", data_input$identifier_column))
  mapping_info_1 = merge(mapping_info_1, annot[,..tmp], by = "fastq_name")
}

mask = mapping_info_1[[data_input$identifier_column]] %in% annot[[data_input$identifier_column]]

mapping_info_1_filtered = mapping_info_1[mask,]

mapping_info_samples = merge(mapping_info_1, annot, by = data_input$identifier_column)

reads_aligned = c()
reads_processed = c()
for (sample in mapping_info_samples[[data_input$identifier_column]]){
  tmp_total = mapping_info_samples[sample,]$reads_processed / 1000000
  tmp_aligned = mapping_info_samples[sample,]$reads_aligned / 1000000
  reads_processed = append(reads_processed, tmp_total)
  reads_aligned = append(reads_aligned, tmp_aligned)
}

if (subgrouping == "") {
  plot_df = data.frame(sample=mapping_info_samples[[data_input$identifier_column]], aligned=reads_aligned, processed=reads_processed, group=mapping_info_samples[[grouping]], check.names = F)
} else {
  plot_df = data.frame(sample=mapping_info_samples[[data_input$identifier_column]], aligned=reads_aligned, processed=reads_processed, group=mapping_info_samples[[grouping]], subgroup=mapping_info_samples[[subgrouping]], check.names = F)
}

plot_df["not aligned"] = plot_df$processed - plot_df$aligned

# remove the total column (will show as the sum in the bar plot)
plot_df$processed = c()
# melt for ggplot
if (subgrouping == "") {
  plot_df = melt(plot_df, id=c("sample", "group"), variable.name="alignment", value.name="counts")
} else {
  plot_df = melt(plot_df, id=c("sample", "group", "subgroup"), variable.name="alignment", value.name="counts")
}
plot_df$alignment = factor(plot_df$alignment, levels=c("not aligned", "aligned"))

# order decreasing by column "counts" (only values where column alignment == aligned)
# get the numeric order
ord = order(plot_df[plot_df$alignment == "aligned",]$counts, decreasing = TRUE)
# get the names
x_order = plot_df[ord,]$sample

# factor(sample, levels=x_order) sets the order
p = ggplot(plot_df, aes(x=factor(sample, levels=x_order), y=counts, fill=alignment)) + 
  geom_bar(stat="identity") +
  geom_hline(yintercept=as.numeric(number_of_reads)/1000000, linetype="dashed", color = "black", size=0.5) +
  labs(x="", y="Processed reads\nin millions") +
  scale_fill_manual(breaks = c("aligned", "not aligned"), values=colors) +
  scale_y_continuous(expand = c(0, 0)) +
  #ggtitle("Aligned counts") +
  theme_classic() +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        axis.text = element_text(angle = 0, family = plots_props$font_family, size = plots_props$font_size),
        axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
        plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header),
        legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        # legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4),
        legend.margin=margin(t = -8), legend.key.size = unit(0.3, "cm"), legend.spacing.y = unit(0.4, "cm")
  )

# Save results to storage
ggsave(sprintf("%s/mapping_stats_9x6.png", output_folder), p, dpi=plots_props$dpi, width=plots_props$image_width, height=plots_props$image_height, units = plots_props$image_units)
ggsave(sprintf("%s/mapping_stats_9x6.svg", output_folder), p, width = plots_props$image_width, height = plots_props$image_height, unit = plots_props$image_units, dpi = plots_props$dpi)

ggsave(sprintf("%s/mapping_stats_6x6.png", output_folder), p, dpi=plots_props$dpi, width = 2/3 * plots_props$image_width, height=plots_props$image_height, units = plots_props$image_units)
ggsave(sprintf("%s/mapping_stats_6x6.svg", output_folder), p, width = 2/3 * plots_props$image_width, height = plots_props$image_height, unit = plots_props$image_units, dpi = plots_props$dpi)


if (subgrouping != "") {
  # grouping
  reads_aligned_group_median = c()
  reads_not_aligned_group_median = c()
  sample_names = c()
  reads_aligned_group_subgroup_median = c()
  reads_not_aligned_group_subgroup_median = c()
  for (g in unique(annot[[grouping]])) {
    mask = (plot_df$group == g)
    plot_df_group = plot_df[mask,]
    
    reads_aligned = median(plot_df_group[plot_df_group$alignment == "aligned",]$counts)
    reads_not_aligned = median(plot_df_group[plot_df_group$alignment == "not aligned",]$counts)
    
    reads_aligned_group_median = append(reads_aligned_group_median, reads_aligned)
    reads_not_aligned_group_median = append(reads_not_aligned_group_median, reads_not_aligned)
    
    for (sg in unique(annot[[subgrouping]])) {
      mask = (plot_df_group$subgroup == sg)
      plot_df_group_subgroup = plot_df_group[mask,]
      
      reads_aligned = median(plot_df_group_subgroup[plot_df_group_subgroup$alignment == "aligned",]$counts)
      reads_not_aligned = median(plot_df_group_subgroup[plot_df_group_subgroup$alignment == "not aligned",]$counts)
      
      sample_names = append(sample_names, sprintf("%s_%s", g, sg))
      reads_aligned_group_subgroup_median = append(reads_aligned_group_subgroup_median, reads_aligned)
      reads_not_aligned_group_subgroup_median = append(reads_not_aligned_group_subgroup_median, reads_not_aligned)
    }
  }
  
  plot_df_group_median = data.frame(sample=unique(annot[[grouping]]), aligned=reads_aligned_group_median, not_aligned=reads_not_aligned_group_median)
  colnames(plot_df_group_median) = c("sample", "aligned", "not aligned")
  plot_df_group_subgroup_median = data.frame(sample=sample_names, aligned=reads_aligned_group_subgroup_median, not_aligned=reads_not_aligned_group_subgroup_median)
  colnames(plot_df_group_subgroup_median) = c("sample", "aligned", "not aligned")
  
  # melt for ggplot
  plot_df_group_median = melt(plot_df_group_median, variable.name="alignment", value.name="counts")
  plot_df_group_median$alignment = factor(plot_df_group_median$alignment, levels=c("not aligned", "aligned"))
  
  # grouping
  # order decreasing by column "counts" (only values where column alignment == aligned)
  # get the numeric order
  ord = order(plot_df_group_median[plot_df_group_median$alignment == "aligned",]$counts, decreasing = TRUE)
  # get the names
  x_order = plot_df_group_median[ord,]$sample
  
  p_orig = ggplot(plot_df_group_median, aes(x=factor(sample, levels=x_order), y=counts, fill=alignment)) + 
    geom_bar(stat="identity") +
    geom_hline(yintercept=as.numeric(number_of_reads)/1000000, linetype="dashed", color = "black", size=0.5) +
    labs(x="", y="Median processed reads\nin millions") +
    scale_fill_manual(breaks = c("aligned", "not aligned"), values=colors) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic() +
    theme(legend.position="bottom",
          legend.title=element_blank(),
          text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          axis.text = element_text(angle = 0, family = plots_props$font_family, size = plots_props$font_size),
          axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
          plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header),
          legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          # legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.margin=margin(t = -8), legend.key.size = unit(0.3, "cm"), legend.spacing.y = unit(0.4, "cm")
  )
  
  # Save results to storage
  ggsave(sprintf("%s/mapping_stats_grouping_%s_9x6.png", output_folder, grouping), p_orig, dpi=plots_props$dpi, width=plots_props$image_width, height=plots_props$image_height, units = plots_props$image_units)
  ggsave(sprintf("%s/mapping_stats_grouping_%s_9x6.svg", output_folder, grouping), p_orig, width = plots_props$image_width, height = plots_props$image_height, unit = plots_props$image_units, dpi = plots_props$dpi)
  
  ggsave(sprintf("%s/mapping_stats_grouping_%s_6x6.png", output_folder, grouping), p_orig, dpi=plots_props$dpi, width = 2/3 * plots_props$image_width, height=plots_props$image_height, units = plots_props$image_units)
  ggsave(sprintf("%s/mapping_stats_grouping_%s_6x6.svg", output_folder, grouping), p_orig, width = 2/3 * plots_props$image_width, height = plots_props$image_height, unit = plots_props$image_units, dpi = plots_props$dpi)
  
  # grouping and subgrouping
  # melt for ggplot
  plot_df_group_subgroup_median = melt(plot_df_group_subgroup_median, variable.name="alignment", value.name="counts")
  plot_df_group_subgroup_median$alignment = factor(plot_df_group_subgroup_median$alignment, levels=c("not aligned", "aligned"))
  
  # order decreasing by column "counts" (only values where column alignment == aligned)
  # get the numeric order
  ord = order(plot_df_group_subgroup_median[plot_df_group_subgroup_median$alignment == "aligned",]$counts, decreasing = TRUE)
  
  tmp = str_split_fixed(plot_df_group_subgroup_median$sample, "_", 2)
  grouping_readable = unname(unlist(xticks_names[[grouping]][tmp[,1]]))
  subgrouping_readable = unname(unlist(xticks_names[[subgrouping]][tmp[,2]]))
  plot_df_group_subgroup_median$sample_readable = sprintf("%s (%s)", grouping_readable, subgrouping_readable) 
  
  # get the names
  x_order = plot_df_group_subgroup_median[ord,]$sample_readable
  
  p_orig = ggplot(plot_df_group_subgroup_median, aes(x=factor(sample_readable, levels=x_order), y=counts, fill=alignment)) + 
    geom_bar(stat="identity") +
    geom_hline(yintercept=as.numeric(number_of_reads)/1000000, linetype="dashed", color = "black", size=0.5) +
    labs(x="", y="Median processed reads\nin millions") +
    scale_fill_manual(breaks = c("Aligned", "Not aligned"), values=colors) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic() +
    theme(legend.position="bottom",
          legend.title=element_blank(),
          text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          axis.text = element_text(angle = 0, family = plots_props$font_family, size = plots_props$font_size),
          axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
          plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header),
          legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          # legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 3),
          #legend.margin=margin(t = 20), legend.key.size = unit(0.3, "cm"), legend.spacing.y = unit(0.4, "cm")
    )
  p_orig
  # Save results to storage
  ggsave(sprintf("%s/mapping_stats_grouping_%s_%s_9x6.png", output_folder, grouping, subgrouping), p_orig, dpi=plots_props$dpi, width=plots_props$image_width, height=plots_props$image_height, units = plots_props$image_units)
  ggsave(sprintf("%s/mapping_stats_grouping_%s_%s_9x6.svg", output_folder, grouping, subgrouping), p_orig, width = plots_props$image_width, height = plots_props$image_height, unit = plots_props$image_units, dpi = plots_props$dpi)

  ggsave(sprintf("%s/mapping_stats_grouping_%s_%s_6x6.png", output_folder, grouping, subgrouping), p_orig, dpi=plots_props$dpi, width= 2/3 * plots_props$image_width, height=plots_props$image_height, units = plots_props$image_units)
  ggsave(sprintf("%s/mapping_stats_grouping_%s_%s_6x6.svg", output_folder, grouping, subgrouping), p_orig, width = 2/3 * plots_props$image_width, height = plots_props$image_height, unit = plots_props$image_units, dpi = plots_props$dpi)
}

# 
#ba10 = mapping_info[mapping_info$fastq_name == "BA10_CA1_cp_26_1" & mapping_info$Mismatches == 1,]
#ba9 = mapping_info[mapping_info$fastq_name == "BA9_CA1_cp_26_1" & mapping_info$Mismatches == 1,]
#ba9_ba10 = rbind(ba9,ba10)
#fwrite(ba9_ba10, "results/report/sample_CA1_cp_26.csv")
##print(ba9_ba10)


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])

