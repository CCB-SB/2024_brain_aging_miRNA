suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
#suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(stringr))


#snakemake = readRDS("snakemake_overlap_sig_correlated_paper_bar_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("overlap_sig_correlated_paper_bar_plot")


#---------------------------------- Functions ----------------------------------
modify_facet_appearance <- function(plot = NULL, strip.background.x.fill = NULL, strip.background.y.fill = NULL, strip.background.x.col = NULL,
                                    strip.background.y.col = NULL, strip.text.x.col = NULL, strip.text.y.col = NULL){
  
  if(is.null(plot)){stop("A ggplot (gg class) needs to be provided!")}
  
  # Generate gtable object to modify the facet strips:
  g <- ggplot_gtable(ggplot_build(plot))
  
  # Get the locations of the right and top facets in g:
  stripy <- which(grepl('strip-r|strip-l', g$layout$name)) # account for when strip positions are switched r-l and/or t-b in facet_grid(switch = )
  stripx <- which(grepl('strip-t|strip-b', g$layout$name))
  
  # Check that the provided value arrays have the same length as strips the plot has:
  lx <- c(length(strip.background.x.fill), length(strip.background.x.col), length(strip.text.x.col))
  if(!all(lx==length(stripx) | lx==0)){stop("x The provided vectors with values need to have the same length and the number of facets in the plot!")}
  ly <- c(length(strip.background.y.fill), length(strip.background.y.col), length(strip.text.y.col))
  if(!all(ly==length(stripy) | ly==0)){stop("y The provided vectors with values need to have the same length and the number of facets in the plot!")}
  
  # Change the strips on the y axis:
  for (i in seq_along(stripy)){ # if no strips in the right, the loop will not be executed as seq_along(stripy) will be integer(0)
    
    # Change strip fill and (border) colour :
    j1 <- which(grepl('strip.background.y', g$grobs[[stripy[i]]]$grobs[[1]]$childrenOrder))
    if(!is.null(strip.background.y.fill[i])){g$grobs[[stripy[i]]]$grobs[[1]]$children[[j1]]$gp$fill <- strip.background.y.fill[i]} # fill
    if(!is.null(strip.background.y.col[i])){g$grobs[[stripy[i]]]$grobs[[1]]$children[[j1]]$gp$col <- strip.background.y.col[i]} # border colour
    
    # Change color of text:
    j2 <- which(grepl('strip.text.y', g$grobs[[stripy[i]]]$grobs[[1]]$childrenOrder))
    if(!is.null(strip.text.y.col[i])){g$grobs[[stripy[i]]]$grobs[[1]]$children[[j2]]$children[[1]]$gp$col <- strip.text.y.col[i]}
    
  }
  
  # Same but for the x axis:
  for (i in seq_along(stripx)){
    
    # Change strip fill and (border) colour :
    j1 <- which(grepl('strip.background.x', g$grobs[[stripx[i]]]$grobs[[1]]$childrenOrder))
    if(!is.null(strip.background.x.fill[i])){g$grobs[[stripx[i]]]$grobs[[1]]$children[[j1]]$gp$fill <- strip.background.x.fill[i]} # fill
    if(!is.null(strip.background.x.col[i])){g$grobs[[stripx[i]]]$grobs[[1]]$children[[j1]]$gp$col <- strip.background.x.col[i]} # border colour
    
    # Change color of text:
    j2 <- which(grepl('strip.text.x', g$grobs[[stripx[i]]]$grobs[[1]]$childrenOrder))
    if(!is.null(strip.text.x.col[i])){g$grobs[[stripx[i]]]$grobs[[1]]$children[[j2]]$children[[1]]$gp$col <- strip.text.x.col[i]}
    
  }
  
  return(g)
}
# Note that it returns a gtable object. This can be ploted with plot() or grid::draw.grid().
# patchwork can handle the addition of such gtable to a layout with other ggplot objects,
# but be sure to use patchwork::wrap_ggplot_grob(g) for proper alignment of plots!
# See: https://patchwork.data-imaginist.com/reference/wrap_ggplot_grob.html


#---------------------------------- Load data ---------------------------------- 
data_input = snakemake@params$data_input
#tissues = snakemake@params$tissues 
plot_props = snakemake@params$plots_props
raw_data_path = snakemake@params$raw_data_path
colors = snakemake@params$colors
xticks_names = snakemake@params$xticks_names
corr_list = snakemake@params$thresholds$corr_th
sig_lvl = snakemake@params$thresholds$adj_p_value
method = snakemake@params$method
plot_size = snakemake@params$plot_size
props = snakemake@params$props
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_overlap_sig_correlated_paper_bar_plot.rds")
#stop()

# input folder
input_folder_tab = sprintf("%s/%s_%s/results_%s/matrices/aggregation_corr", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)

# output folder
output_folder_fig = sprintf("%s/%s_%s/results_%s/figures/corr_plots/overlap_bar_plots",results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_fig, recursive=TRUE)
output_folder_tab = sprintf("%s/%s_%s/results_%s/matrices/corr_plots/overlap_bar_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_tab, recursive=TRUE)

annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 

external_miRNA_list = fread(raw_data_path, sep='\t', header=FALSE) 


#------------------------------------ Script ----------------------------------- 
for (m in method) {
  for (prop in props) {
    
    tissue = prop[1]
    time = prop[2]
     
    # load data
    corr = fread(sprintf("%s/corr_method=%s_%ss_with_%s_per_%s.csv", input_folder_tab, m, data_input$rna_class, time, tissue), sep='\t')
    corr$V1 = c()
    pvalue = fread(sprintf("%s/corr_method=%s_adj_pvalues_%ss_with_%s_per_%s.csv", input_folder_tab, m, data_input$rna_class, time, tissue), sep='\t')
    pvalue$V1 = c()  
    
    for (th in corr_list){
      #stop()
      pos_corr_rnas = c()
      neg_corr_rnas = c()
      sig_pos_corr_rnas = c()
      sig_neg_corr_rnas = c()
      overlap_pos_corr_rnas = c()
      overlap_neg_corr_rnas = c()
      overlap_sig_pos_corr_rnas = c()
      overlap_sig_neg_corr_rnas = c()
      overlap_pos_corr_rnas_number = c()
      overlap_neg_corr_rnas_number = c()
      overlap_sig_pos_corr_rnas_number = c()
      overlap_sig_neg_corr_rnas_number = c()
      #info_table_tissue_tmp = c()
      for (ti in unique(annot[[tissue]])) {
        info_table_tissue = c()
        for (miR in corr[[data_input$rna_class]]) {
          info_table = data.table(feature=miR,
                                  tissue=ti,
                                  ttest_adj_pval=pvalue[pvalue[[data_input$rna_class]] == miR,][[ti]], 
                                  corr=corr[corr[[data_input$rna_class]] == miR,][[ti]]
          )
          
          info_table$corr_cat_pos = rep(FALSE, dim(info_table)[1])
          info_table$corr_cat_neg = rep(FALSE, dim(info_table)[1])
          info_table$ttest_adj_cat_corr_pos = rep(FALSE, dim(info_table)[1])
          info_table$ttest_adj_cat_corr_neg = rep(FALSE, dim(info_table)[1])
          
          #info_table[ttest_adj_pval >= sig_lvl, ttest_adj_cat:="not_sig"]
          #info_table[ttest_adj_pval < sig_lvl, ttest_adj_cat:="sig"]
          info_table[corr >= th, corr_cat_pos:=TRUE]
          info_table[corr <= -th, corr_cat_neg:=TRUE]
          info_table[corr >= th & ttest_adj_pval < sig_lvl, ttest_adj_cat_corr_pos:=TRUE]
          info_table[corr <= -th & ttest_adj_pval < sig_lvl, ttest_adj_cat_corr_neg:=TRUE]
          
          info_table_tissue = rbind(info_table_tissue, info_table)
        }
        
        #info_table_tissue_tmp = rbind(info_table_tissue_tmp, info_table_tissue[info_table_tissue$corr == max(info_table_tissue$corr),])
        
        pos_corr = data.frame("feature" = info_table_tissue[info_table_tissue$corr_cat_pos,]$feature)
        neg_corr = data.frame("feature" = info_table_tissue[info_table_tissue$corr_cat_neg,]$feature)
        sig_pos_corr = data.frame("feature" = info_table_tissue[info_table_tissue$ttest_adj_cat_corr_pos,]$feature)
        sig_neg_corr = data.frame("feature" = info_table_tissue[info_table_tissue$ttest_adj_cat_corr_neg,]$feature)
        
        pos_corr_contain = data.table("feature" = unique(unlist(pos_corr, use.names = FALSE)))
        neg_corr_contain = data.table("feature" = unique(unlist(neg_corr, use.names = FALSE)))
        sig_pos_corr_contain = data.table("feature" = unique(unlist(sig_pos_corr, use.names = FALSE)))
        sig_neg_corr_contain = data.table("feature" = unique(unlist(sig_neg_corr, use.names = FALSE)))
        
        pos_corr_rnas[[ti]] = pos_corr_contain$feature
        neg_corr_rnas[[ti]] = neg_corr_contain$feature
        sig_pos_corr_rnas[[ti]] = sig_pos_corr_contain$feature
        sig_neg_corr_rnas[[ti]] = sig_neg_corr_contain$feature
        
        overlap_pos_corr_rnas[[ti]] = intersect(pos_corr_rnas[[ti]], external_miRNA_list$V1)
        overlap_neg_corr_rnas[[ti]] = intersect(neg_corr_rnas[[ti]], external_miRNA_list$V1)
        overlap_sig_pos_corr_rnas[[ti]] = intersect(sig_pos_corr_rnas[[ti]], external_miRNA_list$V1)
        overlap_sig_neg_corr_rnas[[ti]] = intersect(sig_neg_corr_rnas[[ti]], external_miRNA_list$V1)
        
        overlap_pos_corr_rnas_number[[ti]] = length(overlap_pos_corr_rnas[[ti]])
        overlap_neg_corr_rnas_number[[ti]] = length(overlap_neg_corr_rnas[[ti]])
        overlap_sig_pos_corr_rnas_number[[ti]] = length(overlap_sig_pos_corr_rnas[[ti]])
        overlap_sig_neg_corr_rnas_number[[ti]] = length(overlap_sig_neg_corr_rnas[[ti]])
      }
      
      overlap_pos_corr_rnas_number_df = data.frame(melt(overlap_pos_corr_rnas_number))
      colnames(overlap_pos_corr_rnas_number_df) = c("pos", tissue)
      overlap_neg_corr_rnas_number_df = data.frame(melt(overlap_neg_corr_rnas_number))
      colnames(overlap_neg_corr_rnas_number_df) = c("neg", tissue)
      overlap_sig_pos_corr_rnas_number_df = data.frame(melt(overlap_sig_pos_corr_rnas_number))
      colnames(overlap_sig_pos_corr_rnas_number_df) = c("sig_pos", tissue)
      overlap_sig_neg_corr_rnas_number_df = data.frame(melt(overlap_sig_neg_corr_rnas_number))
      colnames(overlap_sig_neg_corr_rnas_number_df) = c("sig_neg", tissue)
  
      merged_df = Reduce(function(x, y) merge(x, y, by = tissue, all = TRUE), list(overlap_pos_corr_rnas_number_df, overlap_neg_corr_rnas_number_df))
      merged_df_sig = Reduce(function(x, y) merge(x, y, by = tissue, all = TRUE), list(overlap_sig_pos_corr_rnas_number_df, overlap_sig_neg_corr_rnas_number_df))
      
      # Convert data to long format for stacking
      df_long = reshape2::melt(merged_df, id.vars = tissue, variable.name = "direction", value.name = "value")
      df_long_sig = reshape2::melt(merged_df_sig, id.vars = tissue, variable.name = "direction", value.name = "value")
      
      #df_long$brain_region = factor(df_long$brain_region, levels = rev(names(colors[[tissue]])))
      #df_long_sig$brain_region = factor(df_long_sig$brain_region, levels = rev(names(colors[[tissue]])))
      
      df_long$plot_names_split_group = unname(xticks_names[[tissue]])
      df_long$plot_names_split_group = factor(df_long$plot_names_split_group, levels=unname(xticks_names[[tissue]]))
      
      df_long_sig$plot_names_split_group = unname(xticks_names[[tissue]])
      df_long_sig$plot_names_split_group = factor(df_long_sig$plot_names_split_group, levels=unname(xticks_names[[tissue]]))
      
      # Create stacked bar plot
      # bar_plot_corr = ggplot(df_long, aes(x = brain_region, y = Count, fill = Correlation)) +
      #   geom_col() +  # Stacked bar plot
      #   scale_fill_manual(values = unlist(colors$direction)) +  # Set colors
      #   coord_flip() +  # Flip the plot
      #   facet_grid(rows = vars(plot_names_split_group), scales = "free_y", space = "free", switch = "y") +  # Facets on y-axis
      #   labs(x = "", y = "Number of overlaps") +
      #   theme_classic() +
      #   theme(legend.position="none",
      #         legend.title=element_blank(),
      #         text = element_text(family = plot_props$font_family, size = plot_props$font_size),
      #         axis.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
      #         axis.title = element_text(family = plot_props$font_family, size = plot_props$font_size),
      #         plot.title = element_text(family = plot_props$font_family, size = plot_props$font_size_header),
      #         legend.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
      #         strip.text.y.left = element_text(family = plot_props$font_family, angle = 0, size=plot_props$font_size),
      #         # make x invisible
      #         #axis.text.x=element_blank(),
      #         panel.spacing = unit(2, "lines"),  # spacing between groups
      #         strip.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
      #   )
      
      bar_plot_corr = ggplot(df_long, aes(x=factor(brain_region, levels=unname(xticks_names[[tissue]])),
                                              y=value, fill=direction)) + 
        geom_col() +
        facet_grid(rows=vars(plot_names_split_group), scales = "free_y", space = "free", switch = "y") +
        ##theme_cowplot(11) + 
        scale_y_continuous(labels = scales::percent_format(accuracy = 1, scale=1, suffix=""), expand = expansion(mult=c(0, 0.05))) +
        coord_flip() + 
        xlab("") +
        ylab("Number of overlaps") +
        #scale_fill_manual(values=colours_rna_class_filtered,
        #                  labels=function(x) str_trim(gsub("_", " ", gsub("ensembl", "", x))), 
        #                  name="Class",
        #                  guide=guide_legend(ncol=ncol_legend, title.position = "top")
        #) +
        scale_fill_manual(values = unlist(colors$direction)) +  # Set colors
        theme_classic() +
        theme(legend.position = "none", 
              text = element_text(family = plot_props$font_family, size = plot_props$font_size),
              plot.title = element_text(family = plot_props$font_family, size = plot_props$font_size_header),
              legend.title = element_blank(), #element_text(family = plot_props$font_family, size = plot_props$font_size),
              legend.text = element_text(family = plot_props$font_family, size = plot_props$font_size, margin = margin(0, 0, 0, 0, "pt")),
              axis.text = element_text(angle = 0, family = plot_props$font_family, size = plot_props$font_size),
              axis.title = element_text(family = plot_props$font_family, size = plot_props$font_size),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(), 
              axis.line.y = element_blank(),
              strip.text.y.left = element_text(family = plot_props$font_family, angle = 0, size=plot_props$font_size),
              #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4),
              panel.spacing = unit(0, "pt"),
              plot.margin = margin(5.5, 10, 5.5, 0, "pt"),
              legend.margin = margin(-15, 0, 0, 0, "pt"), legend.key.size = unit(0.3, "cm"), legend.spacing.y = unit(0.4, "cm"), legend.spacing.x = unit(0, "cm"),
        )
      #bar_plot_corr
      
      bar_plot_corr_mod = modify_facet_appearance(bar_plot_corr,
                                                  strip.background.y.fill = unlist(colors[[tissue]]),
                                                  strip.background.y.col = unlist(colors[[tissue]]), #rep("white", length(cellline_colors)),
                                                  strip.text.y.col = unlist(colors[[sprintf("%s_text", tissue)]]) #rep("black", length(cellline_colors))
      )
      #plot(bar_plot_corr_mod)
      
      width = plot_size$image_width
      height = plot_size$image_height
      ggsave(sprintf("%s/corr_method=%s_corr_th=%s_grouped_by_%s_%sx%s.png", output_folder_fig, m, th, paste(prop, collapse = "_"), height, width), bar_plot_corr_mod, dpi=plot_props$dpi, width=width, height=height, units = plot_props$image_units)
      ggsave(sprintf("%s/corr_method=%s_corr_th=%s_grouped_by_%s_%sx%s.svg", output_folder_fig, m, th, paste(prop, collapse = "_"), height, width), bar_plot_corr_mod, width = width, height = height, unit = plot_props$image_units, dpi = plot_props$dpi)
      
      # sig
      tmp = str_split_fixed(df_long_sig$direction, "_", 2)
      df_long_sig$direction_colour = tmp[,2]
        
      bar_plot_corr_sig = ggplot(df_long_sig, aes(x=factor(brain_region, levels=unname(xticks_names[[tissue]])),
                                          y=value, fill=direction_colour)) + 
        geom_col() +
        facet_grid(rows=vars(plot_names_split_group), scales = "free_y", space = "free", switch = "y") +
        ##theme_cowplot(11) + 
        scale_y_continuous(labels = scales::percent_format(accuracy = 1, scale=1, suffix=""), expand = expansion(mult=c(0, 0.05))) +
        coord_flip() + 
        xlab("") +
        ylab("Number of overlaps") +
        #scale_fill_manual(values=colours_rna_class_filtered,
        #                  labels=function(x) str_trim(gsub("_", " ", gsub("ensembl", "", x))), 
        #                  name="Class",
        #                  guide=guide_legend(ncol=ncol_legend, title.position = "top")
        #) +
        scale_fill_manual(values = unlist(colors$direction)) +  # Set colors
        theme_classic() +
        theme(legend.position = "none", 
              text = element_text(family = plot_props$font_family, size = plot_props$font_size),
              plot.title = element_text(family = plot_props$font_family, size = plot_props$font_size_header),
              legend.title = element_blank(), #element_text(family = plot_props$font_family, size = plot_props$font_size),
              legend.text = element_text(family = plot_props$font_family, size = plot_props$font_size, margin = margin(0, 0, 0, 0, "pt")),
              axis.text = element_text(angle = 0, family = plot_props$font_family, size = plot_props$font_size),
              axis.title = element_text(family = plot_props$font_family, size = plot_props$font_size),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(), 
              axis.line.y = element_blank(),
              strip.text.y.left = element_text(family = plot_props$font_family, angle = 0, size=plot_props$font_size),
              #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4),
              panel.spacing = unit(0, "pt"),
              plot.margin = margin(5.5, 10, 5.5, 0, "pt"),
              legend.margin = margin(-15, 0, 0, 0, "pt"), legend.key.size = unit(0.3, "cm"), legend.spacing.y = unit(0.4, "cm"), legend.spacing.x = unit(0, "cm"),
        )
      #bar_plot_corr_sig
      
      bar_plot_corr_sig_mod = modify_facet_appearance(bar_plot_corr_sig,
                                                  strip.background.y.fill = unlist(colors[[tissue]]),
                                                  strip.background.y.col = unlist(colors[[tissue]]), #rep("white", length(cellline_colors)),
                                                  strip.text.y.col = unlist(colors[[sprintf("%s_text", tissue)]]) #rep("black", length(cellline_colors))
      )
      #plot(bar_plot_corr_sig_mod)
      
      width = plot_size$image_width
      height = plot_size$image_height
      ggsave(sprintf("%s/corr_method=%s_sig_corr_th=%s_grouped_by_%s_%sx%s.png", output_folder_fig, m, th, paste(prop, collapse = "_"), height, width), bar_plot_corr_sig_mod, dpi=plot_props$dpi, width=width, height=height, units = plot_props$image_units)
      ggsave(sprintf("%s/corr_method=%s_sig_corr_th=%s_grouped_by_%s_%sx%s.svg", output_folder_fig, m, th, paste(prop, collapse = "_"), height, width), bar_plot_corr_sig_mod, width = width, height = height, unit = plot_props$image_units, dpi = plot_props$dpi)
      
      
      # save tables
      overlap_pos_corr_rnas_df = melt(overlap_pos_corr_rnas)
      colnames(overlap_pos_corr_rnas_df) = c("feature", tissue)
      overlap_pos_corr_rnas_df$direction = "pos"
      overlap_neg_corr_rnas_df = melt(overlap_neg_corr_rnas)
      colnames(overlap_neg_corr_rnas_df) = c("feature", tissue)
      overlap_neg_corr_rnas_df$direction = "neg"
      overlap_sig_pos_corr_rnas_df = melt(overlap_sig_pos_corr_rnas)
      colnames(overlap_sig_pos_corr_rnas_df) = c("feature", tissue)
      overlap_sig_pos_corr_rnas_df$direction = "sig_pos"
      overlap_sig_neg_corr_rnas_df = melt(overlap_sig_neg_corr_rnas)
      colnames(overlap_sig_neg_corr_rnas_df) = c("feature", tissue)
      overlap_sig_neg_corr_rnas_df$direction = "sig_neg"
      
      info = rbind(overlap_pos_corr_rnas_df, overlap_neg_corr_rnas_df)
      info_sig = rbind(overlap_sig_pos_corr_rnas_df, overlap_sig_neg_corr_rnas_df)
      tbl = rbind(info, info_sig)
      
      file_name = sprintf("corr_method=%s_corr_th=%s", m, th)
      #if (class(pos_rnas_global_table) == "data.frame") {
      fwrite(tbl, sprintf("%s/%s.csv", output_folder_tab, file_name), sep = "\t", row.names = TRUE)
      write.xlsx(tbl, sprintf("%s/%s.xlsx", output_folder_tab, file_name), colNames = TRUE, rowNames = TRUE, append = FALSE)
      #} else {
      #  print(sprintf("%s is empty", file_name))
      #}
    }
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
