suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(cowplot))
suppressMessages(library(viridisLite))
suppressMessages(library(ggforce))

source("scripts/iwanthue.R")

#snakemake = readRDS("snakemake_umap.rds")

data.table::setDTthreads(1)
set.seed(snakemake@params$seed)

#---------------------------------- Functions ----------------------------------
umap_plot = function(plot_df, annot, props, data_input, plot_props, color_list, output_path_fig) {
  # this is more or less the function from the original script
  # loop over all properties
  prop_mapping = make.names(props)
  names(prop_mapping) = props
  #print(prop_mapping)

  for(a in props){
    plot_df[[prop_mapping[a]]] = as.factor(unlist(lapply(plot_df$ID, function(n) annot[annot$ID == n, a, with=F][1])))
  }
  major_prop = props[1]
  
  colors_all = unlist(color_list[[major_prop]])
  colors = colors_all[unique(plot_df[[major_prop]])]
  if(length(colors) == 1 && colors == "auto") {
    colors = unname(iwanthue(nlevels(plot_df[[prop_mapping[major_prop]]])))
  } else if(length(colors) == 1 && colors == "viridis_c"){
    # make continuous instead of discrete
    plot_df[[prop_mapping[major_prop]]] = unlist(lapply(plot_df$ID, function(n) annot[annot$ID == n, props, with=F][1]))
  } else if(length(colors) > 1 && length(colors) < length(unique(prop_mapping[major_prop]))) {
    stop(sprintf("Not enough colors provided for the property %s! Expected %d, got %d (%s)", major_prop, length(unique(prop_mapping[major_prop])), length(color), paste(names(colors), sep=', ')))
  }
  
  # legend group names in plot
  

  # set font-family for all plots
  theme_set(theme_cowplot(font_family=plot_props$font_family))

  if(is.vector(props) && length(props) == 2){
    p = ggplot(plot_df, aes_string(x="D1",y="D2", color=prop_mapping[major_prop], shape=prop_mapping[props[2]], text="ID"))
  } else {
    p = ggplot(plot_df, aes_string(x="D1",y="D2", color=prop_mapping[major_prop], text="ID"))
  }
  
  if (length(unique(plot_df[[major_prop]])) < 20) {
    p + geom_mark_hull(aes(fill = major_prop, label = major_prop), concavity = 10, con.cap = 0, label.fontsize = plot_props$font_size, label.family = plot_props$font_family, label.colour = colors)
  }

  #p = p + geom_point(size = 1) 
  # for phd thesis change point size
  p = p + geom_point(size = 1.5)
  
  if(length(colors) == 1 && colors == "viridis_c"){
    p = p + scale_colour_viridis_c(name=major_prop)
  } else if(length(colors) == 1 && colors == "viridis_d") {
    p = p + scale_colour_viridis_d(name=major_prop)
  } else {
    p = p + scale_color_manual(name=major_prop, values=colors)
  }

  if(is.vector(props) && length(props) == 2){
    p = p + #scale_shape(name=props[2]) 
      #scale_shape_manual(values=c(16, 17, 15, 3, 8, 7, 10))
      scale_shape_manual(values=c(1, 3, 0, 2, 8, 15, 17))
  }

  if(nrow(plot_df) < 25){
    p = p + geom_text_repel(aes(label = plot_df$ID), show.legend = FALSE)
  }

  #p = p + theme_nothing()

  p = p + theme_classic() +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
    xlab(sprintf("UMAP1")) +
    ylab(sprintf("UMAP2"))
  
  if(length(props) > 1) {
    p = p + theme(legend.position="bottom")
  } else {
    p = p + theme(legend.position="none")
  }
  
  p = p + theme(
                text = element_text(family = plot_props$font_family, size = plot_props$font_size),
                axis.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
                axis.title = element_text(family = plot_props$font_family, size = plot_props$font_size),
                plot.title = element_text(family = plot_props$font_family, size = plot_props$font_size_header),
                legend.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
                legend.title = element_text(family = plot_props$font_family, size = plot_props$font_size, face = "bold"),
                legend.margin=margin(t = -8), legend.key.size = unit(0.3, "cm"), legend.spacing.y = unit(0.4, "cm"),
  ) 
  #p = p + guides(col="none", shape=guide_legend(nrow=2, title="Timepoint"))
    
  needed_rows = 1
  if(is.vector(props) && length(props) == 2){
    #total_len = sum(nchar(unique(as.character(plot_df[[prop_mapping[major_prop]]])))) + sum(nchar(unique(as.character(plot_df[[prop_mapping[props[2]]]]))))
    total_len = sum(nchar(unique(as.character(plot_df[[prop_mapping[props[2]]]]))))
    if(nlevels(plot_df[[prop_mapping[major_prop]]]) + nlevels(plot_df[[prop_mapping[props[2]]]]) > 4 || total_len > 50){
      needed_rows = ceiling(total_len / 15)
      #p = p + theme(legend.box = "vertical")
      #p = p + guides(col=guide_legend(nrow=needed_rows), fill=guide_legend(nrow=needed_rows), shape=guide_legend(nrow=needed_rows))
      p = p + guides(color = "none", shape=guide_legend(nrow=needed_rows, title=xticks_names$categories[sprintf("%s_capital",props[2])]))
    }
  } else if(length(unique(as.character(plot_df[[prop_mapping[major_prop]]]))) == 1) {
    total_len = sum(nchar(unique(as.character(plot_df[[prop_mapping[major_prop]]]))), na.rm=T)
    if(nlevels(plot_df[[prop_mapping[major_prop]]]) > 4 || total_len > 50){
      
      # Tobias setting
      #needed_rows = ceiling(total_len / 50)
      
      # width=9 and height=6
      needed_rows = ceiling(total_len / 15)
      
      # width=4 and height=4
      #needed_rows = ceiling(total_len / 5)
      
      p = p + guides(col=guide_legend(nrow=needed_rows), fill=guide_legend(nrow=needed_rows), shape=guide_legend(nrow=needed_rows))
    }
  }
  
  # output file anpassen
  #save_plot(output_path_fig, p, base_aspect_ratio = 1.1, base_height=6.5*(needed_rows^0.2), dpi=600)
  # width=9 and height=6
  #ggsave(file = sprintf("%s.svg", snakemake@output$png), plot = p, width = 9, height = 6, unit = "cm", dpi = 600)
 
  ggsave(file = sprintf("%s_6x9.png", output_path_fig), plot = p, width = plot_props$image_width, height = plot_props$image_height, unit = plot_props$image_units, dpi = plot_props$dpi)
  ggsave(file = sprintf("%s_6x9.svg", output_path_fig), plot = p, width = plot_props$image_width, height = plot_props$image_height, unit = plot_props$image_units, dpi = plot_props$dpi, limitsize = FALSE, scale = 1)
  
  ggsave(file = sprintf("%s_6x6.png", output_path_fig), plot = p, width = 2/3 * plot_props$image_width, height = plot_props$image_height, unit = plot_props$image_units, dpi = plot_props$dpi)
  ggsave(file = sprintf("%s_6x6.svg", output_path_fig), plot = p, width = 2/3 * plot_props$image_width, height = plot_props$image_height, unit = plot_props$image_units, dpi = plot_props$dpi, limitsize = FALSE, scale = 1)
  
  ggsave(file = sprintf("%s_9x9.png", output_path_fig), plot = p, width = 9, height = 9, unit = plot_props$image_units, dpi = plot_props$dpi)
  ggsave(file = sprintf("%s_9x9.svg", output_path_fig), plot = p, width = 9, height = 9, unit = plot_props$image_units, dpi = plot_props$dpi, limitsize = FALSE, scale = 1)
  
  ggsave(file = sprintf("%s_3x3.png", output_path_fig), plot = p, width = 1/3 * plot_props$image_width, height = plot_props$image_height / 2, unit = plot_props$image_units, dpi = plot_props$dpi)
  ggsave(file = sprintf("%s_3x3.svg", output_path_fig), plot = p, width = 1/3 * plot_props$image_width, height = plot_props$image_height / 2, unit = plot_props$image_units, dpi = plot_props$dpi, limitsize = FALSE, scale = 1)
  
  # width=4 and height=4
  #ggsave(file = snakemake@output$png, plot = p, width = 4, height = 4, unit = "cm", dpi = 600)
}


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
props_list = snakemake@params$props
top_list = snakemake@params$top_list
n_neighbors_s = snakemake@params$n_neighbors
min_dist_s = snakemake@params$min_dist
init_s = snakemake@params$init
metric_s = snakemake@params$metric
norm_s = snakemake@params$norm
seed_s = snakemake@params$seed
#print(n_neighbors_s)
#print(min_dist_s)
#print(init_s)
#print(metric_s)
#print(norm_s)
#print(seed_s)
plot_props = snakemake@params$plots_props
xticks_names = snakemake@params$xticks_names
color_list = snakemake@params$colors
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_umap.rds")

input_folder_tab = sprintf("%s/%s_%s/results_%s/matrices/umap/coords", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 


output_folder_fig = sprintf("%s/%s_%s/results_%s/figures/umap_plots/scatter_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_fig, recursive=TRUE)


#------------------------------------ Script ----------------------------------- 
for (top_features in top_list) {
  for (n_neighbors in n_neighbors_s) {
    for (min_dist in min_dist_s) {
      for (init in init_s) {
        for (metric in metric_s) {
          for (norm in norm_s) {
            for (seed in seed_s) {
              # load the umap results
              # python code used for saving: f"embedding_neigh{n_neighbors}_mdist{min_dist}_{metric}_init{init}_seed{seed}_norm{rna_norm}.csv"
              if (top_features != "all") {
                input_file_path = sprintf("%s/umap_top_%s_expressed_%ss_neigh%s_mdist%s_%s_init%s_seed%s_norm_%s.csv", input_folder_tab, top_features, data_input$rna_class, n_neighbors, min_dist, metric, init, seed, norm)
              } else {
                input_file_path = sprintf("%s/umap_neigh%s_mdist%s_%s_init%s_seed%s_norm_%s.csv", input_folder_tab, n_neighbors, min_dist, metric, init, seed, norm)
              }
              tbl = fread(input_file_path, sep='\t')
              plot_df = as.data.frame(tbl)
              
              for (props in props_list) {
                if (length(props) == 1) {
                  # loop over all properties
                  # not so nice
                  # we just put props[1] in there as a folder name
                  # path with norm and prop + create folder
                  if (top_features != "all") {
                    output_path_fig = sprintf("%s/%s/top_%s_expressed_%ss/%s", output_folder_fig, props[1], top_features, data_input$rna_class, norm)
                  } else {
                    output_path_fig = sprintf("%s/%s/%s/%s", output_folder_fig, props[1], top_features, norm)
                  }
                  dir.create(output_path_fig, recursive=TRUE)
                } else if (length(props) == 2) {
                  if (top_features != "all") {
                    output_path_fig = sprintf("%s/%s/top_%s_expressed_%ss/%s", output_folder_fig, paste(props, collapse="_"), top_features, data_input$rna_class, norm)
                  } else {
                    output_path_fig = sprintf("%s/%s/%s/%s", output_folder_fig, paste(props, collapse="_"), top_features, norm)
                  }
                  dir.create(output_path_fig, recursive=TRUE)
                } else {
                  print(sprintf("Too many parameters provided: %s", paste(props, collapse="_")))
                }
                # make figure filename with all the stuff from the for loops
                output_file_fig = sprintf("%s/umap_neigh%s_mdist%s_%s_init%s_seed%s_norm_%s", output_path_fig, n_neighbors, min_dist, metric, init, seed, norm)
                umap_plot(plot_df, annot, props, data_input, plot_props, color_list, output_file_fig)
              }
            }
          }
        }
      }
    }
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
