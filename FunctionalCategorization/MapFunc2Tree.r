library(tidyverse)
library(ape)
library(dplyr)
library(phytools)
library(ggtree)
library(gridExtra)
library(ggplotify)
library(grid)
library(cowplot)
library(gapminder)
library(patchwork)

args = commandArgs(trailingOnly=TRUE)

setwd('/Users/katzlab/Desktop/Auden-XGT/XGT_Pipeline/Functional_Categorization/Dev/General_Scripts')

og_dirs <- list.files(path='../FunctionSearching')
#og_dirs <- c('OG5_129315')

tree_y <-  function(ggtree, data){
  if(!inherits(ggtree, "ggtree"))
    stop("not a ggtree object")
  left_join(select(data, label), select(ggtree$data, label, y)) %>%
    pull(y)
}


scale_y_tree <- function(expand=expand_scale(0, 0.6), ...){
  scale_y_continuous(expand=expand, ...)
}


tree_ylim <- function(ggtree){
  if(!inherits(ggtree, "ggtree"))
    stop("not a ggtree object")
  range(ggtree$data$y)
}


ggtreeplot <- function(ggtree, data = NULL, mapping = aes(), flip=FALSE,
                       expand_limits=expand_scale(0,.6), ...){
  
  if(!inherits(ggtree, "ggtree"))
    stop("not a ggtree object")
  
  # match the tree limits
  limits <- tree_ylim(ggtree)
  limits[1] <- limits[1] + (limits[1] * expand_limits[1]) - expand_limits[2]
  limits[2] <- limits[2] + (limits[2] * expand_limits[3]) + expand_limits[4]
  
  if(flip){
    mapping <- modifyList(aes_(x=~x), mapping)
    data <- mutate(data, x=tree_y(ggtree, data))
    gg <- ggplot(data=data, mapping = mapping, ...) +
      scale_x_continuous(limits=limits, expand=c(0,0))
  }else{
    mapping <- modifyList(aes_(y=~y), mapping)
    data <- mutate(data, y=tree_y(ggtree, data))
    gg <- ggplot(data=data, mapping = mapping, ...) +
      scale_y_continuous(limits=limits, expand=c(0,0))
  }
  gg
}

i = 0
for(og_dir in og_dirs){
  i = i + 1
  dataframes <- list.files(path=paste('../FunctionSearching/', og_dir, 
                                      '/tree_mapping_dataframes', sep=''))
  
  print(paste('../FunctionSearching/', og_dir, 
              '/', substr(og_dir, 1, 11), '.tre', sep=''))
  
  tree_plot <- ggtree(read.tree(paste('../FunctionSearching/', og_dir, 
                            '/', substr(og_dir, 1, 11), '.tre', sep='')), aes(color=substr(label, 1, 2))) +
 #   geom_treescale(x = 0, y = 100) +
    labs(color='') +
    #geom_tiplab() +
    theme(legend.position = c(.2, .8)) +
    scale_color_discrete(breaks=c('Ba', 'Za', 'Am', 'Sr', 'Op', 'Ex', 'EE', 'Pl'))
    
  for(df in dataframes){
    term_data <- data.frame(
      read.csv(
        paste('../FunctionSearching/', og_dir, '/tree_mapping_dataframes/', df, sep='')
      ))
    
    if(nrow(term_data) > 1){
    
      colnames(term_data)[1] <- 'label'
      
      term_data <- term_data %>%
        filter(label %in% tree_plot$data$label)
      
      heat_map <- ggtreeplot(tree_plot, term_data, aes(x=Term, alpha=as.numeric(Val))) +
        geom_tile(aes(fill = factor(Term))) +
        theme_classic() +
        scale_alpha_continuous(range = c(0, 1)) +
        theme(
          legend.position = "none",
          axis.line = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), 'pt'),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank()
        )
      
      if(strsplit(df, split = '.csv')[1] != 'go' & 
         strsplit(df, split = '.csv')[1] != 'combined_gos' &
         strsplit(df, split = '.csv')[1] != 'gos_from_pfam'){
        heat_map <- heat_map + theme(axis.text.x = element_text(angle=90))
      }else{
        heat_map <- heat_map + theme(axis.text.x = element_blank())
      }
      
      final_plot <- tree_plot + heat_map + 
        plot_layout(nrow = 1, ncol = 2, byrow = FALSE, widths = c(5, 1))
      
      out_handle <- paste('../FunctionSearching/', og_dir, '/mapped_figures/', 
                          strsplit(df, split = '.csv')[1], '.png', sep = '')
      
      ggsave(out_handle, final_plot)
    }
  }
}
















