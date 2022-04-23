library(tidyverse)
library(FactoMineR)
library(factoextra)
library(data.table)
library(ggrepel)
library(grid)
library(gridExtra)
library(gt)
library(dplyr)


taxon = commandArgs(trailingOnly=TRUE)
#taxon = 'Sr_st_Espi'

print(taxon)

setwd('/Users/katzlab/Desktop/Auden-XGT/XGT_Pipeline/Composition/Plastids')

gc3 <- data.frame(read.csv('gc3_master.csv')) %>%
   filter(substr(Seq, 1, 10) == taxon | substr(Seq, 3, 12) == taxon) %>%
   mutate(Type = ifelse(Type == 'egts', 'EGT', 'Conserved'))# %>%
  # filter(Type == 'Conserved')

gc3$GC3.Degen <- as.numeric(gc3$GC3.Degen)
gc3$ObsWrightENc_6Fold <- as.numeric(gc3$ObsWrightENc_6Fold)

enc_null <- data.frame(read_tsv('singletons_comp.ENc.Null.tsv'))

gc3_plot <- ggplot(gc3, aes(as.numeric(GC3.Degen), as.numeric(ObsWrightENc_6Fold))) +
  geom_point(aes(color = Type)) +
 # geom_text_repel(aes(label = label_name, color = type)) + 
  geom_line(data = enc_null, aes(GC3, ENc)) +
  theme_classic() +
  scale_color_manual(values = c('orange', 'darkgreen')) +
  labs(x = 'GC3 Degen', y = 'ObsWrightENc_6Fold') +
  theme(
    legend.position = 'none'
  )

 rscu <- data.frame(read.csv('rscu_master.csv')) %>%
   filter(substr(Seq, 1, 10) == taxon | substr(Seq, 3, 12) == taxon) %>%
   mutate(SeqID = paste(Type, '_', Seq, sep = '')) %>%
   remove_rownames() %>%
   column_to_rownames(var = 'SeqID') %>%
   select(!Seq & !Type)
 
 res.ca <- CA(rscu, ncp = 5, graph = TRUE)
 
 ca_labels.x = fviz_ca_row(res.ca, geom = "point")$labels$x
 ca_labels.y = fviz_ca_row(res.ca, geom = "point")$labels$y
 
 ca_plot_data <- data.table(fviz_ca_row(res.ca, geom = "point")$data) %>%
   mutate(reg_type = ifelse(substr(name, 1, 3) == 'egt', 'EGT', 'Conserved')) %>%
   mutate(Seq = name)
 
 rscu_plot <- ggplot(ca_plot_data, aes(x, y, color = reg_type)) +
   geom_point() +
   #geom_text_repel(aes(label = label_name)) + 
   labs(x = ca_labels.x, y = ca_labels.y, color = 'Type') +
   theme_classic() +
   scale_color_manual(values = c('orange', 'darkgreen')) +
   theme(
    # legend.position = 'none'
   )
 

#comp_plots <- grid.arrange(gc3_plot, rscu_plot, ncol = 2)
#dir.create(paste('Figures/', taxon, sep = ''))
#ggsave(filename=paste('Figures/', taxon, '/comp_plots.png', sep = ''), plot=comp_plots, width = 11, height = 6.5)

 
 
 
#------ MAHALONOBIS DISTANCE TO DETECT OUTLIERS --------#
 
mahal_thresh = .001

#----------------------- RSCU --------------------------#
 
 
ca_plot_data <- ca_plot_data %>%
   select(name | x | y)

names(ca_plot_data) <- c('Seq', 'X.Coord', 'Y.Coord')

egt_seqs <- filter(ca_plot_data, substr(Seq, 1, 3) == 'egt')$Seq

ca_mahal_data <- data.frame(matrix(nrow = 0, ncol = 2))
names(ca_mahal_data) <- c('Seq', 'mahal')

rscu_conserved_mahal <- ca_plot_data %>%
   filter(substr(Seq, 1, 3) == 'con')

rscu_conserved_mahal$mahal <- mahalanobis(rscu_conserved_mahal[, c('X.Coord', 'Y.Coord')], 
                                         colMeans(rscu_conserved_mahal[, c('X.Coord', 'Y.Coord')]), 
                                         cov(rscu_conserved_mahal[, c('X.Coord', 'Y.Coord')]))  
rscu_conserved_mahal$p_value <- pchisq(rscu_conserved_mahal$mahal, df=3, lower.tail=FALSE)

for(egt_seq in egt_seqs){
   
   curr_ca_data <- filter(ca_plot_data, substr(Seq, 1, 3) == 'con' | Seq == egt_seq)
   curr_ca_data$mahal <- mahalanobis(curr_ca_data[, c('X.Coord', 'Y.Coord')], 
                               colMeans(curr_ca_data[, c('X.Coord', 'Y.Coord')]), 
                               cov(curr_ca_data[, c('X.Coord', 'Y.Coord')]))  
   curr_ca_data$p_value <- pchisq(curr_ca_data$mahal, df=3, lower.tail=FALSE)
   ca_mahal_data <- rbind(
      ca_mahal_data, 
      filter(curr_ca_data, Seq == egt_seq)
   )
}

ca_mahal_data <- rbind(ca_mahal_data, rscu_conserved_mahal)

ca_mahal_data <- ca_mahal_data %>%
   mutate(outlying = ifelse(p_value < mahal_thresh & substr(Seq, 1, 3) == 'egt', 'Contamination', 
                            ifelse(p_value >= mahal_thresh & substr(Seq, 1, 3) == 'egt', 'LTGs',
                            'Conserved Genes')))

rscu_mahal_plot <- ggplot(ca_mahal_data, aes(X.Coord, Y.Coord, color = outlying)) +
   geom_point() +
   #geom_text_repel(aes(label = label_name)) + 
   labs(x = ca_labels.x, y = ca_labels.y, color = '') +
   theme_classic() +
   scale_color_manual(values = c('orange', 'red', 'darkblue')) +
   theme(
      # legend.position = 'none'
   )

#----------------------- GC3-Degen --------------------------#

gc3_data <- gc3 %>%
   select(Seq | Type | GC3.Degen | ObsWrightENc_6Fold)

gc3_data$GC3.Degen <- as.numeric(gc3_data$GC3.Degen)
gc3_data$ObsWrightENc_6Fold <- as.numeric(gc3_data$ObsWrightENc_6Fold)

egt_seqs <- filter(gc3_data, Type == 'EGT')$Seq

gc3_mahal_data <- data.frame(matrix(nrow = 0, ncol = 2))
names(gc3_mahal_data) <- c('Seq', 'mahal')

gc3_conserved_mahal <- gc3_data %>%
   filter(Type == 'Conserved')

temp <- gc3_conserved_mahal[, c('GC3.Degen', 'ObsWrightENc_6Fold')]

gc3_conserved_mahal$mahal <- mahalanobis(gc3_conserved_mahal[, c('GC3.Degen', 'ObsWrightENc_6Fold')], 
                                   colMeans(gc3_conserved_mahal[, c('GC3.Degen', 'ObsWrightENc_6Fold')]), 
                                   cov(gc3_conserved_mahal[, c('GC3.Degen', 'ObsWrightENc_6Fold')]))  
gc3_conserved_mahal$p_value <- pchisq(gc3_conserved_mahal$mahal, df=3, lower.tail=FALSE)

for(egt_seq in egt_seqs){
   
   curr_gc3_data <- filter(gc3_data, Type == 'Conserved' | Seq == egt_seq)
   curr_gc3_data$mahal <- mahalanobis(curr_gc3_data[, c('GC3.Degen', 'ObsWrightENc_6Fold')], 
                                     colMeans(curr_gc3_data[, c('GC3.Degen', 'ObsWrightENc_6Fold')]), 
                                     cov(curr_gc3_data[, c('GC3.Degen', 'ObsWrightENc_6Fold')]))  
   curr_gc3_data$p_value <- pchisq(curr_gc3_data$mahal, df=3, lower.tail=FALSE)
   gc3_mahal_data <- rbind(
      gc3_mahal_data, 
      filter(curr_gc3_data, Seq == egt_seq)
   )

}

gc3_mahal_data <- rbind(gc3_mahal_data, gc3_conserved_mahal) %>%
   mutate(outlying = ifelse(p_value < mahal_thresh & Type == 'EGT', 'Contamination', 
                            ifelse(p_value >= mahal_thresh & Type == 'EGT', 'LTGs',
                                   'Conserved Genes')))

gc3_conserved_violin <- ggplot(gc3_conserved_mahal, aes(x = 'Conserved', y = as.numeric(GC3.Degen))) +
   geom_violin(fill = "lightblue", alpha = .5) +
   theme_classic() +
   scale_y_continuous(limits = c(0, 100)) +
   theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.line = element_line(linetype = 'dashed')
   )

gc3_mahal_plot <- ggplot(gc3_mahal_data, aes(GC3.Degen, ObsWrightENc_6Fold)) +
   geom_point(aes(color = outlying)) +
   geom_line(data = enc_null, aes(GC3, ENc)) +
   theme_classic() +
   scale_color_manual(values = c('orange', 'red', 'darkblue')) +
   labs(x = 'GC3 Degen', y = 'ObsWrightENc_6Fold') +
   theme(
      legend.position = 'none'
   ) + 
   annotation_custom(ggplotGrob(gc3_conserved_violin), 
                       xmin = 70, 
                       ymin = 53)

#----------------------- Synthesis --------------------------#


mahal_plots <- grid.arrange(gc3_mahal_plot, rscu_mahal_plot, ncol = 2)


ggsave(filename=paste('Figures/', taxon, '/mahal.png', sep = ''), 
       plot=mahal_plots, width = 11, height = 6.5)

#------------------- Recording Outliers ---------------------#

rscu_outliers <- filter(ca_mahal_data, outlying == 'Plastid')$Seq
gc3_outliers <- filter(gc3_mahal_data, outlying == 'Plastid')$Seq

outlier_df <- data.frame(read.csv('outliers.csv'))
names(outlier_df) <- c('Seq')

for(outlier in rscu_outliers){
   outlier_df <- rbind(outlier_df, outlier)
}

for(outlier in gc3_outliers){
   if(paste('egts_', outlier, sep = '') %in% rscu_outliers){}
   else{
    outlier_df <- rbind(outlier_df, outlier)
   }
}

write.csv(outlier_df, 'outliers.csv', row.names = FALSE)




