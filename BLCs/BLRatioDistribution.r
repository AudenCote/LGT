library(tidyverse)
library(scales)

setwd('')

'%notin%' <- Negate('%in%')

bl_data <- data.frame(read.csv('branch_length_ratios_all_anaerobes.csv')) %>%
  mutate(LGT = ifelse(Type == 'LGT', paste(substr(OG, 1, 10), '_LTG.tre', sep=''), LGT)) %>%
  mutate(LGT = paste(substr(LGT, 1, 10), '_', Cat, sep = '')) %>%
  filter(Ratio < 1000)# %>%
  #filter(LGT == 'OG5_135734')

recips <- data.frame(read.csv('recipients.csv'))

bl_data <- bl_data %>%
  merge(recips, by = 'LGT')

meds <- bl_data %>%
  group_by(LGT) %>%
  summarize(med = median(Ratio))

bl_data <- bl_data %>%
  merge(meds, by = "LGT")

vtgs <- bl_data %>%
  filter(Type == 'Conserved')

ltg <- bl_data %>%
  filter(Type == 'LGT') %>%
  mutate(lgt_rat = Ratio)

rel_col <- ltg %>%
  select(LGT | lgt_rat)

vtgs <- vtgs %>%
  merge(rel_col, by = 'LGT')

vtgs$Recip <- factor(vtgs$Recip, levels = c('Am', 'Sr', 'Ex', 'Intra'))
ltg$Recip <- factor(ltg$Recip, levels = c('Am', 'Sr', 'Ex', 'Intra'))

dist_plot <- ggplot() +
  geom_boxplot(vtgs, mapping = aes(x = reorder(LGT, lgt_rat - med), y = as.numeric(Ratio), fill = Recip)) +
  geom_point(ltg, mapping = aes(x = LGT, y = Ratio), color = '#EE4B2B', size = 2) +
  # coord_cartesian(xlim = c(0, 30), ylim = c(0, .1), expand = c(0, 0)) +
  coord_flip() +
  #scale_y_log10() +
  #scale_y_continuous(trans = 'ln') +
  theme_classic() +
  scale_fill_manual(values = c('#c1d0e8', '#f59597', '#f3ec96', 'white')) +
  labs(y = 'Branch Length Ratio', fill = "") +
  theme(
   # axis.text.y = element_blank(),
    axis.text.y = element_text(size = 7),
    #axis.text.x = element_text(angle = 90)
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  ) +
  facet_grid(rows = vars(Recip), space = "free", scales = 'free')
 

dist_plot
