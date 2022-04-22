library(tidyverse)

setwd('')

`%notin%` <- Negate(`%in%`)

selected_terms <- data.frame(read.csv('ogs_to_use.txt'))

all_data <- data.frame(read.csv('go_by_og.csv')) %>%
  filter(Gene.Family %in% selected_terms$OG)

go_terms_fung_sum <- all_data %>% 
  group_by(GO) %>%
  summarise(fung_sum = sum(N.Fungi))

go_terms_met_sum <- all_data %>%
  group_by(GO) %>%
  summarise(met_sum = sum(N.Met))

all_data <- all_data %>%
  merge(go_terms_fung_sum, by = 'GO') %>%
  merge(go_terms_met_sum, by = 'GO')

go_term_sorting <- all_data %>% 
  group_by(GO) %>%
  summarise(go_sort = sum(ifelse(fung_sum == 0 & met_sum > 0 & N.Met > 0, -1,
                                 ifelse(met_sum == 0 & fung_sum > 0 & N.Fungi > 0, 1,
                                        0))))

all_data <- all_data %>%
  merge(go_term_sorting, by = 'GO') %>%
  filter(GO != '-') %>%
  mutate(clr = ifelse(N.Fungi == 0, 'darkblue', ifelse(N.Met == 0, 'orange', 'purple'))) %>%
  mutate(a = ifelse(N.Fungi + N.Met == 0, 0, 
                    ifelse(N.Fungi + N.Met < 10, 1, 
                           ifelse(N.Fungi + N.Met < 50, 1, 1))))
  #mutate(a = ifelse(N.Met == 0 & N.Fungi == 0, 1, 0))

heat_map <- ggplot(all_data, aes(fct_rev(reorder(GO, as.numeric(go_sort))), 
                                 reorder(Gene.Family, as.numeric(N.Fungi)), alpha = as.factor(a))) +
  geom_tile(aes(fill = clr)) +
  scale_fill_manual(values = c('blue', 'red', 'purple')) +
  theme_classic() +
  scale_alpha_manual(values = c(0, 1, 1, 1)) +
  labs(x = 'GO Term', y = 'Gene Family') +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    # axis.line = element_blank(),
    axis.line = element_line(color = 'darkgrey', size = .2),
    legend.position = 'none',
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10))
  )

fung_counts <- all_data %>%
  group_by(GO) %>%
  summarise(fun_trees = sum(ifelse(N.Fungi > 0, 1, 0)))

tree_type_counts <- all_data %>%
  group_by(GO) %>%
  summarise(met_trees = sum(ifelse(N.Met > 0, 1, 0))) %>%
  merge(fung_counts, by = 'GO') %>%
  merge(go_term_sorting, by = 'GO')

tree_type_curves <- ggplot(tree_type_counts) +
  geom_line(aes(x = reorder(GO, go_sort), y = fun_trees, group = 1), color = "blue") +
  geom_line(aes(x = GO, y = met_trees, group = 1), color = "red") +
  theme_classic()


tree_type_curves <- ggplot(tree_type_counts) +
  geom_density(aes(x = fun_trees), fill = 'red', alpha = .25) +
  geom_density(aes(x = met_trees), fill = 'blue', alpha = .25) +
  labs(x = "Number of Trees (Blue = Metazoa, Red = Fungi)", y = 'Density (GO Terms Shared in Fungi or Metazoa Trees)') +
  coord_cartesian(expand = c(0, 0)) +
  theme_classic()
  
heat_map

write.csv(tree_type_counts, 'tree_type_counts.csv')















