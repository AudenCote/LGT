#library(plyr)
library(tidyverse)
setwd('')


#Reading in the dataframe
presence_data <- read.csv("../Presence_Absence_DataFrames/minor_clade_presence.csv")
presence_data <- data.frame(presence_data)

#The six eukaryotic major clades: Amoebozoa, SAR, Opisthokonta, Excavata, Archaeplastida, and Everything Else
euk_major_clades <- c('am', 'sr', 'op', 'ex', 'pl', 'ee')

#Renaming the first column
names(presence_data)[1] <- "OG"

#Initializing the "counts" dataframe. This will hold the clade-pairing data for our chart
counts <- data.frame(presence_data$OG)
counts <- cbind(counts, host = 0, donor = 0, euks = 0, bacs = 0, am_op = 0)

#For every column in the master minor clade dataframe
for(var in names(presence_data)){
  print(var)
  #If the column is a a target host (e.g. Archaeplastida) minor clade and not bacterial bin data, add its values to the "host" column
  if(substr(var, 0, 7) == "sr_ci_s"){#&& (substr(var, 5, 5) == 'b' | substr(var, 5, 5) == 'u')){
    counts$host <- counts$host + presence_data[, var]
  }
  #If the column is a target donor (e.g. Cyanobacteria) minor clade, add its values to the "donor" column
  if(substr(var, 0, 5) == "ba_pg"){
    counts$donor <- counts$donor + presence_data[, var]
  }
  #If this column is a bacterial minor clade, add its values to the "bacs" column
  if(substr(var, 0, 2) == "ba" | substr(var, 0, 2) == "za"){
    counts$bacs <- counts$bacs + presence_data[, var]
  }
  #If the column is a eukaryotic clade (and not bacterial bin data), add its values to the "euks" column
  if(substr(var, 0, 2) %in% euk_major_clades && substr(var, 5, 5) != 'b'){#  && substr(var, 0, 2) != 'am'){
    # && !(substr(var, 0, 2) == 'am' && substr(var, 0, 5) != 'am_tu')){
    counts$euks <- counts$euks + presence_data[, var]
  }
  #If the column is an amoebozoan or opisthokont and not bacterial bin data, add its values to the "am_op" column
  if((substr(var, 0, 2) == "am" | substr(var, 0, 2) == "op") && substr(var, 5, 5) != 'b'){
    counts$am_op <- counts$am_op + presence_data[, var]
  }
}

#Filtering the "counts" dataframe -- these parameters can be tweaked and are pretty self-explanatory
counts <- counts %>%
  mutate(host_ratio = host / euks) %>%
  mutate(donor_ratio = donor / bacs) %>%
  mutate(am_op_ratio = am_op / euks) %>%
  filter(host >= 5) %>%
  filter(donor >= 5) #%>%
 # filter(hos/N >= .15)
#filter(am_op < 5) %>%
#filter(am_op_ratio < .125)

#Getting rid of any NA values
counts[is.na(counts)] <- 0

#Building our output ggplot
presence_plot <- ggplot(counts, aes(x = bacs, y = host_ratio, color = host/55)) + 
  #Don't change any of these
  geom_point() +
  theme_classic() + 
  coord_cartesian(ylim = c(0, 1)) + #, xlim = c(0, 1)) + 
  scale_y_continuous(expand = c(.01, .01)) +
  scale_x_continuous(expand = c(.01, .01)) +
  #Change color scheme here
  scale_color_continuous(low = '#3cb371', high = 'black', limits = c(0, 1), breaks = c(0, .5, 1)) + 
  #Change the title, subtitle, and axis and legend titles here
  labs(title = "XGT | Spirostomum (N = 55) and Bacteria", subtitle = "Colored by the Proportion of Spirostomum represented", x = "Bacteria + Archaea", y = "Proportion of Eukaryotes that are Spirostomum", color = '') + 
  #Don't change these
  theme(
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
    plot.subtitle = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0))
  ) + 
  geom_hline(yintercept = .5, linetype = "dashed") + 
  geom_hline(yintercept = .75, linetype = "dashed") + 
  geom_hline(yintercept = .25, linetype = "dashed") + 
  #Drawing semi-arbitrary cutoff lines for OG selection
  geom_vline(xintercept = median(counts$donor_ratio), linetype = "solid", color = "red", size = 1) + 
  geom_hline(yintercept = mean(counts$host_ratio) + sd(counts$host_ratio), linetype = "solid", color = "red", size = 1)

#Showing the plot
#presence_plot

#Generating a final OG list based on the semi-arbitrary cutoff lines above
final_og_list <- counts %>%
  filter(bacs >= median(donor_ratio)) %>%
  filter(host_ratio > mean(host_ratio) + sd(host_ratio))# | donor_ratio > .9)

#Writing the final OG list and its supporting information to a .csv file
#write.csv(final_og_list, "../final_og_list.csv")

ogs_run <- read.csv('../Case_studies/all_case_studies_ogs_to_run.csv')
ogs_run <- data.frame(ogs_run)

`%notin%` <- Negate(`%in%`)

counts <- counts %>%
  mutate(run = presence_data.OG %in% ogs_run$Gene.Family & host_ratio > mean(host_ratio) + sd(host_ratio))


#Building our output ggplot
presence_plot <- ggplot(counts, aes(x = bacs, y = host_ratio, color = run)) + 
  #Don't change any of these
  geom_point() +
  theme_classic() + 
  coord_cartesian(ylim = c(0, 1)) + #, xlim = c(0, 1)) + 
  scale_y_continuous(expand = c(.01, .01)) +
  scale_x_continuous(expand = c(.01, .01)) +
  #Change color scheme here
  scale_color_manual(values = c("grey", "darkblue")) +
  #scale_color_continuous(low = '#3cb371', high = 'black', limits = c(0, 1), breaks = c(0, .5, 1)) + 
  #Change the title, subtitle, and axis and legend titles here
  labs(title = "XGT | Spirostomum and Bacteria", x = "Bacteria + Archaea", y = "Proportion of Eukaryotes that are Spirostomum", color = '') + 
  #Don't change these
  theme(
    legend.position = 'none',
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
    plot.subtitle = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0))
  ) + 
  geom_hline(yintercept = .5, linetype = "dashed") + 
  geom_hline(yintercept = .75, linetype = "dashed") + 
  geom_hline(yintercept = .25, linetype = "dashed") + 
  #Drawing semi-arbitrary cutoff lines for OG selection
  geom_vline(xintercept = median(counts$donor_ratio), linetype = "solid", color = "red", size = 1) + 
  geom_hline(yintercept = mean(counts$host_ratio) + sd(counts$host_ratio), linetype = "solid", color = "red", size = 1)


presence_plot















