#---------------------------------------------------------------------------
# Script: Diatom stratigraphic plots (supplementary figures S1-4)
# Paper: Ecological resilience of tropical Andean lakes: a paleolimnological perspective
# Author: Benito, X.
# e-mail: xavier.benito.granell@gmail.com
#---------------------------------------------------------------------------

##loading libraries for functions used
library(analogue) #to make stratiplots
library(rioja) #constrained hierarchical clustering
library(ggplot2) #to make nice plots
library(tidyverse) #allow to summarise variables and manipulate multiple dataframes
library(cluster)


#read diatom data
mergedCores <- read.csv("data/mergedCores_diatomcounts.csv")[,-1] 
diatoms_save <- mergedCores #save dataframe
  
#Gather
spp_thin <- diatoms_save %>% 
  gather(key = taxa, value = count, -depth, -upper_age, -lower_age,  -lake)#don't gather depths, ages and lake variables

#read taxonomic harmonisation diatom names list
changes <- read.csv("data/old_new_nms_cores_counts.csv", stringsAsFactors = FALSE)
#column_new1: ecological groups (ecological guilds, as Passy 2007)
#column_new2: harmonized taxonomic names

#spread--> wide format
spp_wide <- spp_thin %>%
  mutate(taxa = plyr::mapvalues(taxa, from = changes$old, to = changes$new_2)) %>%
  group_by(depth, lake, upper_age, taxa) %>%
  summarise(count = sum(count)) %>%
  spread(key = taxa, value = count)
  
#no spread --> long format
spp_long <- spp_thin %>%
  mutate(taxa = plyr::mapvalues(taxa, from = changes$old, to = changes$new_2)) %>%
  group_by(depth, lake, upper_age, taxa) %>%
  summarise(count = sum(count))
  
#filter cores
  select <- c("llaviucu")
  select <- c("pinan")
  select <- c("yahuarcocha")
  select <- c("fondococha")
      
core_lake <- spp_long %>%
  filter(str_detect(lake, select)) %>% #select lake
  filter(!upper_age==0.0) %>%
  group_by(taxa) %>%
  filter(count > 0) %>% #remove species with 0 counts
  ungroup() %>%
  group_by(depth) %>%
  mutate(relative_abundance_percent = count / sum(count) * 100) %>% #calculate relative abundance
  ungroup()
      
# filter more abundant taxa; format need to be on long-wide format-->no spreaded 
core_common_taxa <- core_lake %>%
   group_by(taxa) %>%
   summarise(max_rel_abund = max(relative_abundance_percent)) %>%
   filter(max_rel_abund >= 5) %>%
   arrange(max_rel_abund) %>%
   pull(taxa)
        
# select from initial species matrix
core_counts_common <- core_lake %>%
  filter(taxa %in% core_common_taxa) %>%
  mutate(taxa = factor(taxa, levels = core_common_taxa)) %>%
  arrange(taxa)
      
#make it wide
  core_counts_wide <- core_counts_common %>%
    dplyr::select(depth, lake, upper_age, taxa, relative_abundance_percent) %>%
    spread(key = taxa, value = relative_abundance_percent)

#check number of species to plot
length(core_common_taxa)
      
#do coniss to add statistically significant stratigraphic zones
  core_counts_wide[is.na(core_counts_wide)] <- 0
      
  diatHel <- decostand(core_counts_wide[,4:ncol(core_counts_wide)], method="hellinger")
  diss <- vegdist(diatHel, method="bray")
  clust <- chclust(diss, method="coniss")
  bstick(clust)
      
  zones <- cutree(clust, k=4)
  locate <- cumsum(rle(zones)$lengths)+1
  zones <- core_counts_wide[locate, ][,3]
  zones <- zones$upper_age
      
## OPTIONAL plot diatom functional groups shifts
#this is to calculate planktic:benthic ratios
  new <- diatoms_save %>% 
      gather(key = taxa, value = count, -depth, -upper_age, -lower_age, -lake) %>%
      mutate(taxa = plyr::mapvalues(taxa, from = changes$old, to = changes$new_1)) %>% #ecological grouping
      group_by(depth, taxa, lake, upper_age, lower_age) %>%
      summarise(count = sum(count)) %>%
      filter(!count == "0" ) %>% #this is to remove empty samples (rows)
      filter(!upper_age==0.0) %>% #this is to drop extraneous ages
      ungroup() %>%
      group_by(depth, lake) %>%
      mutate(relative_abundance_percent = count / sum(count) * 100) %>%
      ungroup()
      

lake <- new %>% filter(lake=="fondococha") 

library(viridis)
seq_palette <- viridis(5)

gg <- ggplot(lake, aes(x=as.numeric(as.character(upper_age)), y=relative_abundance_percent))
plt <- gg + geom_area(aes(colour=taxa, fill=taxa)) +
  ylab("% total assemblage") +
  xlab("Years cal BP")+
  theme_bw()
plt
     
##############################
## 
## plot using analogue Stratiplot (need data input in wide format --> spreaded)
png("Stratiplot.png", width = 11, height = 8, res = 300, units = "in")
    
  Stratiplot(
      core_counts_wide %>% select(-depth, -lake, -upper_age),
      core_counts_wide$upper_age,
      ylab = "Cal yr BP", 
      xlab = "Relative abundance (%)",
      # adds padding to the top of the plot
      # to fix cut-off taxa names
      topPad = 10, 
      # make the plot type a "bar" plot
      type = "h", 
      #sort = "wa",
      # add stratigraphic zones from cluster analyses (regime shifts R file)
      #zones = zones,
      # make the bar colour black
      col = "black")
    
dev.off()  
  
## Using tidypaleo R package (https://fishandwhistle.net/post/2018/stratigraphic-diagrams-with-tidypaleo-ggplot2/)
library(tidypaleo) #remotes::install_github("paleolimbot/tidypaleo")
library(patchwork)
theme_set(theme_bw(9))

diat_plot <- ggplot(core_counts_common, aes(x = relative_abundance_percent, y = upper_age)) +
  geom_col_segsh() +
  scale_y_reverse() +
  facet_abundanceh(vars(taxa)) +
  labs(x = "Relative abundance (%)", y = "Cal yr BP") +
  #add tephra layers Piñan
  # annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=764, ymax=783, alpha=0.7, fill="grey") +
  # annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=964, ymax=1024, alpha=0.7, fill="grey") +
  #add tephra layers for Fondococha
  # annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=2082, ymax=2117, alpha=0.7, fill="grey") +
  # annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=1996, ymax=2063, alpha=0.7, fill="grey") +
  # annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=2255, ymax=2345, alpha=0.7, fill="grey") +
  # add CONISS zones
  geom_hline(yintercept = zones, col = "blue", lty = 1, alpha = 0.7) 
diat_plot




