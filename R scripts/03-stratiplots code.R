##set WD
setwd("/Volumes/xbenitogranell-data/0_project/data/cores")

##loading libraries for functions used
library(analogue) #to join diatom datasets on their common spp
library(rioja) #to merge diatom datasets on their common spp
library(plyr) #allow to join dataframes by common column
library(dplyr) #allow to summarise variables and manipulate multiple dataframes
library(ggplot2) #to make nice plots
library(tidyverse)
library(cluster)


#read diatom data
  #mergedCores <- read.csv("mergedCores_counts3.csv") #read dataframe with diatom absolute counts including Fondococha
  mergedCores <- read.csv("mergedCores_counts4.csv")[,-1] #with new Fondococha agedepth model
  
  diatoms_save <- mergedCores #save dataframe
  
  #Gather
  spp_thin <- diatoms_save %>% 
    gather(key = taxa, value = count, -depth, -upper_age, -lower_age,  -lake)#don't gather depths, ages and lake variables


  #import dataframe wiht old and new names to group
  changes <- read.csv("old_new_nms_cores_counts.csv", stringsAsFactors = FALSE)
    #new1: ecological groups
    #new2: harmonized taxonomic names
    
  
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
      
      select <- c("titicaca")
      select <- c("umayo")
      select <- c("triumfo")
      select <- c("lagunillas")
      

      core_lake <- spp_long %>%
        filter(str_detect(lake, select)) %>% #select lake
        filter(!upper_age==0.0) %>%
        group_by(taxa) %>%
        filter(count > 0) %>% #remove species with 0 counts
        ungroup() %>%
        group_by(depth) %>%
        mutate(relative_abundance_percent = count / sum(count) * 100) %>% #calculate RA
        ungroup()
      
      # filter more abundant taxa; format need to be on long-wide format-->no spreaded 
      core_common_taxa <- core_lake %>%
        group_by(taxa) %>%
        summarise(max_rel_abund = max(relative_abundance_percent)) %>%
        filter(max_rel_abund >= 5) %>%
        arrange(max_rel_abund) %>%
        pull(taxa)
      
      # select from initial table
      core_counts_common <- core_lake %>%
        filter(taxa %in% core_common_taxa) %>%
        mutate(taxa = factor(taxa, levels = core_common_taxa)) %>%
        arrange(taxa)
      
      #make it wide
      core_counts_wide <- core_counts_common %>%
        dplyr::select(depth, lake, upper_age, taxa, relative_abundance_percent) %>%
        spread(key = taxa, value = relative_abundance_percent)

      length(core_common_taxa)
      
      #do coniss
      core_counts_wide[is.na(core_counts_wide)] <- 0
      
      diatHel <- decostand(core_counts_wide[,4:ncol(core_counts_wide)], method="hellinger")
      diss <- vegdist(diatHel, method="bray")
      clust <- chclust(diss, method="coniss")
      bstick(clust)
      
      zones <- cutree(clust, k=4)
      locate <- cumsum(rle(zones)$lengths)+1
      zones <- core_counts_wide[locate, ][,3]
      zones <- zones$upper_age
      

## diatom functional groups shifts
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
     
############################## END

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
      zones = zones,
      # make the bar colour black
      col = "black")
    
dev.off()  
  

#plot using ggplot function
stratPlot <- ggplot(core_counts_common, 
                    aes(y = upper_age, x = relative_abundance_percent)) +
  # draw horizontal lines of the appropriate length for each depth
  geom_segment(aes(xend = 0, yend = upper_age), lwd = 1) +
  # facet by taxon, keeping distance on the x axis comparable between facets
  facet_grid(~taxa, scales = "free_x", space = "free_x") +
  # reverse the y axis for depth
  scale_y_reverse() +
  # make all facets use the same break values
  # (helps with cluttered breaks)
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50)) +
  # set the x and y labels
  labs(x = "Relative Abundance (%)", y = "Age (cal yr BP)", title="Pinan diatoms") +
  # customize the appearance
  theme(
    # rotate the facet labels
    strip.text.x = element_text(angle = 40, hjust = 0, vjust = 0), 
    # turn off the label background
    strip.background = element_blank()
  )+
theme_bw()

stratPlot

# voodoo that makes it so that facet labels can overlap
# https://stackoverflow.com/questions/49740215/ggplot-facet-grid-label-cut-off
stratPlot_grob <- ggplotGrob(stratPlot)
for(i in which(grepl("strip-t", stratPlot_grob$layout$name))){
  stratPlot_grob$grobs[[i]]$layout$clip <- "off"
}

# needed to draw the modified plot_grob
grid::grid.draw(stratPlot_grob)

#save plots
pdf("Pinan.pdf", 
    width = 6.5, 
    height = 4
)
grid::grid.draw(stratPlot_grob)
dev.off()  

## Using tidypaleo R package (https://fishandwhistle.net/post/2018/stratigraphic-diagrams-with-tidypaleo-ggplot2/)
library(tidypaleo)
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
  #add tephra layers Fondococha
  annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=2082, ymax=2117, alpha=0.7, fill="grey") +
  annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=1996, ymax=2063, alpha=0.7, fill="grey") +
  annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=2255, ymax=2345, alpha=0.7, fill="grey") +
  # add CONISS zones
  geom_hline(yintercept = zones, col = "blue", lty = 1, alpha = 0.7) 
diat_plot

## prepare geochemical data 
# FIRST go to lines 264
geochem_plot <- ggplot(geochem_data, aes(x = value, y = upper_age)) +
  geom_lineh() +
  geom_point() +
  scale_y_reverse() +
  facet_geochem_gridh(vars(variable)) +
  labs(x = NULL) 

geochem_plot <- geochem_plot +
  geom_lineh_exaggerate(exaggerate_x = 4, col = "grey70", lty = 2)

strat_plt <- wrap_plots(
  diat_plot + 
    theme(strip.background = element_blank(), strip.text.y = element_blank()),
  geochem_plot +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  nrow = 1,
  widths = c(4, 1)
)
strat_plt
ggsave("fondococha_stratplot.png", strat_plt, height = 6, width = 10)


## add CONISS dendograme
diat_coniss <- core_counts_common %>%
  nested_data(qualifiers = c(upper_age,depth), key = taxa, value = relative_abundance_percent) %>%
  nested_chclust_coniss()

diat_plot +
  layer_dendrogram(diat_coniss, aes(y = depth), param = "CONISS") 


wrap_plots(
  diat_plot + 
    theme(strip.background = element_blank(), strip.text.y = element_blank()),
  geochem_plot +
    layer_dendrogram(diat_coniss, component = "CONISS", aes(y = upper_age)) +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  nrow = 1,
  widths = c(2, 1)
)

## merge species and non-species data  
#Use a left-join to add non-species data
llaviucu_xrf <- read.csv("llaviucu_xrf.csv")
yahuarcocha_geochem <- read.csv("yahuarcocha_geochem_newbatch.csv")
umayo_geochem <- read.csv("umayo_isotopes.csv")
pinan_geochem <- read.csv("pinan_geochem.csv")
triumfo_geochem <- read.csv("triumfo_geochem.csv")
fondococha_xrf <- read.csv("fondococha_xrf.csv", row.names = 1)

#read Llaviucu pollen
llaviucu_pollen_ratios <- read.csv("llaviucu_pollen.csv")
llaviucu_pollen <- read.csv("llaviucu_pollen_raw.csv")

#make it long for joining with non-species data
core_counts_long <- core_counts_common %>%
  select(depth, lake, upper_age, taxa, relative_abundance_percent) 


#llaviucu
    core_diat_geochem <- core_counts_long %>%
      left_join(llaviucu_xrf, by = "depth") %>%
      mutate(Fe_Mn = Fe/Mn) %>%
      mutate(Si_Ti = Si/Ti) %>%
      mutate(K_Ti = K/Ti) %>%
      select(taxa, depth, everything())
    
    #table for tidypaleo stratiplot
    geochem_data <- core_diat_geochem %>% select(Fe_Mn,Si_Ti,K_Ti,upper_age) %>%
        gather(key=variable,value=value,-upper_age) 
      
    core_diat_geochem_wide <- core_diat_geochem %>%
      select(Fe_Mn, Si_Ti, K_Ti, depth, upper_age, taxa, relative_abundance_percent) %>%
      spread(key = taxa, value = relative_abundance_percent)
    
    write.csv(core_diat_geochem_wide, "llaviucu_diat_proxy.csv")
    
        # chunk to join with covariates (code in fit-GAM-PrC)
        #core_diat_covariates <- core_counts_long %>%
          #left_join(sitesData, by = "depth") %>%
          #select(taxa, depth, everything())
        
        #core_diat_covariates_wide <- core_diat_covariates %>%
          #select(depth, upper_age, taxa, relative_abundance_percent, d13C, spor, Fe_Mn, Si_Ti) %>%
          #filter(upper_age <= 2000) %>% 
          #spread(key = taxa, value = relative_abundance_percent)

#Yahuarcocha
    core_diat_geochem <- core_counts_long %>%
      left_join(yahuarcocha_geochem, by = "depth") %>%
      mutate(C_N = orgC / N) %>% #calculate C/N ratio
      select(taxa, depth, everything())
    
    #table for tidypaleo stratiplot
    geochem_data <- core_diat_geochem %>% select(CaCO3,d13C,d18O,C_N,upper_age) %>%
      gather(key=variable,value=value,-upper_age) 
    
    core_diat_geochem_wide <- core_diat_geochem %>%
      select(taxa, relative_abundance_percent, upper_age, d13C, C_N, CaCO3, d18O, depth) %>%
      spread(key = taxa, value = relative_abundance_percent)

    write.csv(core_diat_geochem_wide, "yahuarcocha_diat_proxy.csv")
    
#Umayo
    core_diat_geochem <- core_counts_long %>%
      left_join(umayo_geochem, by = "depth") %>%
      select(taxa, depth, everything())
    
    core_diat_geochem_wide <- core_diat_geochem %>%
      select(taxa, relative_abundance_percent, upper_age, d18O, depth) %>%
      spread(key = taxa, value = relative_abundance_percent)
    
    write.csv(core_diat_geochem_wide, "umayo_diat_proxy.csv")
    
    
#pinan
    core_diat_geochem <- core_counts_long %>%
      left_join(pinan_geochem, by = "depth") %>%
      mutate(C_N = C/N) %>%
      select(taxa, depth, everything())
    
        geochem_data <- core_diat_geochem %>% select(d13C, C_N,upper_age) %>%
          gather(key=variable,value=value,-upper_age) 
        
    core_diat_geochem_wide <- core_diat_geochem %>%
      select(taxa, relative_abundance_percent, upper_age, d13C, C_N, depth) %>%
      spread(key = taxa, value = relative_abundance_percent)
    
    write.csv(core_diat_geochem_wide, "pinan_diat_proxy_2RA.csv")

    
#Fondococha
    core_diat_geochem <- core_counts_long %>%
      left_join(fondococha_xrf, by = "depth") %>%
      mutate(Fe_Mn = Fe/Mn) %>%
      mutate(Si_Ti = Si/Ti) %>%
      mutate(K_Ti = K/Ti) %>%
      select(taxa, depth, everything())
    
    #table for tidypaleo stratiplot
    geochem_data <- core_diat_geochem %>% select(Fe_Mn,Si_Ti,K_Ti,upper_age) %>%
      gather(key=variable,value=value,-upper_age) 
    
    core_diat_geochem_wide <- core_diat_geochem %>%
      select(Fe_Mn, Si_Ti, K_Ti, depth, upper_age, taxa, relative_abundance_percent) %>%
      spread(key = taxa, value = relative_abundance_percent)
    
    write.csv(core_diat_geochem_wide, "fondococha_diat_proxy_3RA.csv")
    
        
#triumfo
    core_diat_geochem <- core_counts_long %>%
      filter(upper_age>0.0) %>%
      left_join(triumfo_geochem, by = "depth") %>%
      #mutate(C_N = C/N) %>%
      select(taxa, depth, everything())
    
    core_diat_geochem_wide <- core_diat_geochem %>%
      select(taxa, relative_abundance_percent, upper_age, d13C, d15N, depth) %>%
      spread(key = taxa, value = relative_abundance_percent)
    

    
png("Stratiplot.png", width = 11, height = 8, res = 300, units = "in")
  
#varTypes argument needs to specify that the non-species variables should have independently sized axes. 
#This should be a vector with the same number of elements as variables in the plot (I’ve use rep() to repeat “relative” and “absolute” the correct number of times.

# code for diat spp and geochemistry
Stratiplot(
  core_diat_geochem_wide %>% select(-depth, -upper_age),
  core_diat_geochem_wide$upper_age, 
  varTypes = c(rep("absolute", 3), rep("absolute", length(core_common_taxa))), 
  ylab = "Age (cal yr BP)", 
  xlab = "Relative abundance (%)",
  topPad = 10, 
  type = c("h"),
  col = "black", 
  zones = zones
)

dev.off()


Stratiplot(
  core_diat_covariates_wide %>% select(-depth, -upper_age),
  core_diat_covariates_wide$upper_age, 
  varTypes = c(rep("relative", 4), rep("absolute", length(core_common_taxa))), 
  ylab = "Age (cal yr BP)", 
  xlab = "Relative abundance (%)",
  topPad = 10, 
  type = "h", 
  col = "black", 
  zones = zones
)

dev.off()

#Llaviucu
core_counts_geochem <- cbind(
  core_counts_wide,
  xrf %>% select(-depth, -upper_depth, -lower_depth, -upper_age, -lower_age) %>%
    mutate(Fe_Mn = Fe/Mn) %>%
    mutate(Si_Ti = log10((Si/Ti)+0.25))
)










## Plot non-species data

data_thin <- xrf %>%
  mutate(Fe_Mn = Fe/Mn) %>%
  mutate(Si_Ti = log10((Si/Ti)+0.25)) %>%
  gather(key = element, value = value, -depth, -upper_age, -lower_age, -upper_depth, -lower_depth)


data_thin <- data %>%
  mutate(C_N = orgC / N) %>% #calculate C/N ratio
  gather(key = element, value = value, -depth, -upper_age, -lower_age)


data_thin %>%
  ggplot(aes(y = depth, x = value)) +
  geom_path() +
  geom_point() +
  facet_wrap(~element, scales = "free_x") +
  scale_y_reverse() +
  labs(x = NULL, y = "Depth (cm)") +
  theme_bw()

llaviucu_age_depth <- data_thin %>%
  select(depth, upper_age)

data_thin %>%
  ggplot(aes(y = depth, x = value)) +
  geom_path() +
  geom_point() +
  facet_wrap(~element, scales = "free_x") +
  scale_y_reverse(
    sec.axis=sec_axis(
      trans = ~approx(llaviucu_age_depth, xout = .)$y,
      name = "Age (cal Yr BP)",
      #breaks = c(2000, 1950, 1900, 1850, 1800, 1750)
    ),
    expand=c(0,0)
  ) +
  labs(x = NULL, y = "Depth (cm)") +
  theme_bw()


yahuarcocha_age_depth <- data %>%
  select(depth, upper_age)

data_thin %>%
  mutate(facet_label=fct_recode(
    element,
    "C(%)" = "C",
    "CaCO3(%)" = "CaCO3", 
    #"delta ^ 13 * C" = "d13C",
    #"delta ^ 13 * C (acidifed)" = "d13C.1",
    #"delta ^ 18 * O" = "d18O",
    "N (%)" = "N", 
    "org C (%)" = "orgC", 
    "TIC (%)" = "TIC"
  )) %>%
  ggplot(aes(y = depth, x = value)) +
  geom_path() +
  geom_point() +
  facet_wrap(~facet_label, scales = "free_x") +
  scale_y_reverse(
    sec.axis=sec_axis(
      trans = ~approx(yahuarcocha_age_depth, xout = .)$y,
      name = "Age (cal Yr BP)",
      #breaks = c(2000, 1950, 1900, 1850, 1800, 1750)
    ),
    expand=c(0,0)
  ) +
  labs(x = NULL, y = "Depth (cm)") +
  theme_bw()






##
diat.groups <- combined %>% 
  mutate(count1 = combined %>% 
           select(starts_with("Aulaco")) %>% 
           rowSums(),
         count2 = diat_red %>% 
           select(starts_with("Fra")) %>% 
           rowSums())


diat_save <- diat_red

diat.groups <- diat_red %>% 
  mutate(Cymbella.cymbiformis = diat %>% 
           select(starts_with("Cymbella")) %>% 
           rowSums(),
         Gomphonema.vibrio = diat %>% 
           select(starts_with("Gomphonema")) %>% 
           rowSums())






