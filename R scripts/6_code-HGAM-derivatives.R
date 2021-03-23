#----------------------------------------------------------------------------------------
# Script: Rate-of-change Generalized Additive Modeling 
# Paper: Ecological resilience of tropical Andean lakes: a paleolimnological perspective
# Author: Benito, X.
# e-mail: xavier.benito.granell@gmail.com
#----------------------------------------------------------------------------------------

##loading libraries for functions used
library(mgcv) #allow to perform GAM analyses
library(ggplot2) #allow to make fancy graphs
library(tidyverse) #allow to subset dataframes
library(gratia)
library(cowplot)

#Function for calculating first derivatives of time series given a before,
#after, and delta step size
calc_1st_deriv = function(fit_before, fit_after,delta) {
  (fit_after-fit_before)/(2*delta)
}

#read dataframe with diatom absolute counts
mergedCores <- read.csv("Data/mergedCores_counts.csv")[-1] #read dataframe with diatom absolute counts including Fondococha

#import dataframe wiht old and new names to group
changes <- read.csv("Data/old_new_nms_cores_counts.csv", stringsAsFactors = FALSE)
#new1: ecological groups
#new2: harmonized taxonomic names

#filter cores
select <- c("llaviucu")
select <- c("yahuarcocha")
select <- c("pinan")
select <- c("fondococha")

# this is to select the target lake 
core <- mergedCores %>% 
  filter(str_detect(lake, select))

agedepth <- core[, names(core) %in% c("depth", "upper_age", "lower_age", "lake")]
diat <- core[, !names(core) %in% c("depth", "upper_age", "lower_age", "lake")]
diat[is.na(diat)] <- 0

#Select most common species 
criteria <- 0.5 #% of the total samples

n.occur <- apply(diat>0, 2, sum)
diat_red <- diat[, n.occur > (dim(diat)[1])*criteria] #
diat <- cbind(agedepth, diat_red)


##This creates the dataset containing the most common species for a single lake core
#Gather
spp_thin <- diat %>% 
  gather(key = taxa, value = count, -depth, -upper_age, -lower_age,  -lake) #don't gather depths, ages and lake variables

diat_data <- spp_thin %>%
  mutate(taxa = plyr::mapvalues(taxa, from = changes$old, to = changes$new_2)) %>%
  group_by(depth) %>%
  mutate(total_sample = sum(count)) %>% 
  filter(!total_sample == "0") %>% #this is to remove empty samples
  filter(!upper_age == 0) %>% #this is to remove ages == 0 (triumfo and fondodocha record)
  mutate(log_total_counts = log10(total_sample+1)) %>%
  mutate(relative_abundance_percent = count / sum(count) * 100) %>%
  mutate(negAge = -upper_age) %>%
  mutate(AgeCE = upper_age*(-1)+1950) %>%
  mutate(elapsedTime = abs(upper_age - lower_age)) %>%
  ungroup() %>%
  group_by(lake) %>%
  mutate(spp = factor(taxa)) 

#check how many species are included
levels(diat_data$spp)

#model S HGAM : similar smootheness between groups (spp) without global smooth 
set.seed(10) #set a seed so this is repeatable

diatom_gam_S <- gam(count ~ s(negAge, spp, k=20, bs="fs") + offset(log_total_counts),
                      weights = elapsedTime / mean(elapsedTime),
                      data=diat_data, family = poisson,
                      method = "REML")

gam.check(diatom_gam_S)
draw(diatom_gam_S)

#model I HGAM: different smootheness for each taxa without global smooth
diatom_gam_I<- gam(count ~ s(negAge, by=spp, k=20, bs="fs") +
                   s(spp, bs="re") + offset(log_total_counts),
                   weights = elapsedTime / mean(elapsedTime),
                   data=diat_data, family = nb,
                   method = "REML")

gam.check(diatom_gam_I)
draw(diatom_gam_I)

#Compare different model fits using AIC
AIC_table <- AIC(diatom_gam_S, diatom_gam_I)%>%
  rownames_to_column(var= "Model")%>%
  mutate(data_source = rep(c("diatom_data")))%>%
  group_by(data_source)%>%
  mutate(deltaAIC = AIC - min(AIC))%>%
  ungroup()%>%
  dplyr::select(-data_source)%>%
  mutate_at(.vars = vars(df,AIC, deltaAIC), 
            .funs = funs(round,.args = list(digits=0)))

#Create synthetic data to predict over a range of ages
diat_plot_data <- with(diat_data, as_tibble(expand.grid(negAge = seq(min(diat_data$negAge), max(diat_data$negAge)),
                                                        spp = factor(levels(diat_data$spp)),
                                                        log_total_counts = mean(log_total_counts))))

diat_modS_fit <- predict(diatom_gam_S, 
                         newdata = diat_plot_data,
                         se.fit = TRUE)

diat_modI_fit <- predict(diatom_gam_I,
                         newdata = diat_plot_data,
                         se.fit = TRUE)

#non-shared trends
diat_plot_data$modS_fit <- as.numeric(diat_modS_fit$fit)
diat_plot_data$modI_fit <- as.numeric(diat_modI_fit$fit)

# comparing non-shared trends
diat_plot_data <- gather(diat_plot_data, key=model, value=fit, modS_fit, modI_fit)
diat_plot_data <- mutate(diat_plot_data, se= c(as.numeric(diat_modS_fit$se.fit),
                                               as.numeric(diat_modI_fit$se.fit)),
                         upper = exp(fit + (2 * se)),
                         lower = exp(fit - (2 * se)),
                         fit   = exp(fit))

#Plot the model output for non-shared trends, with means plus standard deviations for each model.
diat_plot_model_labels <- paste("Model", c("S", "I"))
diat_plot_model_labels <- factor(diat_plot_model_labels, levels = diat_plot_model_labels)

#non-shared trends
theme_set(theme_bw())
theme_update(panel.grid = element_blank())

diat_plot <- ggplot(diat_plot_data) +
  facet_wrap(~spp, nrow = 4,scales = "free_y")+
  geom_ribbon(aes(x=negAge,
                  ymin = lower,
                  ymax = upper,
                  fill = model),
              alpha=0.2)+
  geom_point(data= diat_data, aes(x = negAge, y = count), size=0.06) +
  geom_line(aes(x = negAge, y = fit, color = model))+
  labs(y = "Absolute counts", x = "Age (cal yr BP)") +
  scale_fill_brewer(name = "", palette = "Dark2",
                    labels = diat_plot_model_labels) +
  scale_colour_brewer(name = "",
                      palette = "Dark2", labels = diat_plot_model_labels)+
  theme(legend.position = "top",
        strip.text = element_text(size=10))

diat_plot

#save plot (it takes the last chart created)
# ggsave("HGAM_diat_Fondococha.png", plot = ggplot2::last_plot(), path = "Outputs/",
#        dpi = 300)

##Derivatives and posterior distribution simulation
set.seed(10) #set a seed so this is repeatable
n_sims = 250

years <- seq(min(diat_plot_data$negAge),
             max(diat_plot_data$negAge),
             length.out = 40)

model <- diatom_gam_S
pred <- diat_modS_fit

# Generate multivariate normal simulations of model coefficients
random_coefs <- t(rmvn(n_sims, mu = coef(model),V = vcov(model)))

confint_sims <- crossing(spp=unique(diat_plot_data$spp),
                         negAge = seq(min(diat_plot_data$negAge),
                                      max(diat_plot_data$negAge),
                                      length.out = 40),
                         log_total_counts=0)

map_pred_sims <- predict(model,
                         confint_sims,
                         type = "lpmatrix") %*% random_coefs %>%
  as_data_frame() %>%
  bind_cols(confint_sims)%>%
  gather(key = simulation, value = pred, -negAge, -log_total_counts,-spp)


#specifying the step size for numerical derivative calculations
delta = 0.01

#calculating the predicted value for the current year plus delta
step_ahead_fits = confint_sims %>%
  mutate(negAge = negAge+delta)%>%
  predict(model, 
          ., type = "lpmatrix") %*% random_coefs 


#calculating the predicted value for the current year minus delta
step_behind_fits = confint_sims %>%
  mutate(negAge = negAge-delta)%>%
  predict(model,
          ., type = "lpmatrix") %*% random_coefs 


#using the predicted values for year plus and minus delta to calculate
#derivatives for each species for each simulation
derivs <- calc_1st_deriv(step_behind_fits,step_ahead_fits,delta = delta)%>%
  as_data_frame()%>%
  bind_cols(confint_sims)%>%
  gather(key = simulation,value = deriv, -spp,-negAge, -log_total_counts)

#Creating summaries of derivatives for each simulation for each year
deriv_summaries <- derivs %>%
  group_by(negAge,simulation)%>%
  summarize(deriv_mean = mean(deriv),
            deriv_sd = sd(deriv))%>%
  group_by(negAge)%>% #turning derivative summaries into 95% confidence intervals
  select(-simulation)%>%
  summarize_all(.funs = list(lower = ~quantile(.,probs = 0.025),
                             upper = ~quantile(.,probs = 0.975),
                             med   = ~quantile(.,probs = 0.5)))


#C#reating list of derivatives for each lake
#
llaviucu_deriv <- deriv_summaries
yahuarcocha_deriv <- deriv_summaries
pinan_deriv <- deriv_summaries
fondococha_deriv <- deriv_summaries
  
derivAll <- rbind(pinan_deriv, fondococha_deriv, llaviucu_deriv, yahuarcocha_deriv)
coresnms <- c(rep("Piñan", 40), rep("Fondococha", 40), 
                  rep("Llaviucu", 40), rep("yahuarcocha", 40))
    
derivAll$lake <- coresnms #this is all gam S models
#write.csv(derivAll, "Outputs/mean_sd_deriv_lakes.csv")

# 
derivAll <- derivAll
  mutate(age=-negAge) 
  
derivAll$lake <- factor(derivAll$lake, levels = c("Piñan", "Yahuarcocha", "Fondococha", "Llaviucu"))
    
#Plotting mean rate of change plus the 95% CI
mean_plot <- derivAll %>%
  ggplot(aes(x = age, 
             y = deriv_mean_med, 
             ymin = deriv_mean_lower,
             ymax = deriv_mean_upper))+
  geom_ribbon(fill="grey")+
  geom_line()+
  geom_hline(yintercept = 0, linetype=2) +
  #scale_y_continuous("average rate of change of log-abundance ")+
  ylab(expression(paste ("average rate of change of log-abundance", " (year"^-1, ")")))+
  xlim(-70, 2000)+
  xlab("Years Before Present") +
  theme_bw() +
  facet_wrap(~lake, ncol = 1,  scales = "free")
mean_plot 

ggsave("Outputs/mean_deriv.png",
       plot = mean_plot ,
       width = 5,
       height=6,
       units="in",
       dpi = 400)

#Plotting standard deviation of rate of change plus the 95% CI
sd_plot <- derivAll %>%
  ggplot(aes(x = age, 
             y = deriv_sd_med, 
             ymin=deriv_sd_lower,
             ymax=deriv_sd_upper))+
  geom_ribbon(fill="grey")+
  geom_line()+
  geom_hline(yintercept = 0, linetype=2) +
  #scale_y_continuous("st. dev of rate of change of log-abundance")+
  xlab("Years Before Present") +
  ylab(expression(paste ("st. dev of rate of change of log-abundance", " (year"^-1, ")")))+
  xlim(-70, 2000)+
  facet_wrap(~lake, ncol=1, scales = "free")+
  theme_bw()
sd_plot

ggsave("Outputs/sd_deriv.png",
       plot = sd_plot ,
       width = 5,
       height=6,
       units="in",
       dpi = 400)

# Plot Figure 5 of the paper
mean_sd_pl <- plot_grid(mean_plot, sd_plot, align = "v", ncol = 2, labels = "auto")

ggsave("figures/HGAM_mean_sd_deriv.png",
       plot = mean_sd_pl ,
       width = 10,
       height=8,
       units="in",
       dpi = 400)



