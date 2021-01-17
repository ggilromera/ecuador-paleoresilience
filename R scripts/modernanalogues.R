##Read in the diatom core list and trainingset (absolute counts)
df <- readRDS("data/coresList.rds")

#Remove empty spp resulting from merging dataframes (and drop year & depths vars)
remove <- function(i, cores, ...) {
  core <- cores[[i]]
  core <- core[, -which(names(core) %in% c("depth", "upper_age", "lower_age", "lake", "AgeCE"))] # drop year & depths vars
  # comment the next line when predicting core trajectories in the timetrack analysis
  core <- core[, colSums(core) > 0] #select only present species
  #core <- tran(core, "hellinger")
  return(core)
}

cores <- lapply(seq_along(df), remove, cores=df)

#name list elements
names(cores) <- c("Fondococha", "Lagunillas", "Llaviucu", "Pinan", "Titicaca", "Triumfo", "Umayo", "Yahuarcocha", "trainingset")

#Extract trainingset from the list
trainingset <- cores$trainingset

## Read in environmental lake dataset
lake_env <- read.csv("data/lake_env.csv", row.names = 1)

## Modern analogues
env <- lake_env$Water.T
#change one of c("Llaviucu", "Fondococha", "Yahuarcocha", "Pinan")
lake <- "Fondococha"
lakedepth <- "fondococha"

mod <- rioja::MAT(trainingset, env, dist.method="bray")

pred <- predict(mod,cores[[lake]]/100)
ages <- as.numeric(coresList[[lakedepth]]$upper_age)
# plot(ages, pred$dist.n[,1], ylab="Squared chord distance", xlab="Age cal yr BP")
# goodpoorbad <- quantile(pred$dist.n[,1], probs = c(0.75, 0.95))
# abline(h=goodpoorbad, col=c("orange", "red"))

#code from Richard Telford (https://github.com/richardjtelford/Zabinskie)
dist_to_analogues <- as.data.frame(ages) %>% 
  mutate(
    dist_to_analogues = pred$dist.n[, 1],
    quality  = cut(dist_to_analogues, breaks = c(0, goodpoorbad, Inf), labels = c("good", "poor", "bad"))
  )
attr(dist_to_analogues, which = "goodpoorbad") <- goodpoorbad

dist_to_analogues_plot <- ggpalaeo:::plot_diagnostics(x = dist_to_analogues, x_axis = "ages", y_axis = "dist_to_analogues", 
                                                      goodpoorbad = attr(dist_to_analogues, "goodpoorbad"), fill = c("salmon", "lightyellow", "skyblue"), categories = c("Good", "Fair", "None")) + 
  labs(x = "Cal yr BP", y = "Squared chord distance", fill = "Analogue quality") +
  ggtitle("LLaviucu")

#ggsave("dist_to_analogues_plot_fondococha.png", dist_to_analogues_plot, height = 8, width = 10)

#
dist_to_analogues_pinan <- dist_to_analogues
dist_to_analogues_yahuarcocha <- dist_to_analogues
dist_to_analogues_fondococha <- dist_to_analogues
dist_to_analogues_llaviucu <- dist_to_analogues

dist_to_analogues_all <- bind_rows(dist_to_analogues_pinan, dist_to_analogues_yahuarcocha,
                                   dist_to_analogues_fondococha, dist_to_analogues_llaviucu)

#rename Piñan
dist_to_analogues_all <- dist_to_analogues_all %>% mutate(lake=str_replace(lake,"Pinan", "Piñan"))

#arrange dataframe by latitude of lakes
dist_to_analogues_all$lake <- factor(dist_to_analogues_all$lake, levels = c("Piñan", "Yahuarcocha", "Fondococha", "Llaviucu"))

write.csv(dist_to_analogues_all, "dist_to_analogues_all_lakes.csv")

dist_to_analogues_all_plt <- ggplot(data=dist_to_analogues_all, aes(x=ages, y=dist_to_analogues, col=quality)) +
  geom_point() +
  scale_color_manual(values=c("skyblue", "orange", "red"))+
  facet_wrap(~lake, scales = "free")+
  xlab("Cal yr BP") +
  ylab("Squared chord distance") +
  labs(col = "Quality")+
  theme_bw()

ggsave("dist_to_analogues_plot_all.png", dist_to_analogues_all_plt, height = 8, width = 10)

