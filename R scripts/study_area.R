
library(maps)
library(raster)
library(cowplot)
library(ggrepel)
library(cividis)
library(ggsn)

world <- map_data("world")

southamerica <- ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="lightgrey") +
  theme(legend.position = "right")+
  coord_map("albers", parameters = c(-100, -100),  ylim=c(-30,15), xlim=c(-82,-50)) +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw()


#Plot elevation base map (raster library)
DEM <- raster("data/DEM/dem2.bil")
ext<-extent(-82,-50,-30,15)
ext_ecuador<-extent(-81,-77,-4,2)
altmod<-crop(DEM,ext)
altmod_ecuador <- crop(DEM,ext_ecuador)

#convert the raster to points for plotting
map.p <- rasterToPoints(altmod)
map.p_ecuador <- rasterToPoints(altmod_ecuador)

#Make the points a dataframe for ggplot
df <- data.frame(map.p)
df_ecuador <- data.frame(map.p_ecuador)

#Make appropriate column headings
colnames(df) <- c("long", "lat", "Elevation")
colnames(df_ecuador) <- c("long", "lat", "Elevation")

## select Andes climatic records
andes_clim_records <- read.csv("data/andes_climatic_records.csv") %>%
  filter(record %in% c("Ubaque","Pumacocha", "Palestina", "Triunfo"))

## read Ecuador core geographic coordinates
ecuador_cores <- read.csv("data/ecuador_cores.csv",
                          row.names = 1)[1:4,]
  
## composite plot without dem just for testing
# plt <- southamerica +
#   geom_point(data=env_data_lakes, aes(x=long, y=lat,shape=""),size=2)+
#   geom_point(data=andes_clim_records, aes(x=long, y=lat, col=record),
#              shape=18, size=4) +
#   labs(color = "Hydroclimatic records")+
#   scale_shape_manual(name="",labels = c("Modern Lakes"),
#                      values = c(16))+
#   theme(text= element_text(size=12),
#         legend.direction="vertical",
#         legend.position = c(0.75,0.5),
#         legend.box.just = c("top"), 
#         #legend.title = element_blank(),
#         legend.background =element_rect(fill=alpha('white', 0.4)))+
#   annotate(geom = "rect", ymax = -4, ymin = 2, xmax = -81, xmin = -77, colour = "brown", fill = NA)
# plt  
# 
# ecuador <- ggplot() +
#   geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="lightgrey") +
#   theme(legend.position = "right")+
#   coord_map("albers", parameters = c(-100, -100),  ylim=c(-4,2), xlim=c(-81,-77)) +
#   xlab("Longitude") + ylab("Latitude") +
#   theme_bw()
# 
# plt_ecuador <- ecuador +
#   geom_point(data=env_data_lakes, aes(x=long, y=lat),colour="grey",size=2)+
#   geom_point(data=ecuador_cores, aes(x=Longitude, y=Latitude),shape=19,bg=16,size=3)+
#   geom_text_repel(data=ecuador_cores, aes(Longitude, Latitude, label = lake)) +
#   theme(legend.position = "bottom",
#         legend.title = element_blank())
# plt_ecuador  
# 
# plt_grid <- plot_grid(plt,plt_ecuador, align="hv", axis="t", ncol = 2,
#                  rel_widths = c(1,1),labels = "auto")
# plt_grid
###

## composite map with dem
plot_sa <- ggplot(data=df, aes(y=lat, x=long)) +
  geom_raster(aes(fill=Elevation)) +
  #scale_fill_distiller(palette = "RdYlBu") +
  scale_fill_cividis()+
  #scale_colour_gradient(high = "red") +
  geom_point(data=env_data_lakes, aes(x=long, y=lat, shape=""), colour="grey" ,size=3) +
  geom_point(data=env_data_lakes, aes(x=long, y=lat), shape=1,colour="black" ,size=3) +
  # geom_point(data=env_data_lakes, aes(x=long, y=lat, shape=""), size=3, colour="grey") +
  geom_point(data=andes_clim_records, aes(x=long, y=lat,col=record),shape=18,size=3)+
  labs(color = "Paleoclimatic records")+
  scale_shape_manual(name="",labels = c("Modern Lakes"),
                     values = c(16))+
  guides(fill = guide_colourbar(title="Elevation (m)")) +
  theme(text= element_text(size=12),
        #legend.direction="vertical",
        legend.box.just = c("top"),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.background=element_rect(fill=alpha('white', 0.4)))+
  annotate(geom = "rect", ymax = -4, ymin = 2, xmax = -81, xmin = -77, colour = "red", fill = NA)+
  xlab("Longitude (deg)")+
  ylab("Latitude (deg)")+
  theme_bw() +
  coord_equal()

# select modern lakes within the Ecuador window latitudes
lakes_for_plt <- env_data_lakes %>% filter(between(lat, -4,2))

seq_palette <- viridis(4)

plt_ecuador <- ggplot(data=df_ecuador, aes(y=lat, x=long))+
  geom_raster(aes(fill=Elevation), show.legend = "FALSE") +
  #scale_fill_distiller(palette = "RdYlBu") +
  scale_fill_cividis()+
  geom_point(data=lakes_for_plt, aes(x=long, y=lat, shape=""),colour="grey",size=3)+
  geom_point(data=lakes_for_plt, aes(x=long, y=lat),shape=1,colour="black",size=3)+
  geom_point(data=ecuador_cores, aes(Longitude, Latitude, col=lake), shape=16, size=3)+
  #scale_colour_viridis_d()+
  scale_colour_manual(values=c(seq_palette[3], seq_palette[4], seq_palette[1], seq_palette[2]))+
  geom_text_repel(data=ecuador_cores, aes(Longitude, Latitude, label = lake),col="white") +
  xlab("")+
  ylab("")+
  theme_bw() +
  theme(plot.background = element_rect(colour="red",size=2),
        # panel.background = element_blank(),
        panel.background = element_rect(color = "grey", fill="white"),
        panel.grid.major = element_line(size=0.1, linetype = "solid", color = "grey"),
        panel.grid.minor = element_line(size=0.1, linetype = "solid", color = "grey"),
        legend.position = "none")+
  north(df_ecuador, symbol = 4) + #add north arrow
  scalebar(df_ecuador, dist = 20, dist_unit = "km",
           transform = TRUE, model = "WGS84", st.dist = 0.02, st.size = 2)+
  coord_equal()
plt_ecuador  

plot_composite <- plot_grid(plot_sa,plt_ecuador,align="hv", axis="t", vjust=c(6,6),
                            ncol = 2,  rel_widths = c(1,0.7),labels = "auto")
plot_composite
ggsave("figures/map_figure1.png", plot_composite, height = 8, width = 10)

##
# plot1 <- ggplot(data=df, aes(y=lat, x=long)) +
#   geom_raster(aes(fill=Elevation)) +
#   scale_fill_distiller(palette = "RdYlBu") +
#   #scale_colour_gradient(high = "red") +
#   geom_point(data=env_data_lakes, aes(x=long, y=lat,shape=""),size=2) +
#   geom_point(data=andes_clim_records, aes(x=long, y=lat,col=record),shape=18,size=3)+
#   labs(color = "Records")+
#   scale_shape_manual(name="",labels = c("Modern Lakes"),
#                      values = c(16))+
#   xlab("Longitude")+
#   ylab("Latitude")+
#   theme(legend.title = element_blank())+
#   theme_bw() +
#   coord_equal()

ggsave("map_figure1.png", plot1, height = 8, width = 10)


