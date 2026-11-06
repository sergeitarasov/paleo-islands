# Run from paleo-islands. Highlight only records smaller than Timor.
library(sf)
library(ggplot2)
library(maps)
islands <- readRDS('data/richness/maps/selected-island-polygons.rds')
timor.area <- islands$area_km2[islands$island=='Timor']
islands <- islands[islands$area_km2 < timor.area, ]
world <- st_as_sf(map('world',plot=FALSE,fill=TRUE))
anchors <- list()
for(i in seq_len(nrow(islands))) {
 pieces <- st_cast(st_transform(st_geometry(islands[i, ]),3857),'POLYGON')
 largest <- pieces[which.max(st_area(pieces))]
 anchors[[i]] <- st_transform(st_point_on_surface(largest),4326)[[1]]
}
points <- st_sf(island=islands$island,geometry=st_sfc(anchors,crs=4326))
labels <- cbind(st_drop_geometry(points),st_coordinates(points))
# Explicit offshore placements keep every leader line visible.
positions <- read.table(text='island lon lat
Adonara 125.1 -10.1
Ambon 126 -5
Aru_Islands 137 -7
Bacan 130.6 0.5
Bali 111.5 -7
Basilan 119.2 5.5
Bohol 127.5 9
Buru 123.9 -2.7
Cebu 127 11.7
Flores 122.5 -6
Halmahera 131 2.7
Kei_Islands 135 -3.5
Komodo 118 -11.5
Langkawi 97.5 7.8
Larat 134 -9.3
Leyte 128.4 12.8
Lombok 113 -10.8
Mangoli 122.8 -0.3
Marinduque 122 16
Mindoro 117 14
Negros 119.5 8.2
Palawan 114 10
Panay 117.8 11.8
Rinca 120.5 -12.5
Samar 129 14.4
Seram 131.5 -1.6
Sibuyan 125 15
Sumba 118 -13.2
Sumbawa 116.2 -6.3
Ternate 124.5 3
Tioman 106.5 4
Wetar 128 -6.2',header=TRUE)
positions$island <- gsub('_',' ',positions$island)
labels$island[labels$island=='Mangoli (Mangole)'] <- 'Mangoli'
j <- match(labels$island,positions$island)
stopifnot(nrow(islands)==32,!anyNA(j))
labels$lon <- positions$lon[j]; labels$lat <- positions$lat[j]
p <- ggplot() +
 geom_sf(data=world,fill='#E6E7E4',colour='#C2C6C4',linewidth=.2) +
 geom_sf(data=islands,fill='#D73027',colour='#A02520',linewidth=.22) +
 geom_segment(data=labels,aes(lon,lat,xend=X,yend=Y),colour='#505B60',linewidth=.32) +
 geom_point(data=labels,aes(X,Y),shape=21,fill='#D73027',colour='#505B60',size=1.1,stroke=.25) +
 geom_label(data=labels,aes(lon,lat,label=island),size=3.6,fill='white',colour='#303030',
            label.size=0,label.padding=grid::unit(.15,'lines')) +
 coord_sf(xlim=c(96,139),ylim=c(-14.5,17),expand=FALSE) +
 scale_x_continuous(breaks=seq(100,135,5)) + scale_y_continuous(breaks=seq(-10,15,5)) +
 labs(x=NULL,y=NULL,caption=paste0('Red: ',nrow(islands),' table records smaller than Timor (',format(round(timor.area),big.mark=','),' km²). Larger islands remain grey.\nLeader lines point to the island, or the largest component of an island group. Polygons: Global Islands v3; background: maps.')) +
 theme_minimal(base_size=12) +
 theme(panel.background=element_rect(fill='#F2F8FA',colour=NA),
       panel.grid.major=element_line(colour='#DFE8EA',linewidth=.25),
       axis.text=element_text(colour='#626C70'),plot.caption=element_text(hjust=0,size=10),
       plot.margin=margin(12,12,12,12))
ggsave('data/richness/maps/malay-islands-smaller-than-timor.png',p,width=16,height=12,dpi=220)
ggsave('data/richness/maps/malay-islands-smaller-than-timor.pdf',p,width=16,height=12)
write.csv(st_drop_geometry(islands),'data/richness/maps/islands-smaller-than-timor.csv',row.names=FALSE)
