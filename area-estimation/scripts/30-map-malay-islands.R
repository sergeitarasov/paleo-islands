# Run from paleo-islands. Polygons match the existing island-area table.
library(sf)
library(ggplot2)
library(ggrepel)
library(maps)
islands <- readRDS('data/richness/maps/selected-island-polygons.rds')
world <- st_as_sf(map('world',plot=FALSE,fill=TRUE))
# Label the largest component of each selected record.
anchors <- list()
for(i in seq_len(nrow(islands))) {
 pieces <- st_cast(st_transform(st_geometry(islands[i, ]),3857),'POLYGON')
 largest <- pieces[which.max(st_area(pieces))]
 anchors[[i]] <- st_transform(st_point_on_surface(largest),4326)[[1]]
}
points <- st_sf(island=islands$island,geometry=st_sfc(anchors,crs=4326))
labels <- cbind(st_drop_geometry(points),st_coordinates(points))
labels$island[labels$island=='Mangoli (Mangole)'] <- 'Mangoli'
# Small-island labels can move into adjacent sea, with leader lines.
p <- ggplot() +
 geom_sf(data=world,fill='#E5E5E2',colour='#B9BDBB',linewidth=.18) +
 geom_sf(data=islands,fill='#D73027',colour='#9C211B',linewidth=.18) +
 geom_label_repel(data=labels,aes(X,Y,label=island),seed=27,size=3.4,
  fill=scales::alpha('white',.92),colour='#303030',label.size=.12,
  box.padding=.30,point.padding=.1,min.segment.length=0,
  segment.colour='#606060',segment.size=.28,max.overlaps=Inf,
  max.iter=100000,max.time=30,force=2,force_pull=.5) +
 coord_sf(xlim=c(94,151.5),ylim=c(-12,21.5),expand=FALSE) +
 scale_x_continuous(breaks=seq(100,150,10)) + scale_y_continuous(breaks=seq(-10,20,10)) +
 labs(x=NULL,y=NULL,caption='Red: 40 records in the richness table. Malay Archipelago and New Guinea.\nSelected polygons: Global Islands v3; background coastlines: maps. Present-day geography.') +
 theme_minimal(base_size=12) +
 theme(panel.background=element_rect(fill='#F2F8FA',colour=NA),
  panel.grid.major=element_line(colour='#DDE6E8',linewidth=.25),
  axis.text=element_text(colour='#626C70'),plot.caption=element_text(hjust=0,size=10),
  plot.margin=margin(12,12,12,12))
ggsave('data/richness/maps/malay-archipelago-selected-islands.png',p,width=17,height=10.7,dpi=220)
ggsave('data/richness/maps/malay-archipelago-selected-islands.pdf',p,width=17,height=10.7)
write.csv(st_drop_geometry(islands),'data/richness/maps/mapped-islands.csv',row.names=FALSE)
