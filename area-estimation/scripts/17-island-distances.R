# Run from paleo-islands. Recalculate distances from the archived coastline snapshot.
library(sf)
sf_use_s2(TRUE)

# Global Islands v3 coastline geometries, prepared at 20 m Mercator tolerance.
# The snapshot preserves all 61 focal geometries and every >=10,000 km2 target.
geo <- readRDS('area-estimation/results/isolation/geography/coastlines-20.rds')
ids <- read.csv('area-estimation/results/isolation/geography/focal-polygon-ids.csv')
source <- geo$sources
targets <- geo$targets
mainlands <- targets[grepl('^3:',targets$id), ]

# s2 calculates minimum great-circle distances between complete coastline lines.
main_distance <- as.matrix(units::drop_units(st_distance(source,mainlands)))/1000
large_distance <- as.matrix(units::drop_units(st_distance(source,targets)))/1000
for (i in seq_len(nrow(source))) {
  own <- ids[ids$island==source$island[i], ]
  own_ids <- paste0(own$layer,':',own$objectid_1)
  large_distance[i,targets$id %in% own_ids] <- Inf
}
nearest_main <- max.col(-main_distance,ties.method='first')
nearest_large <- max.col(-large_distance,ties.method='first')
distances <- data.frame(island=source$island,
  distance_mainland_km=main_distance[cbind(seq_len(nrow(source)),nearest_main)],
  nearest_mainland=mainlands$name[nearest_main],
  distance_large_landmass_km=large_distance[cbind(seq_len(nrow(source)),nearest_large)],
  nearest_large_landmass=targets$name[nearest_large],
  nearest_large_landmass_id=targets$id[nearest_large])
stopifnot(all(distances$distance_mainland_km>0),
  all(distances$distance_large_landmass_km>0),
  all(distances$distance_large_landmass_km<=distances$distance_mainland_km+1e-8))
write.csv(distances,'area-estimation/results/isolation/distances.csv',row.names=FALSE)
print(distances)
# The Excel Island isolation tab contains these values, matched by island name.
