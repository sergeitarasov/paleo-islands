# Island isolation: definitions and provenance

Calculated 26 October 2026 for all 61 records, without exclusions.

- **Mainland distance:** minimum great-circle coastline distance to Australia, Africa, South America, North America or Eurasia. Offshore islands are excluded from this mainland definition.
- **Large-landmass distance:** minimum distance to any other island with published `area_geode >= 10000` km², or a continental mainland. The query returned 75 qualifying island polygons. A landmass need not be larger than the focal island. All polygons belonging to the focal record are excluded, including all constituents of a focal archipelago.
- **Names:** verified focal-island names replace inconsistent database placenames for matching polygon IDs; other targets retain their database names. Frederik Hendrik is also called Dolak/Yos Sudarso. The database contains six continental geometries: the sixth is Eurasia's Chukotka segment beyond the antimeridian (its area attribute is zero), and it is retained and labelled Eurasia.
- **Archipelagos:** all polygons used for the original area are retained. Distances are from the nearest constituent, not the main island or centroid. This particularly matters for New Caledonia, whose geographic scope includes remote dependencies. This definition should be revisited if richness covers only the main island.

## Data and computation

Source: Global Islands v3, USGS/Esri/UNEP-WCMC. Service:
https://data-gis.unep-wcmc.org/server/rest/services/Global_Islands/MapServer
Island layer: /2; continental layer: /3. Focal archipelago components also include layers /0 and /1. The original 1,078 selected polygon IDs, geodesic areas and bounds are archived in `results/isolation/geography/focal-polygon-ids.csv`. Their area sums match all current workbook areas to within 0.001 km².

Focal polygons use the same original full-resolution geometry selections as the area compilation. Target queries use `where=area_geode>=10000` (layer 2), `where=1=1` (layer 3), `returnGeometry=true`, `outSR=3857`, `maxAllowableOffset=20` or `100`. Source geometry was made valid and simplified at matching tolerances in EPSG:3857. Complete boundary lines were transformed to longitude/latitude. `sf::st_distance()` with `sf_use_s2(TRUE)` calculates minimum spherical great-circle distances, converted from metres to kilometres. These are boundary-to-boundary distances, not centroid or vertex-only distances, and not routes constrained to stay over water. They are spherical rather than WGS84-ellipsoidal distances; spreadsheet precision is not geodetic accuracy.

The 20 m result is used. A full repeat at 100 m changed mainland distances by at most 0.0441 km and large-landmass distances by at most 0.0836 km, with no nearest-large-landmass identity changes. This checks simplification sensitivity, not accuracy of the original coastline. Check values are in `resolution-check.csv`.

The self-contained compressed `geography/coastlines-20.rds` archives all source lines and all eligible target lines, with names and IDs. `scripts/17-island-distances.R` recalculates the distance CSV from this snapshot. It does not redownload or silently update geometry. The workbook's Island isolation tab stores the same numeric distances, matched by island name, displayed to 0.1 km. Its original eight columns and the original Island richness tab are unchanged.

## Regression

Run `Rscript area-estimation/scripts/18-isolation-regression.R` from the project root. See `reports/isolation-regression.md`. All distances are contemporary: applying this relationship to Sunken Islands at 9 Ma requires a palaeogeographic isolation estimate, not a modern distance substituted without justification.
