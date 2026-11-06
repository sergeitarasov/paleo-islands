# Run from this directory: Rscript make-mcmctree.R
library(ape)
library(MCMCtreeR)

# Keep the original labels so the tree matches the archived MCMCtree input.
full_tree <- read.tree("sunken_run7_70p.tre")
tips_to_keep <- c(
  "Helictopleurus_furcicornis",
  "Helictopleurus_neoamplicollis",
  "Euoniticellus_fulvus",
  "Onthophagus_pipitzi",
  "Onthophagus_taurus",
  "Digitonthophagus_gazella",
  "Pycnopanelus_krikkeni",
  "Dwesasilvasedis_medinae",
  "Metacatharsius_sp1",
  "Epilissus_sp_morpho_1",
  "Arachnodes_sp6",
  "Apterepilissus_ovalis",
  "Copris_fidius",
  "Pedaria_sp2",
  "Temnoplectron_lewisensis",
  "Heliocopris_sp1",
  "Heliocopris_sp2",
  "Synapsis_sp",
  "Nanos_morpho_2",
  "Nanos_sp4",
  "Scarabaeus_n_sp1",
  "Garreta_unicolor",
  "Paragymnopleurus_sp_2",
  "Tanzanolus_sp1",
  "Nesovinsionia_vinsoni",
  "Janssensantus_mabuensis_n_sp1",
  "Panelus_rhodesiensis",
  "S_Neosisyphus_spinipes",
  "Sisyphus_australis",
  "Nesosisyphus_pygmaeus",
  "Epirinus_sp1",
  "Epactoides_major",
  "Grebennikovius_lupanganus_sp_n1",
  "Sarrophorus_costatus",
  "Odontoloma_planatum",
  "Aphodius_immundus",
  "Chrysina_resplendens",
  "Protaetia_brevitarsis",
  "Lepidiota_sp",
  "Hybosoridae_sp1_Madagascar",
  "Lucanus_capreolus",
  "Lamprima_sp",
  "Nicrophorus_vespilloides",
  "Onitis_caffer",
  "Thyregis_relictus",
  "Circellium_bacchus",
  "Chalconotus_convexus",
  "Panelus_sp3",
  "Nesosisyphus_rotundatus",
  "Nesosisyphus_regnardi",
  "Nesosisyphus_vicinus",
  "Bohepilissus_sp1",
  "Ochichanton_sp1"
)
stopifnot(all(tips_to_keep %in% full_tree$tip.label))
reduced_tree <- keep.tip(full_tree, tips_to_keep)
stopifnot(Ntip(reduced_tree) == 53)

# Initialize a label for every node before MCMCtreeR inserts the priors.
reduced_tree$node.labels <- rep(NA_character_, reduced_tree$Nnode)

source("calib_s3_all-points.R")
groups <- tipDes(reduced_tree, calibrations$node)
names(groups) <- calibrations$name

fit <- MCMCtreePhy(
  phy = reduced_tree,
  monoGroups = groups,
  minAges = calibrations$min_age,
  maxAges = calibrations$max_age,
  method = calibrations$method,
  offset = calibrations$offset,
  minProb = calibrations$min_prob,
  estimateScale = TRUE,
  plot = FALSE,
  writeMCMCtree = TRUE,
  MCMCtreeName = "output-all-points.tre"
)
