# Five fossils and six secondary calibrations. Node IDs refer to the reduced tree.
# Time unit: 100 Ma. min_age and max_age are targets for fitting the priors.
calibrations <- read.table(header = TRUE, stringsAsFactors = FALSE, text = "
node  name                   method  min_age  max_age  offset  min_prob
57    Protohybosorus         cauchy  1.49     2.00     0.1     1e-8
60    Antiqusolidus          cauchy  1.25     1.50     0.1     1e-8
56    Lamprima-Lucanus       gamma   1.25     1.60     0.2     0.025
101   O_statzi               skewT   0.25     0.50     0.1     1e-8
81    Gymnopleurus           skewT   0.135    0.50     0.1     1e-8
88    Heliocopris            skewT   0.187    0.50     0.1     1e-8
54    root                   gamma   1.904    2.40     0.2     0.025
55    Scarabaeoidea          gamma   1.52     2.08     0.2     0.025
58    Scarabaeidae           gamma   1.11     1.59     0.2     0.025
61    Aphodius_Scarabaeinae  gamma   0.59     1.25     0.2     0.025
90    Copris_Onthophagus     gamma   0.37     0.80     0.2     0.025
")
