# Run from paleo-islands. AICc-weighted mixture of full area distributions.
library(ggplot2)
weights <- read.csv('area-estimation/results/selected-27/models.csv')
weights$delta <- weights$AICc - min(weights$AICc)
weights$weight <- exp(-weights$delta/2)
weights$weight <- weights$weight / sum(weights$weight)
weights$draw_model <- c('simple_prior1', 'intercepts_prior1', 'separate_prior1')
curves <- read.csv('area-estimation/results/bayesian-selected27-slopes-10ma/area-distributions.csv')
intercepts <- read.csv('area-estimation/results/bayesian-selected27-intercepts-10ma/area-distributions.csv')
curves <- rbind(curves, intercepts)
dir.create('area-estimation/results/aicc-selected27-10ma', recursive = TRUE, showWarnings = FALSE)
dir.create('area-estimation/plots/aicc-selected27-10ma', recursive = TRUE, showWarnings = FALSE)
summary <- data.frame()
mixtures <- data.frame()
for(prior in unique(curves$prior)) {
 selected <- curves[curves$prior == prior, ]
 area <- selected$area_km2[selected$model == 'simple_prior1']
 mass <- rep(0, length(area))
 for(i in 1:nrow(weights)) {
  component <- selected[selected$model == weights$draw_model[i], ]
  stopifnot(identical(component$area_km2, area))
  component_mass <- component$density_log10_area * .01
  stopifnot(abs(sum(component_mass)-1) < 1e-8)
  mass <- mass + weights$weight[i] * component_mass
 }
 stopifnot(abs(sum(mass)-1) < 1e-8)
 z <- log10(area)
 q <- 10^approx(cumsum(mass), z, xout = c(.025, .5, .975), ties = 'ordered')$y
 summary <- rbind(summary, data.frame(prior, mean = sum(area*mass), median = q[2], lower = q[1], upper = q[3]))
 mixtures <- rbind(mixtures, data.frame(prior, area_km2 = area, mass, density_log10_area = mass/.01))
}
write.csv(weights, 'area-estimation/results/aicc-selected27-10ma/model-weights.csv', row.names = FALSE)
write.csv(summary, 'area-estimation/results/aicc-selected27-10ma/area-summary.csv', row.names = FALSE)
write.csv(mixtures, 'area-estimation/results/aicc-selected27-10ma/area-distributions.csv', row.names = FALSE)
show <- mixtures[mixtures$prior == 'log_uniform_10_1million' & mixtures$area_km2 >= 10 & mixtures$area_km2 <= 1e6, ]
endpoints <- show$area_km2 %in% c(10, 1e6)
show$density_log10_area[endpoints] <- 2 * show$density_log10_area[endpoints]
base <- summary[summary$prior == 'log_uniform_10_1million', ]
p <- ggplot(show, aes(area_km2, density_log10_area)) +
 geom_area(fill = '#0072B2', alpha = .18) + geom_line(colour = '#0072B2', linewidth = .8) +
 geom_vline(xintercept = base$median, linetype = 2) + scale_x_log10() + theme_classic(base_size = 12) +
 labs(x = expression('Effective area (km'^2*')'), y = 'Density per log10 area')
ggsave('area-estimation/plots/aicc-selected27-10ma/area-distribution.png', p, width = 8, height = 5, dpi = 160)
ggsave('area-estimation/plots/aicc-selected27-10ma/area-distribution.pdf', p, width = 8, height = 5)
print(weights)
print(summary)
