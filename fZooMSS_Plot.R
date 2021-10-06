
# Plot PPMR
fZooMSS_Plot_PPMR <- function(dat){

  out <- PPMR_plot(dat)

  gg <- ggplot() +
    geom_line(data = out[[2]], mapping = aes(x = Betas, y = y, colour = Species), size = 1) +
    geom_line(data = out[[1]], mapping = aes(x = x, y = y), size = 1.2) +
    theme_bw() +
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
    labs(x = expression('log' [10] * PPMR),
         y = "Zoop. Biomass Proportion", subtitle = "PPMR") +
    geom_vline(data = out[[1]], mapping = aes(xintercept = mn_beta), colour = 'black') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_colour_manual(values = c("Flagellates" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="Flagellates"],
                                   "Ciliates" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="Ciliates"],
                                   "Larvaceans" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="Larvaceans"],
                                   "Salps" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="Salps"],
                                   "Jellyfish" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="Jellyfish"],
                                   "CarnCopepods" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="CarnCopepods"],
                                   "Chaetognaths" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="Chaetognaths"],
                                   "Euphausiids" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="Euphausiids"],
                                   "OmniCopepods" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="OmniCopepods"]))
}

# Plot Size Spectra
fZooMSS_Plot_SizeSpectra <- function(dat) {
  species <- dat$abundances

  rownames(species) <- dat$model$param$Groups$Species
  species <- as_tibble(t(species))

  species <- species %>%
    add_column("Weight" = dat$model$param$w) %>%
    pivot_longer(-Weight, names_to = "Species", values_to = "Abundance") %>%
    filter(Abundance > 0) %>%
    mutate(Species = factor(Species, levels = dat$model$param$Groups$Species))

  gg <- ggplot(data = species, mapping = aes(x = log10(Weight), y = log10(Abundance), colour = Species)) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = dat$model$param$Groups$PlotColour) +
    theme_bw() +
    labs(subtitle = "Abundance Spectrum")
}

# Plot abundance by time
fZooMSS_Plot_AbundTimeSeries <- function(dat){
  tspecies <- rowSums(dat$model$N, dims = 2)
  colnames(tspecies) <- dat$model$param$Groups$Species
  tspecies <- as_tibble(tspecies)
  tspecies$Time <- seq(dat$model$param$dt * dat$model$param$isave,
                       dat$model$param$tmax,
                       dat$model$param$dt * dat$model$param$isave)
  tspecies <- tspecies %>%
    pivot_longer(-Time, names_to = "Species", values_to = "Abundance") %>%
    filter(Abundance > 0) %>%
    mutate(Species = factor(Species, levels = dat$model$param$Groups$Species))

  gg <- ggplot(data = tspecies, mapping = aes(x = Time, y = log10(Abundance), colour = Species)) +
    geom_line(size = 1) +
    geom_point(size = 1.2) +
    scale_color_manual(values = dat$model$param$Groups$PlotColour) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(subtitle = "Abundance") +
    xlab("Time (Years)")
}

# Plot growth by time
fZooMSS_Plot_GrowthTimeSeries <- function(dat){
  gr <- rowSums(dat$model$gg, dims = 2) / length(dat$model$param$w)
  colnames(gr) <- dat$model$param$Groups$Species
  gr <- as_tibble(gr)
  gr$Time <- seq(dat$model$param$dt * dat$model$param$isave,
                 dat$model$param$tmax,
                 dat$model$param$dt * dat$model$param$isave)
  gr <- gr %>%
    pivot_longer(-Time, names_to = "Species", values_to = "Growth") %>%
    filter(Growth > 0) %>%
    mutate(Species = factor(Species, levels = dat$model$param$Groups$Species))

  gg <- ggplot(data = gr, mapping = aes(x = Time, y = log10(Growth), colour = Species)) +
    geom_line(size = 1) +
    geom_point(size = 1.2) +
    scale_color_manual(values = dat$model$param$Groups$PlotColour) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(subtitle = "Growth Rate") +
    xlab("Time (Years)")
}


# Plot predation by time
fZooMSS_Plot_PredTimeSeries <- function(dat){

  Z <- rowSums(dat$model$Z,dims = 2) / length(dat$model$param$w)
  colnames(Z) <- dat$model$param$Groups$Species
  Z <- as_tibble(Z)
  Z$Time <- seq(dat$model$param$dt * dat$model$param$isave,
                dat$model$param$tmax,
                dat$model$param$dt * dat$model$param$isave)
  Z <- Z %>%
    pivot_longer(-Time, names_to = "Species", values_to = "Mortality") %>%
    filter(Mortality > 0) %>%
    mutate(Species = factor(Species, levels = dat$model$param$Groups$Species))

  gg <- ggplot(data = Z, mapping = aes(x = Time, y = Mortality, colour = Species)) +
    geom_line(size = 1) +
    geom_point(size = 1.2) +
    scale_color_manual(values = dat$model$param$Groups$PlotColour) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(subtitle = "Mortality Rate") +
    xlab("Time (Years)")
}
