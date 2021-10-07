
# Plot PPMR
fZooMSS_Plot_PPMR <- function(in_dat){

  ppmr <- PPMR_plot(in_dat)

  gg <- ggplot() +
    geom_line(data = ppmr[[2]], mapping = aes(x = Betas, y = y, colour = Species), size = 1) +
    geom_line(data = ppmr[[1]], mapping = aes(x = x, y = y), size = 1.2) +
    theme_bw() +
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
    labs(x = expression('log' [10] * PPMR),
         y = "Zoop. Biomass Proportion", subtitle = "PPMR") +
    geom_vline(data = ppmr[[1]], mapping = aes(xintercept = mn_beta), colour = 'black') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_colour_manual(values = c("Flagellates" = in_dat$model$param$Groups$PlotColour[in_dat$model$param$Groups$Species=="Flagellates"],
                                   "Ciliates" = in_dat$model$param$Groups$PlotColour[in_dat$model$param$Groups$Species=="Ciliates"],
                                   "Larvaceans" = in_dat$model$param$Groups$PlotColour[in_dat$model$param$Groups$Species=="Larvaceans"],
                                   "Salps" = in_dat$model$param$Groups$PlotColour[in_dat$model$param$Groups$Species=="Salps"],
                                   "Jellyfish" = in_dat$model$param$Groups$PlotColour[in_dat$model$param$Groups$Species=="Jellyfish"],
                                   "CarnCopepods" = in_dat$model$param$Groups$PlotColour[in_dat$model$param$Groups$Species=="CarnCopepods"],
                                   "Chaetognaths" = in_dat$model$param$Groups$PlotColour[in_dat$model$param$Groups$Species=="Chaetognaths"],
                                   "Euphausiids" = in_dat$model$param$Groups$PlotColour[in_dat$model$param$Groups$Species=="Euphausiids"],
                                   "OmniCopepods" = in_dat$model$param$Groups$PlotColour[in_dat$model$param$Groups$Species=="OmniCopepods"],
                                   "Omnivores" = in_dat$model$param$Groups$PlotColour[in_dat$model$param$Groups$Species=="Omnivores"],
                                   "Carnivores" = in_dat$model$param$Groups$PlotColour[in_dat$model$param$Groups$Species=="Carnivores"],
                                   "FilterFeeders" = in_dat$model$param$Groups$PlotColour[in_dat$model$param$Groups$Species=="FilterFeeders"]))
}

# Plot Abundance Size Spectra
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

# Plot Biomass Size Spectra
fZooMSS_Plot_BiomassSizeSpectra <- function(dat) {
  species <- sweep(dat$abundances, 2, dat$model$param$w, "*")
  
  rownames(species) <- dat$model$param$Groups$Species
  species <- as_tibble(t(species))
  
  species <- species %>%
    add_column("Weight" = dat$model$param$w) %>%
    pivot_longer(-Weight, names_to = "Species", values_to = "Biomass") %>%
    filter(Biomass > 0) %>%
    mutate(Species = factor(Species, levels = dat$model$param$Groups$Species))
  
  gg <- ggplot(data = species, mapping = aes(x = log10(Weight), y = log10(Biomass), colour = Species)) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = dat$model$param$Groups$PlotColour) +
    theme_bw() +
    labs(subtitle = "Biomass Spectrum")
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

# Plot biomass by time
fZooMSS_Plot_BiomassTimeSeries <- function(dat){
  tspecies <- rowSums(sweep(dat$model$N, 3, dat$model$param$w, "*"), dims = 2)
  colnames(tspecies) <- dat$model$param$Groups$Species
  tspecies <- as_tibble(tspecies)
  tspecies$Time <- seq(dat$model$param$dt * dat$model$param$isave,
                       dat$model$param$tmax,
                       dat$model$param$dt * dat$model$param$isave)
  tspecies <- tspecies %>%
    pivot_longer(-Time, names_to = "Species", values_to = "Biomass") %>%
    filter(Biomass > 0) %>%
    mutate(Species = factor(Species, levels = dat$model$param$Groups$Species))

  gg <- ggplot(data = tspecies, mapping = aes(x = Time, y = log10(Biomass), colour = Species)) +
    geom_line(size = 1) +
    geom_point(size = 1.2) +
    scale_color_manual(values = dat$model$param$Groups$PlotColour) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(subtitle = "Biomass") +
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


# Written by Jason Everett (UQ/CSIRO/UNSW) and Yiwen Lu (UQ)
# Last Updated 7th October 2021
fZooMSS_PlotPredPrey <- function(Groups){

  Groups <- Groups %>%
    mutate(W0_prey = case_when(is.na(PPMR) ~ W0 - log10(fZooMSS_betas(fZooMSS_g2esd(10^W0), PPMRscale)), # Add prey size
                               !is.na(PPMR) ~ W0 - log10(PPMR)), # Add prey size
           WMax_prey = case_when(is.na(PPMR) ~ Wmax - log10(fZooMSS_betas(fZooMSS_g2esd(10^Wmax), PPMRscale)),
                                 !is.na(PPMR) ~ Wmax - log10(PPMR))) # Add prey size

  gg <- ggplot(data = Groups) +
    geom_segment(aes(x = Species, xend = Species, y = W0, yend = Wmax, colour = Species), size = 7, alpha = 0.6, show.legend = FALSE) +
    geom_segment(aes(x = Species, xend = Species, y = W0_prey, yend = WMax_prey, colour = Species), size = 3, alpha = 1, show.legend = FALSE) +
    coord_flip() +
    scale_colour_manual(values = Groups$PlotColour, breaks = Groups$Species) +
    scale_x_discrete(limits = as.character(Groups$Species))+
    ylab(expression(paste("log"[10],"Body Size (g)"))) +
    theme_bw() +
    theme(axis.title.y = element_blank())

}
