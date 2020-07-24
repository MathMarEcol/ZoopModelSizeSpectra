
## Now lets have a look at how the output changes in 50 year increments

library(tidyverse)
library(patchwork)
fi <- list.files("RawOutput", full.names = TRUE)

int <- 50
tm <- seq(int,1000,int)

# Create empty df to store data
out <- tibble("SST" = numeric(),
              "Chl" = numeric(),
              "Run" = numeric(),
              "Time" = numeric(),
              "Abund" = numeric(),
              "Growth" = numeric(),
              "Mort" = numeric())

for (i in 1:length(fi)){
  dat <- read_rds(fi[i])
  for (j in 1:length(tm)){

    sp <- dat$model$param$Groups$Species

    df <- tibble("SST" = rep(NA, length(sp)),
                 "Chl" = rep(NA, length(sp)),
                 "Run" = rep(NA, length(sp)),
                 "Time" = rep(NA, length(sp)),
                 "Abund" = rep(NA, length(sp)),
                 "Growth" = rep(NA, length(sp)),
                 "Mort" = rep(NA, length(sp)))

    for (k in 1:length(sp)){

      df$Run[k] <- i
      df$Time[k] <- tm[j]
      df$SST[k] <- dat$model$param$sst
      df$Chl[k] <- dat$model$param$chlo
      df$Species[k] <- sp[k]
      df$Abund[k] <- sum(colMeans(dat$model$N[(ceiling(0.5*tm[j]):tm[j]),k,], dims = 1))
      df$Growth[k] <- sum(colMeans(dat$model$gg[(ceiling(0.5*tm[j]):tm[j]),k,], dims = 1))
      df$Mort[k] <- sum(colMeans(dat$model$Z[(ceiling(0.5*tm[j]):tm[j]),k,], dims = 1))


      # df$Abund[k] <- sum(colMeans(dat$model$N[(ceiling(0.5*tm[j]):tm[j]),k,], dims = 1))
      # df$Growth[k] <- sum(colMeans(dat$model$gg[(ceiling(0.5*tm[j]):tm[j]),k,], dims = 1))
      # df$Mort[k] <- sum(colMeans(dat$model$Z[(ceiling(0.5*tm[j]):tm[j]),k,], dims = 1))
    }

    out <- bind_rows(out, df)
    rm(df)

  }
}


x11(width = 12, height = 6)
out %>%
  # filter(Species == "OmniCopepods") %>%
  mutate(Run = as.factor(Run),
         Species = factor(Species, levels=sp),
         ChlSST = factor(paste0(Chl, "Chl-", SST, "SST"))) %>%
  ggplot(aes(x = Time, y = log10(Abund))) +
  # geom_point(size  = 0.5) +
  geom_line(aes(colour = Run), size = 0.1) +
  facet_wrap(~ Species, ncol = 4, scales = "free") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

ggsave("Figures/SpeciesAccumulation.png", dpi = 200)


