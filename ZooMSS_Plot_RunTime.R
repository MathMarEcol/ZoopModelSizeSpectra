library(tidyverse)
library(patchwork)

# Load all files
fi <- list.files("RawOutput/", full.names = TRUE)

# Create empty df to store data
df <- tibble("No" = rep(NA, length(fi)),
             "SST" = rep(NA, length(fi)),
             "Chl" = rep(NA, length(fi)),
             "RunTime" = rep(NA, length(fi)))

for (i in 1:length(fi)){

  dat <- read_rds(fi[i])

  df$No[i] <- i
  df$SST[i] <- dat$model$param$sst
  df$Chl[i] <- dat$model$param$chlo
  df$RunTime[i] <- as.numeric(dat$model$model_runtime[3])/60
  rm(dat)

}
# Plot
gg1 <- ggplot(data = df, aes(x = log10(Chl), y = RunTime)) +
  geom_point() +
  geom_smooth(method = "loess") +
  theme_bw() +
  ylab("Model Run Time (mins)")

gg2 <- ggplot(data = df, aes(x = SST, y = RunTime)) +
  geom_point() +
  geom_smooth(method = "loess") +
  theme_bw() +
  ylab("Model Run Time (mins)")

x11(width = 12, height = 6)
gg1 / gg2
ggsave("Figures/Model_RunTime.png", dpi = 300)

