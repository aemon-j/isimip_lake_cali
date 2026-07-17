setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# clean up
rm(list=ls())
graphics.off()
cat("\14")

library(tidyverse)
library(ggpubr)
library(ggridges)
library(lemon)
library(reshape2)
library(grid)

thm <- theme_pubr(base_size = 17) + grids()

dat_lok <- read_csv("../raw_data/lake_characteristics.csv")
dat_glob <- read_csv("../raw_data/coord_area_depth.csv")


p1 <- dat_glob |> ggplot() + geom_histogram(aes(x = depth), bins = 200) + 
  scale_x_continuous(transform = "log10", limits = c(1, 1500)) +
  geom_vline(aes(xintercept = median(dat_glob$depth)), col = 2, lty = 2, size = 1) +
  xlab("Maximum depth (m)") + thm
p2 <- dat_lok |> ggplot() + geom_histogram(aes(x = max.depth.m)) +
  scale_x_continuous(transform = "log10", limits = c(1, 1500))  +
  geom_vline(aes(xintercept = median(dat_lok$max.depth.m)), col = 2, lty = 2, size = 1) +
  xlab("Maximum depth (m)") + thm
pa1 <- ggarrange(p1, p2, ncol = 1, align = "v", labels = c("A", "B"))

p3 <- dat_glob |> ggplot() + geom_histogram(aes(x = area), bins = 200) +
  scale_x_continuous(transform = "log10", limits = c(0.01, 2.5e4)) +
  geom_vline(aes(xintercept = median(dat_glob$area)), col = 2, lty = 2, size = 1) +
  xlab("Surface area (km²)") + thm
p4 <- dat_lok |> ggplot() + geom_histogram(aes(x = lake.area.sqkm)) +
  scale_x_continuous(transform = "log10", limits = c(0.01, 2.5e4)) +
  geom_vline(aes(xintercept = median(dat_lok$lake.area.sqkm)), col = 2, lty = 2, size = 1) +
  xlab("Surface area (km²)") + thm
pa2 <- ggarrange(p3, p4, ncol = 1, align = "v", labels = c("C", "D"))

ggarrange(pa1, pa2, align = "hv")
ggsave("../Output/dist_depth_area.pdf", width = 9, height = 7)
