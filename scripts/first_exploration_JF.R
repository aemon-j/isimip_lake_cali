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

thm <- theme_pubr(base_size = 17) + grids()
##------------- read in data --------------------##

# list files
files <- list.files(file.path("..", "raw_data"), recursive = TRUE) |>
  data.frame() |> setNames("file") |> filter(file != "README.md") |>
  filter(!grepl(pattern = "lake_characteristics.*", x = file))

# read in files and combine to a single data.frame
dat <- lapply(files$file, function(f) {
  read.csv(file.path("..", "raw_data", f), header = TRUE)}) |>
  reshape2::melt(id.vars = 1:8) |> select(-L1)


meta <- read.csv("../raw_data/lake_characteristics.csv")

##-------- calculations/trends ----------------------------------##

## trends

dat_trend <- dat |> pivot_wider() |>
  group_by(lake, model, cali, scenario) |>
  reframe(sl_surftemp_mean = coefficients(lm(surftemp_mean ~ year))[2],
          ic_surftemp_mean = coefficients(lm(surftemp_mean ~ year))[1],
          m_surftemp_mean = predict(lm(surftemp_mean ~ year), data.frame(year = median(year))), # mean value at middle of time period
          sl_bottemp_mean = coefficients(lm(bottemp_mean ~ year))[2],
          ic_bottemp_mean = coefficients(lm(bottemp_mean ~ year))[1],
          m_t_bottemp_mean = predict(lm(bottemp_mean ~ year), data.frame(year = median(year))), # mean value at middle of time period
          sl_strat_mean = coefficients(lm(strat_sum ~ year))[2],
          ic_strat_mean = coefficients(lm(strat_sum ~ year))[1],
          m_strat_mean = predict(lm(strat_sum ~ year), data.frame(year = median(year))), # mean value at middle of time period
          sl_latentheatf_mean = coefficients(lm(latentheatf_mean ~ year))[2],
          ic_latentheatf_mean = coefficients(lm(latentheatf_mean ~ year))[1],
          m_latentheatf_mean = predict(lm(latentheatf_mean ~ year), data.frame(year = median(year))),
          sl_sensheatf_mean = coefficients(lm(sensheatf_mean ~ year))[2],
          ic_sensheatf_mean = coefficients(lm(sensheatf_mean ~ year))[1],
          m_sensheatf_mean = predict(lm(sensheatf_mean ~ year), data.frame(year = median(year))),
          sl_heat_mean = coefficients(lm(heat_mean ~ year))[2],
          ic_heat_mean = coefficients(lm(heat_mean ~ year))[1],
          m_heat_mean = predict(lm(heat_mean ~ year), data.frame(year = median(year))))



dat_trends_diff <- dat_trend |> pivot_longer(cols = 5:ncol(dat_trend)) |>
  pivot_wider(names_from = cali, values_from = value,
              id_cols = c(model, scenario, lake, name)) |>
  mutate(diff = `uncalibrated` - calibrated) |>
  select(- calibrated, -`uncalibrated`)

##

dat_metr <- dat |> pivot_wider(names_from = cali, values_from = value) |>
  group_by(model, scenario, lake, name) |>
  reframe(R = cor(calibrated, uncalibrated),
          bias = mean(uncalibrated - calibrated, na.rm = TRUE),
          var = var(uncalibrated - calibrated, na.rm = TRUE)) |>
  pivot_longer(5:7, names_to = "metr")


## variance decomposition
var_frac <- function(value, model, gcm, scenario, lake) {
  an <- anova(lm(value ~ model * gcm * scenario * lake))
  totsst <- sum(an$`Sum Sq`)
  sep <- an$`Sum Sq`
  sep <- c(sep[1:4], sum(sep[5:length(sep)]))/totsst
  return(sep)
}
tst <- dat |> filter(name == "surftemp_mean") |>
  pivot_wider(names_from = cali, values_from = value) |>
  group_by(model, scenario, lake, gcm, name) |>
  reframe(R = cor(calibrated, uncalibrated),
          bias = mean(uncalibrated - calibrated, na.rm = TRUE),
          var = var(uncalibrated - calibrated, na.rm = TRUE)) |>
  pivot_longer(6:8, names_to = "metr") |>
  mutate(model = as.factor(model),
         gcm = as.factor(gcm),
         scenario = as.factor(scenario),
         lake = as.factor(lake)) |>
  group_by(metr) |>
  reframe(f_model = var_frac(value, model, gcm, scenario, lake)[1],
          f_gcm = var_frac(value, model, gcm, scenario, lake)[2],
          f_scen = var_frac(value, model, gcm, scenario, lake)[3],
          f_lake = var_frac(value, model, gcm, scenario, lake)[4],
          f_inter = var_frac(value, model, gcm, scenario, lake)[5])


##-------- first plots ----------------------------------##

## temporat developement of difference
for(v in unique(dat$name)) {
  p <- dat |> left_join(meta, by = c(lake = "Lake.Short.Name")) |>
    pivot_wider(names_from = cali, values_from = value) |>
    mutate(diff = uncalibrated - calibrated) |> filter(name == v) |>
    group_by(kmcluster, scenario, name, year) |>
    reframe(med = median(diff, na.rm = TRUE),
            q05 = quantile(diff, 0.05, na.rm = TRUE),
            q25 = quantile(diff, 0.25, na.rm = TRUE),
            q75 = quantile(diff, 0.75, na.rm = TRUE),
            q95 = quantile(diff, 0.95, na.rm = TRUE)) |>
    ggplot() +
    geom_ribbon(aes(x = year, ymin = q25, ymax = q75, fill = kmcluster), alpha = 0.5) +
    geom_line(aes(x = year, y = med, col = kmcluster), lwd = 1.23) +
    facet_wrap(.~scenario, scales = "free_x") + ggtitle(v) +
    scale_fill_viridis_d() + scale_color_viridis_d() + thm
  
  ggsave(file.path("..", "Output", paste0("diff_over_time_", v, ".png")),
         p, width = 13, height = 13)
  
}

## boxplots for metrics comparing cali and uncali

for(v in unique(dat$name)) {
p <- dat_metr |> left_join(meta, by = c(lake = "Lake.Short.Name")) |>
  filter(name == v) |> ggplot() +
  geom_boxplot(aes(x = kmcluster, y = value, fill = kmcluster)) + 
  facet_grid(metr~scenario, scales = "free") + scale_fill_viridis_d() + thm

ggsave(file.path("..", "Output", paste0("boxplot_", v, ".png")),
       p, width = 13, height = 13)
}

## scatter plots

for(v in unique(dat$name)) {
  p <- dat |> filter(name == v) |>
    pivot_wider(names_from = cali, values_from = value) |>
    left_join(meta, by = c(lake = "Lake.Short.Name")) |>
    ggplot() + geom_point(aes(x = calibrated, y = uncalibrated, col = model)) +
    geom_abline(aes(slope = 1, intercept = 0), lty = "dashed") +
    facet_wrap(scenario~kmcluster, scales = "free") +
    scale_color_viridis_d() + thm + ggtitle(v)
  
  ggsave(filename = file.path("..", "Output", paste0("scatter_", v, ".png")),
         plot = p, width = 17, height = 17)
}

# dist trend

for(v in unique(dat_trends_diff$name)) {
  p <- dat_trends_diff |> filter(name == v) |>
    mutate(scenario = case_match(scenario,
                                 "picontrol" ~ "Picontrol",
                                 "historical" ~ "Historical",
                                 "ssp126" ~ "SSP1-2.6",
                                 "ssp370" ~ "SSP3-7.0",
                                 "ssp585" ~ "SSP5-8.5")) |>
    mutate(scenario = factor(scenario,
                             levels = c("Picontrol", "Historical", "SSP1-2.6",
                                        "SSP3-7.0", "SSP5-8.5"))) |>
    left_join(meta, by = c(lake = "Lake.Short.Name")) |>
    ggplot() + geom_density_ridges(aes(x = diff, y = model, fill = model)) +
    facet_grid(scenario ~ kmcluster, scale = "free") +
    theme_pubr(base_size = 16) + grids() +
    xlab("temp slope difference (°C/a)") +
    ggtitle(paste0("Difference in estimated ", v),
            subtitle = "between uncalibirated and claibriated models (uncalibrated − calibrated)") +
    scale_fill_viridis_d("Model", option = "C", end = 0.9) + ylab("") +
    thm
  
  ggsave(filename = file.path("..", "Output", paste0("dist_diff_", v, ".png")),
         plot = p, width = 15, height = 15)
}


