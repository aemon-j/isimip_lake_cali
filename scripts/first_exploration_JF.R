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
  filter(!str_detect(file, "desktop.ini")) |>
  filter(!grepl(pattern = "lake_characteristics.*", x = file)) |>
  filter(!grepl(pattern = "performance", x = file))

# read in files and combine to a single data.frame
dat <- lapply(files$file, function(f) {
  read.csv(file.path("..", "raw_data", f), header = TRUE)}) |>
  reshape2::melt(id.vars = 1:8) |> select(-L1)

dat <- dat |> mutate(scenario = case_match(scenario,
                                        "picontrol" ~ "Picontrol",
                                        "historical" ~ "Historical",
                                        "ssp126" ~ "SSP1-2.6",
                                        "ssp370" ~ "SSP3-7.0",
                                        "ssp585" ~ "SSP5-8.5")) |>
  mutate(scenario = factor(scenario,
                           levels = c("Picontrol", "Historical", "SSP1-2.6",
                                      "SSP3-7.0", "SSP5-8.5")))

meta <- read.csv("../raw_data/lake_characteristics.csv") |>
  mutate(kmcluster = factor(kmcluster,
                            levels = c("Deep lakes",
                                       "Medium temperate lakes",
                                       "Small temperate lakes",
                                       "Large shallow lakes",
                                       "Warm lakes"),
                            labels = c("Deep lakes",
                                       "Medium temperate lakes",
                                       "Small temperate lakes",
                                       "Large shallow lakes",
                                       "Warm lakes")))

##-------- calculations/trends ----------------------------------##

## trends

dat_trend <- dat |> pivot_wider() |>
  group_by(lake, model, cali, scenario) |>
  reframe(sl_surftemp_mean = coefficients(lm(surftemp_mean ~ year))[2],
          ic_surftemp_mean = coefficients(lm(surftemp_mean ~ year))[1],
          m_surftemp_mean = predict(lm(surftemp_mean ~ year), data.frame(year = median(year))), # mean value at middle of time period
          sl_bottemp_mean = tryCatch(coefficients(lm(bottemp_mean ~ year))[2], error = function(e) NA),
          ic_bottemp_mean = tryCatch(coefficients(lm(bottemp_mean ~ year))[1], error = function(e) NA),
          m_t_bottemp_mean = tryCatch(predict(lm(bottemp_mean ~ year), data.frame(year = median(year))), error = function(e) NA), # mean value at middle of time period
          sl_strat_mean = coefficients(lm(strat_sum ~ year))[2],
          ic_strat_mean = coefficients(lm(strat_sum ~ year))[1],
          m_strat_mean = predict(lm(strat_sum ~ year), data.frame(year = median(year))), # mean value at middle of time period
          sl_latentheatf_mean = tryCatch(coefficients(lm(latentheatf_mean ~ year))[2], error = function(e) NA),
          ic_latentheatf_mean = tryCatch(coefficients(lm(latentheatf_mean ~ year))[1], error = function(e) NA),
          m_latentheatf_mean = tryCatch(predict(lm(latentheatf_mean ~ year), data.frame(year = median(year))), error = function(e) NA),
          sl_sensheatf_mean = coefficients(lm(sensheatf_mean ~ year))[2],
          ic_sensheatf_mean = coefficients(lm(sensheatf_mean ~ year))[1],
          m_sensheatf_mean = predict(lm(sensheatf_mean ~ year), data.frame(year = median(year))),
          sl_heat_mean = tryCatch(coefficients(lm(heat_mean ~ year))[2], error = function(e) NA),
          ic_heat_mean = tryCatch(coefficients(lm(heat_mean ~ year))[1], error = function(e) NA),
          m_heat_mean = tryCatch(predict(lm(heat_mean ~ year), data.frame(year = median(year))), error = function(e) NA))


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
  dat <- data.frame(model = model, gcm = gcm,
                    scenario = scenario, lake = lake,
                    value = value) 
  an <- anova(lm(value ~ model * gcm * scenario * lake, data = dat))
  tmp <- print(an)
  totsst <- sum(an$`Sum Sq`)
  sep <- an$`Sum Sq`
  #sep <- c(sep[1:4], sum(sep[5:length(sep)]))/totsst
  #sep <- data.frame(group = c(rownames(tmp)[1:4], "interactions"),
  #                  frac = sep)
  sep <- sep/totsst
  sep <- data.frame(group = c(rownames(tmp)),
                    frac = sep)
  return(sep)
}

# # this is very slow due to the lm with all interactions :(
# # so I saved the outccome in a RDS file and commented the part to calculate it out
# frac_temp_mean <- dat |>
#   filter(name %in% c("surftemp_mean", "bottemp_mean",
#                      "sensheatf_mean", "latentheatf_mean",
#                      "strat_sum", "heat_mean")) |>
#   pivot_wider(names_from = cali, values_from = value) |>
#   group_by(model, scenario, lake, gcm, name) |>
#   reframe(R = cor(calibrated, uncalibrated),
#           bias = mean(uncalibrated - calibrated, na.rm = TRUE)) |>
#   pivot_longer(6:7, names_to = "metr") |>
#   mutate(model = as.factor(model),
#          gcm = as.factor(gcm),
#          scenario = as.factor(scenario),
#          lake = as.factor(lake)) |>
#   group_by(name, metr) |>
#   reframe(fracs = var_frac(value, model, gcm, scenario, lake)) |>
#   unpack(fracs)
# 
# saveRDS(frac_temp_mean, file.path("..", "derived_data", "var_decomp.RDS"))

# load pre-calculated variance decomposition
var_dec <- readRDS(file.path("..", "derived_data", "var_decomp.RDS"))


##-------- first plots just plot everything (this is messy) -----------------------##

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
  
  ggsave(file.path("..", "Output", "first_try",
                   paste0("diff_over_time_", v, ".png")),
         p, width = 13, height = 13)
  
}

## boxplots for metrics comparing cali and uncali

for(v in unique(dat$name)) {
p <- dat_metr |> left_join(meta, by = c(lake = "Lake.Short.Name")) |>
  filter(name == v) |> ggplot() +
  geom_boxplot(aes(x = kmcluster, y = value, fill = kmcluster)) + 
  facet_grid(metr~scenario, scales = "free") + scale_fill_viridis_d() + thm

ggsave(file.path("..", "Output", "first_try",
                 paste0("boxplot_", v, ".png")),
       p, width = 13, height = 13)
}

## scatter plots

# for all variables
for(v in unique(dat$name)) {
  p <- dat |> filter(name == v) |>
    pivot_wider(names_from = cali, values_from = value) |>
    left_join(meta, by = c(lake = "Lake.Short.Name")) |>
    ggplot() + geom_point(aes(x = calibrated, y = uncalibrated, col = model)) +
    geom_abline(aes(slope = 1, intercept = 0), lty = "dashed") +
    facet_wrap(scenario~kmcluster, scales = "free") +
    scale_color_viridis_d() + thm + ggtitle(v)
  
  ggsave(filename = file.path("..", "Output", "first_try",
                              paste0("scatter_", v, ".png")),
         plot = p, width = 17, height = 17)
}


# dist trend

for(v in unique(dat_trends_diff$name)) {
  p <- dat_trends_diff |> filter(name == v) |>
    left_join(meta, by = c(lake = "Lake.Short.Name")) |>
    ggplot() + geom_density_ridges(aes(x = diff, y = model, fill = model)) +
    facet_grid(scenario ~ kmcluster, scale = "free") +
    theme_pubr(base_size = 16) + grids() +
    xlab("temp slope difference (°C/a)") +
    ggtitle(paste0("Difference in estimated ", v),
            subtitle = "between uncalibirated and claibriated models (uncalibrated − calibrated)") +
    scale_fill_viridis_d("Model", option = "C", end = 0.9) + ylab("") +
    thm
  
  ggsave(filename = file.path("..", "Output", "first_try",
                              paste0("dist_diff_", v, ".png")),
         plot = p, width = 15, height = 15)
}


##----------- better plots ----------------------------------

## over time
dat |> filter(name %in% c("surftemp_mean", "bottemp_mean",
                          "sensheatf_mean", "latentheatf_mean",
                          "strat_sum", "heat_mean")) |>
  group_by(year, scenario, cali, name) |>
  reframe(mean = mean(value, na.rm = TRUE),
          median = median(value, na.rm = TRUE),
          q5 = quantile(value, 0.05, na.rm = TRUE),
          q25 = quantile(value, 0.25, na.rm = TRUE),
          q75 = quantile(value, 0.75, na.rm = TRUE),
          q95 = quantile(value, 0.95, na.rm = TRUE)) |>
  ggplot() +
  geom_ribbon(aes(x = year, ymin = q25, ymax = q75, fill = cali), alpha = 0.5) +
  geom_line(aes(x = year, y = median, col = cali), lwd = 1.23) +
  facet_grid(name~scenario, scales = "free") +
  scale_fill_viridis_d() + scale_color_viridis_d() + thm +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(file.path("..", "Output", "ts_all_vars.png"), width = 19, height = 17)

## scatter plots

# just for mean surf & bot temp
dat |> filter(name %in% c("surftemp_mean", "bottemp_mean")) |>
  pivot_wider(names_from = cali, values_from = value) |>
  left_join(meta, by = c(lake = "Lake.Short.Name")) |>
  ggplot() +
  geom_point(aes(x = calibrated, y = uncalibrated, col = kmcluster), alpha = 0.5) +
  geom_abline(aes(slope = 1, intercept = 0), lty = "dashed") +
  facet_grid(scenario~name, scales = "free") +
  scale_color_viridis_d() + thm

ggsave(file.path("..", "Output", "scatter_temp.png"), width = 13, height = 13)

# just for sensible and latent heat flux
dat |> filter(name %in% c("sensheatf_mean", "latentheatf_mean")) |>
  pivot_wider(names_from = cali, values_from = value) |>
  left_join(meta, by = c(lake = "Lake.Short.Name")) |>
  ggplot() +
  geom_point(aes(x = calibrated, y = uncalibrated, col = kmcluster), alpha = 0.5) +
  geom_abline(aes(slope = 1, intercept = 0), lty = "dashed") +
  facet_wrap(scenario~name, scales = "free") +
  scale_color_viridis_d() + thm

ggsave(file.path("..", "Output", "scatter_heatf.png"), width = 13, height = 13)

# just for total heat and strat
dat |> filter(name %in% c("strat_sum", "heat_mean")) |>
  pivot_wider(names_from = cali, values_from = value) |>
  left_join(meta, by = c(lake = "Lake.Short.Name")) |>
  ggplot() +
  geom_point(aes(x = calibrated, y = uncalibrated, col = kmcluster), alpha = 0.5) +
  geom_abline(aes(slope = 1, intercept = 0), lty = "dashed") +
  facet_wrap(scenario~name, scales = "free") +
  scale_color_viridis_d() + thm

ggsave(file.path("..", "Output", "scatter_heat_strat.png"), width = 13, height = 13)

## variance of the variable by lake type
dat |> filter(name %in% c("surftemp_mean", "bottemp_mean",
                          "sensheatf_mean", "latentheatf_mean",
                          "strat_sum", "heat_mean")) |>
  left_join(meta, by = c(lake = "Lake.Short.Name")) |>
  group_by(year, scenario, cali, name, kmcluster) |>
  reframe(range = abs(diff(range(value, na.rm = TRUE))),
          var = var(value, na.rm = TRUE),
          sd = sqrt(var)) |>
  ggplot() +
  geom_line(aes(x = year, y = sd, col = kmcluster, lty = cali), lwd = 1.23) +
  facet_grid(name~scenario, scales = "free") +
  scale_fill_viridis_d() + scale_color_viridis_d() + thm +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(file.path("..", "Output", "var_all_vars.png"), width = 19, height = 17)

## R and bias for all
dat |> filter(name %in% c("surftemp_mean", "bottemp_mean",
                          "sensheatf_mean", "latentheatf_mean",
                          "strat_sum", "heat_mean")) |>
    pivot_wider(names_from = cali, values_from = value) |>
    group_by(model, scenario, lake, gcm, name) |>
    reframe(R = cor(calibrated, uncalibrated),
            bias = mean(uncalibrated - calibrated, na.rm = TRUE)) |>
    pivot_longer(6:7, names_to = "metr") |> ggplot() +
  geom_hline(aes(yintercept = 0), col = "grey42", lty = "dashed") +
  geom_violin(aes(x = name, y = value, fill = name)) +
  facet_wrap(metr~name, scales = "free") + thm +
  scale_fill_viridis_d(option = "H")

ggsave(file.path("..", "Output", "R_bias_all.png"), width = 15, height = 13)

## variance decomposition
var_dec$group <- factor(var_dec$group,
                        levels = unique(var_dec$group),
                        labels = unique(var_dec$group))

cols <- c(viridis::plasma(4), colorspace::desaturate(viridis::mako(12), 0.2))
var_dec |> mutate(interaction = grepl(":", group)) |>
  ggplot() + geom_col(aes(x = metr, y = frac, fill = group)) +
  facet_grid(.~name) + thm + scale_fill_manual("factor", values = cols)

ggsave(file.path("..", "Output", "var_Decomp_R_bias.png"), width = 15, height = 13)

## dist of slope difference

# trend just for surf and bot temp
dat_trends_diff |> filter(name %in% c("sl_surftemp_mean", "sl_bottemp_mean")) |>
  left_join(meta, by = c(lake = "Lake.Short.Name")) |>
  ggplot() + geom_density_ridges(aes(x = diff, y = kmcluster, fill = kmcluster)) +
  facet_grid(scenario ~ name) +
  theme_pubr(base_size = 16) + grids() +
  xlab("temp slope difference (°C/a)") +
  scale_fill_viridis_d("Lake type") + ylab("") +
  thm + theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())

ggsave(file.path("..", "Output", "dist_slope_diff_temp.png"), width = 15, height = 13)

# trend just for sensible and latent heat flux
dat_trends_diff |> filter(name %in% c("sl_sensheatf_mean", "sl_latentheatf_mean")) |>
  left_join(meta, by = c(lake = "Lake.Short.Name")) |>
  ggplot() + geom_density_ridges(aes(x = diff, y = kmcluster, fill = kmcluster)) +
  facet_grid(scenario ~ name) +
  theme_pubr(base_size = 16) + grids() +
  xlab("slope difference (W/m^2/a)") +
  scale_fill_viridis_d("Lake type") + ylab("") +
  thm + theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())

ggsave(file.path("..", "Output", "dist_slope_diff_hflux.png"), width = 15, height = 13)


# trend just for strat dur and total heat
dat_trends_diff |> filter(name %in% c("sl_strat_mean", "sl_heat_mean")) |>
  left_join(meta, by = c(lake = "Lake.Short.Name")) |>
  ggplot() + geom_density_ridges(aes(x = diff, y = kmcluster, fill = kmcluster)) +
  facet_grid(scenario ~ name, scales = "free") +
  theme_pubr(base_size = 16) + grids() +
  xlab("slope difference") +
  scale_fill_viridis_d("Lake type") + ylab("") +
  thm + theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
ggsave(file.path("..", "Output", "dist_slope_diff_strat_heat.png"), width = 15, height = 13)
