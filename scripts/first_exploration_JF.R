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

# log transform the mean internal thermal heat
dat <- dat |> mutate(value = ifelse(grepl("air2water", model ) & name == "heat_mean",
                                    NA, value)) |>
  mutate(value = ifelse(name == "heat_mean", log10(value), value))

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
  group_by(lake, model, gcm, cali, scenario) |>
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
          sl_ice_mean = tryCatch(coefficients(lm(ice_sum ~ year))[2], error = function(e) NA),
          ic_ice_mean = tryCatch(coefficients(lm(ice_sum ~ year))[1], error = function(e) NA),
          m_ice_mean = tryCatch(predict(lm(ice_sum ~ year), data.frame(year = median(year))), error = function(e) NA))


dat_trends_diff <- dat_trend |> pivot_longer(cols = 6:ncol(dat_trend)) |>
  pivot_wider(names_from = cali, values_from = value,
              id_cols = c(model, scenario, lake, gcm, name)) |>
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

## variance decomposition without scenario
var_frac_2 <- function(value, model, gcm, lake) {
  dat <- data.frame(model = model, gcm = gcm,
                    lake = lake,
                    value = value) 
  an <- anova(lm(value ~ model * gcm * lake, data = dat))
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

## variance decomposition without scenario but with calibnrated
var_frac_3 <- function(value, model, gcm, lake, cali) {
  dat <- data.frame(model = model, gcm = gcm,
                    lake = lake, cali = cali,
                    value = value) 
  an <- anova(lm(value ~ model * gcm * lake * cali, data = dat))
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

## variance decomposition of slope of linear model
var_frac_lm <- function(value, model, gcm, lake, scenario, cali) {
  dat <- data.frame(model = model, gcm = gcm,
                    lake = lake, cali = cali,
                    scenario = scenario, value = value) 
  an <- anova(lm(value ~ model * gcm * lake * scenario * cali, data = dat))
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
# # variance partitioning for the R and bias of cali and uncali
# frac_temp_diff <- dat |>
#   filter(name %in% c("surftemp_mean", "bottemp_mean",
#                      "sensheatf_mean", "latentheatf_mean",
#                      "strat_sum", "ice_sum")) |>
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
# saveRDS(frac_temp_diff, file.path("..", "derived_data", "var_decomp_diff.RDS"))
# 
# # variance partitioning for each year for each variable for cali and uncali seperated
# frac_temp_mean <- dat |>
#   filter(name %in% c("surftemp_mean", "bottemp_mean",
#                      "sensheatf_mean", "latentheatf_mean",
#                      "strat_sum", "ice_sum")) |>
#   mutate(model = as.factor(model),
#          gcm = as.factor(gcm),
#          lake = as.factor(lake)) |>
#   group_by(name, cali, year, scenario) |>
#   reframe(fracs = var_frac_2(value, model, gcm, lake)) |>
#   unpack(fracs)
# 
# saveRDS(frac_temp_mean, file.path("..", "derived_data", "var_decomp.RDS"))
# 
# # vaiance partitioning for each year and variable with calibration as a factor
# frac_temp_mean_2 <- dat |>
#   filter(name %in% c("surftemp_mean", "bottemp_mean",
#                      "sensheatf_mean", "latentheatf_mean",
#                      "strat_sum", "ice_sum")) |>
#   mutate(model = as.factor(model),
#          gcm = as.factor(gcm),
#          lake = as.factor(lake),
#          cali = as.factor(cali)) |>
#   group_by(name, year, scenario) |>
#   reframe(fracs = var_frac_3(value, model, gcm, lake, cali)) |>
#   unpack(fracs)
# 
# saveRDS(frac_temp_mean_2, file.path("..", "derived_data", "var_decomp_2.RDS"))
# 
# # vaiance partitioning for slope of linear model
# frac_lm <- dat_trend |> pivot_longer(6:23) |>
#   filter(name %in% c("sl_surftemp_mean", "sl_bottemp_mean",
#                      "sl_sensheatf_mean", "sl_latentheatf_mean",
#                      "sl_strat_mean", "ls_ice_mean")) |>
#   mutate(model = as.factor(model),
#          gcm = as.factor(gcm),
#          lake = as.factor(lake),
#          scenario = as.factor(scenario),
#          cali = as.factor(cali)) |>
#   group_by(name) |>
#   reframe(fracs = var_frac_lm(value, model, gcm, lake, scenario, cali)) |>
#   unpack(fracs)
# 
# saveRDS(frac_lm, file.path("..", "derived_data", "var_decomp_lm.RDS"))

# load pre-calculated variance decomposition
var_dec_diff <- readRDS(file.path("..", "derived_data", "var_decomp_diff.RDS"))
var_dec_diff_i <- var_dec_diff |>
  mutate(group = ifelse(str_detect(group, ":"), "interaction", group)) |>
  group_by(name, metr, group) |> reframe(frac = sum(frac))

var_dec <- readRDS(file.path("..", "derived_data", "var_decomp.RDS"))
var_dec_i <- var_dec |>
  mutate(group = ifelse(str_detect(group, ":"), "interaction", group)) |>
  group_by(name, cali, year, scenario, group) |> reframe(frac = sum(frac))

var_dec_2 <- readRDS(file.path("..", "derived_data", "var_decomp_2.RDS"))
var_dec_2_i <- var_dec_2 |>
  mutate(group = ifelse(str_detect(group, ":"), "interaction", group)) |>
  group_by(name, year, scenario, group) |> reframe(frac = sum(frac))

var_dec_lm <- readRDS(file.path("..", "derived_data", "var_decomp_lm.RDS"))
var_dec_lm_i <- var_dec_lm |>
  mutate(group = ifelse(str_detect(group, ":"), "interaction", group)) |>
  group_by(name, group) |> reframe(frac = sum(frac))


##------------------- first plots just plot everything (this is messy) -----------------------

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
p <- dat |> filter(name %in% c("surftemp_mean", "bottemp_mean",
                          "sensheatf_mean", "latentheatf_mean",
                          "strat_sum", "ice_sum")) |>
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

ggsave(file.path("..", "Output", "ts_all_vars.png"), p, width = 19, height = 17)

## var dist over time for cali and uncali seperated
var_dec$group <- factor(var_dec$group,
                        levels = unique(var_dec$group),
                        labels = unique(var_dec$group))

var_dec_i |> filter(name == "ice_sum") |> ggplot() + geom_area(aes(x = year, y = frac, fill = group),
                                stat="identity", col = 1, lwd = 0.1) +
  facet_grid(cali~scenario, scales = "free_x") + thm + scale_fill_viridis_d() +
  ylab("Fraction of variance (-)") + xlab("Year")

# just look at importance of factor model
var_dec |> filter(group == "model") |> ggplot() +
  geom_line(aes(x = year, y = frac, col = cali)) +
  facet_grid(name~scenario, scales = "free") + thm 

## var dist over time for cali as own factor
var_dec_2_i$group <- factor(var_dec_2_i$group,
                            levels = unique(var_dec_2_i$group),
                            labels = unique(var_dec_2_i$group))

p <- var_dec_2_i |> ggplot() +
  geom_area(aes(x = year, y = frac, fill = group),
            stat="identity", col = 1, lwd = 0.1) +
  facet_grid(name~scenario, scales = "free_x") + thm +
  scale_fill_viridis_d() +
  ylab("Fraction of variance (-)") + xlab("Year")

ggsave(file.path("..", "Output", "var_frac_all.png"), p, width = 19, height = 17)

# just look at importance of factor calibrated
var_dec_2 |> filter(group == "cali") |> ggplot() +
  geom_line(aes(x = year, y = frac, col = scenario), lwd = 1.2) +
  facet_grid(name~., scales = "free") + thm +
  scale_color_viridis_d()


## scatter plots

# just for mean surf & bot temp
p <- dat |> filter(name %in% c("surftemp_mean", "bottemp_mean")) |>
  pivot_wider(names_from = cali, values_from = value) |>
  left_join(meta, by = c(lake = "Lake.Short.Name")) |>
  ggplot() +
  geom_point(aes(x = calibrated, y = uncalibrated, col = model), alpha = 0.5) +
  geom_abline(aes(slope = 1, intercept = 0), lty = "dashed") +
  facet_grid(kmcluster~name, scales = "free") +
  scale_color_viridis_d() + thm

ggsave(file.path("..", "Output", "scatter_temp.png"), p, width = 13, height = 13)

# just for sensible and latent heat flux
p <- dat |> filter(name %in% c("sensheatf_mean", "latentheatf_mean")) |>
  pivot_wider(names_from = cali, values_from = value) |>
  left_join(meta, by = c(lake = "Lake.Short.Name")) |>
  ggplot() +
  geom_point(aes(x = calibrated, y = uncalibrated, col = model), alpha = 0.5) +
  geom_abline(aes(slope = 1, intercept = 0), lty = "dashed") +
  facet_grid(kmcluster~name, scales = "free") +
  scale_color_viridis_d() + thm

ggsave(file.path("..", "Output", "scatter_heatf.png"), p, width = 13, height = 13)

# just for total heat and strat
p <- dat |> filter(name %in% c("strat_sum", "ice_sum")) |>
  pivot_wider(names_from = cali, values_from = value) |>
  left_join(meta, by = c(lake = "Lake.Short.Name")) |>
  ggplot() +
  geom_point(aes(x = calibrated, y = uncalibrated, col = model), alpha = 0.5) +
  geom_abline(aes(slope = 1, intercept = 0), lty = "dashed") +
  facet_grid(kmcluster~name, scales = "free") +
  scale_color_viridis_d() + thm

ggsave(file.path("..", "Output", "scatter_heat_strat.png"), p, width = 13, height = 13)

## variance of the variable by lake type
p <- dat |> filter(name %in% c("surftemp_mean", "bottemp_mean",
                          "sensheatf_mean", "latentheatf_mean",
                          "strat_sum", "ice_sum")) |>
  left_join(meta, by = c(lake = "Lake.Short.Name")) |>
  group_by(year, model, scenario, cali, name, kmcluster) |>
  reframe(range = abs(diff(range(value, na.rm = TRUE))),
          var = var(value, na.rm = TRUE),
          sd = sqrt(var)) |>
  pivot_wider(names_from = cali, values_from = c(var, range, sd)) |>
  ggplot() +
  geom_point(aes(x = year, y = sd_calibrated - sd_uncalibrated, col = model), lwd = 1.23) +
  facet_grid(name~kmcluster, scales = "free") +
  scale_fill_viridis_d() + scale_color_viridis_d() + thm +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(file.path("..", "Output", "var_all_vars.png"), p, width = 19, height = 17)

## R and bias for all
p <- dat |> filter(name %in% c("surftemp_mean", "bottemp_mean",
                          "sensheatf_mean", "latentheatf_mean",
                          "strat_sum", "ice_sum")) |>
    pivot_wider(names_from = cali, values_from = value) |>
    group_by(model, scenario, lake, gcm, name) |>
    reframe(R = cor(calibrated, uncalibrated),
            bias = mean(uncalibrated - calibrated, na.rm = TRUE)) |>
    pivot_longer(6:7, names_to = "metr") |> ggplot() +
  geom_hline(aes(yintercept = 0), col = "grey42", lty = "dashed") +
  geom_violin(aes(x = name, y = value, fill = model)) +
  facet_wrap(metr~name, scales = "free") + thm +
  scale_fill_viridis_d(option = "H")

ggsave(file.path("..", "Output", "R_bias_all.png"), p, width = 15, height = 13)

## variance decomposition
var_dec_diff$group <- factor(var_dec_diff$group,
                             levels = unique(var_dec_diff$group),
                             labels = unique(var_dec_diff$group))

cols <- c(viridis::plasma(4), colorspace::desaturate(viridis::mako(12), 0.2))
p <- var_dec_diff |> mutate(interaction = grepl(":", group)) |>
  ggplot() + geom_col(aes(x = metr, y = frac, fill = group)) +
  facet_grid(.~name) + thm + scale_fill_manual("factor", values = cols)

ggsave(file.path("..", "Output", "var_Decomp_R_bias.png"), p, width = 15, height = 13)

## dist of slope difference

# trend just for surf and bot temp
p <- dat_trends_diff |> filter(name %in% c("sl_surftemp_mean", "sl_bottemp_mean")) |>
  left_join(meta, by = c(lake = "Lake.Short.Name")) |>
  ggplot() + geom_density_ridges(aes(x = diff, y = kmcluster, fill = kmcluster)) +
  facet_grid(scenario ~ name) +
  theme_pubr(base_size = 16) + grids() +
  xlab("temp slope difference (°C/a)") +
  scale_fill_viridis_d("Lake type") + ylab("") +
  thm + theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())

ggsave(file.path("..", "Output", "dist_slope_diff_temp.png"), p, width = 15, height = 13)

# trend just for sensible and latent heat flux
p <- dat_trends_diff |> filter(name %in% c("sl_sensheatf_mean", "sl_latentheatf_mean")) |>
  left_join(meta, by = c(lake = "Lake.Short.Name")) |>
  ggplot() + geom_density_ridges(aes(x = diff, y = kmcluster, fill = kmcluster)) +
  facet_grid(scenario ~ name) +
  theme_pubr(base_size = 16) + grids() +
  xlab("slope difference (W/m^2/a)") +
  scale_fill_viridis_d("Lake type") + ylab("") +
  thm + theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())

ggsave(file.path("..", "Output", "dist_slope_diff_hflux.png"), p, width = 15, height = 13)


# trend just for strat dur and total heat
p <- dat_trends_diff |> filter(name %in% c("sl_strat_mean", "sl_heat_mean")) |>
  left_join(meta, by = c(lake = "Lake.Short.Name")) |>
  ggplot() + geom_density_ridges(aes(x = diff, y = kmcluster, fill = kmcluster)) +
  facet_grid(scenario ~ name, scales = "free") +
  theme_pubr(base_size = 16) + grids() +
  xlab("slope difference") +
  scale_fill_viridis_d("Lake type") + ylab("") +
  thm + theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
ggsave(file.path("..", "Output", "dist_slope_diff_strat_heat.png"), p, width = 15, height = 13)


# variance decompositioning for linear slope
var_dec_lm_i$group <- factor(var_dec_lm_i$group,
                             levels = unique(var_dec_lm_i$group),
                             labels = unique(var_dec_lm_i$group))
p <- var_dec_lm_i |>
  ggplot() + geom_col(aes(x = "", y = frac, fill = group)) +
  facet_grid(.~name) + thm + scale_fill_viridis_d()

ggsave(file.path("..", "Output", "var_Decomp_R_bias.png"), p, width = 15, height = 13)