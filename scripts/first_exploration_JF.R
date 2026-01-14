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

vars_to_evaluate <- c("surftemp_mean", "bottemp_mean",
                      #"sensheatf_mean", "latentheatf_mean",
                      "strat_sum", "ice_sum")

##------------- read in data --------------------

# read in lake meta data
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
# read in meta variable data
vars_meta <- read.csv("../raw_data/variables_description.csv")

# list files
files <- list.files(file.path("..", "raw_data"), recursive = TRUE) |>
  data.frame() |> setNames("file") |> filter(file != "README.md") |>
  filter(!str_detect(file, "desktop.ini")) |>
  filter(!grepl(pattern = "lake_characteristics.*", x = file)) |>
  filter(!grepl(pattern = "variables_description", x = file)) |>
  filter(!grepl(pattern = "performance", x = file))

# read in files and combine to a single data.frame
dat <- lapply(files$file, function(f) {
  read.csv(file.path("..", "raw_data", f), header = TRUE) |>
    filter(name %in% vars_to_evaluate)}) |>
  reshape2::melt(id.vars = 1:8) |> select(-L1)

# rename scenarios
dat <- dat |> mutate(scenario = case_match(scenario,
                                        "picontrol" ~ "Picontrol",
                                        "historical" ~ "Historical",
                                        "ssp126" ~ "SSP1-2.6",
                                        "ssp370" ~ "SSP3-7.0",
                                        "ssp585" ~ "SSP5-8.5")) |>
  mutate(scenario = factor(scenario,
                           levels = c("Picontrol", "Historical", "SSP1-2.6",
                                      "SSP3-7.0", "SSP5-8.5")))


# temporary: filter out data before 1851 and 2021 for ssp scenarios
dat <- dat |> filter(year >= 1851)
# # temporary fix lake names
# dat <- dat |> mutate(lake = case_when(lake == "Mueggelsee" ~ "Muggelsee",
#                                       lake == "NohipaloMustjaerv" ~ "NohipaloMustjarv",
#                                       lake == "NohipaloValgejaerv" ~ "NohipaloValgejarv",
#                                       lake == "Paaijarvi" ~ "Paajarvi",
#                                       lake == "Scharmutzelsee" ~ "Scharmutzel",
#                                       lake == "TroutLake" ~ "Trout",
#                                       lake == "Vortsjaerv" ~ "Vortsjarv",
#                                       .default = lake))


## read in performacnce metrics
perf_files <- list.files(file.path("..", "raw_data", "performance"))

perf <- lapply(perf_files, function(f) {
  read_csv(file.path("..", "raw_data", "performance", f))}) |>
  reshape2::melt(id.vars = 1:5) |> select(-L1)

##-------- calculations/trends ----------------------------------

## improovement in rmse and difference between calibrated and uncalibrated
dat_diff_perf <- dat |>
  pivot_wider(names_from = cali, values_from = value) |>
  mutate(diff = uncalibrated - calibrated) |>
  group_by(model, scenario, lake, name) |>
  reframe(mean = mean(diff),
          median = median(diff))

imp_rmse <- perf |> filter(metric == "rmse") |>
  pivot_wider(names_from = cali, values_from = value) |>
  mutate(impr = uncalibrated - calibrated)

dat_diff_perf <- dat_diff_perf |>
  left_join(imp_rmse, by = c(lake = "lake", model = "model"))

## trends
# calculate slope, intercept and mean value of lm
dat_trend <- dat |>
  group_by(lake, model, gcm, cali, scenario, name) |>
  reframe(slope = tryCatch(coefficients(lm(value ~ year))[2],
                           error = function(e) NA),
          intercept = tryCatch(coefficients(lm(value ~ year))[1],
                               error = function(e) NA),
          mean_h = tryCatch(predict(lm(value ~ year),
                                    data.frame(year = median(year))),
                            error = function(e) NA)) # mean value at middle of time period
      

# calculate difference in slope, intercept and mean value between calibrated and uncalibrated
dat_trends_diff <- dat_trend |>
  pivot_longer(cols = c(slope, intercept, mean_h), names_to = "var_lm") |>
  pivot_wider(names_from = cali, values_from = value,
              id_cols = c(model, scenario, lake, gcm, name, var_lm)) |>
  mutate(diff = `uncalibrated` - calibrated) |>
  select(- calibrated, -`uncalibrated`)


# calculate comparison metrics R, bias, and variance of difference
dat_metr <- dat |> pivot_wider(names_from = cali, values_from = value) |>
  group_by(model, scenario, lake, gcm, name) |>
  reframe(R = cor(calibrated, uncalibrated),
          bias = mean(uncalibrated - calibrated, na.rm = TRUE),
          var = var(uncalibrated - calibrated, na.rm = TRUE)) |>
  pivot_longer(c(R, bias, var), names_to = "metr")

# ## variance decomposition takes some time so its saved as pre-calculated .RDS files
# 
# ## variance decomposition
# var_frac <- function(value, model, gcm, lake) {
#   dat <- data.frame(model = model, gcm = gcm,
#                     lake = lake, value = value)
#   an <- anova(lm(value ~ model * gcm * lake, data = dat))
#   tmp <- print(an)
#   totsst <- sum(an$`Sum Sq`)
#   sep <- an$`Sum Sq`
#   #sep <- c(sep[1:4], sum(sep[5:length(sep)]))/totsst
#   #sep <- data.frame(group = c(rownames(tmp)[1:4], "interactions"),
#   #                  frac = sep)
#   sep <- sep/totsst
#   sep <- data.frame(group = c(rownames(tmp)),
#                     frac = sep)
#   return(sep)
# }
# 
# 
# ## variance decomposition without with calibrated
# var_frac_c <- function(value, model, gcm, lake, cali) {
#   dat <- data.frame(model = model, gcm = gcm,
#                     lake = lake, cali = cali,
#                     value = value)
#   an <- anova(lm(value ~ model * gcm * lake * cali, data = dat))
#   tmp <- print(an)
#   totsst <- sum(an$`Sum Sq`)
#   sep <- an$`Sum Sq`
#   #sep <- c(sep[1:4], sum(sep[5:length(sep)]))/totsst
#   #sep <- data.frame(group = c(rownames(tmp)[1:4], "interactions"),
#   #                  frac = sep)
#   sep <- sep/totsst
#   sep <- data.frame(group = c(rownames(tmp)),
#                     frac = sep)
#   return(sep)
# }
# 
# 
# 
# # variance partitioning for the R and bias of cali and uncali
# frac_temp_diff <- dat_metr |>
#   mutate(model = as.factor(model),
#          gcm = as.factor(gcm),
#          scenario = as.factor(scenario),
#          lake = as.factor(lake)) |>
#   group_by(name, metr, scenario) |>
#   reframe(fracs = var_frac(value, model, gcm, lake)) |>
#   unpack(fracs)
# 
# saveRDS(frac_temp_diff, file.path("..", "derived_data", "var_decomp_metr.RDS"))
# 
# # variance partitioning for each year for each variable for difference between cali and uncali seperated
# frac_temp_mean <- dat |>
#   pivot_wider(names_from = cali, values_from = value) |>
#   mutate(diff = uncalibrated - calibrated) |>
#   mutate(model = as.factor(model),
#          gcm = as.factor(gcm),
#          lake = as.factor(lake)) |>
#   group_by(name, year, scenario) |>
#   reframe(fracs = var_frac(diff, model, gcm, lake)) |>
#   unpack(fracs)
# 
# saveRDS(frac_temp_mean, file.path("..", "derived_data", "var_decomp_diff_ts.RDS"))
# 
# 
# # vaiance partitioning for slope of linear model
# frac_lm <- dat_trend |>
#   mutate(model = as.factor(model),
#          gcm = as.factor(gcm),
#          lake = as.factor(lake)) |>
#   group_by(name, scenario, cali) |>
#   reframe(fracs = var_frac(slope, model, gcm, lake)) |>
#   unpack(fracs)
# 
# saveRDS(frac_lm, file.path("..", "derived_data", "var_decomp_lm.RDS"))

## load pre-calculated variance decomposition
# var decomposition for R and bias
var_dec_diff <- readRDS(file.path("..", "derived_data", "var_decomp_metr.RDS"))
var_dec_diff_i <- var_dec_diff |>
  mutate(group = ifelse(str_detect(group, ":"), "interaction", group)) |>
  group_by(name, metr, group) |> reframe(frac = sum(frac))

# var decomposition for the time series of difference between calibrated and uncalibrated
var_dec <- readRDS(file.path("..", "derived_data", "var_decomp_diff_ts.RDS"))
var_dec_i <- var_dec |>
  mutate(group = ifelse(str_detect(group, ":"), "interaction", group)) |>
  group_by(name, year, scenario, group) |> reframe(frac = sum(frac))

# variance decomposition for slope of the linear model
var_dec_lm <- readRDS(file.path("..", "derived_data", "var_decomp_lm.RDS"))
var_dec_lm_i <- var_dec_lm |>
  mutate(group = ifelse(str_detect(group, ":"), "interaction", group)) |>
  group_by(name, scenario, cali, group) |> reframe(frac = sum(frac))


##-----------  plots ----------------------------------

## performacnce
perf |> ggplot() + geom_density_ridges(aes(x = value, y = cali, fill = cali)) +
  facet_grid(.~metric, scales = "free") + thm + scale_fill_viridis_d("Status") +
  xlab("") + ylab("") + theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank())

perf |> ggplot() + geom_density_ridges(aes(x = value, y = cali, fill = cali)) +
  facet_grid(model~metric, scales = "free") + thm + scale_fill_viridis_d("Status") +
  xlab("") + ylab("") + theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank())

## improovement of performance
perf |> pivot_wider(names_from = cali, values_from = value) |>
  mutate(impr = uncalibrated - calibrated) |> ggplot() +
  geom_density_ridges(aes(x = impr, y = model, fill = model)) +
  facet_grid(metric~., scales = "free") + scale_fill_viridis_d("Model") + thm +
  xlab("Improovement") + ylab("") + theme(axis.text.y = element_blank(),
                                          axis.ticks.y = element_blank()) +
  xlim(-6, 6)

## improovement in rmse vs mean difference calibrated and uncalibrated
dat_diff_perf |> ggplot() + geom_point(aes(x = mean, y = impr)) +
  geom_abline(aes(intercept = 0, slope = 1), col = "grey") +
  facet_grid(scenario~name, scales = "free") + thm


## over time
p <- dat |>
  group_by(year, scenario, cali, name) |>
  reframe(mean = mean(value, na.rm = TRUE),
          median = median(value, na.rm = TRUE),
          q5 = quantile(value, 0.05, na.rm = TRUE),
          q25 = quantile(value, 0.25, na.rm = TRUE),
          q75 = quantile(value, 0.75, na.rm = TRUE),
          q95 = quantile(value, 0.95, na.rm = TRUE),
          min = min(value, na.rm = TRUE),
          max = max(value, na.rm = TRUE)) |>
  left_join(vars_meta[, c(1, 4)], by = c(name = "variable")) |>
  ggplot() +
  #geom_ribbon(aes(x = year, ymin = min, ymax = max, fill = cali), alpha = 0.35) +
  geom_ribbon(aes(x = year, ymin = q5, ymax = q95, fill = cali), alpha = 0.35) +
  geom_ribbon(aes(x = year, ymin = q25, ymax = q75, fill = cali), alpha = 0.35) +
  geom_line(aes(x = year, y = median, col = cali), lwd = 1.23) +
  facet_grid(plot_name~scenario, scales = "free", labeller = label_wrap_gen(23)) +
  scale_fill_viridis_d("") + scale_color_viridis_d("") + thm +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.text.y = element_text(size = 11)) + ylab("") + xlab("Year")

ggsave(file.path("..", "Output", "ts_all_vars.pdf"), p, width = 13, height = 9)

# plot difference
p <- dat |>
  pivot_wider(names_from = cali, values_from = value) |>
  mutate(diff = uncalibrated - calibrated) |>
  group_by(year, scenario, name) |>
  reframe(mean = mean(diff, na.rm = TRUE),
          median = median(diff, na.rm = TRUE),
          q5 = quantile(diff, 0.05, na.rm = TRUE),
          q25 = quantile(diff, 0.25, na.rm = TRUE),
          q75 = quantile(diff, 0.75, na.rm = TRUE),
          q95 = quantile(diff, 0.95, na.rm = TRUE),
          min = min(diff, na.rm = TRUE),
          max = max(diff, na.rm = TRUE)) |>
  left_join(vars_meta[, c(1, 4)], by = c(name = "variable")) |>
  ggplot() +
  #geom_ribbon(aes(x = year, ymin = min, ymax = max), alpha = 0.35) +
  geom_ribbon(aes(x = year, ymin = q5, ymax = q95), alpha = 0.35) +
  geom_ribbon(aes(x = year, ymin = q25, ymax = q75), alpha = 0.35) +
  geom_line(aes(x = year, y = median), lwd = 1.23) +
  facet_grid(plot_name~scenario, scales = "free", labeller = label_wrap_gen(23)) +
  scale_fill_viridis_d("") + scale_color_viridis_d() + thm +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.text.y = element_text(size = 11)) + ylab("") + xlab("Year")

ggsave(file.path("..", "Output", "ts_diff_all_vars.pdf"), p,  width = 13, height = 9)

# plot relative difference
p <- dat |>
  pivot_wider(names_from = cali, values_from = value) |>
  mutate(diff = uncalibrated - calibrated) |>
  group_by(year, scenario, name) |>
  reframe(mean = mean(diff/calibrated, na.rm = TRUE),
          median = median(diff/calibrated, na.rm = TRUE),
          q5 = quantile(diff/calibrated, 0.05, na.rm = TRUE),
          q25 = quantile(diff/calibrated, 0.25, na.rm = TRUE),
          q75 = quantile(diff/calibrated, 0.75, na.rm = TRUE),
          q95 = quantile(diff/calibrated, 0.95, na.rm = TRUE),
          min = min(diff/calibrated, na.rm = TRUE),
          max = max(diff/calibrated, na.rm = TRUE)) |>
  left_join(vars_meta[, c(1, 4)], by = c(name = "variable")) |>
  ggplot() +
  #geom_ribbon(aes(x = year, ymin = min, ymax = max), alpha = 0.35) +
  geom_ribbon(aes(x = year, ymin = q5, ymax = q95), alpha = 0.35) +
  geom_ribbon(aes(x = year, ymin = q25, ymax = q75), alpha = 0.35) +
  geom_line(aes(x = year, y = median), lwd = 1.23) +
  facet_grid(plot_name~scenario, scales = "free", labeller = label_wrap_gen(23)) +
  scale_fill_viridis_d("") + scale_color_viridis_d() + thm +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.text.y = element_text(size = 11)) + ylab("") + xlab("Year")

ggsave(file.path("..", "Output", "ts_diff_rel_all_vars.pdf"), p,  width = 13, height = 9)

# distribution of difference
dat |>
  pivot_wider(names_from = cali, values_from = value) |>
  mutate(diff = uncalibrated - calibrated) |>
  left_join(meta, by = c(lake = "Lake.Short.Name")) |>
  left_join(vars_meta[, c(1, 4)], by = c(name = "variable")) |>
  ggplot() + geom_violin(aes(y = diff, x = "", fill = model)) +
  facet_grid(plot_name~., scale = "free", labeller = label_wrap_gen(23)) +
  scale_fill_viridis_d() + thm +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  

## var diff over time for cali and uncali seperated
cols <- c(viridis::plasma(3), colorspace::desaturate(viridis::mako(5), 0.2))

var_dec$group <- factor(var_dec$group,
                        levels = unique(var_dec$group),
                        labels = unique(var_dec$group))

var_dec_i |> left_join(vars_meta[, c(1, 4)], by = c(name = "variable")) |>
  ggplot() +
  geom_area(aes(x = year, y = frac, fill = group),
            stat="identity", col = 1, lwd = 0.1) +
  facet_grid(plot_name~scenario, scales = "free_x",
             labeller = label_wrap_gen(23)) + thm +
  scale_fill_viridis_d() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.text.y = element_text(size = 11)) +
  ylab("Fraction of variance (-)") + xlab("Year")

p <- var_dec |> left_join(vars_meta[, c(1, 4)], by = c(name = "variable")) |>
  ggplot() +
  geom_area(aes(x = year, y = frac, fill = group),
                                stat="identity", col = 1, lwd = 0.1) +
  facet_grid(plot_name~scenario, scales = "free_x",
             labeller = label_wrap_gen(23)) + thm +
  scale_fill_manual("Factor", values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.text.y = element_text(size = 11)) +
  ylab("Fraction of variance (-)") + xlab("Year")

ggsave("../Output/var_decom_diff_ts.pdf", p, width = 13, height = 9)

# just look at importance of factor model
var_dec |> filter(group == "model") |> ggplot() +
  geom_line(aes(x = year, y = frac), lwd = 1.2) +
  facet_grid(name~scenario, scales = "free") + thm +
  scale_colour_viridis_d(end = 0.85, direction = -1)



## scatter plots

# just for mean surf & bot temp
dat |> filter(name %in% c("surftemp_mean", "bottemp_mean")) |>
  pivot_wider(names_from = cali, values_from = value) |>
  left_join(meta, by = c(lake = "Lake.Short.Name")) |>
  ggplot() +
  geom_point(aes(x = calibrated, y = uncalibrated, col = model), alpha = 0.5) +
  geom_abline(aes(slope = 1, intercept = 0), lty = "dashed") +
  facet_grid(kmcluster~name, scales = "free") +
  scale_color_viridis_d() + thm

# just for ice and strat
dat |> filter(name %in% c("strat_sum", "ice_sum")) |>
  pivot_wider(names_from = cali, values_from = value) |>
  left_join(meta, by = c(lake = "Lake.Short.Name")) |>
  ggplot() +
  geom_point(aes(x = calibrated, y = uncalibrated, col = model), alpha = 0.5) +
  geom_abline(aes(slope = 1, intercept = 0), lty = "dashed") +
  facet_grid(kmcluster~name, scales = "free") +
  scale_color_viridis_d() + thm

# combined with binned data
p1 <- dat |> filter(name %in% c("surftemp_mean", "bottemp_mean")) |>
  pivot_wider(names_from = cali, values_from = value) |>
  left_join(vars_meta[, c(1, 4)], by = c(name = "variable")) |>
  ggplot() + geom_hex(aes(x = calibrated, y = uncalibrated)) +
  geom_abline(aes(slope = 1, intercept = 0), lty = "dashed", size = 1.5, col = "grey") +
  facet_grid(plot_name~scenario, scales = "free",
             labeller = label_wrap_gen(23)) + thm +
  scale_fill_viridis_c("Count", trans = "log10", breaks = c(1, 100, 10000)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.text.y = element_text(size = 11)) +
  xlab("") + ylab("")

p2 <- dat |> filter(name %in% c("strat_sum", "ice_sum")) |>
  pivot_wider(names_from = cali, values_from = value) |>
  left_join(vars_meta[, c(1, 4)], by = c(name = "variable")) |>
  ggplot() + geom_hex(aes(x = calibrated, y = uncalibrated)) +
  geom_abline(aes(slope = 1, intercept = 0), lty = "dashed", size = 1.5, col = "grey") +
  facet_grid(plot_name~scenario, scales = "free",
             labeller = label_wrap_gen(23)) + thm +
  scale_fill_viridis_c("Count", trans = "log10", breaks = c(1, 100, 10000)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.text.y = element_text(size = 11),
        legend.key.size = unit(1.5,"line")) +
  xlab("") + ylab("")

p <- ggarrange(p1, p2, nrow = 2, common.legend = TRUE)

p <- annotate_figure(p, left = textGrob("Uncalibrated", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("Calibrated", gp = gpar(cex = 1.3)))

ggsave(file.path("..", "Output", "scatter_all_vars.pdf"), p, width = 13,
       height = 11, bg = "white")




## variance of the variable by lake type
dat |>
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


## R and bias for all
p <- dat |>
  pivot_wider(names_from = cali, values_from = value) |>
  group_by(model, scenario, lake, gcm, name) |>
  reframe(R = cor(calibrated, uncalibrated),
          bias = mean(uncalibrated - calibrated, na.rm = TRUE)) |>
  pivot_longer(6:7, names_to = "metr") |>
  left_join(vars_meta[, c(1, 4)], by = c(name = "variable")) |>
  ggplot() +
  geom_hline(aes(yintercept = 0), col = "grey42", lty = "dashed") +
  geom_violin(aes(x = plot_name, y = value, fill = model)) +
  facet_wrap(metr~plot_name, scales = "free") + thm +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_viridis_d(option = "H") + xlab("")

ggsave(file.path("..", "Output", "R_bias_all.pdf"), p, width = 13, height = 9)

## variance decomposition
var_dec_diff$group <- factor(var_dec_diff$group,
                             levels = unique(var_dec_diff$group),
                             labels = unique(var_dec_diff$group))

cols <- c(viridis::plasma(3), colorspace::desaturate(viridis::mako(5), 0.2))
p <- var_dec_diff |> mutate(interaction = grepl(":", group)) |>
  left_join(vars_meta[, c(1, 4)], by = c(name = "variable")) |>
  ggplot() + geom_col(aes(x = metr, y = frac, fill = group)) +
  facet_grid(plot_name~scenario, labeller = label_wrap_gen(23)) +
  thm + scale_fill_manual("factor", values = cols) + xlab("") +
  theme(strip.text.y = element_text(size = 11)) +
  ylab("Fraction of variance (-)")

ggsave(file.path("..", "Output", "var_Decomp_R_bias.pdf"), p, width = 13, height = 9)

## slope of linear model
dat_trend |> ggplot() + geom_density_ridges(aes(y = cali, x = slope, fill = cali)) +
  facet_grid(scenario~name, scale = "free") + thm + scale_fill_viridis_d()

## mean value of linear model
dat_trend |> ggplot() + geom_density_ridges(aes(y = cali, x = mean_h, fill = cali)) +
  facet_grid(scenario~name, scale = "free") + thm + scale_fill_viridis_d()

## dist of slope difference
# trend just for surf and bot temp
p1 <- dat_trends_diff |> filter(var_lm == "slope",
                               name %in% c("bottemp_mean", "surftemp_mean")) |>
  left_join(meta, by = c(lake = "Lake.Short.Name")) |>
  left_join(vars_meta[, c(1, 4)], by = c(name = "variable")) |>
  ggplot() + geom_density_ridges(aes(x = diff, y = kmcluster, fill = kmcluster)) +
  facet_grid(scenario ~ plot_name, labeller = label_wrap_gen(23)) +
  theme_pubr(base_size = 16) + grids() + xlim(-0.03, 0.03) +
  xlab("Temp. slope difference (°C/a)") +
  scale_fill_viridis_d("Lake type") + ylab("") +
  thm + theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              strip.text.y = element_text(size = 11))

#ggsave(file.path("..", "Output", "dist_slope_diff_temp.pdf"), p, width = 13, height = 9)




# trend just for strat dur and total heat
p2 <- dat_trends_diff |> filter(var_lm == "slope",
                               name %in% c("ice_sum", "strat_sum")) |>
  left_join(meta, by = c(lake = "Lake.Short.Name")) |>
  left_join(vars_meta[, c(1, 4)], by = c(name = "variable")) |>
  ggplot() + geom_density_ridges(aes(x = diff, y = kmcluster, fill = kmcluster)) +
  facet_grid(scenario ~ plot_name, scales = "free",
             labeller = label_wrap_gen(23)) +
  theme_pubr(base_size = 16) + grids() + xlim(-1.5, 1) +
  xlab("Dur. slope difference (d/a)") +
  scale_fill_viridis_d("Lake type") + ylab("") +
  thm + theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              strip.text.y = element_text(size = 11))

#ggsave(file.path("..", "Output", "dist_slope_diff_strat_heat.pdf"), p, width = 13, height = 9)

p <- ggarrange(p1 + theme(strip.background.y = element_blank(),
                          strip.text.y = element_blank(),
                          plot.margin = margin(t = 10, l = 10, b = 10, r = 0)),
               p2 + theme(plot.margin = margin(t = 10, l = 0, b = 10, r = 10)),
               ncol = 2, common.legend = TRUE)

ggsave("../Output/diff_slope_dist.pdf", p, width = 13, height = 9)


## dist of mean value difference
# trend just for surf and bot temp
p1 <- dat_trends_diff |> filter(var_lm == "mean_h",
                                name %in% c("bottemp_mean", "surftemp_mean")) |>
  left_join(meta, by = c(lake = "Lake.Short.Name")) |>
  left_join(vars_meta[, c(1, 4)], by = c(name = "variable")) |>
  ggplot() + geom_density_ridges(aes(x = diff, y = kmcluster, fill = kmcluster)) +
  facet_grid(scenario ~ plot_name, labeller = label_wrap_gen(23)) +
  theme_pubr(base_size = 16) + grids() + xlim(-10, 7.5) +
  xlab("Temp. mean value difference (°C)") +
  scale_fill_viridis_d("Lake type") + ylab("") +
  thm + theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              strip.text.y = element_text(size = 11))


# trend just for strat dur and total heat
p2 <- dat_trends_diff |> filter(var_lm == "mean_h",
                                name %in% c("strat_sum")) |>
  left_join(meta, by = c(lake = "Lake.Short.Name")) |>
  left_join(vars_meta[, c(1, 4)], by = c(name = "variable")) |>
  ggplot() + geom_density_ridges(aes(x = diff, y = kmcluster, fill = kmcluster)) +
  facet_grid(scenario ~ plot_name, scales = "free",
             labeller = label_wrap_gen(23)) +
  theme_pubr(base_size = 16) + grids() + xlim(-175, 175) +
  xlab("Dur. mean value difference (d)") +
  scale_fill_viridis_d("Lake type") + ylab("") +
  thm + theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              strip.text.y = element_text(size = 11))

# trend just for strat dur and total heat
p3 <- dat_trends_diff |> filter(var_lm == "mean_h",
                                name %in% c("ice_sum")) |>
  left_join(meta, by = c(lake = "Lake.Short.Name")) |>
  left_join(vars_meta[, c(1, 4)], by = c(name = "variable")) |>
  ggplot() + geom_density_ridges(aes(x = diff, y = kmcluster, fill = kmcluster)) +
  facet_grid(scenario ~ plot_name, scales = "free",
             labeller = label_wrap_gen(23)) +
  theme_pubr(base_size = 16) + grids() + xlim(-15, 15) +
  xlab("Dur. mean value difference (d)") +
  scale_fill_viridis_d("Lake type") + ylab("") +
  thm + theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              strip.text.y = element_text(size = 11))


p <- ggarrange(p1 + theme(strip.background.y = element_blank(),
                          strip.text.y = element_blank(),
                          plot.margin = margin(t = 10, l = 10, b = 10, r = 0)),
               p2 + theme(strip.background.y = element_blank(),
                          strip.text.y = element_blank(),
                          plot.margin = margin(t = 10, l = 10, b = 10, r = 0)),
               p3 + theme(plot.margin = margin(t = 10, l = 0, b = 10, r = 10)),
               ncol = 3, common.legend = TRUE, widths = c(2, 1, 1))

ggsave("../Output/diff_mean_h_dist.pdf", p, width = 13, height = 9)



# variance decompositioning for linear slope
var_dec_lm$group <- factor(var_dec_lm$group,
                             levels = unique(var_dec_lm$group),
                             labels = unique(var_dec_lm$group))
p <- var_dec_lm |>
  left_join(vars_meta[, c(1, 4)], by = c(name = "variable")) |>
  ggplot() + geom_col(aes(x = cali, y = frac, fill = group)) +
  facet_grid(scenario~plot_name, labeller = label_wrap_gen(23)) + thm +
  scale_fill_manual("factor", values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 11),
        strip.text.y = element_text(size = 11)) + xlab("") +
  ylab("Fraction of variance (-)")

ggsave(file.path("..", "Output", "var_Decomp_lm.pdf"), p, width = 13, height = 11)
