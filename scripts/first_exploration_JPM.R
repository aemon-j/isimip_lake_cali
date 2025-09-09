# Exploring ISIMIP calibrated vs. uncalibrated results
# Focus on "condensing" the results to understandable tables or graphs

rm(list = ls())
graphics.off()

Sys.setenv(TZ = "UTC")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(data.table)
library(ggplot2)
library(ggtext)
library(stringr)
library(lubridate)
library(LakeEnsemblR)

##### Settings -----
folder_box = "../raw_data/"
folder_out = "../Output/"

vars = list(surftemp_mean = list(longname = "Mean surface temperature",
                                 unit = "&deg;C"),
            strat_sum = list(longname = "Stratification duration",
                             unit = "d"))

capitalisation = list("air2water-4par" = "air2water-4par",
                      "air2water-6par" = "air2water-6par",
                      "flake-ler" = "FLake-LER",
                      "glm" = "GLM",
                      "glm-ler" = "GLM-LER", 
                      "gotm" = "GOTM",
                      "gotm-ler" = "GOTM-LER",
                      "simstrat" = "Simstrat",
                      "simstrat-ler" = "Simstrat-LER")

fig_type = "png"
fig_width = 9
fig_height = 7

# Functions
rmse = function(sim, obs, na.rm = TRUE){
  sqrt(mean((sim - obs)^2, na.rm = na.rm))
}
rsq = function(sim, obs, na.rm = TRUE){
  if(na.rm){
    use_arg = "na.or.complete"
  }else{
    use_arg = "everything"
  }
  cor(sim, obs, use = use_arg)^2
}


##### Loading data -----
models = list.files(folder_box)
models = models[dir.exists(file.path(folder_box, models))]
models = models[models %notin% c("example data", "performance")]

### Load performance data
lst_performance = list()
subfolder = file.path(folder_box, "performance")
for(i in models){
  filename = paste0(subfolder, "/", capitalisation[[i]], "_LHS_performance.csv")
  if(!file.exists(filename)) next
  lst_performance[[length(lst_performance) + 1]] = fread(filename)
}
df_performance = rbindlist(lst_performance)
df_performance = dcast(df_performance, model + lake + metric ~ cali)

### Load simulated data
for(i in names(vars)){
  lst_var = list()
  
  pb = txtProgressBar(min = 0, max = length(models), style = 3)
  progress = 0
  for(j in models){
    subfolders = list.files(file.path(folder_box, j))
    for(k in subfolders){
      the_files = list.files(file.path(folder_box, j, k), pattern = "\\.csv")
      for(l in the_files){
        df = fread(file.path(folder_box, j, k, l))
        lst_var[[length(lst_var) + 1]] = df[name == i]
      }
    }
    
    progress = progress + 1
    setTxtProgressBar(pb, progress)
  }
  
  df_var = rbindlist(lst_var)
  df_var_wide = dcast(df_var, year + model + scenario + gcm + lake ~ cali)
  
  # Plot 1: effect of calibration
  df_plot = df_var_wide[, .(meandiff = mean(calibrated - uncalibrated)),
                        by = .(model, scenario, gcm, lake)]
  
  title = "Mean difference calibrated and uncalibrated"
  subtitle = vars[[i]][["longname"]]
  ylab = paste0("Mean difference (", vars[[i]][["unit"]], ")")
  
  p1 = ggplot(df_plot) +
    geom_violin(aes(scenario, meandiff, fill = gcm)) +
    geom_hline(aes(yintercept = 0.0), linetype = "dashed") +
    facet_wrap(~ model) +
    labs(title = title, subtitle = subtitle, y = ylab) +
    theme_light() +
    theme(axis.title = element_markdown(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(folder_out, "Plot_meandiff_", i, ".", fig_type), plot = p1,
         width = fig_width, height = fig_height)
  
  # Plot 2: plot improvement in performance vs cal-uncal difference
  title = "Mean difference vs. performance"
  xlab = paste0("Difference in RMSE (&deg;C)")
  ylab = paste0("Abs. difference (", vars[[i]][["unit"]], ")")
  
  df_perf_rmse = df_performance[metric == "rmse"]
  df_perf_rmse[, diff_perf := uncalibrated - calibrated]
  df_plot2 = merge(df_plot, df_perf_rmse, by = c("model", "lake"))
  df_plot2[, meandiff := abs(meandiff)]
  df_plot2 = df_plot2[scenario == "historical" & gcm == "gfdl-esm4"] # No need to plot all GCMs and scenarios
  p2 = ggplot(df_plot2) +
    geom_point(aes(diff_perf, meandiff, colour = calibrated)) +
    facet_wrap(~ model) +
    labs(title = title, subtitle = subtitle, x= xlab, y = ylab) +
    theme_light() +
    theme(axis.title = element_markdown())
  ggsave(paste0(folder_out, "Plot_PerformanceVsDiff_", i, ".", fig_type), plot = p2,
         width = fig_width, height = fig_height)
  
  df_rsq = df_plot2[, .(rsq = rsq(meandiff, diff_perf)), by = model]
  
  p2a = ggplot(df_rsq) +
    geom_bar(aes(model, rsq, fill = model), stat = "identity") +
    labs(title = "R-squared - diff performance vs. bias",
         subtitle = subtitle,
         y = "R<sup>2</sup>") +
    theme_light() +
    theme(axis.title = element_markdown())
  ggsave(paste0(folder_out, "Plot_Bar_R2bias_", i, ".", fig_type), plot = p2a,
         width = fig_width, height = fig_height)
  
  # Plot 3: 
  # Fit linear model for each lake/cali combination and assess difference in slope
  df_plot = df_var[, .(slope = coefficients(lm(value ~ year))[2],
                       middle_point = predict(lm(value ~ year), .(year = median(year)))),
                   by = .(model, scenario, gcm, lake, cali)]
  df_plot_wide1 = dcast(df_plot, model + scenario + gcm + lake ~ cali, value.var = "slope")
  df_plot_wide2 = dcast(df_plot, model + scenario + gcm + lake ~ cali, value.var = "middle_point")
  df_plot_wide1[, meandiff := calibrated - uncalibrated]
  df_plot_wide2[, meandiff := calibrated - uncalibrated]
  
  title = "Different slope calibrated and uncalibrated"
  subtitle = vars[[i]][["longname"]]
  ylab = paste0(subtitle, " (", vars[[i]][["unit"]], ")")
  p3a = ggplot(df_plot_wide1) +
    geom_violin(aes(scenario, meandiff, fill = gcm)) +
    geom_hline(aes(yintercept = 0.0), linetype = "dashed") +
    facet_wrap(~ model) +
    labs(title = title, subtitle = subtitle, y = ylab) +
    theme_light() +
    theme(axis.title = element_markdown(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(folder_out, "Plot_diffslope_", i, ".", fig_type), plot = p3a,
         width = fig_width, height = fig_height)
  
  title = "Different middle point linear model calibrated and uncalibrated"
  p3b = ggplot(df_plot_wide2) +
    geom_violin(aes(scenario, meandiff, fill = gcm)) +
    geom_hline(aes(yintercept = 0.0), linetype = "dashed") +
    facet_wrap(~ model) +
    labs(title = title, subtitle = subtitle, y = ylab) +
    theme_light() +
    theme(axis.title = element_markdown(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(folder_out, "Plot_diffmidpoint_", i, ".", fig_type), plot = p3b,
         width = fig_width, height = fig_height)
  
  # Plot 4:
  # Similar correlation as for plot 2, but now with the slope
  title = "Slope difference vs. performance"
  xlab = paste0("Difference in RMSE (&deg;C)")
  ylab = paste0("Slope difference (", vars[[i]][["unit"]], "/yr)")
  
  df_plot4 = merge(df_plot_wide1, df_perf_rmse, by = c("model", "lake"))
  df_plot4[, meandiff := abs(meandiff)]
  df_plot4 = df_plot4[scenario == "historical" & gcm == "gfdl-esm4"] # No need to plot all GCMs and scenarios
  p4 = ggplot(df_plot4) +
    geom_point(aes(diff_perf, meandiff)) +
    facet_wrap(~ model) +
    labs(title = title, subtitle = subtitle, x= xlab, y = ylab) +
    theme_light() +
    theme(axis.title = element_markdown())
  ggsave(paste0(folder_out, "Plot_PerformanceVsSlopeDiff_", i, ".", fig_type), plot = p4,
         width = fig_width, height = fig_height)
  
  df_rsq = df_plot4[, .(rsq = rsq(meandiff, diff_perf)), by = model]
  
  p4a = ggplot(df_rsq) +
    geom_bar(aes(model, rsq, fill = model), stat = "identity") +
    labs(title = "R-squared - diff performance vs. slope diff",
         subtitle = subtitle,
         y = "R<sup>2</sup>") +
    theme_light() +
    theme(axis.title = element_markdown())
  ggsave(paste0(folder_out, "Plot_Bar_R2slope_", i, ".", fig_type), plot = p4a,
         width = fig_width, height = fig_height)
}
