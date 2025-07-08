# This script loops over all aggregated model output files of the
# "ISIMIP lake model calibration" project and checks for naming consistency
# and unrealistic patterns in the data. It writes all findings to a report file
# and optionally it creates figures.

rm(list = ls())
graphics.off()

Sys.setenv(TZ = "UTC")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(data.table)
library(stringr)
library(ggplot2)

##### Settings -----
folder_box = "../raw_data/" # Point to the local root of the "ISIMIP lake model calibration data" Box folder

write_report = T # If FALSE, the report is instead printed to the console
report_file = "../Output/report_aggregated_output.txt" # If write_report == TRUE

make_plots = F
folder_plot = "../Output/" # if make_plots == TRUE
fig_width = 12
fig_height = 10.5

ignore_folders = "example data" # Folders in folder_box to ignore (automatically ignores empty folders)

# Things to check:
all_scen = c("historical", "picontrol", "ssp126", "ssp370", "ssp585")
lst_years = list(historical = c(1850, 2014),
                 picontrol = c(1850, 2100),
                 future = c(2015, 2100))
all_gcms = c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll")
all_lakes = c("Allequash", "Alqueva", "Annie", "Arendsee", "Argyle", "Biel", 
              "BigMuskellunge", "BlackOak", "Bosumtwi", "Bryrup", "BurleyGriffin", 
              "Chao", "Crystal", "CrystalBog", "Delavan", "Dickie", "Eagle", 
              "Ekoln", "Erken", "EsthwaiteWater", "FallingCreek", "Feeagh", 
              "Fish", "GreatPond", "Green", "Harp", "Hassel", "Hulun", "Kilpisjarvi", 
              "Kinneret", "Kivu", "Klicava", "Kuivajarvi", "Langtjern", "Laramie", 
              "LowerLakeZurich", "Mendota", "Monona", "Mozhaysk", "MtBold", 
              "Muggelsee", "Murten", "Neuchatel", "Ngoring", "NohipaloMustjarv", 
              "NohipaloValgejarv", "Okauchee", "Paajarvi", "Rappbode", "Rappbodep", 
              "Rimov", "Rotorua", "Sammamish", "Sau", "Scharmutzel", "Sparkling", 
              "Stechlin", "Sunapee", "Tahoe", "Taihu", "Tarawera", "Thun", 
              "Toolik", "Trout", "TroutBog", "TwoSisters", "Vendyurskoe", "Vortsjarv", 
              "Washington", "Windermere", "Wingra", "Zlutice", "Zurich")
all_metrics = c("bottemp_max", "bottemp_mean", "bottemp_min", "surftemp_max", 
                "surftemp_mean", "surftemp_min", "heat_mean", "ice_start", "ice_end", 
                "ice_sum", "latentheatf_mean", "mixeddepth_mean", "sensheatf_mean", 
                "strat_start", "strat_end", "strat_sum", "stratstrength_mean")
lst_minmax = list(surftemp = list(min = -0.5, max = 50, nas_allowed = F),
                  bottemp = list(min = -0.5, max = 50, nas_allowed = F),
                  strat = list(min = 0, max = 366, nas_allowed = T),
                  mixeddepth = list(min = 0, max = 501, nas_allowed = T), # 501 m is the maximum depth in the ISIMIP lakes
                  ice = list(min = 0, max = 366, nas_allowed = T),
                  heat = list(min = 0, max = Inf, nas_allowed = F),
                  latentheatf = list(min = -Inf, max = Inf, nas_allowed = F),
                  sensheatf = list(min = -Inf, max = Inf, nas_allowed = F),
                  stratstrength = list(min = -0.5, max = 50, nas_allowed = F))
col_names = c("year", "model", "scenario", "gcm", "lake", "cali", "name", "value")

##### Start loop -----
txt_report = paste0("Starting check aggregated output: ", Sys.time())
txt_report = append(txt_report, paste0(" "))

the_folders = list.files(folder_box, include.dirs = T)
the_folders = the_folders[dir.exists(file.path(folder_box, the_folders))]
the_folders = the_folders[!(the_folders %in% ignore_folders)]

plot_vars = sapply(str_split(all_metrics, "_"), function(x) x[-length(x)]) |>
  sapply(function(x) paste0(x, collapse = "_")) |> 
  unique()

if(write_report & !dir.exists(dirname(report_file))){
  dir.create(dirname(report_file))
}
if(make_plots & !dir.exists(folder_plot)){
  dir.create(folder_plot)
}

for(i in the_folders){
  the_subfolders = list.files(file.path(folder_box, i))
  if(length(the_subfolders) == 0L){
    txt_report = append(txt_report, paste0("Folder ", i, " is empty - skipped"))
    txt_report = append(txt_report, paste0(" "))
    next
  }
  
  txt_report = append(txt_report, paste0("=== Started check - Folder ", i, " ==="))
  
  # Check: are all scenarios present?
  if(!(all(all_scen %in% the_subfolders))){
    txt_report = append(txt_report, paste0("Folder ", i, ": missing scenarios: ",
                                           paste0(all_scen[!(all_scen %in% the_subfolders)], collapse = ", ")))
  }
  
  # Check: are there scenarios that are not part of all_scen?
  if(!(any(the_subfolders %in% all_scen))){
    txt_report = append(txt_report, paste0("Folder ", i, ": unknown scenarios: ",
                                           paste0(the_subfolders[!(the_subfolders %in% all_scen)], collapse = ", ")))
    txt_report = append(txt_report, paste0("= Break - Ending check folder ", i, " ="))
    next
  }
  
  # Track progress
  message("\nStarting folder: ", i)
  pb = txtProgressBar(min = 0, max = length(the_subfolders), style = 3)
  progress = 0
  for(j in the_subfolders){
    # In order to reduce the number of plots, we need to collapse calibration status, _mean/_min/_max, and GCMs
    # into one value for minimum, mean, and maximum before plotting. Create empty lists to fill with individual
    # data frames for later merging, aggregating, and plotting
    if(make_plots){
      lst_plot = vector("list", length = length(plot_vars))
      names(lst_plot) = plot_vars
      for(plt in plot_vars){
        lst_plot[[plt]] = list()
      }
    }
    
    the_files = list.files(file.path(folder_box, i, j))
    
    # Check: are there the right number of output files?
    if(length(the_files) != length(all_gcms) * 2){
      txt_report = append(txt_report, paste0("Folder ", i, ", subfolder ", j, ": incorrect number of files: ",
                                             length(the_files), ", should be ", length(all_gcms) * 2))
    }
    
    for(k in the_files){
      parsed_filename = str_split(k, "_")[[1]]
      
      # Check: is filename correct?
      if(length(parsed_filename) != 4L){
        txt_report = append(txt_report, paste0("File ", k, ": incorrect filename: ",
                                               "Should contain three underscores"))
      }
      if(parsed_filename[1] != i){
        txt_report = append(txt_report, paste0("File ", k, ": incorrect filename: ",
                                               "First part should be ", i, ", but is ", parsed_filename[1]))
      }
      if(parsed_filename[2] != j){
        txt_report = append(txt_report, paste0("File ", k, ": incorrect filename: ",
                                               "Second part should be ", j, ", but is ", parsed_filename[2]))
      }
      if(!(parsed_filename[3] %in% all_gcms)){
        txt_report = append(txt_report, paste0("File ", k, ": incorrect filename: ",
                                               "Third part (", parsed_filename[3], ") not part of standard model names"))
      }
      if(!(parsed_filename[4] %in% c("calibrated.csv", "uncalibrated.csv"))){
        txt_report = append(txt_report, paste0("File ", k, ": incorrect filename: ",
                                               "Fourth part (", parsed_filename[4], ") should be 'calibrated.csv' or 'uncalibrated.csv'"))
      }
      
      # Read file
      df = fread(file.path(folder_box, i, j, k))
      
      # Check: are columns correct?
      if(any(!(names(df) %in% col_names))){
        txt_report = append(txt_report, paste0("File ", k, ": ",
                                               "Unknown headers: ", paste0(names(df)[!(names(df) %in% col_names)], collapse = ", ")))
        txt_report = append(txt_report, paste0("= Break - Ending check file ", k, " ="))
        next
      }
      if(any(!(col_names %in% names(df)))){
        txt_report = append(txt_report, paste0("File ", k, ": ",
                                               "Missing headers: ", paste0(col_names[!(col_names %in% names(df))], collapse = ", ")))
        txt_report = append(txt_report, paste0("= Break - Ending check file ", k, " ="))
        next
      }
      
      # Check: are there NA values in columns other than "value"?
      the_cols = names(df)[names(df) != "value"]
      if(any(is.na(df[, ..the_cols]))){
        txt_report = append(txt_report, paste0("File ", k, ": ",
                                               "NA values detected in other columns than 'value'"))
        txt_report = append(txt_report, paste0("= Break - Ending check file ", k, " ="))
        next
      }
      
      
      # Check: are unique-column contents in line with the file name?
      if(length(unique(df$model)) != 1L){
        txt_report = append(txt_report, paste0("File ", k, ": ",
                                               "Multiple unique values in 'model' column"))
      }else if(unique(df$model) != parsed_filename[1]){
        txt_report = append(txt_report, paste0("Folder ", i, ", subfolder ", j, " , file ", k, ": ",
                                               "The 'model' column should be ", parsed_filename[1], " but is ", unique(df$model)))
      }
      if(length(unique(df$scenario)) != 1L){
        txt_report = append(txt_report, paste0("File ", k, ": ",
                                               "Multiple unique values in 'scenario' column"))
      }else if(unique(df$scenario) != parsed_filename[2]){
        txt_report = append(txt_report, paste0("Folder ", i, ", subfolder ", j, " , file ", k, ": ",
                                               "The 'scenario' column should be ", parsed_filename[2], " but is ", unique(df$scenario)))
      }
      if(length(unique(df$gcm)) != 1L){
        txt_report = append(txt_report, paste0("File ", k, ": ",
                                               "Multiple unique values in 'gcm' column"))
      }else if(unique(df$gcm) != parsed_filename[3]){
        txt_report = append(txt_report, paste0("Folder ", i, ", subfolder ", j, " , file ", k, ": ",
                                               "The 'gcm' column should be ", parsed_filename[3], " but is ", unique(df$gcm)))
      }
      calib_status = gsub(".csv", "", parsed_filename[4])
      if(length(unique(df$cali)) != 1L){
        txt_report = append(txt_report, paste0("File ", k, ": ",
                                               "Multiple unique values in 'cali' column"))
      }else if(unique(df$cali) != calib_status){
        txt_report = append(txt_report, paste0("Folder ", i, ", subfolder ", j, " , file ", k, ": ",
                                               "The 'cali' column should be ", calib_status, " but is ", unique(df$cali)))
      }
      
      # Check: Are all lakes present?
      if(any(!(unique(df$lake) %in% all_lakes))){
        txt_report = append(txt_report, paste0("File ", k, ": ",
                                               "Unknown lakes: ", paste0(unique(df$lake)[!(unique(df$lake) %in% all_lakes)], collapse = ", ")))
      }
      if(any(!(all_lakes %in% unique(df$lake)))){
        txt_report = append(txt_report, paste0("File ", k, ": ",
                                               "Missing lakes: ", paste0(all_lakes[!(all_lakes %in% unique(df$lake))], collapse = ", ")))
      }
      
      # Check: is the time period correct?
      if(parsed_filename[2] == "historical"){
        time_period = lst_years[["historical"]]
      }else if(parsed_filename[2] == "picontrol"){
        time_period = lst_years[["picontrol"]]
      }else if(str_detect(parsed_filename[2], "ssp")){
        time_period = lst_years[["future"]]
      }else{
        txt_report = append(txt_report, paste0("File ", k, ": ",
                                               "File name does not point to a valid scenario."))
        txt_report = append(txt_report, paste0("= Break - Ending check file ", k, " ="))
        next
      }
      if(min(df$year) != time_period[1]){
        txt_report = append(txt_report, paste0("File ", k, ": ",
                                               "Wrong starting date: should be ", time_period[1], ", but is ", min(df$year)))
      }
      if(max(df$year) != time_period[2]){
        txt_report = append(txt_report, paste0("File ", k, ": ",
                                               "Wrong ending date: should be ", time_period[2], ", but is ", max(df$year)))
      }
      
      # Check: are all metric names present?
      if(any(!(unique(df$name) %in% all_metrics))){
        txt_report = append(txt_report, paste0("File ", k, ": ",
                                               "Unknown metrics: ", paste0(unique(df$name)[!(unique(df$name) %in% all_metrics)], collapse = ", "), ". Not checking their values in next section."))
      }
      if(any(!(all_metrics %in% unique(df$name)))){
        txt_report = append(txt_report, paste0("File ", k, ": ",
                                               "Missing metrics: ", paste0(all_metrics[!(all_metrics %in% unique(df$name))], collapse = ", ")))
      }
      
      ### From here on, check the values themselves
      the_metrics = unique(df$name)
      for(var in the_metrics){
        df_var = df[name == var, .(year, lake, value)]
        
        var_base = str_split(var, "_")[[1]]
        var_base = paste0(var_base[-length(var_base)], collapse = "_")
        
        # Check: NA values in the values?
        if(any(is.na(df_var$value)) & !lst_minmax[[var_base]][["nas_allowed"]]){
          txt_report = append(txt_report, paste0("File ", k, ", metric ", var, ": ",
                                                 "NA values: ", round(sum(is.na(df_var$value)) / df_var[, .N] * 100, 2L), "%"))
        }
        
        # Check: are minimum and maximum bounds exceeded?
        if(var_base %in% names(lst_minmax)){
          if(min(df_var$value, na.rm = T) < lst_minmax[[var_base]][["min"]]){
            txt_report = append(txt_report, paste0("File ", k, ", metric ", var, ": ",
                                                   "Minimum value (", min(df_var$value, na.rm = T), ") is lower than lower bound (", lst_minmax[[var_base]][["min"]], ")"))
          }
          if(max(df_var$value, na.rm = T) > lst_minmax[[var_base]][["max"]]){
            txt_report = append(txt_report, paste0("File ", k, ", metric ", var, ": ",
                                                   "Maximum value (", max(df_var$value, na.rm = T), ") is higher than higher bound (", lst_minmax[[var_base]][["max"]], ")"))
          }
        }
        
        if(make_plots){
          df_plot = copy(df_var)
          df_plot[, gcm := unique(df$gcm)]
          df_plot[, type := last(str_split(var, "_")[[1]])]
          ind = length(lst_plot[[var_base]]) + 1L
          lst_plot[[var_base]][[ind]] = df_plot
        }
      }
    }
    
    if(make_plots){
      # Here aggregate the data frames and make the plots
      for(plt in plot_vars){
        df_plot_merged = rbindlist(lst_plot[[plt]])
        
        # Different aggregation for mean-min-max and start-end-sum
        if(any(unique(df_plot_merged$type) == "mean")){
          df_plot = df_plot_merged[type != "max" & type != "min", .(centre = mean(value, na.rm = T)),
                                   by = .(year, lake, gcm)]
          df_tmp = df_plot_merged[, .(upper = max(value, na.rm = T)), by = .(year, lake, gcm)]
          df_plot[, upper := df_tmp$upper]
          df_tmp = df_plot_merged[, .(lower = min(value, na.rm = T)), by = .(year, lake, gcm)]
          df_plot[, lower := df_tmp$lower]
        }else{
          if(!(any(unique(df_plot_merged$type) %in% c("start", "end", "sum")))){
            stop("If 'type' column has no max/min/mean, it must have 'start', 'end', and 'sum'!")
          }
          df_plot = df_plot_merged[type == "sum", .(centre = mean(value, na.rm = T)), by = .(year, lake, gcm)]
          df_tmp = df_plot_merged[type == "end", .(upper = mean(value, na.rm = T)), by = .(year, lake, gcm)]
          df_plot[, upper := df_tmp$upper]
          df_tmp = df_plot_merged[type == "start", .(lower = mean(value, na.rm = T)), by = .(year, lake, gcm)]
          df_plot[, lower := df_tmp$lower]
        }
        
        the_name = paste0(i, " - ", j)
        p = ggplot(df_plot) +
          geom_ribbon(aes(year, ymin = lower, ymax = upper, fill = gcm), alpha = 0.3) +
          geom_line(aes(year, centre, colour = gcm)) +
          labs(title = the_name, subtitle = plt) +
          facet_wrap(~ lake) +
          theme_light()
        ggsave(paste0(folder_plot, "Plot_", the_name, "_", plt, ".png"),
               plot = p,
               width = fig_width, height = fig_height)
      }
    }
    progress = progress + 1
    setTxtProgressBar(pb, progress)
  }
  
  txt_report = append(txt_report, paste0("=== Ended check - Folder ", i, " ==="))
  txt_report = append(txt_report, paste0(" "))
}

if(write_report){
  writeLines(txt_report, report_file)
}else{
  message(paste0(txt_report, collapse = "\n"))
}
