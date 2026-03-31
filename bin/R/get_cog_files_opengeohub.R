library(tidyverse)


setwd("/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/OEMC_app_data/COG_SIF_files/")


files_df <- read_csv("/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/OEMC_app_data/mpg_monitor_23.csv")

# for (i in 1:nrow(files_df)) {
#   system(paste("wget", files_df[i,1]))
# }


nrow(files_df)
  