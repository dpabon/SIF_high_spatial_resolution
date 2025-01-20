#!/usr/local/bin/Rscript
# ---------------------------------------------------------------------------- #
# #### analysis_s2het_variablecases.R ####
# ---------------------------------------------------------------------------- #
# Purpose: analyse the effect on s2het metric by changing parameters
# Project: sen4gpp_towercomp
# Authors: Gregory Duveiller
# ---------------------------------------------------------------------------- #


library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# import data
filepath <- '/Volumes/BGI/scratch/zmhamdi/for_greg'
#filepath <- '/Users/gduveiller/work/workspace/sen4gpp_towercomp/data/input_data/new_s2het'

# fig path
figpath <- "figs/s2_match_extended"
dir.create(figpath, showWarnings = F, recursive = T)


# function to get a given set of data for a given radius and process it
get_dat <- function(rad, iSite){
  
  file_name <- list.files(path = filepath, pattern = paste0(iSite, "_", rad, '.csv'), full.names = T)
  
  dat_t <- read_csv(file = file_name)
  
  if("nirv_nonzeros" %in% colnames(dat_t)){ dat_t <- rename(dat_t, nirv_non_zeros = "nirv_nonzeros")}
  if("ndvi_nonzeros" %in% colnames(dat_t)){ dat_t <- rename(dat_t, ndvi_non_zeros = "ndvi_nonzeros")}
  
  dat_t <- dat_t %>% rename_with(~gsub("_non_zeros$", "_n", .), ends_with("_non_zeros")) 
  
  dat_target <- dat_t %>% 
    pivot_longer(cols = matches(".*_(mean|std|n)"),
                 names_to = c("variable", "metric"), 
                 names_pattern = "(.*)_(.*)") %>%
    pivot_wider(names_from = "metric", values_from = "value") %>%
    mutate(mean = ifelse(variable != "ndvi", mean/10000, mean),
           std = ifelse(variable != "ndvi", std/10000, std)) %>%
    mutate(radius_length = as.numeric(rad),
           site = iSite)
  
}

# function to add a given radius to the table and calc S2 metric
add_radius <- function(rad, iSite){
  dat <- left_join(dat_centre, get_dat(rad, iSite), by = c("date", "variable", "site")) %>%
    mutate(F2T_match = 1 - (abs(mean - centre_value) + std)/2) %>%
    #filter(n > 0) %>%
    group_by(variable, radius_length, site) %>%
    mutate(n_max = max(n)) %>%
    ungroup() %>%
    mutate(n_prop = n / n_max) 
}


sitelist <- c("BE-Lcr", "DE-RuS", "DE-Hzd", "DE-Hai", "ES-LM1", "FR-LGt", "IT-BFt")
sitelist <- c("ES-LMa",  "ES-LM2", "ES-LM1")
dat_all <- NULL

for (iSite in sitelist) {

dat_centre <- get_dat("30", iSite) %>% 
  group_by(variable, radius_length) %>%
  mutate(n_max = max(n, na.rm = T)) %>%
  ungroup() %>%
  mutate(n_prop = n / n_max) %>%
  filter(n_prop == 1) %>%
  filter(std/mean < 0.2) %>%
  transmute(site, date, variable, centre_value = mean)


dat_site <- bind_rows(
  add_radius("100", iSite),
  add_radius("250", iSite),
  add_radius("500", iSite),
  add_radius("1000", iSite),
  add_radius("2000", iSite),
  add_radius("3000", iSite)
  )

dat_all <- bind_rows(dat_all, dat_site) 

}

# # TODO: smoothing...
# 
# # dates over where to calculate in time (should cover the period of interest)
# df_dates <- data.frame(date = as.Date(c(seq(as.Date("2018-01-01"), 
#                                             as.Date("2021-12-31"), 
#                                             by = 1)))) %>% tibble()
# 
# # set span parameter for the loess smoothing function
# span_loess_param <- 0.2
# 
# # apply smoother...
# df %>%
# group_by(site, variable) %>% 
#   arrange(site, date) %>% 
#   nest() %>%
#   mutate(smoothed_s2het_index = map(data, function(x){
#     data.frame(df_dates, 
#                s2het_index = loess(s2het_index~as.numeric(date), 
#                                    span = span_loess_param, data = x) %>% 
#                  predict(df_dates) %>% unlist())
#     


variable_order <- c("RED", "VRE_705", "VRE_740", "VRE_783",  "VRE_865", "NIR", "SWIR_1610", "SWIR_2190", "ndvi", "nirv")
dat_all <- dat_all %>% mutate(variable = factor(variable, levels = variable_order))


#### Look at changes in radius for a given variable ####
iSite <- "ES-LM1"; iVariable = "ndvi"
g1 <- ggplot(dat_all %>% filter(variable == iVariable, n_prop > 0.85, site == iSite)) +
  geom_point(aes(x = date, y = F2T_match, colour = F2T_match), show.legend = F) +
  geom_line(aes(x = date, y = F2T_match)) +
  facet_grid(factor(radius_length) ~ .) +
  scale_x_date("") +
  labs(title = paste0("Change in radius for site ", iSite),
       subtitle = paste("F2T_match metric based on", toupper(iVariable)))
    
ggsave(paste0("multirad_", iSite, ".png"), plot = g1, path = figpath, width = 9, height = 6)




#### Look at how the metric changes for different variables ####
iSite <- "ES-LM1"; iRadius <- 1000
var_sel <- c("RED", "VRE_740", "NIR", "SWIR_1610", "ndvi", "nirv")

g2 <- ggplot(dat_all %>% 
               filter(radius_length == iRadius, n_prop > 0.85, site == iSite) %>%
               filter(variable %in% var_sel)) + 
  geom_point(aes(x = date, y = F2T_match, colour = variable), show.legend = FALSE) +
  geom_line(aes(x = date, y = F2T_match), show.legend = FALSE) +
  #scale_y_continuous(limits = c(0.75, 1)) +
  facet_grid(variable~.) +
  scale_x_date("") +
  labs(title = paste0("Variations in input for F2T_match for site ", iSite),
       subtitle = paste("Radius length", iRadius, "m"))

ggsave(paste0("multivar_", iSite, ".png"), plot = g2, path = figpath, width = 9, height = 6)




iVariable <- "ndvi"
site_sel <- c("DE-Hai", "DE-Hzd", "DE-RuS", "ES-LM1", "IT-BFt")


g3 <- ggplot(dat_all %>% 
               filter(variable == iVariable, n_prop > 0.85) %>%
               filter(site %in% site_sel)) + 
  geom_line(aes(x = date, y = F2T_match), show.legend = FALSE) +
  scale_x_date("") + 
  facet_grid(site~radius_length) +
  labs(title = paste0("Various sites for ", toupper(iVariable), " with different radii"))

ggsave(paste0("multisite_", iVariable, ".png"), plot = g3, path = figpath, width = 9, height = 6)



#### Export subset for Richard Nair ####
dat_export <- "data/inter_data/for_export"
dir.create(dat_export)

dat <- dat_all %>% 
  filter(radius_length == 100) %>% 
  select(-n_max, -F2T_match)

write_csv(dat, paste0(dat_export, "/S2-sphet__100m__LasMajadas.csv"))

g_export <- ggplot(dat %>% filter(variable == 'ndvi'), aes(x = date)) + 
  geom_linerange(aes(ymin = mean - std, ymax = mean + std, colour = site)) + 
  geom_point(aes(y = mean, colour = site), size = 0.2) + 
  facet_wrap(~site, nc = 1) 

ggsave("S2-sphet__100m__LasMajadas.png", plot = g_export, path = dat_export, width = 6, height = 8)

