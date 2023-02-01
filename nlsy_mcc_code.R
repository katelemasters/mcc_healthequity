### load packages ###
library(tidyverse)
library(broom)
library(survival)
library(gridExtra)

### set working directory and load data ###
nlsy = read_csv('nlsy_samplecode.csv') 

### create a loop for MCC ###
mcc_weighted = function(data_nested) {
  n0 = sum(data_nested$baseline_wt)
  data_nested %>%
    unnest(cols = c(data)) %>%
    group_by(time) %>%
    summarise(e_j = sum(event*baseline_wt),
              r_jplus1 = sum(lead.death*baseline_wt),
              c_jplus1_mcc = sum(lead.censor_mcc*baseline_wt),
              c_jplus1_km  = sum(lead.censor_km*baseline_wt)) %>%
    mutate(r_j = lag(r_jplus1, default = 0),
           c_j_mcc = lag(c_jplus1_mcc, default = 0), 
           c_j_km  = lag(c_jplus1_km, default = 0),
           n_jminus1_mcc = n0 - cumsum(r_j +c_j_mcc),
           n_jminus1_km  = n0 - cumsum(r_j + c_j_km),
           km = cumprod(1 - r_j/n_jminus1_km),
           mcc = cumsum(e_j/n_jminus1_mcc)*km) %>%
    dplyr::select(time, mcc)
}


### ages of interest ###
times_est = c(20,24,30,33)

### create nested data by individual ###
nlsy_nested = nlsy %>% 
  # group by individual
  group_by(id) %>%
  # event type = 2 is stressful life event
  # remove if not surveyed in waves that stressful life events were asked
  mutate(event = 1*(event_type==2),
         event_measured_this_wave = year %in% c(2002, 2007, 2008, 2009, 2013),
         no_data = !any(event_measured_this_wave)) %>%
  filter(!no_data) %>% 
  mutate(last_event_year = max(year[event_measured_this_wave]),
         #separate censoring variable for the repeated event
         #separate censoring variable for the competing event
         lead.censor_mcc = 1*(year == last_event_year), 
         lead.censor_km  = 1*(year == max(year)),
         raceth = case_when(RaceEthf==0 ~ 'white', RaceEthf==1~'black', RaceEthf==3~'hispanic', RaceEthf==2 ~'other')) %>% 
  dplyr::select(id, baseline_wt, year, time, event, lead.censor_mcc, lead.censor_km, lead.death, raceth) %>%
  group_by(id, baseline_wt, raceth) %>%
  nest()

set.seed(100)

### calculate point estimates ###
mcc_byrace = function(nlsy_nested) {
  nlsy_nested %>% 
    group_by(raceth) %>%
    filter(raceth != 'other') %>%
    nest() %>%
    mutate(mcc = map(data, mcc_weighted)) %>%
    dplyr::select(-data) %>%
    unnest(mcc) %>%
    pivot_wider(names_from = raceth, values_from = mcc) %>%
    mutate(bw_diff = black - white,
           hw_diff = hispanic - white)
}

### estimates overall and by race ###
mcc_estimates_overall = mcc_weighted(nlsy_nested) %>% rename(estimate = mcc)
mcc_estimates_byrace = mcc_byrace(nlsy_nested)
 
  

### Set number of bootstrap samples for confidence intervals ###
n_boots = 500
n_ids = nrow(nlsy_nested)

mcc_bootstraps = tibble(rep         = 1:n_boots,
                        data_nested = replicate(n        = n_boots,
                                                expr     = nlsy_nested[sample(1:n_ids, n_ids, replace=TRUE),],
                                                simplify = FALSE)) %>%
  mutate(mcc_byrace = map(data_nested, mcc_byrace),
         mcc_overall = map(data_nested, mcc_weighted))

### overall data ###
mcc_ses_overall = mcc_bootstraps %>%
  dplyr::select(-data_nested,-mcc_byrace) %>%
  unnest(mcc_overall) %>%
  group_by(time) %>%
  summarise(se = sd(mcc))

mc_ci_overall = mcc_estimates_overall %>%
  left_join(mcc_ses_overall) %>%
  mutate(lwr = estimate - 1.96*se, upr = estimate + 1.96*se)

mcc_ci_overall_100 = mc_ci_overall %>% select(-c(se)) %>%
  mutate(lwr = lwr*100, upr = upr*100, estimate = estimate*100)

### data by race ###
mcc_ses_byrace = mcc_bootstraps %>%
  select(-data_nested) %>%
  unnest(mcc_byrace) %>%
  group_by(time) %>%
  summarise(white = sd(white),
            black = sd(black),
            hispanic = sd(hispanic),
            bw_diff = sd(bw_diff),
            hw_diff = sd(hw_diff))


mcc_ci_byrace = mcc_estimates_byrace %>%
  pivot_longer(-time, names_to = 'parameter', values_to='estimate') %>%
  left_join(mcc_ses_byrace %>% pivot_longer(-time, names_to='parameter', values_to='se')) %>%
  mutate(lwr = estimate - 1.96*se, upr = estimate + 1.96*se)

mcc_ci_byrace_100 = mcc_ci_byrace %>% select(-c(se)) %>%
  mutate(lwr = lwr*100, upr = upr*100, estimate = estimate*100)



