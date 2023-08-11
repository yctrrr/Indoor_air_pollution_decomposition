# PM25 Health Impact Calculate Core
# 
#        - Monte Carlo Uncertainty Analysis
#   WARNING: Monte Carlo Analysis IS Extremely Time-Consuming
# 
#                    By Yifan LIU 2022/5/22
#                    Modified by Yanchuan Shao 2022/12/10

library(tidyverse)
library(furrr)
mode = 'MRBRT'
tell.mode <- function(mode = mode) case_when(
  mode %in% c('5COD', 'NCD+LRI') ~ str_glue('GEMM-{mode}'),
  mode %>% str_detect('IER') ~ mode,
  mode %>% str_detect('MRBRT') ~ mode
)


draw.RR <- function(x) x %>% mutate(
  RR = map2_dbl(MEAN, (UP - LOW) / 3.92, ~ rnorm(n = 1, mean = .x, sd = .y)),
  .keep = 'unused')


draw.PM <- function(x, uncert = .12) x %>% mutate(
  concentration = map_chr(as.numeric(concentration), ~ rnorm(1, .x, .x * uncert/1.96) %>% matchable(1))
)

draw.Time <- function(x, uncert = .12) x %>% mutate(
  Time = map_dbl(Time, ~ rnorm(1, .x, .x * uncert/1.96) )
)

matchable <- function(num, digit = 1) num %>% round(digit) %>% str_c

mortrate_std <- function(x) x %>% mutate(endpoint = tolower(endpoint))

RR_std_for_mtcl <- function(RR_table, mode) {
  
  RR_tbl <- RR_table %>% imap_dfr(
    ~ .x %>% mutate(CI = .y) %>% 
      pivot_longer(
        cols = -c(concentration, CI), values_to = "RR",
        names_to = c("endpoint", "agegroup"),names_sep = '_')
  ) %>% pivot_wider(names_from = 'CI',values_from = 'RR') %>% 
    mutate(endpoint = tolower(endpoint))
  
  RR_reshape <- if (mode == 'NCD+LRI') {
    
    expand_grid(
      concentration = RR_table[[1]] %>% pull(concentration),
      endpoint = c('ncd+lri'),
      agegroup = c('ALL', seq(25, 95, 5) %>% matchable(0))
    ) %>% left_join(
      RR_tbl, by = c('concentration', 'endpoint', 'agegroup')
    ) %>% group_by(concentration, endpoint) %>% fill(MEAN, LOW, UP) %>%
      filter(agegroup != 'ALL') %>% ungroup
    
  } else if (mode == '5COD') {
    
    expand_grid(
      concentration = RR_table[[1]] %>% pull(concentration),
      endpoint = c('copd', 'ihd', 'lc', 'lri', 'stroke'),
      agegroup = c('ALL', seq(25, 95, 5) %>% matchable(0))
    ) %>% left_join(
      RR_tbl, by = c('concentration', 'endpoint', 'agegroup')
    ) %>% group_by(concentration, endpoint) %>% fill(MEAN, LOW, UP) %>%
      filter(agegroup != 'ALL') %>% ungroup
    
  } else if (mode %>% str_detect('IER')) {
    
    expand_grid(
      concentration = RR_table[[1]] %>% pull(concentration),
      endpoint = c('copd', 'ihd', 'lc', 'stroke', 'lri'),
      agegroup = c('ALL', seq(0, 95, 5) %>% matchable(0))
    ) %>% left_join(
      RR_tbl, by = c('concentration', 'endpoint', 'agegroup')
    ) %>% group_by(concentration, endpoint) %>%
      fill(MEAN, LOW, UP) %>% filter(agegroup != 'ALL') %>%
      filter(
        endpoint %in% c('copd', 'ihd', 'lc', 'stroke') &
          as.integer(agegroup) >= 25 | endpoint == 'lri'
      ) %>% ungroup
  } else if (mode == 'MRBRT') {
    
    expand_grid(
      concentration = RR_table[[1]] %>% pull(concentration),
      endpoint = c('copd', 'dm', 'ihd', 'lc', 'lri', 'stroke'),
      agegroup = c('ALL', seq(0, 95, 5) %>% matchable(0))
    ) %>% left_join(
      RR_tbl, by = c('concentration', 'endpoint', 'agegroup')
    ) %>% group_by(concentration, endpoint) %>%
      fill(MEAN, LOW, UP) %>% filter(agegroup != 'ALL') %>% filter(
        endpoint %in% c('copd', 'dm', 'ihd', 'lc', 'stroke') &
          as.integer(agegroup) >= 25 | endpoint == 'lri'
      ) %>% ungroup
  }
}

Mortality_mtcl <- function(RR, PM_r,  ag, mRate, pop, core = 10, n = 10000, mode){
  plan(multisession, workers = core)
  
  MTCL <- 1:n %>% future_map(.options = future_options(seed = FALSE), ~ {
    
    RR_tbl <- RR_std_for_mtcl(RR_table = RR, mode = mode) %>% draw.RR
    
    PM_r %>% 
      left_join(RR_tbl, by = 'concentration') %>% 
      left_join(pop) %>% 
      group_by(endpoint, agegroup) %>%
      mutate(PWRR_real = weighted.mean(RR, Pop)) %>% 
      ungroup %>% 
      rename(RR_cf = RR) %>% 
      left_join(ag, by = 'agegroup') %>%
      left_join(mRate %>% mortrate_std, by  = c('agegroup', 'endpoint')) %>%
      mutate(
        Mort = Pop * AgeStruc * MortRate * (RR_cf - 1) / PWRR_real / 1e5,
        .keep = 'unused'
      ) %>% summarise(Mort = sum(Mort)) %>% deframe
  })
  
  plan(sequential)
  
  return(MTCL)
}

##=====Uncertainty of fuel use fraction======
#' @param cov_step covariates data
#' @param sample_index index indicated which number will be extracted from mcmc chains
Fuel_pop_mtcl <- function(cov_step, plot_years, sample_index, info){
  fuel_names = c('Solid','Biomass','Wood', 'Clean','Coal','Crop')
  area_names = c('Urban','Rural')
  
  plan(multisession, workers = 10)
  options(future.globals.maxSize= 891289600)
  fuel_sample = 1:nrow(info)%>%
    future_map_dfr(function(i){
      fuel_usec = fuel_use_output(plot_id = info$index[i],
                                  plot_index = info$pindex[i], cov_data = cov_step,
                                  plot_years = plot_years, all_samples = model_samples, 
                                  type='survey', output='samples',sampling_bias = TRUE)%>%suppressMessages()
      
      ##  extract sample index
      fuel_usec_mtcl <- fuel_usec[sample_index,,,]
      
      dimnames(fuel_usec_mtcl) <- list(sim = 1:length(sample_index),
                                       year = plot_years,
                                       Area=area_names, Fuel=fuel_names)
      
      melt(fuel_usec_mtcl)%>%
        as_tibble()%>%filter(Fuel%in%c("Biomass","Coal"))%>%
        pivot_wider(names_from = Fuel, values_from = value)%>%
        mutate(index = info$index[i])
    },.options = furrr_options(seed = NULL))
  
  ##  calculate fuel use fraction at provincial level for 1000 times
  fuel_sample%>%
    mutate(Clean = 1 - Biomass - Coal)%>%
    pivot_longer(c("Biomass","Coal","Clean"), names_to = "Fuel", values_to = "Popf")
}

##=====Mortality MTCL=====
mort_mtcl <- function(samp = 1, total = 1000){
  cat("========Uncertainty analysis:", 100*samp/total,"========\n")
  ##  heating fraction
  hd_provf <- match_df%>%
    filter(countyindex > 0)%>%
    group_by(year, pindex)%>%
    summarise(hd_count = mean(hd_count))%>%
    ungroup()%>%suppressMessages()
  
  hd_provsamp <- hd_provf%>%
    mutate(hd_f = hd_count/yday(as.Date(paste0(year,"-12-31"))))%>%
    mutate(hd_f =  map_dbl(hd_f, ~ rnorm(1, .x, .x * 0.1/1.96)))%>%
    suppressMessages()
  
  ##  activity information
  activ_info <- fread("result\\Table\\Activity_information.csv")%>%
    filter(Activity == "Outdoor")%>%
    mutate(Metric = ifelse(Metric == "NHS.Mean", "HS.Mean",
                           ifelse(Metric == "NHS.P5", "HS.P5",
                                  ifelse(Metric == "NHS.P95", "HS.P95", Metric))))%>%
    pivot_wider(names_from = 'Metric',values_from = 'Infl')%>%
    mutate(Infl = map2_dbl(HS.Mean, (HS.P95 - HS.P5)/3.28, ~ rnorm(n = 1, mean = .x, sd = .y)))%>%
    mutate(Condition = ifelse(Condition == "N-H","NH","H"))%>%
    draw.Time(.1)%>%
    dplyr::select(-HS.Mean,-HS.P5,-HS.P95)%>%
    pivot_wider(names_from = "Condition", values_from = c("Time", "Infl"))
  
  
  ## age structure
  ag = AgeGroupt%>%
    mutate(val = map2_dbl(val, (upper - lower)/3.92,
                          ~ rnorm(n = 1, mean = .x, sd = .y)), .keep = 'unused')%>%
    group_by(year)%>% mutate(AgeStruc = val/sum(val))%>% ungroup()%>%
    dplyr::select(year, sex, agegroup, AgeStruc, agediv)
  
  Agesum <- group_by(ag, year, sex, agediv)%>%
    summarise(AgeStruc = sum(AgeStruc))%>%
    ungroup()%>%suppressMessages()
  
  ## mortality rate
  mRate = MortRatet%>%
    mutate(
      MortRate = map2_dbl(MortRate, (MortRateU - MortRateL)/3.92, ~ rnorm(n = 1, mean = .x, sd = .y)),
      .keep = 'unused')%>%
    dplyr::select(year, endpoint, sex, agegroup, MortRate, agediv)
  
  ## relative risk
  RR = RR_std %>% draw.RR
  
  ##  exposure analysis
  HAP_un <- left_join(Pop_provsample, hd_provsamp)%>%
    inner_join(activ_info)%>%
    left_join(Agesum)%>%rowwise()%>%
    # mutate(UB = rnorm(1,223,98/1.96))
    mutate(concentration = ifelse(Area == "Urban"&Fuel == "Biomass",rnorm(1,223,98/1.96)*((1-Time_H)*hd_f+(1-Time_NH)*(1-hd_f)),
                                  ifelse(Area == "Urban"&Fuel == "Coal", rnorm(1,38,10/1.96)*((1-Time_H)*hd_f+(1-Time_NH)*(1-hd_f)),
                                         ifelse(Area == "Rural"&Fuel == "Biomass", rnorm(1,250,70/1.96)*((1-Time_H)*hd_f+(1-Time_NH)*(1-hd_f)),
                                                ifelse(Area == "Rural"&Fuel == "Coal", rnorm(1,117,19/1.96)*((1-Time_H)*hd_f+(1-Time_NH)*(1-hd_f)), 0))))
    )%>%
    group_by(year, Area, Fuel, sex, agediv)%>%
    summarise(concentration_hap = weighted.mean(concentration, Pop*AgeStruc),
              Pop = sum(Pop*AgeStruc)%>%round())%>%
    ungroup()%>%mutate(sample_id = samp)%>%suppressMessages()
  

  AAP_un <- match_df%>%filter(countyindex>0)%>%
    mutate(APM25 =  map_dbl(APM25, ~ rnorm(1, .x, .x * 0.1/1.96)))%>%
    group_by(year, pindex)%>%
    summarise(APM25 = weighted.mean(APM25, Pop_all), Pop_all = sum(Pop_all))%>%
    ungroup()%>% left_join(hd_provsamp)%>%inner_join(activ_info)%>% left_join(Agesum)%>%
    mutate(APM25_H = (APM25*Time_H + APM25*(1 - Time_H)*Infl_H)*hd_f,
           APM25_NH = (APM25*Time_NH + APM25*(1 - Time_NH)*Infl_NH)*(1-hd_f),
           concentration_aap = APM25_H + APM25_NH)%>%
    group_by(year,sex,agediv)%>%
    summarise(concentration_aap = weighted.mean(concentration_aap, Pop_all*AgeStruc),
              Pop_group = sum(Pop_all*AgeStruc)%>%round())%>%
    ungroup()%>%suppressMessages()
  
  Ngroup <- inner_join(HAP_un,AAP_un)%>%
    mutate(concentration = (concentration_hap + concentration_aap)%>%matchable())%>%
    dplyr::select(-concentration_hap,-concentration_aap,-Pop_group)%>%suppressMessages()
  
  Grids <- distinct(Ngroup,Area,Fuel)%>%
    mutate(x = Area, y = Fuel)
  
  PM_wr <- dplyr::select(Ngroup, -Pop)
  pop <- dplyr::select(Ngroup, -concentration)%>%
    group_by(year,Area,Fuel)%>%summarise(Pop = sum(Pop))%>%
    ungroup()%>%suppressMessages()
  
  # Decomposition
  PM25_Shapesumlst <- list(mode = "list", length = length(deco_years)-1)
  for (i in 1:(length(deco_years)-1)) {
    start.y <- deco_years[i]
    end.y <- deco_years[i+1]
    use_CR('MRBRT')
    PM25_Decomptmp <- 1:24%>%
      map(~ Decomposition(serie = .x, start.y = as.character(start.y), 
                          end.y = as.character(end.y),
                          RR = RR, Grids = Grids, Pop = pop, AgeGroup = ag,
                          PM_real = PM_wr, mRate = mRate))
    
    PM25_Shapetmp <- PM25_Decomptmp%>%
      Reduce(function(x,y) left_join(x,y, by = c("x","y","Area","Fuel")), .)%>%
      mutate(Start = rowMeans(select(., starts_with("Start"))),
             GPG = rowMeans(select(., starts_with("PG"))),
             GPA = rowMeans(select(., starts_with("PA"))),
             GORF = rowMeans(select(., starts_with("ORF"))),
             GEXP = rowMeans(select(., starts_with("EXP"))),
             End = rowMeans(select(., starts_with("End"))),
             .keep = "unused")%>%
      dplyr::select(x:y,Fuel,Area,Start:End)
    
    PM25_Shapesumtmp <- PM25_Shapetmp%>%
      # summarise_at(vars(Start:End),sum)%>%
      mutate(year = deco_years[i])%>%
      dplyr::select(-x,-y)%>%
      pivot_longer(cols = Start:End, names_to = "Decomp")
    
    # filter(Decomp%in%c("Start", "GPG", "GPA", "GORF", "GEXP", "End"))
    PM25_Shapesumlst[[i]] <- PM25_Shapesumtmp
  }
  
  PM25_Shapesum_wide <- reduce(PM25_Shapesumlst, full_join)%>%
    mutate(sample_id = samp)%>%suppressMessages()
  
  return(PM25_Shapesum_wide)
}
