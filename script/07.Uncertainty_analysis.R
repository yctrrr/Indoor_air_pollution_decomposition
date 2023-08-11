library(Cairo)
library(terra)
library(tictoc)
library(data.table)
library(dplyr)
library(foreach)
library(furrr)
library(ggpubr)

setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
dat_dir <- "D:\\shaoyanchuan\\data\\"
source('./script/Core_MonteCarlo.R',encoding = 'UTF8')
source('./script/Core.R', encoding = 'UTF8')
source('./script/Decomposition.R', encoding = 'UTF8')
source(paste0(dat_dir,"Population\\Fuel_Use_WHO\\Code\\HAP_Data_Setup.R"))
source(paste0(dat_dir,"Population\\Fuel_Use_WHO\\Code\\HAP_Model_parallel_provincial.R"))
source(paste0(dat_dir,"Population\\Fuel_Use_WHO\\Code\\HAP_Functions_parallel2.R"))

## Gridded variables containing AAP, population and heating degree
match_df = 2000:2020%>%
  purrr::map(~ rast(paste0("result\\Tif\\Match_data\\Matching_covariate_raster_g0.1_",.x,".tif"))%>%
        as.data.frame(match_rast, xy = TRUE, na.rm =FALSE)%>%
        dplyr::select(x,y,APM25,hd_count,Pop_all,Pop_U,Pop_O,pindex,countyindex)%>%
        na.omit()%>%mutate(year = .x))%>%reduce(full_join)

prov_cov <- fread("result\\Table\\Input_data\\covariates_by_province.csv", encoding = "UTF-8")
prov_info <- fread("result\\Table\\Input_data\\province_information.csv", encoding = "UTF-8")%>%
  mutate(index = pindex)

##=======Uncertainty analysis of fuel use population (Provincial level)========
deco_years = c(2000,2007,2013,2020)
METE = c("hd","cd","t2m")
ECO = c("GDP_pc_from_CSMAR","Urban_expenditure_pc","Rural_expenditure_pc")

# cov_data <- fread("result\\Table\\Input_data\\covariate_data_by_province.csv")

## fuel use fraction model
model_samples <- readRDS("result\\Rdata\\Fuel_pop\\Fit_all_samples.Rds")

## random select 1000 results of bayesian output
n_sim = dim(model_samples$processed_samples$beta_cov)[1]
sample_index = sample(1:n_sim, 1000)

## Uncertainty of fuel use fraction and variable contribution
deco_years
tic("Fuel use decomposition for mtcl")
Decomp_fuel = Pop_Decomposition_mtcl(deco_years, prov_cov, prov_info, sample_index)
toc()

## Uncertainty of fuel use fraction for c(2002,2007,2012,2017)
tic("Fuel use decomposition for mtcl")
Decomp_fuel_add = Pop_Decomposition_mtcl(c(2002,2007,2012,2017),
                                     prov_cov, prov_info, sample_index)
toc()

##  calculate SHAPE value (relative contribution) of two factors
fwrite(Decomp_fuel, "result\\Table\\Uncertainty\\Provincial_population_decomposition_uncertainty.csv")
fwrite(Decomp_fuel_add, "result\\Table\\Uncertainty\\Provincial_population_decomposition_uncertainty_add.csv")
Decomp_fuel%>%filter(sim == 1, year == 2013)%>%dplyr::select(year,sim,Area,Fuel,index,Start,End)
Decomp_fuel%>%filter(sim == 1, year == 2007)%>%dplyr::select(year,sim,Area,Fuel,index,Start,End)


Decomp_fuel <- fread("result\\Table\\Uncertainty\\Provincial_population_decomposition_uncertainty.csv", encoding = "UTF-8")

##  derive fuel use population uncertainty for further calculation (1000 samples)
Pop_fuel1 <- Decomp_fuel%>%
  dplyr::select(sim,Area,index,Fuel,year,Start,ECO,MET)%>%
  dplyr::rename(Popf = Start)

Pop_fuel2 <- Decomp_fuel%>%filter(year == 2013)%>%
  dplyr::select(sim,Area,index,Fuel,year,End)%>%
  dplyr::rename(Popf = End)%>%mutate(year = 2020, ECO = 0, MET = 0)

Pop_fuel <- full_join(Pop_fuel1, Pop_fuel2)
fwrite(Pop_fuel, "result\\Table\\Uncertainty\\Provincial_population_uncertainty.csv")

Pop_fuel_add1 <- Decomp_fuel_add%>%
  dplyr::select(sim,Area,index,Fuel,year,Start,ECO,MET)%>%
  dplyr::rename(Popf = Start)

Pop_fuel_add2 <- Decomp_fuel_add%>%filter(year == 2012)%>%
  dplyr::select(sim,Area,index,Fuel,year,End)%>%
  dplyr::rename(Popf = End)%>%mutate(year = 2017, ECO = 0, MET = 0)

Pop_fuel_add <- full_join(Pop_fuel_add1, Pop_fuel_add2)
fwrite(Pop_fuel_add, "result\\Table\\Uncertainty\\Provincial_population_uncertainty_add.csv")


fread("result\\Table\\Uncertainty\\Provincial_population_uncertainty.csv")%>%
  # filter(year == 2013)
  filter(sim == 2)%>%
  group_by(year,Area,Fuel)%>%
  summarise(ECO = sum(ECO), MET = sum(MET))

## Nonlinear relationship 
library(Metrics)
Decomp_fuel <- fread("result\\Table\\Uncertainty\\Provincial_population_decomposition_uncertainty.csv", encoding = "UTF-8")

mae(Decomp_fuel$MET_0, Decomp_fuel$MET_1)
cor(Decomp_fuel$ECO_0, Decomp_fuel$ECO_1)
mae(Decomp_fuel$ECO_0, Decomp_fuel$ECO_1)

##========Uncertainty analysis of PM2.5 death========
##**Mortality are calculated at national level##
##**Population are derived from satellited-based data##

## local function to simulate mortality uncertainty

#'@param HAP_id 0 HAP from Zhao et al., PNAS. 1 replace clean energy from Hu et al., EST. 2 self-defined HAP
#'@param HAP_list A vector indicates Urban_clean, Urabn_coal, Urban_biomass, Rural...
#'@param Pop_fuel simulations from hierarchical bayesian model
mort_mtcl <- function(samp = 1, total = 1000, HAP_id = 0, HAP_list = NULL,
                      deco_years = c(2000,2007,2013,2020), Pop_fuel = NULL){
  cat("========Uncertainty analysis:", 100*samp/total,"========\n")
  ## provincial population
  Pop_prov <- fread("result\\Table\\Input_data\\province_population.csv")%>%
    group_by(Province,Year)%>%
    mutate(urban = ifelse(urban*1.1 > total, urban, rnorm(1, urban, urban * 0.1/1.96)),
           rural = total - urban)%>%ungroup()%>%
    dplyr::rename(Urban = urban, Rural = rural,year = Year)%>%
    dplyr::select(Province, year ,Urban, Rural)%>%
    filter(year%in%deco_years)%>%
    pivot_longer(c(Urban,Rural), names_to = "Area", values_to = "Pop")

  ## provincial population with uncertainty analysis
  Pop_provsample <- Pop_fuel%>%
    dplyr::rename(sample_id = sim)%>%filter(sample_id == samp)%>%
    left_join(prov_info)%>%left_join(Pop_prov)

  ##  heating fraction
  hd_provf <- match_df%>%
    filter(year%in%deco_years, countyindex > 0)%>%
    group_by(year, pindex)%>%
    dplyr::summarise(hd_count = mean(hd_count))%>%
    ungroup()%>%suppressMessages()

  hd_provsamp <- hd_provf%>%
    mutate(hd_f = hd_count/yday(as.Date(paste0(year,"-12-31"))))%>%
    mutate(hd_f =  map_dbl(hd_f, ~ rnorm(1, .x, .x * 0.1/1.96)))%>%
    suppressMessages()
  
  ##  activity information
  activ_info <- fread("result\\Table\\Input_data\\activity_information.csv")%>%
    filter(Activity == "Outdoor")%>%
    mutate(Metric = ifelse(Metric == "NHS.Mean", "HS.Mean",
                           ifelse(Metric == "NHS.P5", "HS.P5",
                                  ifelse(Metric == "NHS.P95", "HS.P95", Metric))))%>%
    pivot_wider(names_from = 'Metric',values_from = 'Infl')%>%
    mutate(Infl = map2_dbl(HS.Mean, (HS.P95 - HS.P5)/3.28, ~ rnorm(n = 1, mean = .x, sd = .y)))%>%
    mutate(Condition = ifelse(Condition == "N-H","NH","H"))%>%
    draw.Time(.1)%>%
    dplyr::select(-HS.Mean,-HS.P5,-HS.P95)%>%
    pivot_wider(names_from = "Condition", values_from = c("Time", "Infl"))%>%
    dplyr::rename(whoname = Province)
  
  ## age structure
  ag = AgeGroupt%>%
    mutate(val = map2_dbl(val, (upper - lower)/3.92,
                          ~ rnorm(n = 1, mean = .x, sd = .y)), .keep = 'unused')%>%
    group_by(year)%>% mutate(AgeStruc = val/sum(val))%>% ungroup()%>%
    dplyr::select(year, sex, agegroup, AgeStruc, agediv)

  Agesum <- group_by(ag, year, sex, agediv)%>%
    dplyr::summarise(AgeStruc = sum(AgeStruc))%>%
    ungroup()%>%suppressMessages()

  ## mortality rate
  mRate = MortRatet%>%
    mutate(
      MortRate = map2_dbl(MortRate, (MortRateU - MortRateL)/3.92, ~ rnorm(n = 1, mean = .x, sd = .y)),
      .keep = 'unused')%>%
    dplyr::select(year, endpoint, sex, agegroup, MortRate, agediv)
  
  ## relative risk
  RR = RRstd %>% draw.RR

  ## exposure analysis
  ## provincial population-weighted household PM2.5 concentration

  if(HAP_id == 1){
    HAP_un <- left_join(Pop_provsample, hd_provsamp)%>%
      inner_join(activ_info)%>%
      left_join(Agesum)%>%rowwise()%>%
      mutate(HAP_I = ifelse(Area == "Urban"&Fuel == "Biomass",rnorm(1,223,98/1.96)*((1-Time_H)*hd_f+(1-Time_NH)*(1-hd_f)),
                     ifelse(Area == "Urban"&Fuel == "Coal", rnorm(1,38,10/1.96)*((1-Time_H)*hd_f+(1-Time_NH)*(1-hd_f)),
                     ifelse(Area == "Rural"&Fuel == "Biomass", rnorm(1,250,70/1.96)*((1-Time_H)*hd_f+(1-Time_NH)*(1-hd_f)),
                     ifelse(Area == "Rural"&Fuel == "Coal", rnorm(1,117,19/1.96)*((1-Time_H)*hd_f+(1-Time_NH)*(1-hd_f)), 
                            28.6*((1-Time_H)*hd_f+(1-Time_NH)*(1-hd_f)))))),
             Pop = round(Pop*AgeStruc))%>%
      dplyr::select(year, Area, Fuel, sex, agediv, pindex, Province, Popf, ECO, MET, Pop, HAP_I)%>%
      suppressMessages()
    
  }else if(HAP_id == 2){
    HAP_un <- left_join(Pop_provsample, hd_provsamp)%>%
      inner_join(activ_info)%>%
      left_join(Agesum)%>%rowwise()%>%
      mutate(HAP_I = ifelse(Area == "Urban"&Fuel == "Biomass",rnorm(1,HAP_list[3],HAP_list[3]*0.1)*((1-Time_H)*hd_f+(1-Time_NH)*(1-hd_f)),
                     ifelse(Area == "Urban"&Fuel == "Coal", rnorm(1,HAP_list[2],HAP_list[2]*0.1)*((1-Time_H)*hd_f+(1-Time_NH)*(1-hd_f)),
                     ifelse(Area == "Rural"&Fuel == "Biomass", rnorm(1,HAP_list[6],HAP_list[6]*0.1)*((1-Time_H)*hd_f+(1-Time_NH)*(1-hd_f)),
                     ifelse(Area == "Rural"&Fuel == "Coal", rnorm(1,HAP_list[5],HAP_list[5]*0.1)*((1-Time_H)*hd_f+(1-Time_NH)*(1-hd_f)), 
                     ifelse(Area == "Rural"&Fuel == "Clean", rnorm(1,HAP_list[4],HAP_list[3]*0.1)*((1-Time_H)*hd_f+(1-Time_NH)*(1-hd_f)),
                            rnorm(1,HAP_list[1],HAP_list[1]*0.1)*((1-Time_H)*hd_f+(1-Time_NH)*(1-hd_f))))))),
             Pop = round(Pop*AgeStruc))%>%
      dplyr::select(year, Area, Fuel, sex, agediv, pindex, Province, Popf, ECO, MET, Pop, HAP_I)%>%
      suppressMessages()
  }else{
    HAP_un <- left_join(Pop_provsample, hd_provsamp)%>%
      inner_join(activ_info)%>%
      left_join(Agesum)%>%rowwise()%>%
      mutate(HAP_I = ifelse(Area == "Urban"&Fuel == "Biomass",rnorm(1,223,98/1.96)*((1-Time_H)*hd_f+(1-Time_NH)*(1-hd_f)),
                            ifelse(Area == "Urban"&Fuel == "Coal", rnorm(1,38,10/1.96)*((1-Time_H)*hd_f+(1-Time_NH)*(1-hd_f)),
                                   ifelse(Area == "Rural"&Fuel == "Biomass", rnorm(1,250,70/1.96)*((1-Time_H)*hd_f+(1-Time_NH)*(1-hd_f)),
                                          ifelse(Area == "Rural"&Fuel == "Coal", rnorm(1,117,19/1.96)*((1-Time_H)*hd_f+(1-Time_NH)*(1-hd_f)), 0)))),
             Pop = round(Pop*AgeStruc))%>%
      dplyr::select(year, Area, Fuel, sex, agediv, pindex, Province, Popf, ECO, MET, Pop, HAP_I)%>%
      suppressMessages()
  }

  ## provincial population-weighted ambient PM2.5 concentration (insensitive to fuel use population)
  AAP_un <- match_df%>%filter(year%in%deco_years, countyindex > 0)%>%
    mutate(APM25 =  map_dbl(APM25, ~ rnorm(1, .x, .x * 0.1/1.96)))%>%
    group_by(year, pindex)%>%
    dplyr::summarise(APM25 = weighted.mean(APM25, Pop_all), Pop_all = sum(Pop_all))%>%
    ungroup()%>% left_join(hd_provsamp)%>%inner_join(activ_info)%>% left_join(Agesum)%>%
    mutate(AAP = APM25*Time_H*hd_f + APM25*Time_NH*(1-hd_f),
           HAP_O =  APM25*(1 - Time_H)*Infl_H*hd_f + APM25*(1 - Time_NH)*Infl_NH*(1-hd_f))%>%
    group_by(year,sex,pindex,agediv)%>%
    dplyr::summarise(AAP = weighted.mean(AAP, Pop_all*AgeStruc),
              HAP_O = weighted.mean(HAP_O, Pop_all*AgeStruc))%>%
    ungroup()%>%suppressMessages()

  ## overall concentration
  Ngroup <- inner_join(HAP_un,AAP_un)%>%
    mutate(concentration = (HAP_O + AAP + HAP_I)%>%matchable())%>%
    dplyr::select(-HAP_O,-AAP,-HAP_I,-Province)%>%suppressMessages()

  Grids <- distinct(Ngroup,Area,Fuel,pindex)%>%
    mutate(x = pindex, y = pindex)%>%dplyr::select(-pindex)
  
  PM_wr <- dplyr::select(Ngroup, -Pop,-ECO,-MET)%>%
    mutate(x = pindex, y = pindex)%>%dplyr::select(-pindex)
  
  ## urban or rural population and fuel use fraction
  Pop_group <- dplyr::select(Ngroup, -concentration)%>%
    group_by(year,Area,Fuel,pindex,Popf,ECO,MET)%>%
    dplyr::summarise(Pop = sum(Pop))%>%
    ungroup()%>%mutate(x = pindex, y = pindex)%>%
    dplyr::select(-pindex)%>%suppressMessages()

  Popf <- dplyr::select(Pop_group, -Pop)
  Pop <- dplyr::select(Pop_group, -Popf,-ECO,-MET)

  # Decomposition
  PM25_Shapesumlst <- vector(mode = "list", length = length(deco_years)-1)
  
  for (i in 1:(length(deco_years)-1)) {
    start.y <- deco_years[i]
    end.y <- deco_years[i+1]
    use_CR('MRBRT')
    
    PM25_Decomptmp <-  Decomposition(start.y = start.y,end.y = end.y,
                          RR = RR, Grids = Grids, Popf = Popf, Pop = Pop, 
                          AgeGroup = ag, PM_real = PM_wr, mRate = mRate)%>%
      mutate(year = deco_years[i])%>%
      suppressMessages()
    
    PM25_Shapesumlst[[i]] <- PM25_Decomptmp
  }

  PM25_Shapesum <- reduce(PM25_Shapesumlst, full_join)%>%
    mutate(sample_id = samp)%>%suppressMessages()

  return(PM25_Shapesum)
}

## Prepare necessary data for Montacaro simulation
deco_years = c(2000,2007,2013,2020)
use_CR('MRBRT')
# Data load
read_files(
  GRID = './result/Table/Integrated_PM2.5_exposure/GRID_information.csv',
  Pop = paste0('./result/Table/Population/Fuel_use_population_g0.1.csv'),
  PM_real = paste0("./result/Table/Integrated_PM2.5_exposure/Integrated_PM2.5_2000.csv"),
  PM_cf = './data/PM_Ctrl.csv', # PM_cf works only in counter-fact scenario
  MortRate = './data/GBD_incidence_China_2000-2019.csv',
  AgeGroup = './data/GBD_agestructure_China_2000-2017.csv'
)

AgeGroupt = fread('./data/IHME_GBD_2019_POP_SYA_2000-2020.csv')%>%
  # dplyr::select(sex, agegroup, val, upper,lower,year)%>%
  mutate(agediv = cut(agegroup, breaks = c(0,5,15,65,96),
                      labels = 1:4,
                      include.lowest = TRUE, right = FALSE)%>%as.numeric())%>%
  mutate(sex = tolower(sex), agegroup = as.character(agegroup))

MortRatet = fread('./data/IHME_GBD_2019_Deathrate_mod_2000-2020.csv')%>%
  filter(metric == "Rate")%>%
  dplyr::select(endpoint, sex, agegroup, MortRate, MortRateU, MortRateL, year)%>%
  mutate(agediv = cut(agegroup, breaks = c(0,5,15,65,96),
                      labels = 1:4, include.lowest = TRUE, right = FALSE)%>%as.numeric())%>%
  mutate(endpoint = tolower(endpoint),sex = tolower(sex), agegroup = as.character(agegroup))

RR = (RR_table$MEAN)%>%RR_std()%>%
  mutate(agediv = cut(as.numeric(agegroup), breaks = c(0,5,15,65,96),
                      labels = 1:4, include.lowest = TRUE, 
                      right = FALSE)%>%as.numeric())

# Relative risk
RRstd <- RR_std_for_mtcl(RR_table = RR_table, mode = mode) %>%
  mutate(agediv = cut(as.numeric(agegroup), breaks = c(0,5,15,65,96),
                      labels = 1:4,
                      include.lowest = TRUE, right = FALSE)%>%as.numeric())
## Uncertainty of fuel use population
Pop_fuel <- fread("result\\Table\\Uncertainty\\Provincial_population_uncertainty.csv")
Pop_fuel_add <- fread("result\\Table\\Uncertainty\\Provincial_population_uncertainty_add.csv")

## simulate Mortality uncertainty for 1000 times
tic("Mort")
plan(multisession, workers = 10)
options(future.globals.maxSize= 891289600)
Shape_sample = 1:1000%>%
  future_map_dfr(function(x){
    mort_mtcl(samp = x, total = 1000, Pop_fuel = Pop_fuel)
  },.options = furrr_options(seed = NULL))
toc()

fwrite(Shape_sample,"result\\Table\\Uncertainty\\Mortality_Decompisition_uncertainty.csv")


mort_mtcl(samp = 1, total = 1000, HAP_id = 1)
tic("Mort for clean fuel HAP")
plan(multisession, workers = 8)
options(future.globals.maxSize= 891289600)
Shape_sample_sens = 1:1000%>%
  future_map_dfr(function(x){
    mort_mtcl(samp = x, total = 1000, HAP_id = 1, Pop_fuel = Pop_fuel)
  },.options = furrr_options(seed = NULL))
toc()
fwrite(Shape_sample_sens,"result\\Table\\Uncertainty\\Mortality_Decompisition_Clean_Accounted_uncertainty.csv")

## PURE-AIR test
tic("Mort for PURE-AIR")
plan(multisession, workers = 5)
options(future.globals.maxSize= 891289600)
Shape_sample_PURE = 1:1000%>%
  future_map_dfr(function(x){
    mort_mtcl(samp = x, total = 1000, Pop_fuel = Pop_fuel,
              HAP_id = 2, HAP_list = c(22.95,42.66,39.72,22.95,42.66,39.72))
  },.options = furrr_options(seed = NULL))
toc()
fwrite(Shape_sample_PURE,"result\\Table\\Uncertainty\\Mortality_Decompisition_PURE_uncertainty.csv")

## PURE-baseline test
tic("Mort for PURE BASELINE")
plan(multisession, workers = 5)
options(future.globals.maxSize= 891289600)
Shape_sample_PUREB = 1:1000%>%
  future_map_dfr(function(x){
    mort_mtcl(samp = x, total = 1000, Pop_fuel = Pop_fuel,
              HAP_id = 2, HAP_list = c(4.86,53.83,51.22,4.86,53.83,51.22))
  },.options = furrr_options(seed = NULL))
toc()
fwrite(Shape_sample_PUREB,"result\\Table\\Uncertainty\\Mortality_Decompisition_PUREB_uncertainty.csv")

## WHO2018 test
tic("Mort for WHO")
plan(multisession, workers = 3)
options(future.globals.maxSize= 891289600)
Shape_sample_WHO = 1:1000%>%
  future_map_dfr(function(x){
    mort_mtcl(samp = x, total = 1000, Pop_fuel = Pop_fuel,
              HAP_id = 2, HAP_list = c(5.91,93.91,125.91,5.91,93.91,125.91))
  },.options = furrr_options(seed = NULL))
toc()
fwrite(Shape_sample_WHO,"result\\Table\\Uncertainty\\Mortality_Decompisition_WHO_uncertainty.csv")

## Comparison with Geng et al., 2021 Nature Geoscience (2002,2007,2012,2017)
tic("Mort for diifferent periods")
plan(multisession, workers = 10)
options(future.globals.maxSize= 891289600)
Shape_sample_aap = 1:1000%>%
  future_map_dfr(function(x){
    mort_mtcl(samp = x, total = 1000, Pop_fuel = Pop_fuel_add,
              deco_years = c(2002,2007,2012,2017))
  },.options = furrr_options(seed = NULL))
toc()
fwrite(Shape_sample_aap,"result\\Table\\Uncertainty\\Mortality_Decompisition_AAP-ECO_uncertainty.csv")

#========Visualization of population decomposition========
Pop_shape_sample <- fread("result\\Table\\Uncertainty\\Provincial_population_decomposition_uncertainty.csv", encoding = "UTF-8")%>%
  dplyr::select(year:sample_id,Start,PG,MET,ECO,End)%>%
  pivot_longer(Start:End, names_to = "Decomp")%>%
  group_by(year,Area,Fuel,sample_id,Decomp)%>%
  dplyr::summarise(value = sum(value))%>%ungroup()

color_palette = c('PG' = '#f9320c','MET' = '#f9c00c','ECO' = '#00b9f1',
                  'dock' = 'white')
deco_num = 3
deco_yname = c("2000",rep("",deco_num),"2007",rep("",deco_num),"2013",rep("",deco_num),"2020")

fuel_table <- expand.grid(Area = c("Urban","Rural"),Fuel = c("Clean","Biomass","Coal"))%>%
  mutate(name = paste0(Area,"_",Fuel))

## visualize the drivers of different fuel use population
1:nrow(fuel_table)%>%
  purrr::map(function(x){
    areaf = fuel_table[x,"Area"]
    fuelf = fuel_table[x,"Fuel"]
  
    Pop_shape_sum = Pop_shape_sample%>%
      filter(Area == areaf, Fuel == fuelf)%>%
      Decomp_to_plot()%>%
      mutate_at(vars(value,dock_upper,dock_lower), ~round(.x/1000))%>%
      mutate_at(vars(number_txt, number_txt1, number_txt2), ~(as.numeric(.x)/1000)%>%round()%>%as.character())

    Pop_shape_sumplot = ggplot(Pop_shape_sum)+
      geom_bar(aes(x = seq, y = value, fill = Decomp),color = "white",
               stat = 'identity', position ="identity", width = .8)+
      geom_text(
        aes(x = seq, y = value, label = var_txt1),
        hjust = 1.1, vjust = .5,angle = 90, size = 5) +
      geom_text(
        aes(x = seq, y = dock_lower, label = var_txt2),
        hjust = 1.1, vjust = .5,angle = 90, size = 5) +
      geom_text(
        aes(x = seq, y = value, label = number_txt1), 
        hjust = .5,vjust = -.5,angle = 0, size = 5)+
      geom_text(
        aes(x = seq, y = dock_upper, label = number_txt2), 
        hjust = .5,vjust = -.5,angle = 0, size = 5)+
      geom_errorbar(mapping = aes(x = seq, ymin = dock_lower, ymax = dock_upper),
                    width = .3, color = "black") +
      labs(y = NULL)+
      scale_fill_manual(values = color_palette)+
      # scale_color_manual(values = color_palette)+
      scale_x_discrete(labels = deco_yname) +
      scale_y_continuous(
        name = "Population (million)",
        expand = c(0, 0), limits = c(-800,800)) +
      theme(
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white",colour = "black", size = .5),
        legend.position = 'none',
        plot.margin = unit(c(.2, .5, .1, .2), "cm"),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = 'black', size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18,hjust = 0.5))
  
    ggsave(paste0("result/Fig/Decomposition/Decomposition_of_Population_",
                  areaf,"_", fuelf,"_Uncertainty.png"),
           Pop_shape_sumplot, width = 8, height = 4.5)
  }
)

#========Visualization of mortality decomposition========
library(dplyr)
fuel_table <- expand.grid(Area = c("Urban","Rural"),Fuel = c("Clean","Biomass","Coal"))%>%
  mutate(name = paste0(Area,"_",Fuel))

Shape_sum <- fread("result\\Table\\Uncertainty\\Mortality_Decompisition_uncertainty.csv")%>%
  mutate(Mort = replace_na(Mort,0))%>%
  group_by(year, Decomp, sample_id)%>%
  dplyr::summarise(value = sum(Mort,na.rm = TRUE))%>%ungroup()%>%
  dplyr::filter(Decomp != "FF")%>%
  Decomp_to_plot()

Shape_sum_AAPECO <- fread("result\\Table\\Uncertainty\\Mortality_Decompisition_AAP-ECO_uncertainty.csv")%>%
  mutate(Mort = replace_na(Mort,0))%>%
  group_by(year, Decomp, sample_id)%>%
  dplyr::summarise(value = sum(Mort,na.rm = TRUE))%>%ungroup()%>%
  dplyr::filter(Decomp != "FF")
  Decomp_to_plot(deco_years = c(2002,2007,2012,2017))

## Uncertainty of ECO to Exposure from AAP
fread("result\\Table\\Uncertainty\\Mortality_Decompisition_AAP-ECO_uncertainty.csv")%>%
  na.omit%>%filter(Decomp %in% c("EXP"))%>%
  mutate(ECO_AAP = ifelse(year == 2002, 2.19*Mort, ifelse(year == 2007, -3.47*Mort, -0.82*Mort)))%>%
  group_by(sample_id)%>%
  dplyr::summarise(value = sum(ECO_AAP))%>%ungroup()%>%
  # group_by(year)%>%
  dplyr::summarise(lower = quantile(value, 0.025)%>%round, mean = quantile(value, 0.5)%>%round,
                   upper = quantile(value, 0.975)%>%round)


Shape_sum_sens <- fread("result\\Table\\Uncertainty\\Mortality_Decompisition_Clean_Accounted_uncertainty.csv")%>%
  group_by(year, Decomp, sample_id)%>%
  dplyr::summarise(value = sum(Mort,na.rm = TRUE))%>%ungroup()%>%
  dplyr::filter(Decomp != "FF")%>%
  Decomp_to_plot()

Shape_sum_PURE <- fread("result\\Table\\Uncertainty\\Mortality_Decompisition_PURE_uncertainty.csv")%>%
  group_by(year, Decomp, sample_id)%>%
  dplyr::summarise(value = sum(Mort,na.rm = TRUE))%>%ungroup()%>%
  dplyr::filter(Decomp != "FF")%>%
  Decomp_to_plot()

Shape_sum_PUREB <- fread("result\\Table\\Uncertainty\\Mortality_Decompisition_PUREB_uncertainty.csv")%>%
  group_by(year, Decomp, sample_id)%>%
  dplyr::summarise(value = sum(Mort,na.rm = TRUE))%>%ungroup()%>%
  dplyr::filter(Decomp != "FF")%>%
  Decomp_to_plot()

Shape_sum_WHO <- fread("result\\Table\\Uncertainty\\Mortality_Decompisition_WHO_uncertainty.csv")%>%
  group_by(year, Decomp, sample_id)%>%
  dplyr::summarise(value = sum(Mort,na.rm = TRUE))%>%ungroup()%>%
  dplyr::filter(Decomp != "FF")%>%
  Decomp_to_plot()

## Table output
c("WHO","PURE","PUREB")%>%
  map_dfr(function(x){
    Shape_sum_un <- get(paste0("Shape_sum_",x))%>%
      mutate(mean = ifelse(Decomp == "Start"| Decomp == "End", dock/1e3, mean),
             lower = lower/1e3, upper = upper/1e3)%>%
      filter(Decomp != "dock")%>%
      mutate(value = str_c(round(mean)," (",round(lower),"-",round(upper),")"))%>%
      dplyr::select(year, Decomp, value)%>%
      mutate_name(factor = FALSE)%>%
      mutate(name = ifelse(name == "",Decomp,name))%>%
      pivot_wider(names_from = year, values_from = value)%>%
      mutate(Study = x)
  })%>%fwrite("result\\Table\\Uncertainty\\Mortality_Decompisition_turb_uncertainty_summary.csv")

## Outdoor PM2.5 validation
valid_AAP <- match_df%>%filter(year%in%deco_years, countyindex > 0)%>%
  group_by(year, pindex)%>%
  dplyr::summarise(APM25 = weighted.mean(APM25, Pop_all), Pop_all = sum(Pop_all))%>%
  ungroup()%>%left_join(hd_provsamp)%>%inner_join(activ_info)%>% left_join(Agesum)%>%
  mutate(AAP = APM25*Time_H*hd_f + APM25*Time_NH*(1-hd_f),
         HAP_O =  APM25*(1 - Time_H)*Infl_H*hd_f + APM25*(1 - Time_NH)*Infl_NH*(1-hd_f),
         Pop = Pop_all*AgeStruc)%>%
  group_by(year,sex,pindex,agediv)%>%
  dplyr::summarise(AAP = weighted.mean(AAP, Pop_all*AgeStruc),
                   HAP_O = weighted.mean(HAP_O, Pop_all*AgeStruc),
                   APM25 = weighted.mean(APM25, Pop_all*AgeStruc),
                   Pop = sum(Pop))%>%
  ungroup()%>%suppressMessages()
fread(paste0("./result/Table/Integrated_PM2.5_exposure/Integrated_PM2.5_",year,".csv"))

## output table
valid_AAP%>%
  mutate(AAP_infl = AAP+HAP_O)%>%
  group_by(year)%>%
  dplyr::summarise(APM_infl = weighted.mean(AAP_infl, Pop),
                   AAP = weighted.mean(AAP, Pop),
                   HAP_O = weighted.mean(HAP_O, Pop),
                   APM25 = weighted.mean(APM25, Pop))%>%
  ungroup()%>%fwrite("result\\Table\\Uncertainty\\Population-weighted_AAP_2002-2017.csv")
  
  
Shape_sum_AAPECO%>%
  mutate(mean = ifelse(Decomp == "Start"| Decomp == "End", dock/1e3, mean),
         lower = lower/1e3, upper = upper/1e3)%>%
  filter(Decomp != "dock")%>%
  mutate(value = str_c(round(mean)," (",round(lower),"-",round(upper),")"))%>%
  dplyr::select(year, Decomp, value)%>%
  mutate_name(factor = FALSE)%>%
  mutate(name = ifelse(name == "",Decomp,name))%>%
  pivot_wider(names_from = year, values_from = value)%>%
  fwrite("result\\Table\\Uncertainty\\Mortality_Decompisition_uncertainty_2002-2017.csv")


## Main decomposition plot
data = fread("result\\Table\\Uncertainty\\Mortality_Decompisition_uncertainty.csv")%>%
  mutate(Mort = replace_na(Mort,0))%>%
  group_by(year, Decomp, sample_id)%>%
  dplyr::summarise(value = sum(Mort,na.rm = TRUE))%>%ungroup()%>%
  dplyr::filter(Decomp != "FF", year%in%c(2000,2013,2020))

Shape_by_pop <- 1:6%>%
  purrr::map_dfr(~ fread("result\\Table\\Uncertainty\\Mortality_Decompisition_uncertainty.csv")%>%
                   mutate(Mort = replace_na(Mort,0))%>%
                   group_by(year, Area, Fuel, Decomp, sample_id)%>%
                   summarise(value = sum(Mort, na.rm = TRUE))%>%ungroup()%>%
                   filter(Area == fuel_table$Area[.x], Fuel == fuel_table$Fuel[.x])%>%
                   dplyr::filter(Decomp != "FF")%>%
                   Decomp_to_plot()%>%
                   mutate(name = .x, Area = fuel_table$Area[.x], Fuel = fuel_table$Fuel[.x])
  )

Shape_by_pop_sens <- 1:6%>%
  purrr::map_dfr(~ fread("result\\Table\\Uncertainty\\Mortality_Decompisition_Clean_Accounted_uncertainty.csv")%>%
                   group_by(year, Area, Fuel, Decomp, sample_id)%>%
                   summarise(value = sum(Mort, na.rm = TRUE))%>%ungroup()%>%
                   filter(Area == fuel_table$Area[.x], Fuel == fuel_table$Fuel[.x])%>%
                   dplyr::filter(Decomp != "FF")%>%
                   Decomp_to_plot()%>%
                   mutate(name = .x, Area = fuel_table$Area[.x], Fuel = fuel_table$Fuel[.x])
  )

Shape_raw <- 1:6%>%
  purrr::map_dfr(~ fread("result\\Table\\Uncertainty\\Mortality_Decompisition_uncertainty.csv")%>%
                   group_by(year, Area, Fuel, Decomp, sample_id)%>%
                   summarise(value = sum(Mort, na.rm = TRUE))%>%ungroup()%>%
                   filter(Area == fuel_table$Area[.x], Fuel == fuel_table$Fuel[.x])%>%
                   dplyr::filter(Decomp != "FF")
  )

##  ready to plot
deco_num = 6
deco_yname = c("2000",rep("",deco_num),"2007",rep("",deco_num),"2013",rep("",deco_num),"2020")
color_palette = c('FF' = '#d95f0e' ,'ECO' = "#f9320c", 'MET' = "#9ecae1",
  'PC' = '#bcbddc','PA' = '#f9c00c','EXP' = '#00b9f1',
  'ORF' = '#5CAB7D','dock' = 'white')

breakna <- function(x) {
  x <- breaks_pretty()(x)
  x[x<0] <- NA
  x
}

findnear <- function(x, value){
  map_dbl(value, ~ x[which.min(abs(.x - x))])
}

Shape_plot <- function(data, limits = NULL){
  rang <- limits[2] - limits[1]
  
  brack_df = data%>%
    filter(Decomp %in%c("MET","ECO"))%>%
    mutate(start = findnear(seq(2,21, by = 7), as.numeric(seq)))%>%
    dplyr::select(Decomp,start,dock_upper)%>%
    pivot_wider(names_from = "Decomp", values_from = "dock_upper")%>%
    mutate(label = "Fuel Use Fraction", 
           label.pos = pmax(MET,ECO) + rang*0.095, end = start +1)

  ggplot(data)+
    geom_bar(aes(x = seq, y = value, fill = Decomp),color = "white",
             stat = 'identity', position ="identity", width = .8)+
    geom_text(
      aes(x = seq, y = value, label = var_txt1),
      hjust = 1.1, vjust = .5,angle = 90, size = 7.5)+
    geom_text(
      aes(x = seq, y = dock_lower, label = var_txt2),
      hjust = 1.1, vjust = .5,angle = 90, size = 7.5)+
    geom_text(
      aes(x = seq, y = value, label = number_txt1),
      hjust = .5,vjust = -.5,angle = 0, size = 7.5)+
    geom_text(
      aes(x = seq, y = dock_upper, label = number_txt2), 
      hjust = .5,vjust = -.5,angle = 0, size = 7.5)+
    geom_errorbar(mapping = aes(x = seq, ymin = dock_lower, ymax = dock_upper),
                  width = .3, color = "black")+
    geom_bracket(aes(xmin = start, xmax = end, y.position = label.pos, label = label),
                 data = brack_df, size = .6, label.size = 7, vjust = -0.5)+
    labs(y = NULL)+
    scale_fill_manual(values = color_palette)+
    # scale_color_manual(values = color_palette)+
    scale_x_discrete(labels = deco_yname)+
    scale_y_continuous(
      name = expression(PM['2.5'] *  ~ Health ~ Burden ~ '(' * thousand * ')'),
      expand = c(0, 0), limits = limits, breaks = breakna)+
    theme(
      strip.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white",colour = "black", linewidth= .5),
      legend.position = 'none',
      plot.margin = unit(c(.2, .5, .1, .2), "cm"),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_line(colour = "black"),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(color = 'black', size = 25),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 28,hjust = 0.5))
  
}

##  Overall plot
Shape_sumplot <- Shape_plot(Shape_sum, limits = c(0, 2250))

CairoPNG(paste0("result/Fig/Decomposition/Decomposition_of_Weighted_PM2.5_Death_Uncertainty.png"),
         width = 3000,height = 1300, res = 180)
print(Shape_sumplot)
dev.off()

##  Sensitivity analysis of HAP for claen fuel users
Shape_sum_sensplot <- Shape_plot(Shape_sum_sens, limits = c(0, 2800))
CairoPNG(paste0("result/Fig/Decomposition_sensitivity/Decomposition_of_Weighted_PM2.5_Death_Clean_Accounted_Uncertainty.png"),
         width = 3000,height = 1300, res = 180)
print(Shape_sum_sensplot)
dev.off()

##  Decomposition by population groups
##  Seperate Death
c(1:6)%>%
  purrr::map(~ {
    Shape_pop <- Shape_by_pop%>%
      filter(name == .x)
    if(.x %in%c(3)){
      upv = 1.5
      dov = .9
    }else{
      upv = 1.45
      dov = 1.35
    }
    limitu <- max(Shape_pop$dock_upper%>%na.omit)*upv
    limitl <- 0- max(Shape_pop$dock_upper%>%na.omit)*dov
    
    CairoPNG(paste0("result/Fig/Decomposition/Decomposition_of_Weighted_PM2.5_Death_",
                    fuel_table$name[.x],"_Uncertainty.png"),
             width = 3000,height = 1300, res = 180)
    Shape_plot(Shape_pop, limits = c(limitl, limitu))%>%print()
    dev.off()
  })

1:6%>%
  purrr::map(~ {
    Shape_pop <- Shape_by_pop_sens%>%
      filter(name == .x)
    
    limitu <- max(Shape_pop$dock_upper%>%na.omit)*1.4
    limitl <- 0- max(Shape_pop$dock_upper%>%na.omit)*1
    
    CairoPNG(paste0("result/Fig/Decomposition_sensitivity/Decomposition_of_Weighted_PM2.5_Death_",
                    fuel_table$name[.x],"_Clean_Accounted_Uncertainty.png"),
             width = 3000,height = 1300, res = 180)
    Shape_plot(Shape_pop, limits = c(limitl, limitu))%>%print()
    dev.off()
  })

Shape_aggreplot <- Shape_by_pop%>%mutate(var_txt1 = "", var_txt2 = "")%>%
  Shape_plot()+
  facet_grid(Fuel ~ Area, labeller = labeller(groupwrap = label_wrap_gen(10)), scales = "free")+
  theme(legend.position = 'right')

ggsave(paste0("result/Fig/Decomposition/Decomposition_of_Weighted_PM2.5_Death_Group_Uncertainty.png"),
       Shape_aggreplot, width = 20, height = 14.5)
