

setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
dat_dir <- "D:\\shaoyanchuan\\data\\"
library(data.table)
library(dplyr)
library(readr)
library(aomisc)
library(minpack.lm)
library(tidyr)
library(stargazer)
library(purrr)

prov_cov <- fread("result\\Table\\Input_data\\covariates_by_province.csv", encoding = "UTF-8")%>%
  dplyr::rename(hdd = hd, cdd = cd, year = Year)%>%
  mutate_at(vars(GDP_pc_from_CSMAR,Urban_expenditure_pc,Rural_expenditure_pc),
            function(x){x/1e3})
input_data <- fread("result/Table/Input_data/input_survey_data.csv", encoding = 'UTF-8')
prov_info <- fread("result/Table/Input_data/province_information.csv", encoding = 'UTF-8')%>%
  mutate(Province = whoname)
IV = c("GDP_pc_from_CSMAR","Urban_expenditure_pc","Rural_expenditure_pc",
       "hdd","cdd","t2m")
fuel_usepop_lookup <- fread("result/Table/Population/Fuel_use_population_by_province.csv")

##=====Comparison with asynptotic regression=======
## Compare the usage of clean fuel with Tao et al.,(Nature Energy 2018)
area_indx = c("urban","rural","total")
for(areai in 1:2){
  survey_data_train <- input_data%>%
    left_join(prov_info)%>%
    dplyr::rename(index = pindex)%>%
    filter(area == areai)%>%
    left_join(prov_cov)%>%
    dplyr::select(Province, year, clean, all_of(IV))%>%
    mutate(clean = clean*100)

  prenon <- drm(clean ~ GDP_pc_from_CSMAR,
                data = survey_data_train,
                fct = DRC.asymReg())


  parm <- prenon$parmMat%>%as.vector()
  NLS_model1 <- nlsLM(clean ~ plateau*(1-exp(-m*GDP_pc_from_CSMAR)) + cov*hdd + int,
                  data = survey_data_train,
                  start = list(plateau = parm[3], m = parm[2], cov = 0, int = parm[1]),
                  control = list(maxiter = 1000))
  summary(NLS_model1)
  
  NLS_model2 <- nlsLM(clean ~ plateau*(1-exp(-m*GDP_pc_from_CSMAR-p*Urban_expenditure_pc-q*Rural_expenditure_pc)) + cov1*hdd + cov2*cdd +cov3*t2m+ int,
                  data = survey_data_train,
                  start = list(plateau = parm[3], m = parm[2], p = parm[2], q = parm[2], cov1 = 0, cov2 = 0, cov3 = 0,int = parm[1]),
                  control = list(maxiter = 1000))

  summary(NLS_model2)
  assign(paste0("NLS_model_",area_indx[areai],"1"),NLS_model1)
  assign(paste0("NLS_model_",area_indx[areai],"2"),NLS_model2)
}

NLS_cotable <- 1:2%>%purrr::map_dfr(function(areai){
  1:2%>%map_dfr(function(mi){
    NLS_model <- get(paste0("NLS_model_",area_indx[areai],mi))
    df_model = summary(NLS_model)$parameters%>%
      as.data.frame()
    colnames(df_model) <- c("Coef.", "Std. error", "t-stat.", "p")
    summary(NLS_model)
    NLS_LM <- lm(clean ~ predict(NLS_model), data = survey_data_train)%>%broom::glance()
    mat <- data.frame(df_model)%>%
      mutate(Area = area_indx[areai], Model = mi, Parm = row.names(.),
             p_star = ifelse(p < 0.01,"***", ifelse(p < 0.05, "**", ifelse(p < 0.1, "*",""))),
             RS = NLS_LM$r.squared, RS_adj = NLS_LM$adj.r.squared)
  })
})

fwrite(NLS_cotable,"result/Table/Decomposition/Model_comparison_NLS_model.csv")

NLS_fcon <- 1:2%>%map_dfr(function(areai){
  1:2%>%map_dfr(function(mi){
    NLS_model <- get(paste0("NLS_model_",area_indx[areai],mi))
    pre <- predict(NLS_model, prov_cov)
    pre1 <- predict(NLS_model, prov_cov%>%mutate_at(vars(hdd,cdd,t2m),function(x){0}))
    pre2 <- predict(NLS_model, prov_cov%>%mutate_at(vars(GDP_pc_from_CSMAR,Urban_expenditure_pc,Rural_expenditure_pc),function(x){10e3}))
    
    fuel_scen <- prov_cov%>%dplyr::select(Province, year)%>%
      left_join(fuel_usepop_lookup%>%dplyr::filter(Fuel == "Clean"))%>%
      mutate(clean_nls = pre, clean_nls1 =  pre1, clean_nls2 = pre2)%>%
      dplyr::rename(Popa = area_indx[areai])%>%group_by(year)%>%
      dplyr::summarise(f = sum(clean_nls*Popa)/sum(Popa),
                       f1 = sum(clean_nls1*Popa)/sum(Popa),
                       f2 = sum(clean_nls2*Popa)/sum(Popa))%>%ungroup()
    
    fuel_scenp <- fuel_scen%>%
      dplyr::filter(year%in%c(2000,2007,2013,2020))%>%
      mutate(Area = area_indx[areai], Model = mi,
             ECO_diff = f1 - lag(f1), MET_diff = f2 - lag(f2))%>%
      .[-1,]
    
    fuel_scen%>%
      dplyr::filter(year%in%c(2000,2020))%>%
      mutate(Area = area_indx[areai], Model = mi, year = 20002020,
             ECO_diff = f1 - lag(f1), MET_diff = f2 - lag(f2))%>%
      .[-1,]%>%full_join(fuel_scenp)
  })
})
fwrite(NLS_fcon, "result/Table/Decomposition/Model_comparison_NLS_contribution.csv")


##=====Comparison with OLS model=======
area_indx = c("urban","rural","total")
for(areai in 1:2){
  survey_data_train <- input_data%>%
    left_join(prov_info)%>%
    dplyr::rename(index = pindex)%>%
    filter(area == areai)%>%
    left_join(prov_cov)%>%
    dplyr::select(Province, year, clean, all_of(IV))%>%
    mutate(clean = clean*100)

  fm1 <- paste0("clean ~ ", paste(c("GDP_pc_from_CSMAR","hdd"), collapse = " + "))
  fm2 <- paste0("clean ~ ", paste(IV, collapse = " + "))

  LM_model1 <- lm(fm1,survey_data_train)
  LM_model2 <- lm(fm2,survey_data_train)
  assign(paste0("LM_model_",area_indx[areai],"1"),LM_model1)
  assign(paste0("LM_model_",area_indx[areai],"2"),LM_model2)

}
stargazer(LM_model_urban1, LM_model_urban2, LM_model_rural1, LM_model_rural2, order = IV,
          type = "html", single.row = TRUE, out = "result/Table/Decomposition/Model_comparison_LM_model.html")


LM_fcon <- 1:2%>%map_dfr(function(areai){
  1:2%>%map_dfr(function(mi){
    LM_model <- get(paste0("LM_model_",area_indx[areai],mi))
    pre <- predict(LM_model, prov_cov)
    pre1 <- predict(LM_model, prov_cov%>%mutate_at(vars(hdd,cdd,t2m),function(x){0}))
    pre2 <- predict(LM_model, prov_cov%>%mutate_at(vars(GDP_pc_from_CSMAR,Urban_expenditure_pc,Rural_expenditure_pc),function(x){0}))
    
    fuel_scen <- prov_cov%>%dplyr::select(Province, year)%>%
      left_join(fuel_usepop_lookup%>%dplyr::filter(Fuel == "Clean"))%>%
      mutate(clean_lm = pre, clean_lm1 =  pre1, clean_lm2 = pre2)%>%
      dplyr::rename(Popa = area_indx[areai])%>%group_by(year)%>%
      dplyr::summarise(f = sum(clean_lm*Popa)/sum(Popa),
                       f1 = sum(clean_lm1*Popa)/sum(Popa),
                       f2 = sum(clean_lm2*Popa)/sum(Popa))%>%ungroup()
    
    fuel_scenp <- fuel_scen%>%
      dplyr::filter(year%in%c(2000,2007,2013,2020))%>%
      mutate(Area = area_indx[areai], Model = mi,
             All_dff = f - lag(f), ECO_diff = f1 - lag(f1), MET_diff = f2 - lag(f2))%>%
      .[-1,]
    
    fuel_scen%>%
      dplyr::filter(year%in%c(2000,2020))%>%
      mutate(Area = area_indx[areai], Model = mi, year = 20002020,
             ECO_diff = f1 - lag(f1), MET_diff = f2 - lag(f2))%>%
      .[-1,]%>%full_join(fuel_scenp)
  })
})
fwrite(LM_fcon, "result/Table/Decomposition/Model_comparison_LM_contribution.csv")

##========HAP exposure comparison=========
library(terra)
library(stringr)
library(data.table)
library(dplyr)
library(tidyverse)
library(purrr)
not_any_na <- function(x) all(!(is.na(x)| x == ""))
prov_info <- fread("result/Table/Input_data/province_information.csv", encoding = 'UTF-8')%>%
  dplyr::rename(Region = Province)%>%
  dplyr::select(pindex, Region)

WHO_HAP <- fread("data/Global_hap_measurement/WHO2018_annual_HAP_China.csv")

PURE_HAP <- fread("data/Global_hap_measurement/PURE_HAP_China.csv")%>%
  # filter(Study == "PURE-AIR")%>%
  left_join(prov_info)%>%mutate(pindex = ifelse(is.na(pindex),0,pindex))%>%
  filter(Primary_Cooking_Fuel != "Animal dung")%>%
  mutate(Fuel = ifelse(Primary_Cooking_Fuel == "Gas"|Primary_Cooking_Fuel == "Electricity","Clean",
                       ifelse(Primary_Cooking_Fuel == "Coal","Coal","Biomass")))%>%
  dplyr::group_by(Study,Fuel,pindex)%>%
  dplyr::summarise(Concentration = weighted.mean(Annual_Concentration,N),
                   N = sum(N))%>%ungroup()

PURE_HAP%>%filter(Study == "PURE-AIR")
infl_info <- fread("result\\Table\\Input_data\\activity_information.csv")%>%
  filter(Activity == "Outdoor", Metric%in%c("HS.Mean","NHS.Mean"))%>%
  distinct(Province, pindex, Condition, Infl)

Pop_PURE <- fread(paste0('./result/Table/Population/Fuel_use_population_g0.1.csv'))%>%
  filter(year %in% c(2017,2018,2019))%>%
  group_by(x,y,year)%>%
  dplyr::summarise(Popall = sum(Pop*Popf))%>%ungroup()%>%
  mutate_at(c("x","y"),.funs = function(num, dgt = 2) num %>% round(dgt) %>% str_c)

Pop_PUREB <- fread(paste0('./result/Table/Population/Fuel_use_population_g0.1.csv'))%>%
  filter(year %in% 2005:2010)%>%
  group_by(x,y,year)%>%
  dplyr::summarise(Popall = sum(Pop*Popf))%>%ungroup()%>%
  mutate_at(c("x","y"),.funs = function(num, dgt = 2) num %>% round(dgt) %>% str_c)

Pop_WHO <- fread(paste0('./result/Table/Population/Fuel_use_population_g0.1.csv'))%>%
  filter(year %in% 2008:2016)%>%
  group_by(x,y,year)%>%
  dplyr::summarise(Popall = sum(Pop*Popf))%>%ungroup()%>%
  mutate_at(c("x","y"),.funs = function(num, dgt = 2) num %>% round(dgt) %>% str_c)

APM_PURE <- 2017:2019%>%
  map_dfr(~ rast(paste0("result\\Tif\\Match_data\\Matching_covariate_raster_g0.1_",.x,".tif"))%>%
            as.data.frame(xy = TRUE, na.rm =FALSE)%>%mutate(year = .x))%>%
  dplyr::select(x,y,APM25,hd_count,pindex,countyindex,year)%>%na.omit()%>%
  mutate_at(c("x","y"),.funs = function(num, dgt = 2) num %>% round(dgt) %>% str_c)%>%
  mutate(hd_f = hd_count/365)%>%
  inner_join(infl_info)%>%mutate(HAP_O = Infl*APM25)%>%
  dplyr::group_by(x,y,year)%>%
  dplyr::summarise_at(vars(APM25,HAP_O),
                      ~ .x[Condition == "H"]*hd_f[Condition == "H"] + 
                        .x[Condition == "N-H"]*(1-hd_f[Condition == "H"]))%>%ungroup()

APM_PUREB <- 2005:2010%>%
  map_dfr(~ rast(paste0("result\\Tif\\Match_data\\Matching_covariate_raster_g0.1_",.x,".tif"))%>%
            as.data.frame(xy = TRUE, na.rm =FALSE)%>%mutate(year = .x))%>%
  dplyr::select(x,y,APM25,hd_count,pindex,countyindex,year)%>%na.omit()%>%
  mutate_at(c("x","y"),.funs = function(num, dgt = 2) num %>% round(dgt) %>% str_c)%>%
  mutate(hd_f = hd_count/365)%>%
  inner_join(infl_info)%>%mutate(HAP_O = Infl*APM25)%>%
  dplyr::group_by(x,y,pindex,year)%>%
  dplyr::summarise_at(vars(APM25,HAP_O),
                      ~ .x[Condition == "H"]*hd_f[Condition == "H"] + 
                        .x[Condition == "N-H"]*(1-hd_f[Condition == "H"]))%>%ungroup()

APM_WHO <- 2008:2016%>%
  map_dfr(~ rast(paste0("result\\Tif\\Match_data\\Matching_covariate_raster_g0.1_",.x,".tif"))%>%
            as.data.frame(xy = TRUE, na.rm =FALSE)%>%mutate(year = .x))%>%
  dplyr::select(x,y,APM25,hd_count,pindex,countyindex,year)%>%na.omit()%>%
  mutate_at(c("x","y"),.funs = function(num, dgt = 2) num %>% round(dgt) %>% str_c)%>%
  mutate(hd_f = hd_count/365)%>%
  inner_join(infl_info)%>%mutate(HAP_O = Infl*APM25)%>%
  dplyr::group_by(x,y,year)%>%
  dplyr::summarise_at(vars(APM25,HAP_O),
                      ~ .x[Condition == "H"]*hd_f[Condition == "H"] + 
                        .x[Condition == "N-H"]*(1-hd_f[Condition == "H"]))%>%ungroup()

HAP_PURE <- left_join(APM_PURE,Pop_PURE)%>%
  dplyr::summarise_at(vars(APM25,HAP_O),
                      ~ weighted.mean(.x, Popall))%>%.[,"HAP_O"]

## remove the contributions from outdoor sources
HAP_PUREB <- left_join(APM_PUREB,Pop_PUREB)%>%
  group_by(pindex)%>%
  dplyr::summarise_at(vars(APM25,HAP_O),
                      ~ weighted.mean(.x, Popall))%>%ungroup()

HAP_WHO <- left_join(APM_WHO, Pop_WHO)%>%
  dplyr::summarise_at(vars(APM25,HAP_O),
                      ~ weighted.mean(.x, Popall))%>%.[,"HAP_O"]
PURE_HAP_I <- PURE_HAP%>%
  filter(Study == "PURE-AIR")%>%
  dplyr::select(Study,Fuel,Concentration)%>%
  mutate(Concentration = Concentration - as.numeric(HAP_PURE))

PUREB_HAP_I <- PURE_HAP%>%
  filter(Study == "PURE")%>%left_join(HAP_PUREB)%>%
  mutate(Concentration = Concentration - HAP_O)%>%
  group_by(Study,Fuel)%>%
  dplyr::summarise(Concentration = weighted.mean(Concentration, N))%>%ungroup()

WHO_HAP_I <- WHO_HAP[1,]%>%
  dplyr::rename(Biomass = Wood, Clean = Gas)%>%
  dplyr::select(Biomass, Clean, Coal)%>%
  mutate(Study = "WHO2018")%>%
  pivot_longer(cols = Biomass:Coal, names_to = "Fuel", values_to = "Concentration")%>%
  mutate(Concentration = Concentration - as.numeric(HAP_WHO))

list(WHO_HAP_I,PURE_HAP_I,PUREB_HAP_I)%>%reduce(full_join)%>%
  fwrite("data/Global_hap_measurement/Public_annual_HAP_China.csv")
