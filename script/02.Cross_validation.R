#==========Survey data==========
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
dat_dir <- "D:\\shaoyanchuan\\data\\"
library(data.table)
library(dplyr)
library(readr)
source(paste0(dat_dir,"Population\\Fuel_Use_WHO\\Code\\HAP_Data_Setup.R"))
source(paste0(dat_dir,"Population\\Fuel_Use_WHO\\Code\\HAP_Provincial_Model.R"))

survey_data <- readRDS(paste0(dat_dir, "Population/Chinese_Census/Fuel_type_survey_province_long_format.Rds"))


input_data <- fread("result/Table/Input_data/input_survey_data.csv", encoding = 'UTF-8')
prov_info <- fread("result/Table/Input_data/province_information.csv", encoding = 'UTF-8')%>%
  mutate(Province = whoname)
survey_data_lf <- fread("result/Table/Input_data/survey_data_lf.csv", encoding = 'UTF-8')

##  calculate urban population for each province
population <-  readRDS(paste0("data/panel_data/Chinese_provincial_panel_data.Rds"))%>%
  dplyr::select(Province, Year, Pop_all, Pop_U, Pop_O)%>%
  dplyr::rename(total = Pop_all, rural = Pop_O, urban = Pop_U)%>%
  mutate(urban_proportion = urban/total, County = Province)


## Prepare other covariates
IV = c("GDP_pc_from_CSMAR","hd","cd","t2m",
       "Urban_expenditure_pc","Rural_expenditure_pc")
prov_cov <- readRDS(paste0("data/panel_data/Chinese_provincial_panel_data.Rds"))%>%
  dplyr::select(Province, Year, all_of(IV))%>%
  left_join(prov_info,.)%>%
  mutate(loc_id = pindex)

fwrite(prov_cov, "result\\Table\\Input_data\\covariates_by_province.csv")


#==========Leave-one-year-out Cross validation==========
fuel_names = c('Solid','Biomass','Wood', 'Clean','Coal','Crop')
area_names = c('Urban','Rural')
yearseq = c(2000,2005,2010,2020)
plot_years = c(2000,2005,2010,2020)
for (y in yearseq) {
  cat("=====Year:",y,"=====\n")
  survey_data_train <- input_data%>%
    left_join(prov_info)%>%
    dplyr::rename(index = pindex)%>%
    filter(!year%in%y)

  survey_data_test <- input_data%>%
    left_join(prov_info)%>%
    dplyr::rename(index = pindex)%>%
    filter(year%in%y)

  aggreg_input <- data_aggreg(survey_data_train)

  model_data = aggreg_input%>%
    map(~ prov_cov%>%
          rename(country = whoname, year = Year)%>%
          dplyr::select(all_of(IV),country,year)%>%
          left_join(.x)%>%filter(!is.na(unique_country_index)))


  # all_samples$processed_samples$phi
  all_samples = run_model(model_data, IV,niter=4000,nburnin=1000,nchains=1,thin=1)


  province_samples = abind(foreach(i=1:nrow(prov_info))%do%{
    return(fuel_use_output(plot_id = prov_info$pindex[i],
                           plot_index = prov_info$pindex[i], cov_data = prov_cov,
                           plot_years = plot_years, all_samples = all_samples,
                           type='survey', output='samples',sampling_bias = TRUE))
  },along=5)
  
  n_sim = dim(all_samples$processed_samples$beta)[1]
  
  # Compute 2.5%, 50% and 97.5% quantiles.
  province_quantiles <- apply(province_samples, c(2,3,4,5),quantile,c(0.025,0.5,0.975),na.rm = TRUE)
  
  # province_quantiles[,11,1,1,1]
  dimnames(province_quantiles) <- list(quantile=c('lower','median','upper'),
                                       year=as.character(plot_years),
                                       Area=area_names, Fuel=fuel_names,
                                       Province = prov_info$whoname)

  data_joined <- province_quantiles%>%
    melt(value.name = "pred")%>%
    left_join(survey_data_lf,.)%>%
    na.omit()%>%filter(Year == y)%>%
    spread(quantile,pred)

  if(y == yearseq[1]){
    data_joined_tmp <- data_joined
  }else{
    data_joined_tmp <- full_join(data_joined_tmp, data_joined)
  }
}
data_joined_full <- data_joined_tmp

saveRDS(data_joined_full, "result\\Rdata\\Fuel_pop\\model_CV_fuel_use.Rds")

#==========Model Fitting==========
plot_years = 2000:2020
fuel_names = c('Solid','Biomass','Wood', 'Clean','Coal','Crop')
area_names = c('Urban','Rural')
IV = c("GDP_pc_from_CSMAR","hd","cd","t2m",
       "Urban_expenditure_pc","Rural_expenditure_pc")

survey_data_train <- input_data%>%
  left_join(prov_info)%>%
  rename(index = pindex)

aggreg_input <- data_aggreg(survey_data_train)

fit_model_data <- aggreg_input%>%
  purrr::map(~ prov_cov%>%
        rename(country = whoname, year = Year)%>%
        dplyr::select(all_of(IV),country,year)%>%
        left_join(.x)%>%filter(!is.na(unique_country_index)))


fit_all_samples <- run_model(fit_model_data, IV, niter=8000, nburnin=1000,nchains=2,thin=1)
saveRDS(fit_all_samples, "result\\Rdata\\Fuel_pop\\Fit_all_samples.Rds")

fit_all_samples <- readRDS("result\\Rdata\\Fuel_pop\\Fit_all_samples.Rds")
fit_province_samples=abind(foreach(i=1:nrow(prov_info))%do%{
  return(fuel_use_output(plot_id = prov_info$pindex[i],
                         plot_index = prov_info$pindex[i], cov_data = prov_cov,
                         plot_years = plot_years, all_samples = fit_all_samples,
                         type='survey', output='samples',sampling_bias = TRUE))
},along=5)

# Compute 2.5%, 50% and 97.5% quantiles.
fit_province_quantiles <- apply(fit_province_samples, c(2,3,4,5),quantile,c(0.025,0.5,0.975),na.rm = TRUE)

dimnames(fit_province_quantiles) <- list(quantile=c('lower','median','upper'),
                                    year=as.character(plot_years),
                                    Area=area_names, Fuel=fuel_names,
                                    Province = prov_info$whoname)
## Provincial predictions
fit_province_quantiles_df <- fit_province_quantiles%>%
  melt(value.name = "pred")%>%
  spread(quantile,pred)

saveRDS(fit_province_quantiles_df, "result\\Rdata\\Fuel_pop\\Fuel_use_province_prediction.Rds")

## Modeling fitting validation 
fit_province_quantiles_df
data_joined_fit <- left_join(survey_data_lf,fit_province_quantiles_df)%>%na.omit()
saveRDS(data_joined_fit, "result\\Rdata\\Fuel_pop\\model_Fit_fuel_use.Rds")

#==========CV and Model fitting results==========
##  read in CV and MF dataframe
data_joined_1 <- read_rds("result\\Rdata\\Fuel_pop\\model_CV_fuel_use.Rds")
data_joined_2 <- read_rds("result\\Rdata\\Fuel_pop\\model_fit_fuel_use.Rds")

save_name <- c("CV","Fit")
fuel_indx <- c("Clean","Biomass","Coal")
for (num in 1:2) {
  data_joined <- get(paste0("data_joined_", num))%>%
      group_by(year, Area, Province)%>%
      mutate(median = ifelse(Fuel == "Clean", 1 - median[Fuel == "Biomass"] - median[Fuel == "Coal"],median))%>%
      ungroup()
  for(m in 1:length(area_names)){
    for(p in 1:length(fuel_indx)){
      data_joined_df <- data_joined%>%
        filter(Area == area_names[m], Fuel == fuel_indx[p])%>%
        dplyr::select(value, median)

      corr_value <- (cor(data_joined_df)[1,2])^2
      
      cat("####Area:",area_names[m]," Fuel:",fuel_indx[p]," R2 =",corr_value,save_name[num],"#####\n")
      corr_df <- data.frame(Area = area_names[m], Fuel = fuel_indx[p], corr = corr_value)
      if(m==1&p==1){
        corr_df_tmp <- corr_df
      }else{
        corr_df_tmp <- full_join(corr_df_tmp, corr_df)%>%suppressMessages()
      }
    }
  }
  corr_df_all  <- corr_df_tmp
  saveRDS(corr_df_all, paste0("result\\Rdata\\Fuel_pop\\",save_name[num],"_fuel_use_cor_by_all.Rds"))
  for (y in yearseq) {
    for(m in 1:length(area_names)){
      for(p in 1:length(fuel_names)){
        data_joined_df <- data_joined%>%
          filter(Area == area_names[m], Fuel == fuel_names[p], year == y)%>%
          dplyr::select(value, median)
        
        corr_value <- (cor(data_joined_df)[1,2])^2
        corr_df <- data.frame(Area = area_names[m], Fuel = fuel_names[p], corr = corr_value, year = y)
        if(m==1&p==1&y==yearseq[1]){
          corr_df_tmp <- corr_df
        }else{
          corr_df_tmp <- full_join(corr_df_tmp, corr_df)%>%suppressMessages()
        }
      }
    }
  }
  corr_df_by_year  <- corr_df_tmp
  saveRDS(corr_df_by_year, paste0("result\\Rdata\\Fuel_pop\\",save_name[num],"_fuel_use_cor_by_year.Rds"))
}

#==========Visualization of scatter plot==========
library(egg)
library(tagger)
library(Cairo)
library(grid)

fuel_use_scatter <- function(data, cor){
  xlim = ylim = c(0,1)
  cv_plot <- ggplot(data_to_plotf, aes(y, ypred, color = Area))+
    geom_point(size = 1,  alpha = 0.5,shape = 20)+
    facet_grid(Area ~ Fuel, labeller = labeller(groupwrap = label_wrap_gen(10)))+
    theme(strip.background = element_blank())+ ##change the facetwrap background
    # theme(strip.background = element_rect(color="black", fill="white", size = 0.5))+ ##change the facetwrap background
    theme(strip.text.x = element_text(size = 28),
          strip.text.y = element_text(size = 28, angle = 0))+
    coord_cartesian(xlim = xlim, ylim = ylim) +
    tag_facets(tag_pool = (corr_to_plotf$tag), tag_suffix = "", 
               position = list(x = .005, y = .825, hjust = 0))+
    theme(tagger.panel.tag.text = element_text(color = "black", size = 26),
          tagger.panel.tag.background = element_rect(color = "black", size = 2))+
    geom_abline(slope = 1,linetype = "dashed",
                intercept = 0, colour = "black",size = 0.1)+
    theme(panel.background = element_rect(linewidth = 0.5, fill = NA, colour = 'black'),
          panel.border = element_rect(linewidth = 1,colour = "black", fill = NA))+
    labs(x = expression(atop("Observed Fuel Use")),
         y = expression("Predicted Fuel Use"))+
    guides(color = "none")+
    scale_color_manual(values = c("#CE161F","#0A4990"))+ ## filled color of different areas
    scale_x_continuous(limits = c(0, 1))+
    scale_y_continuous(limits = c(0, 1))+
    theme(axis.text = element_text(size = 20),
          axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm"), size = 30),
          axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm"), size = 30),
          panel.spacing.x = unit(8, "mm") , panel.spacing.y = unit(8, "mm"))

  cv_plot
}

fuel_indx <- c("Clean","Biomass","Coal")
save_name <- c("CV","Fit")

for (num in 1:2) {
  data_to_plot <- read_rds(paste0("result\\Rdata\\Fuel_pop\\model_",save_name[num],"_fuel_use.Rds"))%>%
    group_by(year, Area, Province)%>%
    mutate(median = ifelse(Fuel == "Clean", 1 - median[Fuel == "Biomass"] - median[Fuel == "Coal"],median))%>%
    ungroup()%>%filter(Area%in%c("Urban","Rural"), Fuel%in%fuel_indx)%>%
    rename(y = value, ypred = median)%>%
    mutate(Area = factor(Area, levels = c("Urban","Rural")),
           Fuel = factor(Fuel, levels = c("Clean","Biomass","Coal")))%>%
    as.data.frame()

  corr_to_plot <- read_rds(paste0("result\\Rdata\\Fuel_pop\\",save_name[num],"_fuel_use_cor_by_all.Rds"))%>%
    filter(Area%in%c("Urban","Rural"), Fuel%in%fuel_indx)%>%
    mutate(corr = formatC(sqrt(corr), digits = 2, width = 2, format = "f"),
           Area = factor(Area, levels = c("Urban","Rural")),
           Fuel = factor(Fuel, levels = c("Clean","Biomass","Coal")) )%>%
    arrange(Fuel,Area)%>%
    mutate(number = 1:nrow(.), alpha = LETTERS[number], 
           tag = paste0("(",alpha,")","<br><br>R = ",corr))
  
  
  data_to_plotf <- filter(data_to_plot, Fuel%in%fuel_indx)
  corr_to_plotf <- filter(corr_to_plot, Fuel%in%fuel_indx)
  
  cv_plot <- fuel_use_scatter(data_to_plotf, corr_to_plotf)

  CairoPNG(paste0("result/Fig/Fuel_use/Fuel_use_",save_name[num],"_validation.png"),width = 2500,height = 1900,res = 200)
  print(cv_plot)
  dev.off()
}

#==========Calculate and visualize overall fuel use fraction during 2000-2020==========
area_indx = c("urban","rural","total")
fuel_indx = c("Clean","Biomass","Coal")

prov_cov <- readRDS(paste0("data/panel_data/Chinese_provincial_panel_data.Rds"))%>%
  dplyr::select(Province, Year, all_of(IV))%>%
  left_join(prov_info,.)%>%
  mutate(loc_id = pindex)

fuel_use_prov <- read_rds("result\\Rdata\\Fuel_pop\\Fuel_use_province_prediction.Rds")%>%
  group_by(year, Area, Province)%>%
  mutate(median = ifelse(Fuel == "Clean", 1 - median[Fuel == "Biomass"] - median[Fuel == "Coal"],median))%>%
  ungroup()%>%filter(Fuel%in%fuel_indx)

population_prov <- fread("result\\Table\\Input_data\\province_population.csv")%>%
  dplyr::select(-Province)%>%dplyr::rename(Province = whoname,year = Year)

## group_by 31 provinces
population_overall <- population_prov%>%
  filter(!Province%in%c("台湾","澳门","香港"))%>%
  group_by(year)%>%dplyr::summarise(total = sum(total),
  urban = sum(urban), rural = sum(rural))%>%ungroup()

fuel_usepop_prov <- fuel_use_prov%>%
  dplyr::select(-lower,-upper)%>%
  pivot_wider(names_from = c(Area), values_from = median)%>%
  left_join(population_prov, by = c("year","Province"))

fuel_usepop_lookup <- fuel_usepop_prov%>%
  left_join(prov_cov%>%dplyr::rename(year = Year))

vname <- c("GDP_pc_from_CSMAR")
fuel_usepop_lookup%>%
  filter(Fuel == "Clean", year > 2017)%>%
  dplyr::select(year,Province,all_of(vname))%>%
  pivot_wider(names_from = year, values_from = all_of(vname))%>%
  filter(`2019`-`2018`<0)

write.csv(fuel_usepop_lookup, "result/Table/Population/Fuel_use_population_by_province.csv",
          fileEncoding = "GBK", row.names = FALSE)


##========Uncertainty analysis========
source('./script/Core_MonteCarlo.R',encoding = 'UTF8')
## fuel use fraction model
model_samples <- readRDS("result\\Rdata\\Fuel_pop\\Fit_all_samples.Rds")
prov_cov <- fread("result\\Table\\Input_data\\covariates_by_province.csv", encoding = "UTF-8")
prov_info <- fread("result\\Table\\Input_data\\province_information.csv", encoding = "UTF-8")%>%
  mutate(index = pindex)

## total simulations generated from MCMC chain
n_sim = dim(model_samples$processed_samples$beta_cov)[1]
sample_index = sample(1:n_sim, 1000)
## Uncertainty of fuel use fraction and variable contribution
tic("Fuel use fraction simulation")
Fuelf_simu = Fuel_pop_mtcl(prov_cov, 2000:2020, sample_index, prov_info)
toc()

fuel_use_mtcl <- function(samp = 1, total = 1000, clean = FALSE){
  ## provincial population

  Pop_survey <- survey_data%>%dplyr::rename(whoname = Province,year = Year,Area = Population_type)%>%
    left_join(prov_info)%>%filter(Area !="all")%>%
    dplyr::select(Province,year,Area,sum)%>%
    mutate(Area = ifelse(Area == "urban","Urban",ifelse(Area == "rural", "Rural","Overall")))
  Pop_prov <- fread("result\\Table\\Input_data\\province_population.csv")%>%
    group_by(Province,Year)%>%
    dplyr::rename(Urban = urban, Rural = rural,year = Year)%>%
    dplyr::select(Province, year ,Urban, Rural)%>%
    pivot_longer(c(Urban,Rural), names_to = "Area", values_to = "Pop")%>%
    left_join(Pop_survey)%>%
    mutate(Pop = ifelse(is.na(sum),Pop,sum))
  
  ## national population with uncertainty analysis
  Fuel_nsample <- Fuelf_simu%>%
    dplyr::rename(sample_id = sim, pindex = index)%>%filter(sample_id == samp)%>%
    left_join(prov_info)%>%left_join(Pop_prov)%>%
    group_by(year, Area, Fuel)%>%
    dplyr::summarise(Popf = weighted.mean(Popf, Pop), Pop = sum(Pop))%>%ungroup()
  
  Fuel_nsample%>%
    group_by(year, Fuel)%>%
    dplyr::summarise(Popf = weighted.mean(Popf, Pop), Pop = sum(Pop))%>%ungroup()%>%
    mutate(Area = "Overall")%>%full_join(Fuel_nsample,.)
}

tic("Sim")
plan(multisession, workers = 10)
options(future.globals.maxSize= 891289600)
Fuel_nsim = 1:1000%>%
  future_map_dfr(function(x){
    fuel_use_mtcl(samp = x, total = 1000)
  },.options = furrr_options(seed = NULL))
toc()
fwrite(Fuel_nsim%>%dplyr::select(-Pop), "result/Table/Population/Fuel_use_fraction_uncertainty(national).csv")


#==========Fuel use fraction predictions and uncertainties during 2000-2020==========
fuel_usepop_lookup <- fread("result/Table/Population/Fuel_use_population_by_province.csv")
## national fuel use fraction
fuel_usepop_lookup%>%
  group_by(year,Fuel)%>%
  dplyr::summarise(Total = sum(urban*Urban+rural*Rural)/sum(total))%>%ungroup()%>%
  pivot_wider(names_from = Fuel, values_from = Total)%>%
  filter(year%in%c(2000,2005,2010,2020))%>%as.data.frame()

fuel_use_overall <- fuel_usepop_lookup%>%
  group_by(year, Fuel)%>%
  dplyr::summarise(total_fpop = sum(Urban*urban+Rural*rural),
            urban_fpop = sum(Urban*urban),
            rural_fpop = sum(Rural*rural))%>%
  ungroup()%>%left_join(population_overall, by = c("year"))%>%
  mutate(total_fpopf = total_fpop/total,urban_fpopf = urban_fpop/urban,
         rural_fpopf = rural_fpop/rural,Fuel = factor(Fuel, levels = fuel_indx))

fuel_use_lf <- fuel_use_overall%>%
  dplyr::select(year, Fuel, paste0(area_indx,"_fpopf"))%>%
  setnames(c("year", "Fuel", "Urban", "Rural", "Overall"))%>%
  pivot_longer(cols = Urban:Overall, names_to = "Area")%>%
  mutate(Area = factor(Area, levels = c("Urban","Rural","Overall")),
         Index = ifelse(Area == "Urban","a",ifelse(Area == "Rural","b","c")))

fit_pp <- ggplot(fuel_use_lf)+
  facet_wrap("Index", nrow = 1, labeller = labeller(groupwrap = label_wrap_gen(20)))+
  geom_density(aes(x = year, y = value, fill = Fuel),
               stat="identity", position="stack", color = "white")+
  labs(x= NULL, y = "Fuel  Use  Proportion")+
  scale_fill_brewer(name = NULL, palette = "Spectral")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(size = 1.2,fill='transparent'),
        panel.border = element_rect(size = 1.2,colour = "black", fill=NA),
        panel.spacing.x = unit(8, "mm") , panel.spacing.y = unit(8, "mm"))+
  theme(strip.background = element_blank(), ##change the facetwrap background
        strip.text.x = element_text(size = 35, color = "black", face = "bold",
                                    margin = margin(0, 0, .5, 0, "cm")))+
  theme(axis.text = element_text(colour = "black",size = 25), ##axis text size
        axis.title = element_text(colour = "black",size = 30))+ ##axis title size
  theme(legend.position = "bottom",
        legend.text = element_text(colour = "black", size = 25), ##legend text size
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "white"),
        legend.key.width=unit(1, "cm"),
        legend.key.height=unit(1, "cm"))

CairoPNG(paste0("result\\Fig\\Fuel_use\\Fit_prediction_multiyear.png"),height = 1000,width = 3400, res = 150)
print(fit_pp)
dev.off()

Fuel_uc <- fread("result/Table/Population/Fuel_use_fraction_uncertainty(national).csv")%>%
  group_by(year, Area, Fuel)%>%
  dplyr::summarise(lower = quantile(Popf, 0.025), 
            mean = mean(Popf),
            upper = quantile(Popf, 0.975))%>%ungroup()
Survey_obs <- fread("result/Table/Population/National_survey_fuel_use_fraction.csv")%>%
  dplyr::select(Year, Area, Fuel, Survey)%>%dplyr::rename(year = Year)
Clean_uc <- Fuel_uc%>%
  filter(Fuel == "Clean")%>%
  left_join(Survey_obs)%>%
  mutate(Area = factor(Area, levels = c("Urban","Rural","Overall")))

Clean_uc_plot <- ggplot(Clean_uc, aes(x = year)) +
  geom_point(aes(y = Survey/100, color = Area), size = 3) +
  geom_line(aes(y = mean, color = Area))+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Area),alpha = 0.4)+
  scale_fill_manual(name = NULL,values=c("red3", "royalblue3", "palegreen3")) +
  scale_color_manual(name = NULL, values=c("red3", "royalblue3", "palegreen3")) +
  labs(x = NULL, y = "Clean Fuel Use Proportion") +
  theme_bw()+
  theme(
        panel.grid.major = element_line(size = .8,colour = "grey", linetype = 2),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(size = 1.2,fill='transparent'),
        panel.border = element_rect(size = 1.2,colour = "black", fill=NA))+
  theme(axis.text=element_text(colour = "black",size = 25),
        axis.title=element_text(colour = "black",size = 30, margin = margin(0,20,0,0)),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        legend.text=element_text(colour = "black",size = 25),
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.key.width=unit(1, "cm"),
        legend.key.height=unit(1, "cm"),
        legend.spacing.y = unit(.8, 'cm'),
        legend.spacing.x = unit(.5, 'cm'))+
  guides(color = guide_legend(byrow = TRUE), title.hjust = 1)
ggsave("result/Fig/Fuel_use/Fuel_use_fraction_uncertainty.png",Clean_uc_plot,
       width = 10, height = 7)

##======Explore fuel use change in regions with different economic status======
fuel_usepop_lookup <- fread("result/Table/Population/Fuel_use_population_by_province.csv")
## find the provinces with highest and lowest GDP
Prov_by_ECO <- fuel_usepop_lookup%>%
  group_by(Province, pindex)%>%
  summarise_at(vars(GDP_pc_from_CSMAR, Urban_expenditure_pc, Rural_expenditure_pc), mean)%>%
  ungroup()

library(pinyin)
mypy <- pydic(method = 'toneless', dic = c("pinyin2")) ## convert Chinese character into English

fuel_usepop_sr <- fuel_usepop_lookup%>%
  filter((pindex %in%c(4,7,24,2))&Fuel == "Clean")%>%as.data.frame()%>%
  mutate(Total = (Urban*urban + Rural*rural)/total,
         pcGDP = GDP_pc_from_CSMAR/1e3,
         Province = str_to_title(gsub("_","",py(whoname, dic = mypy))))%>%
  mutate(Province = factor(Province, levels = c("Gansu","Guizhou","Beijing","Shanghai")))

# Modify the code to apply same colour to fill and color arguments
GDP_fuel_use <- ggplot(fuel_usepop_sr, aes(x = pcGDP, y = Total, color = Province,fill = Province)) +
  geom_point() +
  geom_smooth() +
  labs(x = "pcGDP (Thousand Chinese yuan/person)", y = "Clean Fuel Use Proportion") +
  scale_color_viridis_d(name = "") +
  scale_fill_viridis_d(name = "") + #Add this line to apply same color to fill and color arguments
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(size = 1.2,fill='transparent'),
        panel.border = element_rect(size = 1.2,colour = "black", fill=NA))+
  theme(axis.text=element_text(colour = "black",size = 25),
        axis.title=element_text(colour = "black",size = 30, margin = margin(0,20,0,0)),
        legend.text=element_text(colour = "black",size = 25),
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.key.width=unit(1, "cm"),
        legend.key.height=unit(1, "cm"),
        legend.spacing.y = unit(.8, 'cm'),
        legend.spacing.x = unit(.5, 'cm'))+
  guides(color = guide_legend(byrow = TRUE), title.hjust = 1)
ggsave("result/Fig/Fuel_use/GDP_with_fuel_use_trend.png",GDP_fuel_use,
       width = 10, height = 7)
