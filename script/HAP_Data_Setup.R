##############################
## Household Air Pollution: ##
##   A Hierarchical Model   ##
##############################
# setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
dat_dir <- "D:\\shaoyanchuan\\data\\"

# data = survey_data_train
data_aggreg <- function(data){
  filtered_data_1 <- data
  # Divide fuel observations by their sum or by 1-(no response+unlisted fuel+no cooking),
  # whichever is greater. This ensures that the sum of individual fuel values can't exceed one.
  scaling_factor <- apply(cbind(filtered_data_1$total_fuels,1-filtered_data_1$total_others),1,max)
  
  filtered_data_1$solid_scaled=filtered_data_1$solid_2/scaling_factor
  filtered_data_1$biomass_scaled=filtered_data_1$biomass/scaling_factor
  filtered_data_1$total_biomass_scaled=filtered_data_1$total_biomass/scaling_factor
  filtered_data_1$wood_scaled=filtered_data_1$wood/scaling_factor
  filtered_data_1$cropwaste_scaled=filtered_data_1$cropwaste/scaling_factor
  filtered_data_1$dung_scaled=filtered_data_1$dung/scaling_factor
  filtered_data_1$charcoal_scaled=filtered_data_1$charcoal/scaling_factor
  filtered_data_1$coal_scaled=filtered_data_1$coal/scaling_factor
  filtered_data_1$kerosene_scaled=filtered_data_1$kerosene/scaling_factor
  filtered_data_1$electricity_scaled=filtered_data_1$electricity/scaling_factor
  filtered_data_1$gas_scaled=filtered_data_1$gas/scaling_factor
  filtered_data_1$clean_scaled=filtered_data_1$clean/scaling_factor
  
  
  
  filtered_data_1$other_scaled = 1-(filtered_data_1$biomass_scaled+
                                      filtered_data_1$charcoal_scaled+
                                      filtered_data_1$coal_scaled+
                                      filtered_data_1$electricity_scaled+
                                      filtered_data_1$gas_scaled+
                                      filtered_data_1$kerosene_scaled)
  
  # Filter out surveys above the 'incompleteness' threshold.
  filtered_data_2 <- filtered_data_1
  n_id <- length(unique(filtered_data_2$id)) # Number of unique surveys.
  
  # Put the data into a list of equal size data frames for modelling.
  empty <- data.frame(id=rep(NA,n_id),country=rep('NA',n_id),country_index=rep(NA,n_id),unique_country_index=rep(NA,n_id),region=rep('NA',n_id),
                      unique_region_index=rep(NA,n_id),region_index=rep(NA,n_id),super_region=rep('NA',n_id),
                      unique_super_region_index=rep(NA,n_id),super_region_index=rep(NA,n_id),
                      year=rep(NA,n_id),respondents=rep(NA,n_id),
                      solid_scaled=rep(NA,n_id),solid_count=rep(NA,n_id),kerosene_scaled=rep(NA,n_id),kerosene_count=rep(NA,n_id),
                      gas_scaled=rep(NA,n_id),gas_count=rep(NA,n_id),electricity_scaled=rep(NA,n_id),electricity_count=rep(NA,n_id),
                      biomass_scaled=rep(NA,n_id),biomass_count=rep(NA,n_id),charcoal_scaled=rep(NA,n_id),charcoal_count=rep(NA,n_id),
                      wood_scaled=rep(NA,n_id),wood_count=rep(NA,n_id),cropwaste_scaled=rep(NA,n_id),cropwaste_count=rep(NA,n_id),
                      coal_scaled=rep(NA,n_id),coal_count=rep(NA,n_id),dung_scaled=rep(NA,n_id),dung_count=rep(NA,n_id),
                      clean_scaled = rep(NA,n_id),clean_count = rep(NA,n_id))
  
  unique_country_index = unique(filtered_data_2$index)
  
  model_data=list(urban=empty,rural=empty,overall=empty)
  
  model_data[[1]]$id=model_data[[2]]$id=model_data[[3]]$id=unique(filtered_data_2$id)
  
  N = 100000
  
  
  for(j in 1:3){
    for(i in 1:n_id){
      data_row=which(filtered_data_2$id==model_data[[j]]$id[i]&filtered_data_2$area==j)
      if(!is.na(data_row[1])){
        if(length(data_row)>1){
          data_row = data_row[filtered_data_2$`Household or population weighting`[data_row]=='p']
        }
        if(is.na(data_row[1])){
          print('Check: Duplicate survey and unable to select population weighting only.')
          browser()
        }else{
          model_data[[1]]$country[i]=model_data[[2]]$country[i]=
            model_data[[3]]$country[i]=filtered_data_2$whoname[data_row]
          
          model_data[[1]]$country_index[i]=model_data[[2]]$country_index[i]=
            model_data[[3]]$country_index[i]=filtered_data_2$index[data_row]
          
          ## incase inconsistent country number
          ind_sel = which(unique_country_index==filtered_data_2$index[data_row])

          model_data[[1]]$unique_country_index[i]=model_data[[2]]$unique_country_index[i]=
            model_data[[3]]$unique_country_index[i]= unique_country_index[ind_sel]

          
          model_data[[1]]$year[i]=model_data[[2]]$year[i]=
            model_data[[3]]$year[i]=filtered_data_2$year[data_row]
          
          
          model_data[[j]]$electricity_scaled[i]=filtered_data_2$electricity_scaled[data_row]
          model_data[[j]]$electricity_count[i]=floor(filtered_data_2$electricity_scaled[data_row]*N)
          model_data[[j]]$kerosene_scaled[i]=filtered_data_2$kerosene_scaled[data_row]
          model_data[[j]]$kerosene_count[i]=floor(filtered_data_2$kerosene_scaled[data_row]*N)
          model_data[[j]]$solid_scaled[i]=filtered_data_2$solid_scaled[data_row]
          model_data[[j]]$solid_count[i]=floor(filtered_data_2$solid_scaled[data_row]*N)
          model_data[[j]]$gas_scaled[i]=filtered_data_2$gas_scaled[data_row]
          model_data[[j]]$gas_count[i]=floor(filtered_data_2$gas_scaled[data_row]*N)
          model_data[[j]]$coal_scaled[i]=filtered_data_2$coal_scaled[data_row]
          model_data[[j]]$coal_count[i]=floor(filtered_data_2$coal_scaled[data_row]*N)
          model_data[[j]]$biomass_count[i]=floor(filtered_data_2$biomass_scaled[data_row]*N)
          model_data[[j]]$charcoal_scaled[i]=filtered_data_2$charcoal_scaled[data_row]
          if(is.na(model_data[[j]]$charcoal_scaled[i])&!is.na(model_data[[j]]$biomass_scaled[i])&!is.na(model_data[[j]]$coal_scaled[i])){
            model_data[[j]]$charcoal_scaled[i] <- round(model_data[[j]]$solid_scaled[i]-model_data[[j]]$biomass_scaled[i]-model_data[[j]]$coal_scaled[i],8) # Sometimes goes negative due to roundoff.
          }
          model_data[[j]]$charcoal_count[i]=floor(model_data[[j]]$charcoal_scaled[i]*N)
          model_data[[j]]$wood_scaled[i]=filtered_data_2$wood_scaled[data_row]
          model_data[[j]]$wood_count[i]=floor(filtered_data_2$wood_scaled[data_row]*N)
          model_data[[j]]$dung_scaled[i]=filtered_data_2$dung_scaled[data_row]
          model_data[[j]]$dung_count[i]=floor(filtered_data_2$dung_scaled[data_row]*N)
          model_data[[j]]$cropwaste_scaled[i]=filtered_data_2$cropwaste_scaled[data_row]
          if(is.na(model_data[[j]]$cropwaste_scaled[i])&!is.na(model_data[[j]]$dung_scaled[i])&!is.na(model_data[[j]]$wood_scaled[i])){
            model_data[[j]]$cropwaste_scaled[i] <- round(model_data[[j]]$biomass_scaled[i]-model_data[[j]]$dung_scaled[i]-model_data[[j]]$wood_scaled[i],8)
          }
          model_data[[j]]$cropwaste_count[i]=floor(model_data[[j]]$cropwaste_scaled[i]*N)
          
          model_data[[j]]$clean_scaled[i]=filtered_data_2$clean_scaled[data_row]
          model_data[[j]]$clean_count[i] = floor(model_data[[j]]$clean_scaled[i]*N)
        }}
    }
  }
  return(model_data)
}
