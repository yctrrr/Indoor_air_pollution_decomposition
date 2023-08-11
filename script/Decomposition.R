## Multiple decomposition functions

##======Decomposition of overall morality======
library(combinat)
scenario_extract <- function(var, start.y, end.y, Popf, Pop, AgeGroup, PM_real, mRate){
  if("PA"%in%var){
    AgeGroupS <- filter(AgeGroup, year == end.y)
  }else{
    AgeGroupS <- filter(AgeGroup, year == start.y)
  }
  
  if("ORF"%in%var){
    mRateS <- filter(mRate, year == end.y)
    PM_rS <- filter(PM_real, year == end.y)
  }else{
    mRateS <- filter(mRate, year == start.y)
    PM_rS <- filter(PM_real, year == start.y)
  }

  if("EXP"%in%var){
    PM_cS <- filter(PM_real, year == end.y)
  }else{
    PM_cS <- filter(PM_real, year == start.y)
  }
  
  if("FF"%in%var){
    PopfS <- filter(Popf, year == end.y)%>%dplyr::select(-ECO,-MET)
  }else{
    PopfS <- filter(Popf, year == start.y)%>%dplyr::select(-ECO,-MET)
  }
  
  if("PC"%in%var){
    PopS <- filter(Pop, year == end.y)
  }else{
    PopS <- filter(Pop, year == start.y)
  }
  
  
  Slist <- list(PopfS,PopS,AgeGroupS,PM_cS,mRateS,PM_rS)
  names(Slist) <- c("FF",'PC','PA','EXP','ORF',"PMr")
  return(Slist)
}

## calculate contributions for each feature
## num: variable position
posi_contrib <- function(num, vname,MortS,FF_scale){
  ## For specific variable position(comb1-5), all possible permutations given sequence
  seqnum = c(24,6,4,6,24)
  # x=3
  comb <- filter(MortS, comb_id == num)$comb
  comb_former <- filter(MortS, comb_id == num-1)$comb
  Mortlst <- filter(MortS, comb_id == num)$Mortlst
  Mortlst_former <- filter(MortS, comb_id == num-1)$Mortlst

  VShape <- 1:length(comb)%>%
    purrr::map(function(x){
      seq <- comb[[x]]%>%setdiff(vname)
      if(length(seq)==0){seq = ""}
      
      ## extract the former scenario index
      seq_former <- comb_former%>%purrr::map_lgl(~ all(seq%in%.x))
      
      ## if any True, then calculate the contribution
      if(any(seq_former)){
        if(vname != "FF"){
          left_join(Mortlst[[x]]%>%dplyr::rename(MortE = Mort),
                    Mortlst_former[[which(seq_former)]]%>%dplyr::rename(MortS = Mort))%>%
            mutate(Mort = MortE - MortS, Decomp = vname, .keep = "unused")%>%suppressMessages()
        }else{
          left_join(Mortlst[[x]]%>%dplyr::rename(MortE = Mort),
                    Mortlst_former[[which(seq_former)]]%>%dplyr::rename(MortS = Mort))%>%
            left_join(FF_scale)%>%
            mutate(ECO = ECO*MortS, FF = MortE - MortS, MET = FF - ECO, .keep ="unused")%>%
            pivot_longer(c(ECO,MET,FF), names_to = "Decomp", values_to = "Mort")%>%
            suppressMessages()
          
        }
      }
    })
  
  VShape %>% purrr::discard(is.null)%>%imap_dfr(
    ~ .x %>% mutate(Step = .y))%>%group_by(x,y,Area,Fuel,Decomp)%>%
    dplyr::summarise(Mort = seqnum[num]*sum(Mort)/120)%>%ungroup()%>%suppressMessages()
}

##  Mortality decomposition
Decomposition <- function(start.y, end.y, RR, Grids, Popf, Pop,
                          AgeGroup, PM_real, mRate) {
  serie = c("FF",'PC','PA','EXP','ORF')
  ##  scaled value of fuel use fraction
  # FF_scale%>%filter(x==4)
  FF_scale <- filter(Popf, year == start.y)%>%
    mutate(ECO = (ECO/Popf),MET = (MET/Popf))%>%
    dplyr::select(x,y,Area,Fuel,ECO,MET)

  ## calculate mortalities for all scenarios
  MortS <- 0:5%>%purrr::map_dfr(function(i){
    ## all scenarios when selecting i th variable
    comb_all <- if(i == 0){list("")}else{combn(serie, i, simplify = FALSE)}
    
    Mortlst <- comb_all%>%
      purrr::map(function(x){
        
        Slist <- scenario_extract(x, start.y, end.y, Popf, Pop, AgeGroup, PM_real, mRate)
       
        Pop_group <- left_join(Slist$FF%>%dplyr::select(-year),
                               Slist$PC%>%dplyr::select(-year))%>%
          mutate(Pop = Pop*Popf, .keep = "unused")
        
        Mort_Scomb <- Mortality(
          RR = RR, Grids = Grids, pop = Pop_group,
          ag = Slist$PA%>%dplyr::select(-year), PM_c = Slist$EXP%>%dplyr::select(-year),
          PM_r = Slist$PMr%>%dplyr::select(-year), mRate = Slist$ORF%>%dplyr::select(-year))

      })%>%suppressMessages()
    
    tibble(comb = comb_all, Mortlst = Mortlst)%>%
      mutate(comb_id = i)
  })

  ## contribution of each variable
  VShape_mean <- 1:length(serie)%>%
    purrr::map_dfr(function(i){
      ## calculate separate contributions of one variable in 5 different positions
      1:5%>%
        purrr::map_dfr(~ posi_contrib(.x, serie[i],MortS,FF_scale)%>%mutate(posi = .x))%>%
        group_by(x,y,Area,Fuel,Decomp)%>%dplyr::summarise(Mort = sum(Mort))%>%ungroup()%>%
        suppressMessages()
    })
  

  Decompf <- list((filter(MortS, comb_id == 0)$Mortlst[[1]])%>%mutate(Decomp = "Start"),
                  (filter(MortS, comb_id == 5)$Mortlst[[1]])%>%mutate(Decomp = "End"),
                  VShape_mean)%>%reduce(full_join)%>%arrange(x,y,Area,Fuel)
  # Decompf%>%filter(x==4,Decomp%in%c("Start","FF","ECO","MET","End"))
  return(Decompf)
}


##======Decomposition of fuel use fraction change======
Fuel_pop <- function(cov_step, plot_years, info){
  fuel_names = c('Solid','Biomass','Wood', 'Clean','Coal','Crop')
  area_names = c('Urban','Rural')

  ##  calculate fuel use fraction
  plan(multisession, workers = 20)
  options(future.globals.maxSize= 891289600)
  fuel = 1:nrow(info)%>%
    future_map_dfr(function(i){
      fuel_usec = fuel_use_output(plot_id = info$index[i],
                                  plot_index = info$pindex[i], cov_data = cov_step,
                                  plot_years = plot_years, all_samples = model_samples, 
                                  type='survey', output='samples',sampling_bias = TRUE)%>%
        apply(c(2,3,4),quantile,c(0.5),na.rm = TRUE)%>%suppressMessages()

      dimnames(fuel_usec) <- list(year = plot_years, Area = area_names, Fuel = fuel_names)

      reshape2::melt(fuel_usec)%>%
        as_tibble()%>%filter(Fuel%in%c("Biomass","Coal"))%>%
        pivot_wider(names_from = Fuel, values_from = value)%>%
        mutate(index = info$index[i])%>%suppressMessages()

    },.options = furrr_options(seed = NULL))

  ##  long format
  fuel%>%
    mutate(Clean = 1 - Biomass - Coal)%>%
    pivot_longer(c("Biomass","Coal","Clean"), names_to = "Fuel", values_to = "Popf")
  
}

Pop_Decomposition <- function(deco_years, cov, info, Popfy) {
  ##  variable names
  METv = c("hd","cd","t2m")
  ECOv = c("GDP_pc_from_CSMAR","Urban_expenditure_pc","Rural_expenditure_pc")
  decoendy = deco_years[length(deco_years)-1]
  
  ##  create scenario index
  scenario_table = expand.grid(year = deco_years[1:(length(deco_years)-1)], 
                               serie = c(1.1,1.2)%>%as.character())%>%
    mutate(plot_year = paste0(year,"_",serie))
  
  plot_years = scenario_table$plot_year

  cov_step <- 1:(length(deco_years)-1)%>%
    purrr::map_dfr(function(i){
      start.y <- deco_years[i]
      end.y <- deco_years[i+1]
      
      ## MET
      cov_step_1.1 = {
        cov1 = cov%>%filter(Year == end.y)%>%
          dplyr::select(any_of(METv),Year,loc_id)%>%mutate(Year = start.y)
        
        cov2 = cov%>%filter(Year == start.y)%>%
          dplyr::select(any_of(ECOv),Year,loc_id)
        
        left_join(cov1, cov2)%>% mutate(Year = paste0(start.y,"_1.1"))
      }
      ## ECO
      cov_step_1.2 = {
        cov1 = cov%>%filter(Year == end.y)%>%
          dplyr::select(any_of(ECOv),Year,loc_id)%>%mutate(Year = start.y)
        
        cov2 = cov%>%filter(Year == start.y)%>%
          dplyr::select(any_of(METv),Year,loc_id)
        
        left_join(cov1, cov2)%>% mutate(Year = paste0(start.y,"_1.2"))
      }

      ## Join results
      list(cov_step_1.1,cov_step_1.2)%>%
        reduce(full_join)
    })%>%suppressMessages()

  Popf = Fuel_pop(cov_step, plot_years, info)%>%
    mutate(year = as.character(year))

  ## match on grids
  Gridsf <- distinct(Popfy,x,y,Area,Fuel,index)
  Popf_grid <- left_join(Gridsf,Popf)


  Decomp = 1:(length(deco_years)-1)%>%
    purrr::map_dfr(
      function(i){
        start.y = deco_years[i]
        end.y = deco_years[i+1]

        if(i != length(deco_years)-1){
          seq_order = c(paste0(start.y,c("_0","_1.1","_1.2")),paste0(end.y,"_0"))
        }else{
          seq_order = paste0(start.y,c("_0","_1.1","_1.2","_2"))
        }
        
        Popf_start = Popfy%>%filter(year == start.y)%>%
          dplyr::select(x,y,year,Area,Fuel,index,Popf)%>%
          mutate(year = seq_order[1])
        
        Popf_end = Popfy%>%filter(year == end.y)%>%
          dplyr::select(x,y,year,Area,Fuel,index,Popf)%>%
          mutate(year = seq_order[4])
        Popf_grid%>%filter(year%in%seq_order)
        
        
        list(Popf_start,Popf_grid,Popf_end)%>%
          reduce(full_join)%>%
          dplyr::rename(plot_year = year)%>%
          filter(plot_year%in%seq_order)%>%
          pivot_wider(names_from = plot_year, values_from = Popf)%>%
          mutate(year = deco_years[i],
                 MET_0 := !!rlang::sym(seq_order[2]) - !!rlang::sym(seq_order[1]),
                 MET_1 := !!rlang::sym(seq_order[4]) - !!rlang::sym(seq_order[3]),
                 ECO_0 := !!rlang::sym(seq_order[3]) - !!rlang::sym(seq_order[1]),
                 ECO_1 := !!rlang::sym(seq_order[4]) - !!rlang::sym(seq_order[2]),
                 Start := !!rlang::sym(seq_order[1]),
                 MET = (MET_0 + MET_1)/2,
                 ECO = (ECO_0 + ECO_1)/2,
                 End := !!rlang::sym(seq_order[4]),
                 .keep = 'unused'
          )
      })
  Decomp
}




Pop_Decomposition_mtcl <- function(deco_years, cov, info, sample_index) {
  ##  variable names
  METv = c("hd","cd","t2m")
  ECOv = c("GDP_pc_from_CSMAR","Urban_expenditure_pc","Rural_expenditure_pc")
  decoendy = deco_years[length(deco_years)-1]
  
  ##  create scenario index
  scenario_table = expand.grid(year = deco_years[1:(length(deco_years)-1)], 
                               serie = c(0,1.1,1.2)%>%as.character())%>%
    mutate(plot_year = paste0(year,"_",serie))%>%
    full_join(data.frame(year = decoendy, serie = as.character(2), plot_year = paste0(decoendy,"_2")))
  
  plot_years = scenario_table$plot_year
  
  cov_step <- 1:(length(deco_years)-1)%>%
    purrr::map_dfr(function(i){
      start.y <- deco_years[i]
      end.y <- deco_years[i+1]
      ## Start
      cov_step_0 = cov%>%filter(Year == start.y)%>%
        dplyr::select(any_of(ECOv),any_of(METv),Year,loc_id)%>%
        mutate(Year = paste0(start.y,"_0"))
      
      ## MET
      cov_step_1.1 = {
        cov1 = cov%>%filter(Year == end.y)%>%
          dplyr::select(any_of(METv),Year,loc_id)%>%mutate(Year = start.y)
        
        cov2 = cov%>%filter(Year == start.y)%>%
          dplyr::select(any_of(ECOv),Year,loc_id)
        
        left_join(cov1, cov2)%>% mutate(Year = paste0(start.y,"_1.1"))
      }
      ## ECO
      cov_step_1.2 = {
        cov1 = cov%>%filter(Year == end.y)%>%
          dplyr::select(any_of(ECOv),Year,loc_id)%>%mutate(Year = start.y)
        
        cov2 = cov%>%filter(Year == start.y)%>%
          dplyr::select(any_of(METv),Year,loc_id)
        
        left_join(cov1, cov2)%>% mutate(Year = paste0(start.y,"_1.2"))
      }
      ## End
      if(i == length(deco_years)-1){
        cov_step_2 = cov%>%filter(Year == end.y)%>%
          dplyr::select(any_of(ECOv),any_of(METv),Year,loc_id)%>% 
          mutate(Year = paste0(start.y,"_2"))
      }else{
        cov_step_2 = cov_step_0[F,]
      }
      ## Join results
      list(cov_step_0,cov_step_1.1,cov_step_1.2,cov_step_2)%>%
        reduce(full_join)
    })%>%suppressMessages()

  Popf = Fuel_pop_mtcl(cov_step, plot_years, sample_index, info)

  Decomp = 1:(length(deco_years)-1)%>%
    purrr::map_dfr(
      function(i){
        if(i != length(deco_years)-1){
          seq_order = c(paste0(deco_years[i],c("_0","_1.1","_1.2")),paste0(deco_years[i+1],"_0"))
        }else{
          seq_order = paste0(deco_years[i],c("_0","_1.1","_1.2","_2"))
        }
        
        Popf%>%
          dplyr::rename(plot_year = year)%>%
          filter(plot_year%in%seq_order)%>%
          pivot_wider(names_from = plot_year, values_from = Popf)%>%
          mutate(year = deco_years[i],
                 MET_0 := !!rlang::sym(seq_order[2]) - !!rlang::sym(seq_order[1]),
                 MET_1 := !!rlang::sym(seq_order[4]) - !!rlang::sym(seq_order[3]),
                 ECO_0 := !!rlang::sym(seq_order[3]) - !!rlang::sym(seq_order[1]),
                 ECO_1 := !!rlang::sym(seq_order[4]) - !!rlang::sym(seq_order[2]),
                 Start := !!rlang::sym(seq_order[1]),
                 MET = (MET_0 + MET_1)/2,
                 ECO = (ECO_0 + ECO_1)/2,
                 End := !!rlang::sym(seq_order[4]),
                 .keep = 'unused'
          )
      })
  Decomp
}


##======Decomposition plot======
mutate_name <- function(data, factor = TRUE){
  name_seq <- c("Fuel Use Fraction","Temperature Influence",
                "Economic Growth", "Population Change",
                "Population Aging", "Improved Health", "Concentration Change")
  result <- mutate(data,
         name = ifelse(Decomp == "PC", "Population Change",
                ifelse(Decomp == "ECO", "Economic Growth",
                ifelse(Decomp == "MET", "Temperature Influence",
                ifelse(Decomp == "FF", "Fuel Use Fraction",
                ifelse(Decomp == "PA", "Population Aging",
                ifelse(Decomp == "ORF", "Improved Health",
                ifelse(Decomp == "EXP", "Concentration Change",""))))))))
  if(factor == TRUE){
    result%>%mutate(name = factor(name, name_seq))
  }else{
    result
  }
}

Decomp_to_plot <- function(data, deco_years = c(2000,2007,2013,2020)) {
  penu = rev(deco_years)[2]
  last = rev(deco_years)[1]
  first = deco_years[1]

  Shape_sample <- data%>%
    group_by(year,Decomp)%>%
    dplyr::summarise(lower = quantile(value, 0.025)%>%round, 
              mean = mean(value)%>%round,
              upper = quantile(value, 0.975)%>%round)%>%ungroup()%>%
    dplyr::filter(!(year != penu & Decomp == "End"))%>%
    mutate(year = ifelse(year == penu & Decomp == "End", last, year),
           Decomp = factor(Decomp, levels = c("Start","MET","ECO","FF","PC","PA","ORF","EXP","End")))%>%
    arrange(year,Decomp)%>%
    mutate(lower_diff = round((lower - mean)), upper_diff = round((upper - mean)))%>%
    mutate_at(vars(mean), ~ ifelse(year !=first & Decomp== "Start", 0,
                                          ifelse(year ==last & Decomp == "End", 0, round(.x))))%>%suppressMessages()
  
  
  seq_len = nrow(Shape_sample)

  Shape_mean = Shape_sample%>%
    mutate(dock = cumsum(mean), dock_lower = dock + lower_diff, dock_upper = dock + upper_diff,
           basement = ifelse(year == first & Decomp == "Start",dock, dock - mean))%>%
    mutate(shape_top = ifelse(basement > dock, basement, dock),
           shape_bott = ifelse(basement < dock, basement, dock))%>%
    mutate_at(vars(shape_top,shape_bott,dock_lower,dock_upper,mean), ~ round(.x/1000))%>%
    mutate(seq = 1:n(), number_txt = ifelse(Decomp%in%c("Start","End"),"",str_c(mean)))%>%
    mutate_name(factor = F)%>%dplyr::rename(var_txt = name)

  Shape_top <- Shape_mean%>%
    dplyr::rename(value = shape_top)%>%
    mutate(var_txt1 = "", var_txt2 = "", 
           number_txt1 = ifelse(dock_upper<=value, number_txt, ""),
           number_txt2 = ifelse(dock_upper>value, number_txt, ""))%>%
    mutate(seq = factor(seq, levels = 1:seq_len))

  Shape_bott <- Shape_mean%>%
    dplyr::rename(value = shape_bott)%>%
    dplyr::filter(!Decomp %in% c("Start","End"))%>%
    mutate(Decomp = "dock", number_txt1 = "", number_txt = "",
           var_txt1 = ifelse(dock_lower>=value, var_txt, ""),
           var_txt2 = ifelse(dock_lower<value, var_txt, ""))%>%
    mutate(seq = factor(seq, levels = 1:seq_len))

  Shape_sum <- full_join(Shape_top, Shape_bott)%>%
    suppressMessages()
  
}

