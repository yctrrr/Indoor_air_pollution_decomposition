##############################
## Household Air Pollution: ##
##   A Hierarchical Model   ##
##     Modified version     ##
##############################
library(nimble)
library(data.table)
library(ggplot2)
library(tidyverse)
library(coda)
library(gridExtra)
library(grid)
library(viridis)
library(readxl)
library(reshape2)
library(doParallel)
library(mgcv)
library(abind)
library(openxlsx)
library(ggfan)
library(scales)
library(RColorBrewer)
library(rgdal)

# Random number generator seeds (for each chain).
seed1=c(874326,383625,256846,383790,688352,132732,894322,848314) # For generating initial values.
seed2=c(642647,147376,368742,868372,582743,970287,185935,274763) # For MCMC. Markov chain Monte Carlo


inits_function=function(data,seed,N,C,K,Cov_num,n_fuels=3){
  
  set.seed(seed)
  
  inits=list(y=array(NA,dim=c(N,2,n_fuels)),
             phi=array(rgamma(C*2*n_fuels,2,0.05),dim=c(C,2,n_fuels)),
             beta_cov=array(rnorm((Cov_num+1)*2*n_fuels*C,0,1),dim=c(Cov_num+1,2,n_fuels,C)),
             probit_rho=rnorm(N,-4,2))
  
  for(i in 1:N){
    for(j in 1:2){
      for(f in 1:n_fuels){
        if(is.na(data$y[i,j,f])){
          inits$y[i,j,f]=0
        }
      }
    }
  }
  return(inits)
}

run_model=function(model_data, IV,
                   niter=2000, nburnin=1000,nchains=2,thin=1){

  nimbleOptions(oldConjugacyChecking = FALSE)
  nimbleOptions(useNewConfigureMCMC = TRUE)

  # Add the Beta-Binomial distribution to NIMBLE (required functions are in HAP_Functions.R).
  registerDistributions(list(dbetabin=list(
    BUGSdist='dbetabin(mu,phi,size,rho,weight)',discrete=TRUE)))

  # Centre splines over the data period.
  # jags.data$X has five knots
  N = dim(model_data$urban)[1] # Number of surveys.
  C = length(unique(model_data$urban$unique_country_index)) # Number of countries in the data.
  
  Cov_num = length(IV)

  Cov_df = dplyr::select(model_data[[1]],IV)%>%
    scale()

  center <- attributes(Cov_df)[[3]]
  scale <- attributes(Cov_df)[[4]]

  Cov_data = data.frame(intercept = rep(1,N))%>%
    cbind(Cov_df)
  
  Cov_df_bind = Cov_data%>%
    cbind(dplyr::select(model_data[[1]], year, country))

  # The model is written in the BUGS language.
  nimble_code=nimbleCode({
    ## i: number of surveys (total surveys equal to the number of data rows), j: urban, rural or overall, f:fuel types
    for(i in 1:N){
      # Urban proportion.
      # pi[i] <- ilogit(logit(p[i]))
      
      ## nu:relative means;  mu: proportion of corresponding fuels
      for(j in 1:2){ # Urban and rural usage models.
        ##  Way 1

        nu[i,j,1] <- ilogit(
          # sum(X[i,1:K]*beta[1:K,j,1,cindex[i]])+
          # sum(Cov_data[i,1:(Cov_num+1)]*beta_cov[1:(Cov_num+1),j,1])
          sum(Cov_data[i,1:(Cov_num+1)]*beta_cov[1:(Cov_num+1),j,1,cindex[i]])
          # sum(Cov_data[i,4:6]*beta_cov[4:6,j,1,cindex[i]])
        ) # Solid fuels
        
        nu[i,j,2] <- ilogit(
          # sum(X[i,1:K]*beta[1:K,j,2,cindex[i]])+
          # sum(Cov_data[i,1:(Cov_num+1)]*beta_cov[1:(Cov_num+1),j,2])
          sum(Cov_data[i,1:(Cov_num+1)]*beta_cov[1:(Cov_num+1),j,2, cindex[i]])
            # sum(Cov_data[i,7:9]*beta_cov[7:9,j,2, cindex[i]])
        )## Biomass
 
        
        nu[i,j,3] <- ilogit(
          # sum(X[i,1:K]*beta[1:K,j,3,cindex[i]])+
          # sum(Cov_data[i,1:(Cov_num+1)]*beta_cov[1:(Cov_num+1),j,3])
          sum(Cov_data[i,1:(Cov_num+1)]*beta_cov[1:(Cov_num+1),j,3,cindex[i]])
          # sum(Cov_data[i,7:9]*beta_cov[7:9,j,3,cindex[i]])
        )## Wood
        
        mu[i,j,1] <- nu[i,j,1] # Solid fuels
        # Disaggregate solid fuels.
        mu[i,j,2] <- nu[i,j,2]*mu[i,j,1] # Biomass.
        # Disaggregate biomass.
        mu[i,j,3] <- nu[i,j,3]*mu[i,j,2] # Wood.
      }
      
      for(j in 1:2){ # Data conditional likelihoods.
        y[i,j,1] ~ dbetabin(nu[i,j,1],phi[cindex[i],j,1],100000,rho[i],weight[i,j,1]) # Solid Fuels.
        # Disaggregate solid fuels.
        y[i,j,2] ~ dbetabin(nu[i,j,2],phi[cindex[i],j,2],y[i,j,1],rho[i],weight[i,j,2]) # Biomass.
        # Disaggregate biomass.
        y[i,j,3] ~ dbetabin(nu[i,j,3],phi[cindex[i],j,3],y[i,j,2],rho[i],weight[i,j,3]) # Wood.
        
      }

      # Outlier mixing parameters.
      rho[i] <- iprobit(probit_rho[i])
      probit_rho[i] ~ dnorm(-4,2)
    }
    
    for (c in 1:C) {
      for (j in 1:2) {
        for (f in 1:3) {
          beta_cov[1,j,f,c] ~  dnorm(0,4)
          # beta_cov[1,j,f] ~  dflat()
          for (v in c(1:Cov_num)) {
            beta_cov[(v+1),j,f,c] ~ dnorm(0, 4)  # Abundance slopes
          }
        }
      }
    }

    # Variance parameter model.
    for(c in 1:C){
      for(j in 1:2){
        for(f in 1:3){
          phi[c,j,f] ~ dgamma(2,0.05) # Dispersion parameter. 2 and 0.05 are global hyperparameters
          # phi[j,f] ~ dgamma(2,0.05) # Dispersion parameter. 2 and 0.05 are global hyperparameters
        }
      }
    }
  })## End of nimble code


  nimble_data=nimble_constants=nimble_inits=nimble_model=nimble_compiled_model=
    nimble_mcmc_config=nimble_mcmc=nimble_compiled_mcmc=list()

  # Set up nimble for each chain so we can run them in parallel.
  for(k in 1:nchains){
    # NIMBLE takes in constants (things like indices which don't relate
    # to the probabilistic model) separately from the data.

    nimble_constants[[k]]=list(N=N,C=C,Cov_num=Cov_num, cindex=model_data$urban$unique_country_index)

    # Data are things like the observed values and often covariates or structure matrices for splines.
    # blank_jagam: smooth function of time for 21 years (2000-2020 for examples)
    # X: survey data result

    nimble_data[[k]]=list(y=array(NA,dim=c(N,2,3)),
                          weight=array(1,dim=c(N,2,3)),
                          Cov_data=Cov_data)
  
    for(j in 1:2){
      nimble_data[[k]]$y[,j,]=as.matrix(dplyr::select(model_data[[j]],c('solid_count','biomass_count','wood_count')))
    }


    for(i in 1:N){ # Trick to prevent sampling of wholly missing survey areas (urban/rural/overall).
      for(j in 1:2){
        if(sum(is.na(nimble_data[[k]]$y[i,j,]))==5){
          nimble_data[[k]]$weight[i,j,]=0
          nimble_data[[k]]$y[i,j,]=0
        }
      }
    }

    # Initialise model parameters.
    nimble_inits[[k]]=inits_function(nimble_data[[k]],seed=seed1[k],N,C,K,Cov_num,3)
    # Construct the model object (which can be manipulated to obtain prior probabilities etc.)

    nimble_model[[k]]=nimbleModel(nimble_code,nimble_constants[[k]],nimble_data[[k]],nimble_inits[[k]])

    # Compile the model object.
    nimble_compiled_model[[k]]=compileNimble(nimble_model[[k]])

    # Configure the MCMC, telling it whether to look for conjugate (Gibbs) relationships
    # and what parameters to save.
    nimble_mcmc_config[[k]]=configureMCMC(nimble_model[[k]],useConjugacy = FALSE,
                                     # monitors = c('nu','pi','phi','rho','beta_cov','cov_gama'),
                                     monitors = c('nu','phi','rho','beta_cov'),
                                     # monitors = c('nu','pi','phi','rho','beta_cov','beta'),
                                     control=list(adaptInterval=1000,scale=0.1))
    
    # Build the MCMC object.
    nimble_mcmc[[k]]=buildMCMC(nimble_mcmc_config[[k]])

    # Compile the MCMC object.
    nimble_compiled_mcmc[[k]]=compileNimble(nimble_mcmc[[k]])
  }

  # Run the MCMC.
  cat("=======Run MCMC=======\n")
  # nimble_model[[1]]$getNodeNames()
  # nimble_model[[1]]$plotGraph()

  MCMC_list <- vector(mode = "list", length = nchains)

  for (k in 1:nchains) {
    MCMC_result <- runMCMC(nimble_compiled_mcmc[[k]], niter=niter, nburnin=nburnin,
            nchains = 1, thin=thin, inits=nimble_inits[[k]],
            samplesAsCodaMCMC = TRUE,setSeed = seed2[k])
    MCMC_list[[k]] = MCMC_result
  }

  samples = mcmc.list(MCMC_list)
  
  # Combine the chains into one matrix.
  { 
    combined_samples=do.call('rbind',samples)
    n_sim=dim(combined_samples)[1]
    
    output=list(mcmc=samples,processed_samples=list(),index=list())
    
    # To save memory we only save a selection of the model parameters for plotting etc.
    
    output$index$nu=which(dimnames(combined_samples)[[2]]=='nu[1, 1, 1]'):which(dimnames(combined_samples)[[2]]==paste('nu[',N,', 2, 3]',sep=''))
    output$processed_samples$nu=array(combined_samples[,output$index$nu],dim=c(n_sim,N,2,3))
    
    output$index$phi=which(dimnames(combined_samples)[[2]]=='phi[1, 1, 1]'):which(dimnames(combined_samples)[[2]]==paste('phi[',C,', 2, 3]',sep=''))
    output$processed_samples$phi=array(combined_samples[,output$index$phi],dim=c(n_sim,C,2,3))
    
    # output$index$phi=which(dimnames(combined_samples)[[2]]=='phi[1, 1]'):which(dimnames(combined_samples)[[2]]==paste('phi[3, 3]',sep=''))
    # output$processed_samples$phi=array(combined_samples[,output$index$phi],dim=c(n_sim,3,3))

    # output$index$kappa=which(dimnames(combined_samples)[[2]]=='kappa[1, 1]'):which(dimnames(combined_samples)[[2]]==paste('kappa[',K,', ',C,']',sep=''))
    # output$processed_samples$kappa=array(combined_samples[,output$index$kappa],dim=c(n_sim,K,C))
    
    # output$index$beta=which(dimnames(combined_samples)[[2]]=='beta[1, 1, 1, 1]'):which(dimnames(combined_samples)[[2]]==paste('beta[',K,', 2, 3, ',C,']',sep=''))
    # output$processed_samples$beta=array(combined_samples[,output$index$beta],dim=c(n_sim,K,2,3,C))
    

    output$index$rho=which(dimnames(combined_samples)[[2]]=='rho[1]'):which(dimnames(combined_samples)[[2]]==paste('rho[',N,']',sep=''))
    output$processed_samples$rho=array(combined_samples[,output$index$rho],dim=c(n_sim,N))
    
    # output$index$pi=which(dimnames(combined_samples)[[2]]=='pi[1]'):which(dimnames(combined_samples)[[2]]==paste('pi[',N,']',sep=''))
    # output$processed_samples$pi=array(combined_samples[,output$index$pi],dim=c(n_sim,N))
  
    output$index$beta_cov=which(dimnames(combined_samples)[[2]]=='beta_cov[1, 1, 1, 1]'):which(dimnames(combined_samples)[[2]]==paste('beta_cov[',Cov_num+1,', 2, 3, ',C,']',sep=''))
    output$processed_samples$beta_cov=array(combined_samples[,output$index$beta_cov],dim=c(n_sim,Cov_num+1,2,3,C))
    
    # output$replicates=v_replicates/100000
    # output$model_matrix=blank_jagam$jags.data$X
    output$Cov_df_bind =  Cov_df_bind
    output$center = center
    output$scale = scale
  }
  return(output)
}


## Output the (downscaled) fuel use prediction
#' @param plot_id unique location id
#' @param plot_index superior index in all_samples (province ind in orginal paper)
#' @param cov_data covariates corresponding to location id
#' @param all_samples outputs of Bayesian model
fuel_use_output <- function(plot_id, plot_index, cov_data,
                            plot_years = 2000:2020, all_samples,
                            type='underlying', output='samples',sampling_bias=TRUE){
  
  model_cov = all_samples$Cov_df_bind
  processed_samples = all_samples$processed_samples
  model_matrix = all_samples$model_matrix
  center = all_samples$center
  scale = all_samples$scale
  Cov_num = length(center)
  IV = names(center)
  
  cov_data_scaled <- dplyr::select(cov_data, IV)
  
  
  for (v in 1:Cov_num) {
    cov_data_scaled <- cov_data_scaled%>%
      mutate(!!paste0(IV[v]) := (UQ(rlang::sym(IV[v])) - center[IV[v]])/scale[IV[v]])
  }
  
  
  cov_data_scaled_join <- mutate(cov_data_scaled,Year = cov_data$Year,
                                 loc_id = cov_data$loc_id)
  # dim(all_samples$processed_samples$beta_cov)
  
  
  plot_phi=processed_samples$phi[,plot_index,1:2,]
  # plot_phi=processed_samples$phi[,,]
  Y = length(plot_years)
  
  ## urban population proportion
  n_sim = dim(processed_samples$beta)[1]
  plot_nu=plot_mu=plot_omega=array(NA,dim=c(n_sim,Y,2,3))
  plot_v=array(NA,dim=c(n_sim,Y,2,6))

  
  for(i in 1:Y){
    plot_cov_Y <- cov_data_scaled_join%>%
      filter(Year == plot_years[i], loc_id == plot_id)%>%
      dplyr::select(IV)%>%
      cbind(intercept = 1 ,.)
    
    for(j in 1:2){

      plot_nu[,i,j,1] <- expit(
        processed_samples$beta_cov[,1:(Cov_num+1),j,1,plot_index]%*%as.numeric(plot_cov_Y[,1:(Cov_num+1)])
      )# Solid fuels

      plot_nu[,i,j,2] <- expit(
        processed_samples$beta_cov[,1:(Cov_num+1),j,2,plot_index]%*%as.numeric(plot_cov_Y[,1:(Cov_num+1)])
      )# Biomass
      
      plot_nu[,i,j,3] <- expit(
        processed_samples$beta_cov[,1:(Cov_num+1),j,3,plot_index]%*%as.numeric(plot_cov_Y[,1:(Cov_num+1)])
      )# Wood
      
      
      plot_mu[,i,j,1] <- plot_nu[,i,j,1]
      # Disaggregate solid fuels.
      plot_mu[,i,j,2] <- plot_nu[,i,j,2]*plot_nu[,i,j,1] # Biomass.
      
      # Disaggregate biomass.
      plot_mu[,i,j,3] <- plot_nu[,i,j,3]*plot_mu[,i,j,2] # Wood.
     }
  }

  if(type=='underlying'){
    plot_u=array(NA,dim=c(n_sim,Y,2,6))
    plot_u[,,,1:6]=plot_mu
  }
  if(type=='survey'){
    for(i in 1:Y){
      # Truncate nu for numerical stability.
      plot_nu[,i,,]=structure(vapply(plot_nu[,i,,],function(x){min(max(x,expit(-6)),expit(6))}, numeric(1)),dim=dim(plot_nu[,i,,]))
      
      ##  alpha and beta follow beta-binomial, alpha = plot_nu, beta = 1-plot_nu

      plot_omega[,i,,] <- rbeta(n_sim*2*3,plot_nu[,i,,]*plot_phi,(1-plot_nu[,i,,])*plot_phi)
      
      for(j in 1:2){
        plot_v[,i,j,1] <- rbinom(n_sim,100000,plot_omega[,i,j,1])
  
        # Disaggregate solid fuels.
        plot_v[,i,j,2] <- rbinom(n_sim,plot_v[,i,j,1],plot_omega[,i,j,2]) # Biomass.
 
        plot_v[,i,j,3] <- rbinom(n_sim, plot_v[,i,j,2],plot_omega[,i,j,3]) # Wood.
      }
      plot_v[,,,4] <- 100000 - plot_v[,,,1]
      plot_v[,,,5] <- plot_v[,,,1]-plot_v[,,,2] # Coal.
      plot_v[,,,6] <- plot_v[,,,2]-plot_v[,,,3] # Crop waste
    }
    plot_u=array(NA,dim=c(n_sim,Y,2,6))
    plot_u[,,,1:6]=plot_v/100000 ## samples, year, area and fuel type
  }
  
  if(output=='samples'){
    return(plot_u)
  }
  if(output=='quantiles'){
    return(apply(plot_u,c(2,3,4),quantile,c(0.025,0.5,0.975)))
  }
}

