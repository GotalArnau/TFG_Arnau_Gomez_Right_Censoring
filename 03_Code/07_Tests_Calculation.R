
calc_power_plus <- function(folder, alpha = 0.1, h0, n){
  library(GofCens)
  library(epitools)
  library(gridExtra)
  library(fitdistrplus)
  
  results <- data.frame(
    pvalsCvM = numeric(),
    pvalsKS = numeric(),
    pvalsKS_noboot = numeric(),
    pvalsAD = numeric(),
    t_max = numeric(),
    km_min = numeric(),
    SEC_min = numeric()
  )
  
  
  files <- list.files(folder, pattern = "\\.RData$", full.names = TRUE)
  
  for(file in files){ #---------------------------------------OJO--------------------------------------
    
    load(file)
    
    pvalsCvM <- numeric(n)
    pvalsKS <- numeric(n)
    pvalsAD <- numeric(n)
    pvalsKS_noboot <- numeric(n)
    km_min <- numeric(n)
    t_max <- numeric(n)
    SEC_min <- numeric(n)
    km <- NULL
    all3 <- NULL
    SEC <- NULL
    
    for (i in 1:n) {
      print(i)
      km <- NULL
      SEC <- NULL
      dados <- as.data.frame(data[i])
      times <- dados$times
      cens <- dados$cens
      
      # Calculate p-values for each test
      all3 <- try(gofcens(times, cens, distr = h0))
      
      pvalsCvM[i] <- try(all3$pval[2])
      pvalsKS[i] <- try(all3$pval[1])
      pvalsAD[i] <- try(all3$pval[3])
      SEC <- try(KScens(times, cens, distr = h0, boot = FALSE))
      pvalsKS_noboot[i] <- SEC$Test[[2]]
      km <- survfit(Surv(times, cens) ~ 1)
      km_min[i] <- tail(km$surv, 1)
      t_max[i] <- tail(km$time, 1)
      
      SEC_min[i] <- 1 - switch(h0,
                               weibull = pweibull(q = t_max[i], shape = SEC$Estimates[[1]], scale = SEC$Estimates[[2]]),
                               lognormal = plnorm(q = t_max[i], meanlog = SEC$Estimates[[1]], sdlog = SEC$Estimates[[2]]),
                               logistic = plogis(q = t_max[i], location = SEC$Estimates[[1]], scale = SEC$Estimates[[2]]))
    }
    
    results <- data.frame(
      pvalsCvM = pvalsCvM,
      pvalsKS =  pvalsKS,
      pvalsAD = pvalsAD,
      pvalsKS_noboot = pvalsKS_noboot,
      t_max = t_max,
      km_min = km_min,
      SEC_min = SEC_min
    )
    
    output_dir <- paste0(folder,"/df_power_all4_",h0)
    
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    output_file <- paste0(folder,"/df_power_all4_", h0, "/", sub("\\.RData$", "", basename(file)), "h0_", h0, "_Stmax_", ".RData" )
    
    save(results, file = output_file)
  }
  
  
}


#-----------------------------------------------------------------------------------------------------------------------------------------------
#Random Censoring
#-----------------------------------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------------------------------
#X ~ Weibull
#-----------------------------------------------------------------------------------------------------------------------------------------------

calc_power_plus("C:/Users/arnau.gomez/Desktop/GofCensSimulatios-Study/02_Data/Weibull_RandomCensoring", 
                alpha = 0.1, h0 = "lognormal", n = 1000)

calc_power_plus("C:/Users/arnau.gomez/Desktop/GofCensSimulatios-Study/02_Data/Weibull_RandomCensoring", 
                alpha = 0.1, h0 = "logistic", n = 1000)

calc_power_plus("C:/Users/arnau.gomez/Desktop/GofCensSimulatios-Study/02_Data/Weibull_RandomCensoring", 
                alpha = 0.1, h0 = "weibull", n = 1000)



#-----------------------------------------------------------------------------------------------------------------------------------------------
#X ~ Lognormal
#-----------------------------------------------------------------------------------------------------------------------------------------------

calc_power_plus("C:/Users/arnau.gomez/Desktop/GofCensSimulatios-Study/02_Data/Lognormal_RandomCensoring", 
                alpha = 0.1, h0 = "weibull", n = 1000)

calc_power_plus("C:/Users/arnau.gomez/Desktop/GofCensSimulatios-Study/02_Data/Lognormal_RandomCensoring",
                alpha = 0.1, h0 = "lognormal", n = 1000)

calc_power_plus("C:/Users/arnau.gomez/Desktop/GofCensSimulatios-Study/02_Data/Lognormal_RandomCensoring", 
                alpha = 0.1, h0 = "logistic", n = 1000)







#-----------------------------------------------------------------------------------------------------------------------------------------------
#Administrative Censoring
#-----------------------------------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------------------------------
#X ~ Weibull
#-----------------------------------------------------------------------------------------------------------------------------------------------

calc_power_plus("C:/Users/arnau.gomez/Desktop/GofCensSimulatios-Study/02_Data/Weibull_AdminCensoring", 
                alpha = 0.1, h0 = "lognormal", n = 1000)

calc_power_plus("C:/Users/arnau.gomez/Desktop/GofCensSimulatios-Study/02_Data/Weibull_AdminCensoring", 
                alpha = 0.1, h0 = "logistic", n = 1000)

calc_power_plus("C:/Users/arnau.gomez/Desktop/GofCensSimulatios-Study/02_Data/Weibull_AdminCensoring", 
                alpha = 0.1, h0 = "weibull", n = 1000)



#-----------------------------------------------------------------------------------------------------------------------------------------------
#X ~ Lognormal
#-----------------------------------------------------------------------------------------------------------------------------------------------

calc_power_plus("C:/Users/arnau.gomez/Desktop/GofCensSimulatios-Study/02_Data/Lognormal_AdminCensoring", 
                alpha = 0.1, h0 = "weibull", n = 1000)

calc_power_plus("C:/Users/arnau.gomez/Desktop/GofCensSimulatios-Study/02_Data/Lognormal_AdminCensoring", 
                alpha = 0.1, h0 = "logistic", n = 1000)

calc_power_plus("C:/Users/arnau.gomez/Desktop/GofCensSimulatios-Study/02_Data/Lognormal_AdminCensoring",
                alpha = 0.1, h0 = "lognormal", n = 1000)
