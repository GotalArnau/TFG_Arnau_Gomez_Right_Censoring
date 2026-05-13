################################################################################
#-------------------------------------------------------------------------------
# Function to calculate the power of the tests using the .RData files generated with script 07.
#-------------------------------------------------------------------------------
################################################################################

Power_calc <- function(folder, alpha = 0.1){
  
  is_rejected <- function(p_value, alpha) {
    return(as.numeric(p_value < alpha))
  }
  
  files <- list.files(folder, pattern = "\\.RData", full.names = TRUE)
  
  power_CvM <- numeric()
  power_KS <- numeric()
  power_AD <- numeric()
  power_KS_noboot <- numeric()
  km_min_median <- numeric()
  SEC_min_median <- numeric()
  n <- numeric()
  cens_aux <- numeric()
  
  df2 <- data.frame(
    power = numeric(),
    km_min_median = numeric(),
    SEC_min_median = numeric(),
    test = character(),
    cens = character(),
    n = numeric(),
    min_mean = numeric(),
    min_median = numeric()
  )
  
  
  for(file in files){
    
    load(file)
    
    power_CvM <- c(power_CvM, sum(is_rejected(results$pvalsCvM, alpha),
                                  na.rm = TRUE) / length(na.omit(results$pvalsCvM)))
    power_KS <- c(power_KS, sum(is_rejected(results$pvalsKS, alpha),
                                na.rm = TRUE) / length(na.omit(results$pvalsKS)))
    power_AD <- c(power_AD, sum(is_rejected(results$pvalsAD, alpha),
                                na.rm = TRUE) / length(na.omit(results$pvalsAD)))
    power_KS_noboot <- c(power_KS_noboot, sum(is_rejected(results$pvalsKS_noboot, alpha),
                                              na.rm = TRUE) / length(na.omit(results$pvalsKS_noboot)))
    km_min_median <- c(km_min_median, median(results$km_min))
    
    SEC_min_median <- c(SEC_min_median, median(results$SEC_min))
    
    split1 <- strsplit(strsplit(file, "/")[[1]][4],"_")[[1]]
    n <- c(n, as.numeric(split1[7]))
    cens_aux <- c(cens_aux, paste0((1 - as.numeric(split1[5]))*100, "% Censura"))
  }
  
  if(split1[3] == "lognormal"){
    sim <- rep(c("Symmetric", "Asymmetric", "Extremly asymmetric"))
  }else{
    sim <- rep(c( "Asymmetric", "Extremly asymmetric", "Symmetric"))
  }
  
  n_t <- length(km_min_median)
  df2 <- data.frame(
    power = as.numeric(c(power_AD, power_CvM, power_KS, power_KS_noboot)),
    km_min_median = as.numeric(rep(km_min_median, 4)),
    SEC_min_median = as.numeric(rep(SEC_min_median, 4)),
    cens = rep(cens_aux, 4),
    test = rep(c("Anderson-Darling", "Cramér von Mises", "Kolmogorov-Smirnov", "KS No Bootstrap"), each = n_t),
    symmetry = rep(sim, 9), 
    n = n
  )
  
  directorio <- paste0(folder, "/Power")
  
  if (!dir.exists(paste0(folder, "/Power"))){
    dir.create(paste0(folder, "/Power"), recursive = TRUE)
  }
  
  output_file <- paste0(folder, "/Power", "/_Distribution_", split1[3], "_resultados_h0_", split1[11], "_all_.RData")
  save(df2, file = output_file)
  
}



#-----------------------------------------------------------------------------------------------------------------------------------------------
#Random Censoring
#-----------------------------------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------------------------------
#X ~ Weibull
#-----------------------------------------------------------------------------------------------------------------------------------------------
Power_calc("02_Data/Weibull_RandomCensoring/df_power_all4_logistic")

Power_calc("02_Data/Weibull_RandomCensoring/df_power_all4_lognormal")

Power_calc("02_Data/Weibull_RandomCensoring/df_power_all4_weibull")


#-----------------------------------------------------------------------------------------------------------------------------------------------
#X ~ Lognormal
#-----------------------------------------------------------------------------------------------------------------------------------------------

graph_comp("02_Data/Lognormal_RandomCensoring/df_power_all4_logistic")

graph_comp("02_Data/Lognormal_RandomCensoring/df_power_all4_weibull")

graph_comp("02_Data/Lognormal_RandomCensoring/df_power_all4_lognormal")


#-----------------------------------------------------------------------------------------------------------------------------------------------
#Administrative Censoring
#-----------------------------------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------------------------------
#X ~ Weibull
#-----------------------------------------------------------------------------------------------------------------------------------------------
graph_comp("02_Data/Weibull_AdminCensoring/df_power_all4_logistic")

graph_comp("02_Data/Weibull_AdminCensoring/df_power_all4_lognormal")

graph_comp("02_Data/Weibull_AdminCensoring/df_power_all4_weibull")


#-----------------------------------------------------------------------------------------------------------------------------------------------
#X ~ Lognormal
#-----------------------------------------------------------------------------------------------------------------------------------------------

graph_comp("02_Data/Weibull_AdminCensoring/df_power_all4_logistic")

graph_comp("02_Data/Weibull_AdminCensoring/df_power_all4_weibull")

graph_comp("02_Data/Weibull_AdminCensoring/df_power_all4_lognormal")

