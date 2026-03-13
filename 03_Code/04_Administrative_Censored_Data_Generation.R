#Function to generate the data to be used in the simulations

set.seed(28657)
gendata <- function(distr, shape, scale, meanlog, sdlog, cens, n){
  if (distr== "weibull"){
    
    dat <- rweibull(n, shape, scale)
    
  }else if(distr == "lognormal"){
    dat <- rlnorm(n, meanlog, sdlog)
  }
  
  #generate the cesoring 
  quant <- quantile(dat,cens)
  cens <- c()
  cens <- ifelse(dat < quant, 1, 0)
  dat[dat > quant] <- quant
  
  
  return(data.frame(times = dat, cens))
  
}

# set the parameters to create the data 
# three level of censoring 
# three sample sizes
# if weibul use hape and scale, if lognorm use meanlog and sdlog

set.seed(28657)
sets <- list(distr = c("weibull", "lognormal"),
             shape = c(3.5, 1.7, 1),
             scale = 20,
             meanlog = 2,
             sdlog = c(0.01, 0.4, 1.5),
             cens = c(1,0.4,0.7),
             n = c(75,150,300))

leng <- sapply(sets, length)
for (a in 1:leng[1]) {
  distr <- sets$distr[a]
  
  if (distr == "weibull") {
    # Iterate over relevant parameters for Weibull
    for (b in 1:leng[2]) {
      for (c in 1:leng[3]) {
        for (f in 1:leng[6]) {
          for (g in 1:leng[7]) {
            data <- vector("list", 1000) #generate 100 Data sets for each condition
            for (i in 1:1000) {
              data[[i]] <- gendata(
                distr = distr,
                shape = sets$shape[b],
                scale = sets$scale[c],
                cens = sets$cens[f],
                n = sets$n[g]
              )
            }
            comment(data) <- paste(
              "Distribution", distr,
              "Shape", sets$shape[b],
              "Scale", sets$scale[c],
              "Censoring", sets$cens[f],
              "n", sets$n[g]
            )
            file_name <- paste(
              "C:/Users/arnau.gomez/Desktop/GofCensSimulatios-Study/02_Data/Weibull_AdminCensoring/",
              "Distribution", distr,
              "Censoring", sets$cens[f],
              "n", sets$n[g],
              "Shape", sets$shape[b],
              ".RData", sep = '_'
            )
            save(data, file = file_name)
            rm(data, i)
          }
        }
      }
    }
  } else {
    # Iterate over relevant parameters for Lognormal
    for (d in 1:leng[4]) {
      for (e in 1:leng[5]) {
        for (f in 1:leng[6]) {
          for (g in 1:leng[7]) {
            data <- vector("list", 1000)
            for (i in 1:1000) {
              data[[i]] <- gendata(
                distr = distr,
                meanlog = sets$meanlog[d],
                sdlog = sets$sdlog[e],
                cens = sets$cens[f],
                n = sets$n[g]
              )
            }
            comment(data) <- paste(
              "Distribution", distr,
              "Meanlog", sets$meanlog[d],
              "Sdlog", sets$sdlog[e],
              "Censoring", sets$cens[f],
              "n", sets$n[g]
            )
            file_name <- paste(
              "C:/Users/arnau.gomez/Desktop/GofCensSimulatios-Study/02_Data/Lognormal_AdminCensoring/",
              "Distribution", distr,
              "Censoring", sets$cens[f],
              "n", sets$n[g],
              "Sdlog", sets$sdlog[e],
              ".RData", sep = '_'
            )
            save(data, file = file_name)
            rm(data, i)
            
          }
        }
      }
    }
  }
}
