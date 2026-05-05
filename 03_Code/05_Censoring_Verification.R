#################################################################################
#--------------------------------------------------------------------------------
# We compute the censoring on avg for each scenario in the Random censoring case
#  - We have 1000 data sets per scenario
#  - Data sets have been generated using random censoring generating the
#    specific censoring we want. Nonetheless, this censoring is not exact
#    so we should compute the censoring on avg per scenario
#--------------------------------------------------------------------------------
#################################################################################

calculate_censoring_avg <- function(folder) {
  
  # Initialize an empty data frame to store results
  results_df <- data.frame(
    File = character(),
    AvgCens = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Get list of .RData files in the specified folder
  files <- list.files(folder, pattern = "\\.RData$", full.names = TRUE)
  
  for (file in files) {
    # Load data from the file
    load(file)
    #number of datasets to be used per scenario!!!
    n <- 1000
    censor <- numeric(n)
    
    #obs: we are using just n/n0 datasets we have in the .RData
    for (i in 1:n) {
      print(i)
      dados <- data[i]
      dados <- as.data.frame(dados)
      censor[i] <- sum(dados$cens==0)/length(dados$cens)
    }
    
    # Add the result to the data frame
    results_df <- rbind(results_df, data.frame(
      File = sub("\\.RData$", "", basename(file)),  # Remove .RData from filename
      AvgCens = mean(censor)
    ))
    #we print the results just to check
    print(results_df)
  }
  
  # Define output file path for CSV
  output_file <- paste0(folder, "/AVG_censoring_13march2026.csv")
  
  # Write the results data frame to a CSV file
  write.csv2(results_df, file = output_file, row.names = FALSE)
}

###weibull data
calculate_censoring_avg("C:/Users/arnau.gomez/Desktop/GofCensSimulatios-Study/02_Data/Weibull_AdminCensoring")

calculate_censoring_avg("C:/Users/arnau.gomez/Desktop/GofCensSimulatios-Study/02_Data/Weibull_RandomCensoring")

###lognormal data
calculate_censoring_avg("C:/Users/arnau.gomez/Desktop/GofCensSimulatios-Study/02_Data/Lognormal_AdminCensoring")

calculate_censoring_avg("C:/Users/arnau.gomez/Desktop/GofCensSimulatios-Study/02_Data/Lognormal_RandomCensoring")

