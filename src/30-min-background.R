thirty.min.Background <- function(inv_data){

  print("Starting background estimation", quote = FALSE)
  
  # Remove times where source is off
  
  mintime <- min(inv_data$Time) # Get earliest time in the data set
  maxtime <- max(inv_data$Time) # Get latest time in the data set
  
  ## Create an array of cut points (every thirty minutes)
  thirtymin_cuts <- seq(from = mintime, to = maxtime + minutes(30), by = 1800)
  head(thirtymin_cuts)
  
  # Add 30 min cuts
  inv_data$interval <- cut(inv_data$Time, thirtymin_cuts)
  bgconc <- inv_data %>% group_by(interval) %>%
    summarise(bg = mean(sort(Concentration[Concentration != 0])[1:5]))
  
  # join to inv_data
  inv_data_bg <- left_join(inv_data, bgconc, by = "interval")
  inv_data <- inv_data_bg
  
  # subtract off background
  inv_data$bg_conc <- inv_data$Concentration - inv_data$bg
  
  # remove NAs
  inv_data <- filter(inv_data, bg_conc != "NA")
  
  print("Background estimation complete", quote = FALSE)
  
  return(inv_data)
}  
