suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(lubridate))

##############################################################

### Read in functions needed for background estimation

source("upwind-downwind.R")
source("U-D-background.R")
source("30-min-background.R")

##############################################################

### Values you may wish to change

angle <- 10 # angle which determines what is upwind and what is downwind, near the plume
wind_speed_cutoff <- 0.5

##############################################################

### Read in files created from Python's interactive script

print("Processing data", quote = FALSE)

single_values <- read.csv("Rdata_single_vals.csv", sep = ',', header = TRUE)
outfile <- as.character(single_values[1, 1])
bg_type <- as.character(single_values[1, 2])
exp_setup <- as.character(single_values[1, 3])
no_of_insts <- as.numeric(single_values[1, 4])
datafile <- as.character(single_values[1, 5]) # where to get data from
tot_insts <- as.numeric(single_values[1, 6]) # total number of instruments possible to use
true_val_known <- as.character(single_values[1, 7])
true_Q <- as.numeric(single_values[1, 8])

if (no_of_insts == 1){
  inst_numbers <- as.numeric(single_values[1, 9])
} else {
  insts <- read.csv("Rdata_inst_nums.csv", sep = ',', skip = 1, header = FALSE)
  inst_numbers <- as.numeric(insts[1, ])
}

if (true_val_known == "No"){
    source_on_rate <- 0
} else {
    source_on_rate <- (true_Q * 60) - 0.1
} 


################################################################

### Read in data file and begin background estimation

inv_data <- read.csv(datafile, header = TRUE)
inv_data <- filter(inv_data, wind_speed > wind_speed_cutoff)

MyDate <- as.POSIXct(inv_data$Time, format = "%Y-%m-%d %H:%M")
Month <- month(MyDate)
Day <- day(MyDate)
Year <- year(MyDate) + 2000
Hour <- hour(MyDate)
Min <- minute(MyDate)
MyTime <- as.POSIXct(paste0(Year,"/",Month,"/",Day," ",Hour,":",Min))
inv_data$Time <- MyTime

# Filter out times when source is off if using 30-min-averaging BG method
if (bg_type == "30 min averaging"){
  inv_data <- filter(inv_data, release_rate > source_on_rate - 0.1)
}

# Filter out unwanted instruments
if (no_of_insts == 1){
  inv_data <- filter(inv_data, inst_no == inst_numbers)
} else {
  inst_remove <- c()
  for (i in 1:nrow(inv_data)){
    if (length(which(inst_numbers == inv_data$inst_no[i])) > 0){
      inst_remove[i] <- 0
    } else {
      inst_remove[i] <- i
    }
  }
  inst_remove <- inst_remove[inst_remove != 0]
  if (length(inst_remove) > 0){
    inv_data <- inv_data[-inst_remove, ]
  }
}

# Stop running if all observations have been removed, warn if not many observations left
if (nrow(inv_data) == 0){
    stop("No data values!")
    terminate
} else if (nrow(inv_data) <= 50){
    print(paste0("WARNING: low number of observations (", nrow(inv_data), ")"), quote = FALSE)
}

### Background

if (bg_type == "Upwind-downwind"){ # upwind-downwind background
  inv_data_temp <- inv_data
  inv_data <- Upwind.Downwind(inv_data = inv_data_temp)
  if (exp_setup == "April 23 -- June 7" | exp_setup == "June 8 -- June 12"){
      inv_data <- filter(inv_data, (Concentration < 2 & Position == "upwind") | (Concentration >= 2 & Position == "downwind"))
  }
  
  inv_data <- arrange(inv_data, Time)
  dim3 <- as.numeric(nrow(inv_data))
  
  # Get upwind count
  upwind_count <- length(inv_data$Position[inv_data$Position == "upwind"])
  # Add warning if upwind count low, stop running if upwind count is 0
  if (upwind_count < 1){
    stop("No upwind values!")
    terminate
  } else if (upwind_count <= 40){
    print(paste0("WARNING: Upwind count low (", upwind_count, ")"), quote = FALSE)
  }
  
  # Run the background estimation function
  inv_data_temp <- inv_data
  inv_data <- U.D.Background(inv_data = inv_data_temp)
  
  # remove all upwind values
  inv_data <- filter(inv_data, Position == "downwind")
  inv_data$bg_conc <- inv_data$Concentration - inv_data$background
  
  # Now remove times where source is off
  inv_data <- filter(inv_data, release_rate >= source_on_rate - 0.1)
  
  # Get downwind count now
  downwind_count <- length(inv_data$Position[inv_data$Position == "downwind"])
  # Add warning if downwind count low, stop running if downwind count is 0
  if (downwind_count == 0){
    stop("No downwind values!")
    terminate
  } else if (downwind_count <= 30){
    print(paste0("WARNING: Downwind count low (", downwind_count, ")"), quote = FALSE)
  }
} else { # 30 min average background
  # Run the background estimation function
  inv_data_temp <- inv_data
  inv_data <- thirty.min.Background(inv_data = inv_data_temp)
}

######################################################################

#### Create input file for Python script

# Add upwind and downwind counts if applicable
if (bg_type == "Upwind-downwind"){
  upwind_total <- rep(upwind_count, nrow(inv_data))
  downwind_total <- rep(downwind_count, nrow(inv_data))
} else {
  upwind_total <- rep(NA, nrow(inv_data))
  downwind_total <- rep(NA, nrow(inv_data))
}

# create empty data frame to be python input
python_data <- data.frame(inv_data[, c("air_temp", "air_pressure", 
                                       "wind_speed", "wind_dir", "L", "inst_no", 
                                       "source_x", "source_y", "z", "x1", 
                                       "y1", "x2", "y2", "bg_conc", "release_rate")], 
                                       upwind_total, downwind_total)


# write csv
write.csv(python_data, file = outfile, row.names = FALSE, na = "NaN")
