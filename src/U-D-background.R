U.D.Background <- function(inv_data){

    print("Starting background estimation", quote = FALSE)
  
    inv_data$background <- rep(NA, dim3) # Create empty vector to save BG values in
  
    for (i in 1:dim3){
        # print progress
        if (((i / dim3) <= 0.25) & (((i + 1) / dim3) > 0.25)){
            print("25% done", quote = FALSE)
        }
        if (((i / dim3) <= 0.5) & (((i + 1) / dim3) > 0.5)){
            print("50% done", quote = FALSE)
        }
        if (((i / dim3) <= 0.75) & (((i + 1) / dim3) > 0.75)){
            print("75% done", quote = FALSE)
        }
        # only find bg for downwind values as upwind vals give us no useful information really
        if (inv_data$Position[i] == "downwind"){
            inst <- inv_data$inst_no[i]
            # Set instrument number of the observation under consideration
            k <- 1
            # See if there are 5 or more upwind observations for the same instrument number as the observation under consideration within a 2-hour interval centred on the observation under consideration. If not, extend the interval to 4 hours, then 6, and so on, stopping at 10 hours. If 5 or more upwind observations found, take the mean of the concentrations corresponding to these values. If not, do nothing.
            while (k <= 5 & is.na(inv_data$background[i]) == TRUE){
                ATdata_sub <- filter(inv_data, Time[i] - hours(k) <= Time)
                ATdata_sub <- filter(ATdata_sub, Time <= inv_data$Time[i] + hours(k))
                ATdata_sub <- filter(ATdata_sub, Concentration != 0)
                if (length(ATdata_sub$Concentration[ATdata_sub$inst_no == inst & ATdata_sub$Position == "upwind"]) >= 5){
                    inv_data$background[i] <- mean(ATdata_sub$Concentration[ATdata_sub$inst_no == inst & ATdata_sub$Position == "upwind"])
                } else {
                    k <- k + 1
                }
            }
            j <- 1
            # Find all upwind values within a 2-hour interval centred on the observation under consideration. While the number of upwind observations found is less than 5, continue to expand the interval by 2 hours at a time. Once 5 or more are found, take the average of the concentrations of these values.
            if (is.na(inv_data$background[i]) == TRUE){
                ATdata_sub <- filter(inv_data, Time[i] - hours(j) <= Time)
                ATdata_sub <- filter(ATdata_sub, Time <= inv_data$Time[i] + hours(j))
                ATdata_sub <- filter(ATdata_sub, Concentration != 0)
                while (length(ATdata_sub$Concentration[ATdata_sub$Position == "upwind"]) < 5){
                    j <- j + 1
                    ATdata_sub <- filter(inv_data, Time[i] - hours(j) <= Time)
                    ATdata_sub <- filter(ATdata_sub, Time <= inv_data$Time[i] + hours(j))
                    ATdata_sub <- filter(ATdata_sub, Concentration != 0)
                }
                inv_data$background[i] <- mean(ATdata_sub$Concentration[ATdata_sub$Position == "upwind"])
            }
        }
    }
  
    print("Background estimation complete", quote = FALSE)
  
    return(inv_data)
}
