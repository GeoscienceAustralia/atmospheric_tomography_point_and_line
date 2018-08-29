################################################################################

# Load required package

suppressPackageStartupMessages(library("coda"))

################################################################################

# Read in data from PyMC

traces <- read.csv("traces.csv", skip = 2, header = TRUE)
traces2 <- as.mcmc.list(lapply(traces, mcmc))
vals <- read.csv("traces.csv", header = FALSE)

vals <- vals[c(1, 2, 3, 4), ]
setup <- as.character(vals[1, 1])
file_name <- as.character(vals[1, 2])
true_val_known <- as.character(vals[1, 3])
true_Q <- as.numeric(vals[1, 4])

Q <- traces2[[1]]
tau <- traces2[[2]]

if (setup == "April 23 -- June 7"){
  true_Q <- 0.09667
} else if (setup == "June 8 -- June 12"){
  true_Q <- 0.08333
} else if (true_val_known == "No"){
  true_Q <- mean(Q)
}

################################################################################

# Generate plots

pdf(file = paste0(file_name,"-plots.pdf"))
par(mfrow = c(2, 3))
# Q
hist(Q, main = "Histogram of Q", freq = FALSE, xlab = "Q", col = "skyblue", xlim = c(min(traces$Q, true_Q), max(traces$Q, true_Q)))
if (setup == "April 23 -- June 7" || setup == "June 8 -- June 12" || (setup == "Other" & true_val_known == "Yes")){
    abline(v = true_Q, col = "red", lty = "dashed", lwd = 2)
}
plot(Q, density = FALSE, auto.layout = FALSE, col = "indianred1", main = "Trace of Q", ylab = "Q")
autocorr.plot(Q, auto.layout = FALSE, main = "Autocorrelation of Q")
# tau
hist(tau, main = "Histogram of tau", freq = FALSE, col = "skyblue", xlab = "tau")
plot(tau, density = FALSE, auto.layout = FALSE, col = "indianred1", main = "Trace of tau", ylab = "tau")
autocorr.plot(tau, auto.layout = FALSE, main = "Autocorrelation of tau")
dev.off()

