# Atmospheric Tomography point and line

The Atmospheric Tomography software is a command line tool written in R and Python to estimate the emission rate of a point source from concentration data. It implements a Bayesian inversion method.

A user manual is provided for the software, titled "AT_2018_user_manual", and contains detailed instructions on installing all of the components required for the software to run.

All of the source code is located within the "src" folder. A collection of simulated data has also been provided, to allow the user to test the software. This file is also located within the "src" folder, and the code needed to run the software on this data set is given in the "Quick Start" section below, as well as in Appendix A of the software user manual. 

# Quick Start

To run the software on the simulated data set provided, enter the following inputs. The number of iterations, burnin, thinning variable, number of instruments you are considering (no more than three for this example), instrument numbers (enter each number with a single space between them if you don't wish to consider all three instruments), and background estimation method can be changed if desired, however the other values must stay the same. 

```
will-scarlett:src Laura$ python3 atmospheric-tomography.py sim-data
	
	Enter the number of iterations for the MCMC simulation: 100000
	
	Enter the burn in for the MCMC simulation: 20000
	
	Enter the thining variable for the MCMC simulation: 10
	
	[?] Which setup?: Other
	April 23 -- June 7
	June 8 -- June 12
	> Other
	
	[?] Is the true release rate known?: Yes
	> Yes
	No
	
	Enter the true release rate in grams per second: 0.07
	
	Enter the name of the data file: simulated_data.csv
	
	What is the total number of instruments in the experiment? 3
	
	Enter the height of the gas source in metres: 0.35
	
	Enter the molar mass of the gas in grams: 16.04
	
	How many instruments are you considering? 3
	
	Enter the instrument numbers to consider: all
	
	[?] Which background method?: Upwind-downwind
	30 min averaging
	> Upwind-downwind
```

