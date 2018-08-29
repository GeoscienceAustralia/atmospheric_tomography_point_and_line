import numpy as np
import util
import stats_calcs as sc
import params
import pymc as pm
import math
import rpy2.robjects as rp
import inquirer
import csv
import os
import sys

from itertools import zip_longest
from pymc import deterministic
from pymc.Matplot import plot

############################################################################
# Set MCMC values
############################################################################

runs = int(input("Enter the number of iterations for the MCMC simulation: "))
burnin = int(input("Enter the burn in for the MCMC simulation: "))
thin = int(input("Enter the thining variable for the MCMC simulation: "))

############################################################################
# Set constants
############################################################################

outfile = "python-data.csv"
plume_Q = 1
samples = 100

############################################################################
# Determine what data we will be working with.
############################################################################

# Which month?
setup_q = [
inquirer.List('setup', message = "Which setup?", choices = ['April 23 -- June 7', 'June 8 -- June 12', 'Other'],),
]
exp_setup = inquirer.prompt(setup_q)

if exp_setup["setup"] == "Other":
    true_val_known_q = [
    inquirer.List('true_val_known', message = "Is the true release rate known?", choices = ['Yes', 'No'],),
    ]
    true_val_k = inquirer.prompt(true_val_known_q)
    true_val_known = true_val_k["true_val_known"]
    if true_val_known == "Yes":
        true_Q = input('Enter the true release rate in grams per second: ')
    else:
        true_Q = "Unknown"
    data_file = input("Enter the name of the data file: ")
    tot_insts = int(input("What is the total number of instruments in the experiment? "))
    h_source = float(input("Enter the height of the gas source in metres: "))
    molar_mass = float(input("Enter the molar mass of the gas in grams: "))
elif exp_setup["setup"] == "April 23 -- June 7":
    data_file = "GA-inversion-data-1.csv"
    tot_insts = 25
    h_source = 0.3
    molar_mass = 16.04
    true_val_known = "Yes"
    true_Q = 0.09667
else:
    data_file = "GA-inversion-data-2.csv"
    tot_insts = 25
    h_source = 0.3
    molar_mass = 16.04
    true_val_known = "Yes"
    true_Q = 0.08333

no_of_insts = int(input('How many instruments are you considering? '))
if no_of_insts == 1:
    inst_numbers = input('Enter the instrument number you wish to consider: ')
else:
    insts = input('Enter the instrument numbers to consider: ')
    if insts == "all":
        inst_numbers = np.asarray(range(1, tot_insts + 1, 1))
    else:
        inst_numbers = np.asarray(list(map(int, insts.split(' '))))

# Which background method?
background_q = [
inquirer.List('background', message = "Which background method?", choices = ['30 min averaging', 'Upwind-downwind'],),
]
bg_type = inquirer.prompt(background_q)

############################################################################
# Write this data to csv for R input.
############################################################################

# if only 1 instrument to be used then put the instrument number into the single vals
if no_of_insts == 1:
	R_data_single_vals = np.array([outfile, bg_type["background"], exp_setup["setup"], no_of_insts, data_file, tot_insts, true_val_known, true_Q, inst_numbers])

	with open("Rdata_single_vals.csv", "w") as f:
		wr = csv.writer(f, delimiter = ',')
		wr.writerow(['outfile', 'bg_type', 'exp_month', 'no_of_insts', 'data_file', 'tot_insts', 'true_val_known', 'true_Q', 'inst_number'])
		wr.writerow(R_data_single_vals)

else:
	R_data_single_vals = np.array([outfile, bg_type["background"], exp_setup["setup"], no_of_insts, data_file, tot_insts, true_val_known, true_Q])

	with open("Rdata_single_vals.csv", "w") as f:
		wr = csv.writer(f, delimiter = ',')
		wr.writerow(['outfile', 'bg_type', 'setup', 'no_of_insts', 'data_file', 'tot_insts', 'true_val_known', 'true_Q'])
		wr.writerow(R_data_single_vals)
	# save instrument numbers separately if there are multiple instruments to be considered
	with open("Rdata_inst_nums.csv", "w") as f:
		wr = csv.writer(f, delimiter = ',')
		wr.writerow(['inst_numbers'])
		wr.writerow(inst_numbers)

############################################################################
# Put data in necessary format, calling R.
############################################################################

rp.r.source("background-estimation.R")

############################################################################
# Read off the background corrected data.
############################################################################

infile = outfile

measured = np.loadtxt(infile, delimiter=',', skiprows=1)
temp = measured[:, 0] # in Kelvin
pressure = measured[:, 1]  # in pascals.
w_speed = measured[:, 2] # in metres
wind_dir = measured[:, 3] # Theta in degrees East of North
plume_dir = ((wind_dir - 270) % 360) # Phi in degrees North of East
lvals = measured[:, 4] # stability class
source_x = measured[:, 6] # x-coordinate of source
source_y = measured[:, 7] # y-coordinate of source
z = measured[:, 8] # height of instrument tower
x1 = measured[:, 9] # location of mirror for lasers, location of tower for points
y1 = measured[:, 10] # location of mirror for lasers, location of tower for points
x2 = measured[:, 11] # location of reflector for lasers
y2 = measured[:, 12] # location of reflector for lasers
perturbation = measured[:, 13] # observed measurements
upwind_count = measured[0, 15] # count of upwind values
downwind_count = measured[0, 16] # count of downwind values

############################################################################
# Set up the model
############################################################################

tau = pm.Gamma('tau', alpha = 1.058, beta = 0.621)
Q = pm.HalfNormal('Q', tau = 1 / 2.25)

dims = len(perturbation)

conc = [0] * dims
for i in range(dims):
    if math.isnan(x2[i]) == True:
      rotcor = [0] * 2
      rotcor[:] = util.rotated(np.array([x1[i] - source_x[i], y1[i] - source_y[i]]), plume_dir[i])
      conc[i] = 0 if rotcor[0] < 0 else util.predicted2(rotcor[0], rotcor[1], z[i], plume_Q, lvals[i], w_speed[i], h_source, 'gaussian')
      conc[i] = util.gm_3toppmv(conc[i], temp[i], pressure[i], molar_mass)

    else:
      params = {'z': z[i], 'L': lvals[i], 'U': w_speed[i], 'H': h_source}
      conc[i] = util.line_average([source_x[i], source_y[i]], [x2[i], y2[i]], [x1[i], y1[i]], z[i], samples, plume_Q, h_source, plume_dir[i], temp[i], pressure[i], params, 'gaussian', molar_mass)

@deterministic
def predicted(Q = Q):
    pred = Q * conc
    return pred

res_sim = pm.Normal('residual', mu = predicted, tau = tau, value = perturbation)
residual = pm.Normal('residual', mu = predicted, tau = tau, value = perturbation, observed = True)

############################################################################
# Run the MCMC
############################################################################

S = pm.MCMC([Q, tau, predicted, res_sim, residual])

S.use_step_method(pm.Metropolis, Q)
S.use_step_method(pm.Metropolis, tau)
S.use_step_method(pm.Metropolis, res_sim)

S.sample(runs, burnin, thin)
trpred = S.trace('predicted')[:]
trres = S.trace('residual')[:]
trQ = S.trace('Q')[:]
trtau = S.trace('tau')[:]

tau_stats = sc.post_stats(trtau, "tau")
Q_stats = sc.post_stats(trQ, "Q")
pred_stats = sc.post_stats(trpred, "predicted")
res_stats = sc.post_stats(trres, "residual")
fd = [tau_stats, Q_stats, pred_stats, res_stats]

############################################################################
# Summarise MCMC results and save to csv
############################################################################

with open("summary.csv", 'w') as resultFile:
    wr = csv.writer(resultFile, delimiter = ',')
    wr.writerow(['parameter', 'mean', 'sd', 'var', 'lower', 'upper'])
    for i in fd:
        wr.writerow(i)
    wr.writerow([" ", " ", " ", " ", " ", " "])
    wr.writerow(['Upwind Count:', upwind_count, 'Downwind Count:', downwind_count])

############################################################################
# Save traces to send to R for plots
############################################################################

Trace = np.array([trQ, trtau])
if exp_setup["setup"] == "Other" and true_val_known == "Yes":
    with open("traces.csv", 'w') as Tr:
        wr = csv.writer(Tr, delimiter = ',')
        wr.writerow([exp_setup["setup"], sys.argv[1], true_val_known, true_Q])
        wr.writerow([" ", " ", " ", " "])
        wr.writerow(['Q', 'tau'])
        for i in np.transpose(Trace):
            wr.writerow(i)
else:
    with open("traces.csv", 'w') as Tr:
        wr = csv.writer(Tr, delimiter = ',')
        wr.writerow([exp_setup["setup"], sys.argv[1], "none", "none"])
        wr.writerow([" ", " "])
        wr.writerow(['Q', 'tau'])
        for i in np.transpose(Trace):
            wr.writerow(i)

rp.r.source("trace-histograms.R") # Generate and save plots

############################################################################
# Finally clean up and remove csv files which are not needed
############################################################################

fname = sys.argv[1]
sname = fname + '-summary.csv'
os.rename("summary.csv", sname)

os.remove("Rdata_single_vals.csv")
if no_of_insts > 1:
   os.remove("Rdata_inst_nums.csv")
os.remove("python-data.csv") # remove this line if you want to keep the filtered data set which the inversion is run on
os.remove("traces.csv")


