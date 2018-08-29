import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from numpy import linalg as LA
from scipy.constants import R
import plume_models as pmods
import params


def stability_class(L):
    if L < -1e5 or L > 1e5:
        stability = 'C'
    elif -1e5 < L < -1e2:
        stability = 'B'
    elif -1e2 < L < 0:
        stability = 'A'
    elif 0 < L < 1e2:
        stability = 'E'
    elif 1e2 < L < 1e5:
        stability = 'D'
    return stability


def plume_index_gaussian(z, L):
    zoverl = z / L
    if zoverl > -350 and zoverl < -1.75:
        index = '1'
    elif zoverl > -1.75 and zoverl < -1.25:
        index = '2'
    elif zoverl > -1.25 and zoverl < -0.9:
        index = '3'
    elif zoverl > -0.9 and zoverl < -0.7:
        index = '4'
    elif zoverl > -0.7 and zoverl < -0.5:
        index = '5'
    elif zoverl > -0.5 and zoverl < -0.35:
        index = '6'
    elif zoverl > -0.35 and zoverl < -0.25:
        index = '7'
    elif zoverl > -0.25 and zoverl < -0.175:
        index = '8'
    elif zoverl > -0.175 and zoverl < -0.135:
        index = '9'
    elif zoverl > -0.135 and zoverl < -0.105:
        index = '10'
    elif zoverl > -0.105 and zoverl < -0.075:
        index = '11'
    elif zoverl > -0.075 and zoverl < -0.023:
        index = '12'
    elif zoverl > -0.023 and zoverl < -0.00375:
        index = '13'
    elif zoverl > -0.00375 and zoverl < -2.3e-5:
        index = '14'
    elif zoverl > -2.3e-5 and zoverl < 2.3e-5:
        index = '15'
    elif zoverl > 2.3e-5 and zoverl < 0.00375:
        index = '16'
    elif zoverl > 0.00375 and zoverl < 0.01125:
        index = '17'
    elif zoverl > 0.01125 and zoverl < 0.0225:
        index = '18'
    elif zoverl > 0.0225 and zoverl < 0.045:
        index = '19'
    elif zoverl > 0.045 and zoverl < 0.046:
        index = '20a'
    elif zoverl > 0.046 and zoverl < 0.05:
        index = '20b'
    elif zoverl > 0.05 and zoverl < 0.054:
        index = '20c'
    elif zoverl > 0.054 and zoverl < 0.057:
        index = '20d'
    elif zoverl > 0.057 and zoverl < 0.06:
        index = '20e'
    elif zoverl > 0.06 and zoverl < 0.063:
        index = '20f'
    elif zoverl > 0.063 and zoverl < 0.067:
        index = '20g'
    elif zoverl > 0.067 and zoverl < 0.07:
        index = '20h'
    elif zoverl > 0.07 and zoverl < 0.075:
        index = '20i'
    elif zoverl > 0.075 and zoverl < 0.105:
        index = '21'
    elif zoverl > 0.105 and zoverl < 0.135:
        index = '22'
    elif zoverl > 0.135 and zoverl < 0.175:
        index = '23'
    elif zoverl > 0.175 and zoverl < 0.25:
        index = '24'
    elif zoverl > 0.25 and zoverl < 0.35:
        index = '25'
    elif zoverl > 0.35 and zoverl < 0.5:
        index = '26'
    elif zoverl > 0.5 and zoverl < 0.7:
        index = '27'
    elif zoverl > 0.7 and zoverl < 0.9:
        index = '28'
    elif zoverl > 0.9 and zoverl < 1.25:
        index = '29'
    elif zoverl > 1.25 and zoverl < 1.75:
        index = '30'
    else:  # zoverl > 1.75 and zoverl < 55:
        index = '31'
    # print('Plume index = ',index)
    return index


#################
def plume_index(z, L):
    zoverl = z / L
    if zoverl > -350 and zoverl < -1.75:
        index = '1'
    elif zoverl > -1.75 and zoverl < -1.25:
        index = '2'
    elif zoverl > -1.25 and zoverl < -0.9:
        index = '3'
    elif zoverl > -0.9 and zoverl < -0.7:
        index = '4'
    elif zoverl > -0.7 and zoverl < -0.5:
        index = '5'
    elif zoverl > -0.5 and zoverl < -0.35:
        index = '6'
    elif zoverl > -0.35 and zoverl < -0.25:
        index = '7'
    elif zoverl > -0.25 and zoverl < -0.175:
        index = '8'
    elif zoverl > -0.175 and zoverl < -0.135:
        index = '9'
    elif zoverl > -0.135 and zoverl < -0.105:
        index = '10'
    elif zoverl > -0.105 and zoverl < -0.075:
        index = '11'
    elif zoverl > -0.075 and zoverl < -0.045:
        index = '12'
    elif zoverl > -0.045 and zoverl < -0.0225:
        index = '13'
    elif zoverl > -0.0225 and zoverl < -0.01125:
        index = '14'
    elif zoverl > -0.01125 and zoverl < -0.00375:
        index = '15'
    elif zoverl > -0.00375 and zoverl < 0.00375:
        index = '16'
    elif zoverl > 0.00375 and zoverl < 0.01125:
        index = '17'
    elif zoverl > 0.01125 and zoverl < 0.0225:
        index = '18'
    elif zoverl > 0.0225 and zoverl < 0.045:
        index = '19'
    elif zoverl > 0.045 and zoverl < 0.075:
        index = '20'
    elif zoverl > 0.075 and zoverl < 0.105:
        index = '21'
    elif zoverl > 0.105 and zoverl < 0.135:
        index = '22'
    elif zoverl > 0.135 and zoverl < 0.175:
        index = '23'
    elif zoverl > 0.175 and zoverl < 0.25:
        index = '24'
    elif zoverl > 0.25 and zoverl < 0.35:
        index = '25'
    elif zoverl > 0.35 and zoverl < 0.5:
        index = '26'
    elif zoverl > 0.5 and zoverl < 0.7:
        index = '27'
    elif zoverl > 0.7 and zoverl < 0.9:
        index = '28'
    elif zoverl > 0.9 and zoverl < 1.25:
        index = '29'
    elif zoverl > 1.25 and zoverl < 1.75:
        index = '30'
    else:  # zoverl > 1.75 and zoverl < 55:
        index = '31'
    # print('Plume index = ',index)
    return index


##############################################
# concentartion is in g/m^3
# Pressure is in Pa
# T is in Kelvin
# molar_mass is in g/mol
###########################################################################################################
def gm_3toppmv(concentration, T, P, molar_mass):
    num_moles = concentration / molar_mass
    vol = 1e6 * num_moles * R * T / P  # To get PPM.
    return vol

###########################################################################################################
def powerlaw(x, scale, exponent):
    return scale * (x ** exponent)


###########################################################################################################
def line(t, p0, p1):
    l00, l01 = p0[0], p0[1]
    l10, l11 = p1[0], p1[1]
    return [l00 * (1 - t) + l10 * t, l01 * (1 - t) + l11 * t]


###########################################################################################################
def rotated(vec, theta):
    theta = np.deg2rad(theta)
    rot_mat = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    return rot_mat.dot(vec)


###########################################################################################################
def predicted(x, y, z, Q, dict, method):
    if method == 'gaussian':
        stability = stability_class(dict['L'])
        pars = params.gaussian[stability]
        sigmay = powerlaw(x, pars['c'], pars['d'])
        sigmaz = powerlaw(x, pars['a'], pars['b'])
        return pmods.gaussian_plume(x, y, z, Q, dict['U'], dict['H'], sigmay, sigmaz)
    elif method == 'semi-gaussian':
        index = plume_index(dict['z'], dict['L'])
        pars = params.semigaussian[index]
        return pmods.semigaussian(y, x, Q, pars)  # Apparently this works, do not change!
    elif method == 'gauss_poly':
        index = plume_index(dict['z'], dict['L'])
        pars = params.gauss_poly[index]
        return pmods.gauss_poly(x, y, Q / dict['U'], pars)
    elif method == 'gaussian2':
        index = plume_index_gaussian(dict['z'], dict['L'])
        pars = params.gaussian2[index]
        sigmay = powerlaw(x, pars['a'], pars['b'])
        sigmaz = powerlaw(x, pars['c'], pars['d'])
        return pmods.gaussian_plume(x, y, z, Q, dict['U'], dict['H'], sigmay, sigmaz)


def predicted2(x, y, z, Q, L, U, H, method):
    if method == 'gaussian':
        stability = stability_class(L)
        pars = params.gaussian[stability]
        sigmay = powerlaw(x, pars['c'], pars['d'])
        sigmaz = powerlaw(x, pars['a'], pars['b'])
        return pmods.gaussian_plume(x, y, z, Q, U, H, sigmay, sigmaz)


###########################################################################################################
def samples(x, y, z, Q, params, method):
    s = [0 if xi < 0 else predicted(xi, yi, z, Q, params, method) for xi, yi in zip(x, y)]
    return s


###########################################################################################################
def line_average(s_loc, p0, p1, height, n, Q, H, theta, T, P, dict, method, m_mass):
    vec = rotate_los(np.array(s_loc), np.array(p0), np.array(p1), theta, n)
    line_samples = np.asarray(samples(vec[0], vec[1], height, Q, dict, method))
    if method != 'semi-gaussian':
        line_samples = [gm_3toppmv(ls, T, P, m_mass) for ls in line_samples]
    lavg = np.mean(line_samples)
    return lavg


###########################################################################################################
def rotate_los(s_loc, p0, p1, theta, n):
    t = np.linspace(0, 1, n)
    points = line(t, p0, p1)
    xvals = points[0] - s_loc[0]
    yvals = points[1] - s_loc[1]
    vec = rotated(np.array([xvals, yvals]), theta)
    return vec


###########################################################################################################
# p0 and p1 are the 2-tuples indicating the start and end of the line over which the line integral is wanted.
# p0 and p1 are of type np.array.
# height is the z co-ordinate.
# n is the number of samples over the spline.
# Q is the source strength in grams per second
# method can be one of gaussian, semi-gaussian, gauss_poly
###########################################################################################################
def line_integral(s_loc, p0, p1, height, n, Q, H, theta, T, P, dict, method, m_mass):
    t = np.linspace(0, 1, n)
    vec = rotate_los(s_loc, p0, p1, theta, n)
    line_samples = samples(vec[0], vec[1], height, Q, dict, method)
    if method != 'semi-gaussian':
        line_samples = [gm_3toppmv(ls, T, P, m_mass) for ls in line_samples]

    conc = simps(line_samples, t)
    line_intgrl = 2 * LA.norm(p1 - p0) * conc
    return line_intgrl
