import numpy as np


def gaussian(x, sigma):
    den = sigma * np.sqrt(2 * np.pi)
    ret = (np.exp(-x ** 2 / (2 * sigma ** 2))) / den
    return ret


# x downwind distance in meters
# y crosswind distance in meters
# Q source emission rate in grams per second
# U average wind speed at stack height in meters per second
# H stack height. Since we are not modeling the plume rise, this is the same as the height of the source.
# ret ground level concentration at x,y in grams per cubic meter
def gaussian_plume(x, y, z, Q, U, H, sigmay, sigmaz):
    exp1 = gaussian(y, sigmay)
    exp2 = gaussian(z - H, sigmaz) + gaussian(z + H, sigmaz)
    conc = (Q / U) * exp1 * exp2
    return conc

def semigaussian(x, y, Q, params):
    A = params['a']
    B = params['b']
    C = params['c']
    D = params['d']
    E = params['e']
    F = params['f']
    G = params['g']
    H = params['h']
    exp1 = (B * (x ** 2) * (C + y ** D) ** 2) / (y ** E)
    exp2 = F / (y ** G)
    conc = (Q * A * np.exp(-exp1) * np.exp(-exp2)) / (y ** H)
    return conc


def gauss_poly(x, y, Q, params):
    c0 = params['c0']
    c1 = params['c1']
    c2 = params['c2']
    c3 = params['c3']
    amplitude = params['amplitude']
    exponent = params['exponent']
    height = c0 + c1 * x + c2 * x ** 2 + c3 * x ** 3
    width = amplitude * (x ** exponent)
    conc = Q * height * np.exp(-y ** 2 / (2 * (width ** 2)))
    return conc
