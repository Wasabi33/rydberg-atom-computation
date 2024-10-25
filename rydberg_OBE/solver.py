
# Contains the complex differential equation solver

import matplotlib.pyplot as plt  # Import library for direct plotting functions
import numpy as np  # Import Numerical Python
from IPython.display import display, HTML  # Import HTML for formatting output
from scipy.integrate import ode  # Import ODE solver
import scipy.constants as cts  # Import physical constants
from arc import *
from scipy.integrate import solve_ivp
from numpy.linalg import eig

def complex_de_solver(eqs, y0, time_span, parameters):
    # the complex differential equation solver
    # eqs is the function that defines the differential equations
    # y0 is the initial conditions
    # t is the time array
    # parameters is the parameters array
    # rhos are the arrays to store the results
    r = solve_ivp(eqs, time_span, y0, method = 'RK45', args = parameters, vectorized = True)
    rhos = r.y
    t = r.t
    return t, rhos

def matrix_solver(t,y0,Omega1, Omega2, gamma1, gamma2, delta1, delta2, b1, b2):
    obef3 = -1j*np.array([[-gamma2,0,0,1j/2*Omega2,-1j/2*Omega2, 0, 0, 0, 0], 
                    [b2*gamma2,-gamma1,0, -1j/2*Omega2, 1j/2*Omega2, 1j/2*Omega1, -1j/2*Omega1, 0, 0],
                    [0, b1*gamma1, 0, 0, 0, -1j/2*Omega1, 1j/2*Omega1, 0, 0],
                    [1j/2*Omega2, -1j/2*Omega2, 0, (-(gamma2+gamma1)/2+1j*delta2), 0, 0, 0, 1j/2*Omega1, 0],
                    [-1j/2*Omega2, 1j/2*Omega2, 0, 0, (-(gamma2+gamma1)/2-1j*delta2), 0, 0, 0, -1j/2*Omega1],
                    [0, 1j/2*Omega1, -1j/2*Omega1, 0, 0, -gamma1/2+1j*delta1, 0, -1j/2*Omega2, 0],
                    [0, -1j/2*Omega1, 1j/2*Omega1, 0, 0, 0, -gamma1/2-1j*delta1, 0, 1j/2*Omega2],
                    [0, 0, 0, 1j/2*Omega1, 0, -1j/2*Omega2, 0, -gamma2/2+1j*(delta1+delta2),0],
                    [0, 0, 0, 0, -1j/2*Omega1, 0, 1j/2*Omega2, 0,-gamma2/2-1j*(delta1+delta2)]],dtype=complex)
    
    w, v = eig(obef3)
    c = np.linalg.solve(v, y0)
    
    rho_ee = c[0]*v[0,0]*np.exp(1j*w[0]*t) + c[1]*v[0,1]*np.exp(1j*w[1]*t) + \
    c[2]*v[0,2]*np.exp(1j*w[2]*t) + c[3]*v[0,3]*np.exp(1j*w[3]*t) + \
    c[4]*v[0,4]*np.exp(1j*w[4]*t) + c[5]*v[0,5]*np.exp(1j*w[5]*t) + \
    c[6]*v[0,6]*np.exp(1j*w[6]*t) + c[7]*v[0,7]*np.exp(1j*w[7]*t) + \
    c[8]*v[0,8]*np.exp(1j*w[8]*t)


    return rho_ee



