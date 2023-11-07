import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import scipy.constants as SI
import scipy.optimize as opt

Ee = 7500e9*SI.e
e = SI.e
me = SI.m_e
re = SI.physical_constants['classical electron radius'][0]

n0 = 7e21 #m⁻³
L_plasma = 5 #m

wp = np.sqrt(n0*e**2/me/SI.epsilon_0)
kp = wp/SI.c
E0 = 96e9* np.sqrt(n0/1e24) #V/m
Ez = 0.56*E0 #V/m
tau_r = 2*re/3/SI.c

K = kp/ np.sqrt(2)
A = tau_r * SI.c**2 * K**2
B = SI.c**2 * K**2
C = wp*Ez/E0
D = tau_r * SI.c**2 * K**4

gamma0 = Ee/me/(SI.c**2)

beta_matched = np.sqrt(2*gamma0)/kp
lambda_beta = 2*np.pi*beta_matched
emittance_norm = 10e-6 # m rad

sig_x = np.sqrt(beta_matched*emittance_norm/gamma0)

z0 = 0
x0 = 3*sig_x # m
x_dot0 = 0 # m/s
gamma0 = Ee/me/(SI.c**2)

def oscillator(t, x):
    # x= [x,v_x, gamma]
    v_x = x[1]
    gamma = x[2]
    a_x = -(C/gamma + A)*v_x - B/gamma*x[0]
    d_gamma = C - D*gamma**2*x[0]**2
    return np.array([v_x, a_x, d_gamma])

def evolve_betatron_motion(x_vec, ):
    sysinit = np.array([x0, x_dot0, gamma0, z0, 0])
    solution = solve_ivp(fun = oscillator, y0 = sysinit, method='RK45', t_span = (0,T), t_eval = t)

    x = solution.y[0]
    gamma = solution.y[2]

    return x, gamma
    


