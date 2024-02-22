import scipy.constants as SI
import numpy as np

def dgamma(gamma, x, C, D):
    return C - D*gamma**2*(x[0]**2 + x[1]**2)
    
def acc_func(y, A, B, C, D):
    # y = [x, y, vx, vy, gamma]
    dy = np.zeros(y.size) # [vx, vy, ax, xy, dgamma]

    a = -(A/y[-1] + B)*y[2:4] - C/y[-1]*y[:2]
    dgamma = A - D*y[-1]**2*(y[0]**2 + y[1]**2)
    dy[:2] = y[2:4]
    dy[2:4] = a
    dy[-1] = dgamma
    return dy

def integrator(z, y0, n0):
    # y0 = [x0, vx0, y0, vy0, gamma0]
    dz = z[1]-z[0]
    solution = np.zeros((y0.size, z.size))
    solution[:,0] = y0
    
    re = SI.physical_constants['classical electron radius'][0]
    wp = np.sqrt(n0*SI.e**2/SI.m_e/SI.epsilon_0)
    kp = wp/SI.c
    E0 = SI.m_e*SI.c*wp/SI.e
    Ez = 3.2e9 #V/m
    tau_r = 2.*re/3./SI.c
    
    K = kp/ np.sqrt(2)
    A = kp*Ez/E0
    C = K**2
    B = tau_r * SI.c * C
    D = tau_r * SI.c * C**2

    
        
    for i in range(z.size -1):
        y = solution[:,i]
        k1 = acc_func(y, A, B, C, D)
        
        k2 = acc_func(y + k1 * dz/2, A, B, C, D)
        
        k3 = acc_func(y + k2 * dz/2, A, B, C, D)

        k4 = acc_func(y + k3 * dz, A, B, C, D)
        
        k_av = 1/6*(k1+2*k2+2*k3+k4)
        
        solution[:,i+1] = y + k_av*dz
        if i%(round(z.size/20)) < 0.5:
            print(round(i/z.size * 100), '%', '\n',solution[0,i+1]*1e6)
        if solution[0,i+1] > 100e-6 or solution[2,i+1] > 100e-6:
            print('Divergence')
            return

        """
        gamma = solution[4,i]
        numerator = 1-dz/2 *(C/gamma + A)
        denominator = 1+dz/2 *(C/gamma + A)

        solution[1,i+1] = solution[1,i] * numerator/denominator - B/gamma*solution[0,i]*dz/denominator
        solution[3,i+1] = solution[3,i] * numerator/denominator - B/gamma*solution[2,i]*dz/denominator
        solution[0,i+1] += solution[1,i+1]*dz
        solution[2,i+1] += solution[3,i+1]*dz
        """
    return z, solution