{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "613672b0-7c09-4684-b4cb-82d680e98c2e",
   "metadata": {},
   "outputs": [],
   "source": [
    "import numpy as np\n",
    "import matplotlib.pyplot as plt\n",
    "from scipy.integrate import solve_ivp"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "565e1a13-4fee-4455-b871-8c39c18b617b",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "975728.7474082205\n",
      "6.59993615865274\n"
     ]
    }
   ],
   "source": [
    "#Constants\n",
    "e = 1.6e-19\n",
    "Ee = 500e9*e\n",
    "re = 2.8179403262e-15 \n",
    "c = 3e8\n",
    "me = 9.11e-31\n",
    "n0 = 1.5e22\n",
    "eps0 = 8.85418782e-12\n",
    "wp = np.sqrt(n0*e**2/me/eps0)\n",
    "kp = wp/c\n",
    "E0 = me*c*wp/e\n",
    "Ez = 0.56*E0 #V/m\n",
    "tau_r = 2*re/3/c\n",
    "K = kp**2 / 2\n",
    "A = tau_r * c**2 * K**2\n",
    "B = c**2 * K**2\n",
    "C = wp*Ez/E0\n",
    "D = tau_r * c**2 * K**4\n",
    "\n",
    "x0 = 1e-10 # m\n",
    "x_dot0 = 0 # m/s\n",
    "gamma0 = Ee/me/(c**2)\n",
    "print(gamma0)\n",
    "print(Ez*1e-9)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "385c5139-f17a-410d-bcf3-a5d649ef29a8",
   "metadata": {},
   "outputs": [],
   "source": [
    "def oscillator(t, x):\n",
    "    # x= [x,v_x, gamma]\n",
    "    v_x = x[1]\n",
    "    a_x = -(C/x[2] + A)*v_x - B/x[2]*x[0]\n",
    "    d_gamma = C - D*x[2]**2*x[0]**2\n",
    "    #d_gamma= 1\n",
    "    return np.array([v_x, a_x, d_gamma])\n",
    "\n",
    "T = 1.6e-9\n",
    "n = int(1e5)\n",
    "t = np.linspace(0,T,n)\n",
    "\n",
    "sysinit = np.array([x0, x_dot0, gamma0])\n",
    "solution = solve_ivp(fun = oscillator, y0 = sysinit, method='RK45', t_span = (0,T), t_eval = t)\n",
    "x_ivp = solution.y[0,:]\n",
    "gamma_ivp = solution.y[2,:]\n",
    "t_ivp = solution.t * c\n",
    "energy_ivp = gamma_ivp*me*c**2"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "419c66ac-9e58-46f7-8a28-1f409bd73193",
   "metadata": {},
   "outputs": [],
   "source": [
    "fig, ax = plt.subplots(2,1, figsize=(8,8))\n",
    "ax[0].plot(t_ivp,x_ivp)\n",
    "ax[0].grid()\n",
    "ax[1].plot(t_ivp,energy_ivp)\n",
    "ax[1].grid()"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.9.12"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
