#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

c = 2.99792458e8
S = 1
N_lambda = 15
lambda_0 = 1.55e-6
lambda_1 = lambda_0 / 1.45
dx = lambda_0 / N_lambda
dt = dx * S / c
T = lambda_0 / c

nbx = int((10 * lambda_0 / dx) + 1)
nbt = 1000
xmin = 0
xmax = 20 * lambda_0

x = np.linspace(xmin, xmax, nbx)

unm = np.zeros(nbx)
un = np.zeros(nbx)
unp = np.zeros(nbx)

esp_r = np.ones(nbx)
#mid = int(nbx / 2)
#esp_r[mid-10:mid+10] = 1.45**2
for i in range(nbx):
    if x[i] > 10*lambda_0 - lambda_1/8 and x[i] < 10*lambda_0 + lambda_1/8:
        esp_r[i] = 1.45**2

fig = plt.figure() # initialise la figure
line1, = plt.plot([],[])
line2, = plt.plot([],[])

plt.xlim(xmin, xmax)
plt.ylim(-1.5,1.5)

# fonction à définir quand blit=True
# crée l'arrière de l'animation qui sera présent sur chaque image
def init():
    line1.set_data([],[])
    line2.set_data([],[])
    return line1, line2,

def animate(n): 
    tnp = (n+1) * dt

    for i in range(1, nbx-1):
        unp[i] = S**2 / esp_r[i] * (
            un[i+1] - 2*un[i] + un[i-1]) + 2*un[i] - unm[i]     

    unp[0] = np.cos(2*np.pi * tnp / T) * \
             np.exp(-((tnp - 6*T) / (2*T))**2)

    line1.set_data(x, unp)
    line2.set_data(x, esp_r)

    unm[:] = un[:]
    un[:] = unp[:]

    return line1, line2,
 
ani = animation.FuncAnimation(
    fig, animate, init_func=init, frames=nbt,
    blit=True, interval=20, repeat=False)

plt.show()
