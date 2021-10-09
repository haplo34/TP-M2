#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

c = 1
S = 1
dx = 0.01
dt = dx * S / c
t_0 = 20 * dt
tau = 8 * dt

nbx = 200
nbt = 500
xmin = 0
xmax = xmin + (nbx-1) * dx

x = np.linspace(xmin, xmax, nbx)

unm = np.zeros(nbx)
un = np.zeros(nbx)
unp = np.zeros(nbx)

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
        unp[i] = S**2 * (
            un[i+1] - 2*un[i] + un[i-1]) + 2*un[i] - unm[i]     

    unp[0] = np.exp(-((tnp-t_0) / tau)**2)

    line1.set_data(x, unp)
    line2.set_data(x, -unp)

    unm[:] = un[:]
    un[:] = unp[:]

    return line1, line2,
 
ani = animation.FuncAnimation(
    fig, animate, init_func=init, frames=nbt,
    blit=True, interval=20, repeat=False)

plt.show()