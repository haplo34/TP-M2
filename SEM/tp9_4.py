#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

x = np.linspace(-3, 13, 161)
dt = 4 / 41

fig = plt.figure()
line, = plt.plot([],[])
plt.xlim(-3, 13)
plt.ylim(-.5,1.5)

def init():
    line.set_data([],[])
    return line,

def animate(i):
    t = i * dt
    y = 1 / (1 + (x - 4.8 * t)**2)
    line.set_data(x, y)
    return line,
 
ani = animation.FuncAnimation(
    fig, animate, init_func=init, frames=100,
    blit=True, interval=20, repeat=False)

plt.show()