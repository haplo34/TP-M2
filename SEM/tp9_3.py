#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

def f(z):
    return np.exp(-2j * np.pi * z) * np.exp(-z**2)

x = np.linspace(-3, 3, 101)
y = np.linspace(-3, 3, 101)
X, Y = np.meshgrid(x, y)

Z = X + 1j * Y
C = np.angle(f(Z)) 


plt.pcolormesh(X, Y, C, shading="gouraud")
plt.colorbar()
plt.xlabel("partie reelle")
plt.ylabel("partie imaginaire")
plt.show()