#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return np.sin(x) / x * np.sqrt(x**2 +1)

I = (-6, 6)
N = 121
x = I[0]
dx = (abs(I[1] - I[0])) / N
X = []
Y = []

while x < I[1]:
    if x != 0:
        X.append(x)
        Y.append(f(x))
    x += dx
    
plt.plot(X, Y)
plt.xlabel("x")
plt.ylabel("f(x)")
plt.show()