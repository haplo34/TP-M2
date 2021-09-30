#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from sys import argv
import matplotlib.pyplot as plt

try:
    # Marcheurs
    M = int(argv[1])
    # Nombre de pas
    N = int(argv[2])
except ValueError:
    print("Les arguments doivent être des entiers !!")
    
# Positions
P = np.zeros((M, N))

for i in range(M):
    for j in range(N-1):
        P[i][j+1] = P[i][j] + np.random.choice([-1, 1])
        

plt.subplot(3, 1, 1)
for i in range(M):
    plt.plot(P[i])
    
plt.subplot(3, 1, 2)
for i in range(M):
    plt.plot(P[i]**2)

plt.subplot(3, 1, 3)
plt.plot(np.mean(P**2, axis=0))

plt.show()

