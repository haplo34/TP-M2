#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from sys import argv

coord = [[0, 0], [6, 0], [3, np.sqrt(27)], [0, 0]]

xs, ys = zip(*coord) #create lists of x and y values

points = [[0.7], [1]]

for i in range(int(argv[1])):
    rd = np.random.choice([1, 2, 3])
    if rd==1:
        points[0].append(points[0][i] / 2)
        points[1].append(points[1][i] / 2)
    elif rd==2:
        points[0].append((points[0][i] + 6) / 2)
        points[1].append(points[1][i] / 2)
    elif rd==3:
        points[0].append((points[0][i] + 3) / 2)
        points[1].append((points[1][i] + np.sqrt(27))/ 2)

plt.plot(xs, ys)
plt.scatter(points[0], points[1])
plt.show()

