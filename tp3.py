#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 14:09:18 2021

@author: gromoff-q
"""
import numpy as np
from numpy import linalg as LA


A = np.array([[complex(), complex(0, -1)],
               [complex(0, 1), complex()]])


print(A)
print(np.conj(A.T))
print((np.conj(A.T)==A).all())
print("A = A* <==> A Hermitien.\n")


valp, vecp = LA.eig(A)
print("Valeurs propres et vecteurs propres:")
print(valp)
print(vecp)
print('\n')


a1 = vecp[:, :1]
a2 = vecp[:, 1:2]

print(a1)
print(a2)
print('\n')
print("Produit scalaire des vecteurs propres:")
print(np.dot(a1.T, a2))
print('\n')

p1 = np.dot(a1, np.conj(a1.T))
p2 = np.dot(a2, np.conj(a2.T))

print("Projecteurs P1 et P2:")
print(p1)
print(p2)
print('\n')

p3 = p1 + p2
I = np.identity(2)

print("P1 + P2:")
print(p3)
print(np.allclose(p3, I))
print('\n')

p1p2 = np.dot(p1, p2)
p2p1 = np.dot(p2, p1)
M0 = np.zeros(2)

print("Produits matriciels P1P2 et P2P1:")
print(p1p2)
print(np.allclose(p1p2, M0))
print(p2p1)
print(np.allclose(p2p1, M0))
print('\n')

