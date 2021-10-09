# -*- coding: utf-8 -*-
#!/usr/bin/python3

"""
Résolution équation de Poisson (électrostatique) via Gauss-Seidel
"""


import numpy as np
import matplotlib.pyplot as plt


#Création de la densité de charges
RHO = np.zeros((50, 50), dtype=float)
for i in range(10, 40):
    RHO[i,25] = 1


#Variables et arguments
H = 1                                                                          #Pas de discrétisation
IT_MAX = 100                                                                   #Nombre de boucles
W = 2/(1+np.pi/RHO.shape[0])                   	                               #Paramètre de relaxation
U = np.zeros((RHO.shape[0]+2, RHO.shape[0]+2), dtype=float)                    #Potentiel électrosatique


#Algorithme de Gauss-Seidel
def _gs_(u, w, h, rho):
    for i in range(rho.shape[0]):
        for j in range(rho.shape[0]):
            u[i+1,j+1] = ((1-w)
                *u[i+1,j+1]+(w/4)
                *(u[i,j+1]+u[i+2,j+1]+u[i+1,j]+u[i+1,j+2]-4*h*h*rho[i,j]))
    return u


#Appel de Gauss-Seidel jusqu'à la convergence
def _poisson_(v, it_max):
    conv = []
    for i in range(it_max):
        v = _gs_(v, W, H, RHO)
        conv.append(np.sqrt(np.sum(v**2))/(v.shape[0]-2)**2)
        if len(conv) >= 2:
            if round(conv[-1], 9) == round(conv[-2], 9):
                break
    return v


#plt.matshow(_poisson_(U, IT_MAX), cmap='gray')
plt.contour(_poisson_(U, IT_MAX), 20)
figure = plt.gcf()
plt.savefig("fil_charge_2.png", bbox_inches='tight')