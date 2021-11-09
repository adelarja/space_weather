# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 20:07:03 2021

@author: aguli
"""
import matplotlib.pyplot as plt

import numpy as np

from rotar import calculo_gamma as cg

import scipy.special as sp

B0 = 10
x1 = np.linspace(-2.38, 2.38, 100)
bx_nube = np.zeros(100)
by_nube = sp.jv(1, x1) * 10
bz_nube = sp.jv(0, x1) * 10


x_gse = np.array([1, 0, 0]).reshape(3, 1)
y_gse = np.array([0, 1, 0]).reshape(3, 1)
z_gse = np.array([0, 0, 1]).reshape(3, 1)


theta_deg = 30
phi_deg = 45
# from deg to radians
theta = np.deg2rad(theta_deg)
phi = np.deg2rad(phi_deg)

gamma = cg(phi, theta)
# building  the R matrix
R = np.zeros((3, 3))
R[0][0] = np.cos(gamma) * np.sin(theta) * np.cos(phi) - np.sin(gamma) * np.sin(
    phi
)
R[0][1] = np.cos(gamma) * np.sin(theta) * np.sin(phi) + np.sin(gamma) * np.cos(
    phi
)
R[0][2] = -np.cos(gamma) * np.cos(theta)
R[1][0] = -np.sin(gamma) * np.sin(theta) * np.cos(phi) - np.cos(
    gamma
) * np.sin(phi)
R[1][1] = -np.sin(gamma) * np.sin(theta) * np.sin(phi) + np.cos(
    gamma
) * np.cos(phi)
R[1][2] = np.sin(gamma) * np.cos(theta)
R[2][0] = np.cos(theta) * np.cos(phi)
R[2][1] = np.cos(theta) * np.sin(phi)
R[2][2] = np.sin(theta)


# applying the rotation by hand
x_mc = R[0][0] * x_gse + R[0][1] * y_gse + R[0][2] * z_gse
y_mc = R[1][0] * x_gse + R[1][1] * y_gse + R[1][2] * z_gse
z_mc = R[2][0] * x_gse + R[2][1] * y_gse + R[2][2] * z_gse

RT = np.transpose(R)

x_gse_ = R @ x_mc
y_gse_ = R @ y_mc
z_gse_ = R @ z_mc

bx_gse = bx_nube * x_mc[0] + by_nube * x_mc[1] + bz_nube * x_mc[2]
by_gse = bx_nube * y_mc[0] + by_nube * y_mc[1] + bz_nube * y_mc[2]
bz_gse = bx_nube * z_mc[0] + by_nube * z_mc[1] + bz_nube * z_mc[2]


plt.figure(1)
plt.plot(x1, by_nube, label='By_nube')
plt.plot(x1, bz_nube,label='Bz_nube')
plt.plot(x1, bx_nube,label='Bx_nube')
plt.title('Components of Magnetic Field')
plt.xlabel('time')
plt.ylabel('nT')
plt.legend()
#plt.show()
plt.savefig('Figura_Nube.png')

plt.figure(2)
plt.plot(x1, by_nube, label='By_gse')
plt.plot(x1, bz_nube,label='Bz_gse')
plt.plot(x1, bx_nube,label='Bz_gse')
plt.title('Components of Magnetic Field')
plt.xlabel('time')
plt.ylabel('nT')
plt.legend()
#plt.show()
plt.savefig('Figura_Nube_gse.png')


bx_mc = RT[0][0] * bx_gse + RT[0][1] * by_gse + RT[0][2] * bz_gse
by_mc = RT[1][0] * bx_gse + RT[1][1] * by_gse + RT[1][2] * bz_gse
bz_mc = RT[2][0] * bx_gse + RT[2][1] * by_gse + RT[2][2] * bz_gse

plt.figure(3)
plt.plot(x1, by_nube, label='By_mc')
plt.plot(x1, bz_nube,label='Bz_mc')
plt.plot(x1, bx_nube,label='Bz_mc')
plt.title('Components of Magnetic Field')
plt.xlabel('time')
plt.ylabel('nT')
plt.legend()
#plt.show()
plt.savefig('Figura_Nube_mc.png')

np.save("bx.npy", bx_gse)  # save
np.save("by.npy", by_gse)  # save
np.save("bz.npy", bz_gse)  # save
