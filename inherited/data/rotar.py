# -*- coding: utf-8 -*-
#This file is part of the
#	Solarwindpy Project(https://github.com/adelarja/space_weather).

#Copyright (c) 2021, Adriana Gulisano, Adel Arja, Ricardo Pafundigit, Violeta Bazzano.
#All rights reserved.

#License: BSD 3-Clause License
#	Full Text: https://github.com/adelarja/space_weather/blob/main/LICENSE


import numpy as np
from cycler import cycler

# 'Rotate Magnetic clouds from GSE_Coordinates as input'
# 'Delivers MC_Coordinates as output'
# these should be replace at some point por the rotation functionality of
# numpy or scypy when we understand how to apply it


def norm_vec(mc_cordinates):
    normvec = (
        mc_cordinates[0] ** 2 + mc_cordinates[1] ** 2 + mc_cordinates[2] ** 2
    ) ** 0.5
    return normvec


def calculo_gamma(phi, theta):
    """
    Here we test if x_MC dot x_gse >0
    """
    epsilon = 10e-15
    gamma = np.arctan(-np.tan(phi) / np.sin(theta + epsilon))
    xmc_dot_xgse = np.cos(gamma) * np.sin(theta) * np.cos(phi) - np.sin(
        gamma
    ) * np.sin(phi)
    if xmc_dot_xgse < 0:
        gamma = np.arctan(-np.tan(phi) / np.sin(theta + epsilon)) + np.pi
        xmc_dot_xgse = np.cos(gamma) * np.sin(theta) * np.cos(phi) - np.sin(
            gamma
        ) * np.sin(phi)
        if xmc_dot_xgse < 0:
            print("Error in gamma cuadrant")
    return gamma


def rotacion(x_gse, y_gse, z_gse, theta_deg, phi_deg):
    # from deg to radians
    theta = theta_deg * np.pi / 180.0
    phi = phi_deg * np.pi / 180.0
    gamma = calculo_gamma(phi, theta)
    # building  the R matrix
    R = np.zeros((3, 3))
    R[0][0] = np.cos(gamma) * np.sin(theta) * np.cos(phi) - np.sin(
        gamma
    ) * np.sin(phi)
    R[0][1] = np.cos(gamma) * np.sin(theta) * np.sin(phi) + np.sin(
        gamma
    ) * np.cos(phi)
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

    return x_mc, y_mc, z_mc

def plot(time, bx, by, bz, ax=None, scatter_bx=None,  scatter_by=None, scatter_bz=None):

    """ The aim of this class is to create a plot of the components of Cloud's Magnetic Field
    """
    ax = plt.gca() if ax is None else ax
    scatter_bx.setdefault("linewidth", 5)
    scatter_by.setdefault("linewidth", 5)
    scatter_bz.setdefault("linewidth", 5)

    scatter_bx.setdefault("color", "r")
    scatter_by.setdefault("color", "g")
    scatter_bz.setdefault("color", "b")


  #  ax.scatter = {} if scatter_kws is None else scatter_kws

    ax.plot(time, bx, **scatter_bx)
    ax.plot(time, by, **scatter_by)
    ax.plot(time, bz, **scatter_bz)

    return ax
