# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 15:42:34 2021

@author: aguli
"""
import numpy as np

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
