# -*- coding: utf-8 -*-
# This file is part of the
#   Instituto de Astronomico y Fisica del Espacio (IAFE, CONICET-UBA)
# Copyright (c) 2021, Departamento de Fisica, FCEN (UBA).
# v1.2, by S. Dasso, Dec 20, 2007.
#   Full Text: https://github.com/adelarja/space_weather/blob/main/LICENSE.md

# =====================TEMPORARY (REMOVE AFTERWARDS)===========================
'''
Q: El metodo rotate funciona para nubes y vientos? Debe mantenerse generico?
A:

Q: Cuando podemos tener una muestra del archivo cdf?
A:

A:
'''
# =============================================================================


# =============================================================================
# IMPORTS
# =============================================================================

import sys
import numpy as np
import pandas as pd
import scipy.linalg as la
import cdflib

# ============================================================================
# CLASSES
# ============================================================================

class Wind():

    def __init__(self):
        self.M = 1
        self.epsilon = 1e-4
        pass

    def clean(self):
        # Remove data gaps and produce for bx_gse inside the cloud

        pass

    def ordered_eigen(self, M):
        if M.shape[0] == M.shape[1]:
            eigvals, eigvecs = la.eig(M)
            #do ordering
        else:
            sys.exit('Matrix should be square')

    def get_cloud(self):
        pass

    def plot(self):
        pass

    def rotate(self, x_versor, y_versor, z_versor, bx, by, bz):
        '''
        :param x_versor:
        :param y_versor:
        :param z_versor:
        :param bx:
        :param by:
        :param bz:
        :return: cloud/wind three rotated coordinates
        '''
        bx_rotated = bx * x_versor[1] + by * x_versor[2] + bz * x_versor[3]
        by_rotated = bx * y_versor[1] + by * y_versor[2] + bz * y_versor[3]
        bz_rotated = bx * z_versor[1] + by * z_versor[2] + bz * z_versor[3]
        return bx_rotated, by_rotated, bz_rotated


# ============================================================================
# FUNCTIONS
# ============================================================================


def read_cdf(input_file):
    '''
    Based on library https://github.com/MAVENSDC/cdflib
    :param input_file:
    :return: n-dimensional array that represents a solar wind
    '''
    cdf_file = cdflib.CDF(input_file)
    return cdf_file


cdf_file = read_cdf('./data/ac_h2s_swe_20100401000000_20100406000000.cdf')
cdf_file.cdf_info()