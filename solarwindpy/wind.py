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

Q: COmo definimos una nueva en la clase de solardwind?
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

class SolarWind():

    def __init__(self):
        self.M = 1
        self.epsilon = 1e-4
        pass

    def clean_wind(self):
        # Remove data gaps and produce for bx_gse inside the cloud
        '''
                epsilon=1e-4;
        if gap>0,
           cond_good=find(bx_gse_ext_<gap-epsilon & by_gse_ext_<gap-epsilon & bz_gse_ext_<gap-epsilon);
        elseif gap<0,
           cond_good=find(bx_gse_ext_>gap+epsilon & by_gse_ext_>gap+epsilon & bz_gse_ext_>gap+epsilon);
        else
           'ERROR: gap flag cannot be zero'
        end
        bx_gse_ext=bx_gse_ext_(cond_good);
        by_gse_ext=by_gse_ext_(cond_good);
        bz_gse_ext=bz_gse_ext_(cond_good);
        date_ext=date_ext_(cond_good);
        B=B__(cond_good);
        clear cond_good
                :return:
        '''
        pass


    def get_cloud(self):
        '''
        Tiene que recibir como input las fechas

                epsilon=1e-4;
        if gap>0,
           cond_good=find(bx_gse_ext_<gap-epsilon & by_gse_ext_<gap-epsilon & ...
           bz_gse_ext_<gap-epsilon & date_ext_>=initial_date & date_ext_<=end_date);
        elseif gap<0,
           cond_good=find(bx_gse_ext_>gap+epsilon & by_gse_ext_>gap+epsilon & ...
           bz_gse_ext_>gap+epsilon & date_ext_>=initial_date & date_ext_<=end_date);
        else
           'ERROR: gap flag cannot be zero'
        end
        % Comment: for simplicity in notation define bi as b_gse_i (i=x,y,z)
        bx=bx_gse_ext_(cond_good);
        by=by_gse_ext_(cond_good);
        bz=bz_gse_ext_(cond_good);
        date_=date_ext_(cond_good);
        B_=B__(cond_good);
        clear cond cond_good
                :return:
        '''
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