#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 10:03:31 2021

@author: adriana
"""

import numpy as np
def clean(bx,by,bz, error_flag=1e+20):
    bx[abs(bx) < error_flag] = np.NaN
    by[abs(by) < error_flag] = np.NaN
    bz[abs(bz) < error_flag] = np.NaN
    return bx, by, bz