# -*- coding: utf-8 -*-
# This file is part of the
# 	Solarwindpy Project(https://github.com/adelarja/space_weather).

# Copyright (c) 2021, Adriana Gulisano, Adel Arja, Ricardo Pafundigit,
# Violeta Bazzano.
# All rights reserved.

# License: BSD 3-Clause License
# 	Full Text: https://github.com/adelarja/space_weather/blob/main/LICENSE
"""Module that handles the plots of the project."""

from typing import List

import matplotlib.pyplot as plt

from solarwindpy.data_manager import MagneticField
from solarwindpy.rotation import RotatedWind


def plot_mf(cloud: List[MagneticField], ax=None):
    """Function to plot a no rotated MagneticField."""
    x_components = [magnetic_field.bgse0 for magnetic_field in cloud]
    y_components = [magnetic_field.bgse1 for magnetic_field in cloud]
    z_components = [magnetic_field.bgse2 for magnetic_field in cloud]

    ax = plt.gca() if ax is None else ax

    ax.plot(x_components)
    ax.plot(y_components)
    ax.plot(z_components)

    return ax


def plot_rw(cloud: RotatedWind, ax=None):
    """Function to plot a rotated MagneticField."""
    ax = plt.gca() if ax is None else ax

    ax.plot(cloud.bgse0)
    ax.plot(cloud.bgse1)
    ax.plot(cloud.bgse2)

    return ax


def plot_rw_and_mf(
    cloud: List[MagneticField], rotated_cloud: RotatedWind, axs=None
):
    """Function to plot and compare Magnetic Field and Rotated Cloud."""
    axs = plt.gca() if axs is None else axs

    axs[0].plot([magnetic_field.bgse0 for magnetic_field in cloud])
    axs[0].plot([magnetic_field.bgse1 for magnetic_field in cloud])
    axs[0].plot([magnetic_field.bgse2 for magnetic_field in cloud])

    axs[1].plot(rotated_cloud.bgse0)
    axs[1].plot(rotated_cloud.bgse1)
    axs[1].plot(rotated_cloud.bgse2)

    return axs
