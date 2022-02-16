# -*- coding: utf-8 -*-
# This file is part of the
# 	Solarwindpy Project(https://github.com/adelarja/space_weather).

# Copyright (c) 2021, Adriana Gulisano, Adel Arja, Ricardo Pafundigit,
# Violeta Bazzano.
# All rights reserved.

# License: BSD 3-Clause License
# 	Full Text: https://github.com/adelarja/space_weather/blob/main/LICENSE
"""Module that handles the plots of the project."""

from dataclasses import dataclass
from typing import List, Union

import matplotlib.pyplot as plt

from solarwindpy.data_manager import MagneticField
from solarwindpy.rotation import RotatedWind


@dataclass
class MagneticCloudPlotter:
    """Class used to abstract all the Magnetic Cloud plots."""

    cloud: Union[List[MagneticField], RotatedWind]

    def plot_mf(self, ax=None):
        """Method to plot a no rotated MagneticField."""
        x_components = [magnetic_field.bgse0 for magnetic_field in self.cloud]
        y_components = [magnetic_field.bgse1 for magnetic_field in self.cloud]
        z_components = [magnetic_field.bgse2 for magnetic_field in self.cloud]

        ax = plt.gca() if ax is None else ax

        ax.plot(x_components)
        ax.plot(y_components)
        ax.plot(z_components)

        plt.show()

        return ax

    def plot_rw(self, ax=None):
        """Method to plot a rotated MagneticField."""
        ax = plt.gca() if ax is None else ax

        ax.plot(self.cloud.bgse0)
        ax.plot(self.cloud.bgse1)
        ax.plot(self.cloud.bgse2)

        plt.show()

        return ax
