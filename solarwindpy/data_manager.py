# -*- coding: utf-8 -*-
# This file is part of the
# 	Solarwindpy Project(https://github.com/adelarja/space_weather).

# Copyright (c) 2021, Adriana Gulisano, Adel Arja, Ricardo Pafundigit,
# Violeta Bazzano.
# All rights reserved.

# License: BSD 3-Clause License
# 	Full Text: https://github.com/adelarja/space_weather/blob/main/LICENSE
"""This module helps to obtain Magnetic Fields Data."""

from dataclasses import dataclass
from datetime import datetime

import heliopy.data.wind as wind

import numpy as np


@dataclass
class Period:
    """Dataclass to manage the dates involved in the analysis.

    TODO: Add some properties/setters to obtain relative dates. For example,
        to obtain relative dates from 'date_to', receiving seconds, minutes,
        hours, days, months, etc. The same for date_from.
    """

    date_from: datetime
    date_to: datetime


@dataclass
class MagneticField:
    """A dataclass to represent a Magnetic Field components.

    The aim of this dataclass, is to make it easier to perform the
    calculations when rotating or when doing any process to the data.

    Instead of working over a DataFrame row, we will work over a MagneticField
    (MagneticField has the data corresponding to one row of the dataframe
    obtained using DataManager object).
    """

    time: datetime
    bgse0: float
    bgse1: float
    bgse2: float


def not_nan_neither_inf(magnetic_field: MagneticField):
    """Method used to filter invalid MAgneticField values."""
    invalid_gse0 = np.isnan(magnetic_field.bgse0) or np.isinf(
        magnetic_field.bgse0
    )
    invalid_gse1 = np.isnan(magnetic_field.bgse1) or np.isinf(
        magnetic_field.bgse1
    )
    invalid_gse2 = np.isnan(magnetic_field.bgse2) or np.isinf(
        magnetic_field.bgse2
    )
    return not any([invalid_gse0, invalid_gse1, invalid_gse2])


class DataManager:
    """This class abstracts methods needed to obtain Magnetic Fields Data."""

    @staticmethod
    def get_gse_magnetic_vector(period: Period) -> list[MagneticField]:
        """Obtain the Magnetic Field Vector in GSE.

        GSE is the acronum for Geocentric Solar Ecliptic System.

        Args:
            period (Period): A Period object with the dato_from and date_to
            information.


        Returns:
            A List of MagneticField objects that would be used to perform
            calculations.

        Example:
            >>> date_from = datetime(2021, 1, 1, 0, 0, 0)
            >>> date_to = datetime(2021, 1, 1, 1, 0, 0)
            >>> data_period = Period(date_from, date_to)
            >>> cdf_file_data = (
                    DataManager.get_gse_magnetic_vector(data_period)
                )

        TODO: What happens if we don't retrieve any data? Or if we passed a
           date_from greater than date_to?
        """
        cdf_data = wind.mfi_h0(period.date_from, period.date_to)
        cdf_data = cdf_data.to_dataframe()
        cdf_data = cdf_data.loc[:, ["BGSE_0", "BGSE_1", "BGSE_2"]]

        return [
            MagneticField(
                time,
                components["BGSE_0"],
                components["BGSE_1"],
                components["BGSE_2"],
            )
            for time, components in cdf_data.iterrows()
        ]

    @staticmethod
    def filter_nan_and_inf_values(
        magnetic_fields_data: list[MagneticField],
    ) -> list[MagneticField]:
        """Filter the nan and inf components in the magnetic field elements.

        If any of the components has an inf or nan value, the entire
        magnetic field element is deleted from the list.

        Args:
            magnetic_fields_data (list[MagneticField]): Data about magnetic
                field in gse coordinates.

        Returns:
            A list of MagneticField elements after filtering the
            MagneticField elements with one or more coordinates with
            NaN or INF values.
        """
        return [
            magnetic_field
            for magnetic_field in magnetic_fields_data
            if not_nan_neither_inf(magnetic_field)
        ]
