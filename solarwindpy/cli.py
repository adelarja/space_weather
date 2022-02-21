# -*- coding: utf-8 -*-
# This file is part of the
# 	Solarwindpy Project(https://github.com/adelarja/space_weather).

# Copyright (c) 2021, Adriana Gulisano, Adel Arja, Ricardo Pafundigit,
# Violeta Bazzano.
# All rights reserved.

# License: BSD 3-Clause License
# 	Full Text: https://github.com/adelarja/space_weather/blob/main/LICENSE
"""CLI module for swindpy package usage."""

import csv
from datetime import datetime

import matplotlib.pyplot as plt

from solarwindpy.data_manager import DataManager, Period
from solarwindpy.plotter import plot_mf, plot_rw, plot_rw_and_mf
from solarwindpy.rotation import RotatedWind

import typer

app = typer.Typer()


class InvalidDateException(Exception):
    """Raise this exception when the date is invalid."""

    pass


def valid_dates(date_from: str, date_to: str) -> bool:
    """Check date format and values."""
    try:
        df = datetime.strptime(date_from, "%Y-%m-%d")
        dt = datetime.strptime(date_to, "%Y-%m-%d")
        return (dt - df).days > 0
    except ValueError:
        return False


def get_magnetic_field_data(date_from: str, date_to: str):
    """Function to obtain Magnetic Field information."""
    if not valid_dates(date_from, date_to):
        raise InvalidDateException("Invalid dates.")

    df = datetime.strptime(date_from, "%Y-%m-%d")
    dt = datetime.strptime(date_to, "%Y-%m-%d")
    period = Period(df, dt)

    cdf_data = DataManager.get_gse_magnetic_vector(period)
    filtered_data = DataManager.filter_nan_and_inf_values(cdf_data)

    return filtered_data


@app.command()
def to_csv(date_from: str, date_to: str, filename: str):
    """Magnetic Field data to csv file."""
    try:
        cdf_data = get_magnetic_field_data(date_from, date_to)
    except InvalidDateException:
        return

    with open(filename + ".csv", "w", newline="") as csv_file:
        writer = csv.writer(csv_file)
        for magnetic_field in cdf_data:
            writer.writerow(
                [
                    magnetic_field.time,
                    magnetic_field.bgse0,
                    magnetic_field.bgse1,
                    magnetic_field.bgse2,
                ]
            )


@app.command()
def plot_cloud(date_from: str, date_to: str):
    """Plot a no rotated wind."""
    try:
        cdf_data = get_magnetic_field_data(date_from, date_to)
        plot_mf(cdf_data)
        plt.show()
    except InvalidDateException:
        return


@app.command()
def plot_rotated_cloud(date_from: str, date_to: str):
    """Plot a rotated wind."""
    try:
        filtered_data = get_magnetic_field_data(date_from, date_to)
        rotated = RotatedWind.get_rotated_wind(filtered_data)
        plot_rw(rotated)
        plt.show()
    except InvalidDateException:
        return


@app.command()
def plot_rotated_and_non_rotated(date_from: str, date_to: str):
    """Compare magnetic fields and rotated wind."""
    try:
        fig, axs = plt.subplots(2)
        filtered_data = get_magnetic_field_data(date_from, date_to)
        rotated = RotatedWind.get_rotated_wind(filtered_data)

        plot_rw_and_mf(filtered_data, rotated, axs)
        plt.show()

    except InvalidDateException:
        return


def main():
    """Main cli method."""
    app()
