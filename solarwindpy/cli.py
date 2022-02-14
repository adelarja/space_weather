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

from solarwindpy.data_manager import DataManager, Period

import typer

app = typer.Typer()


def valid_dates(date_from: str, date_to: str) -> bool:
    """Check date format and values."""
    try:
        df = datetime.strptime(date_from, "%Y-%m-%d")
        dt = datetime.strptime(date_to, "%Y-%m-%d")
        return (dt - df).days > 0
    except ValueError:
        return False


@app.command()
def to_csv(date_from: str, date_to: str, filename: str):
    """Magnetic Field data to csv file."""
    if not valid_dates(date_from, date_to):
        return

    df = datetime.strptime(date_from, "%Y-%m-%d")
    dt = datetime.strptime(date_to, "%Y-%m-%d")
    period = Period(df, dt)

    cdf_data = DataManager.get_gse_magnetic_vector(period)

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


def main():
    """Main cli method."""
    app()
