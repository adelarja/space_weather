# This file is part of the
# Solarwindpy Project(https://github.com/adelarja/space_weather).

# Copyright (c) 2021, Adriana Gulisano, Adel Arja, Ricardo Pafundigit,
# Violeta Bazzano.
# All rights reserved.

# License: BSD 3-Clause License
# Full Text: https://github.com/adelarja/space_weather/blob/main/LICENSE

from datetime import datetime
from unittest import mock

import numpy as np

import pandas as pd
from pandas import Timestamp

from solarwindpy.cli import app
from solarwindpy.data_manager import MagneticField

from typer.testing import CliRunner


runner = CliRunner()

FAKE_DATAFRAME = pd.DataFrame(
    {
        "Time": [
            datetime(2021, 1, 1, 0, 0, 0),
            datetime(2021, 1, 1, 0, 30, 0),
            datetime(2021, 1, 1, 1, 0, 0),
        ],
        "BGSE_0": [1.0, -1.0, 0.0],
        "BGSE_1": [2.0, -2.0, 1.0],
        "BGSE_2": [3.0, -3.0, 2.0],
    }
).set_index("Time")

EXPECTED_CSV_RESULT = (
    "2021-01-01 00:00:00,1.0,2.0,3.0\n"
    "2021-01-01 00:30:00,-1.0,-2.0,-3.0\n"
    "2021-01-01 01:00:00,0.0,1.0,2.0\n"
)

FAKE_MAGNETIC_FIELDS = [
    MagneticField(Timestamp(2021, 1, 1, 0, 0, 0), np.inf, 2.0, 3.0),
    MagneticField(Timestamp(2021, 1, 1, 0, 30, 0), np.inf, np.inf, np.inf),
    MagneticField(Timestamp(2021, 1, 1, 0, 30, 0), np.nan, -2.0, -3.0),
    MagneticField(Timestamp(2021, 1, 1, 0, 30, 0), np.nan, np.nan, np.nan),
    MagneticField(Timestamp(2021, 1, 1, 0, 0, 0), 1.0, 2.0, 3.0),
    MagneticField(Timestamp(2021, 1, 1, 0, 30, 0), -1.0, -2.0, -3.0),
    MagneticField(Timestamp(2021, 1, 1, 1, 0, 0), 0.0, 1.0, 2.0),
]


def get_fake_dataframe():
    return FAKE_DATAFRAME


def test_to_csv(tmpdir):
    p = tmpdir.join("test.csv")
    path = str(p)
    fake_generic_series_object = mock.MagicMock()
    fake_generic_series_object.to_dataframe = get_fake_dataframe
    with mock.patch(
        "heliopy.data.wind.mfi_h0", return_value=fake_generic_series_object
    ):
        result = runner.invoke(
            app, ["to-csv", "2021-01-01", "2021-01-02", path[:-4]]
        )
        assert result.exit_code == 0
        assert "".join(p.read()) == EXPECTED_CSV_RESULT
