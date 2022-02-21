# This file is part of the
# Solarwindpy Project(https://github.com/adelarja/space_weather).

# Copyright (c) 2021, Adriana Gulisano, Adel Arja, Ricardo Pafundigit,
# Violeta Bazzano.
# All rights reserved.

# License: BSD 3-Clause License
# Full Text: https://github.com/adelarja/space_weather/blob/main/LICENSE

from datetime import datetime
from unittest import mock

import pandas as pd

import pytest

from solarwindpy.cli import app

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


def test_plot_cloud(mocker):
    fake_generic_series_object = mock.MagicMock()
    fake_generic_series_object.to_dataframe = get_fake_dataframe
    mocker.patch(
        "heliopy.data.wind.mfi_h0", side_effect=fake_generic_series_object
    )
    mocker.patch("matplotlib.pyplot.show", return_value=True)
    mocker.patch("solarwindpy.plotter.plot_mf")

    result = runner.invoke(app, ["plot-cloud", "2021-01-01", "2021-01-02"])
    assert result.exit_code == 0


@pytest.mark.parametrize(
    "command, date_1, date_2",
    [
        ("plot-cloud", "fake_date_1", "fake_date_2"),
        ("plot-cloud", "2021-01-01", "fake_date_2"),
        ("plot-cloud", "2021-01-02", "2021-01-01"),
        ("plot-rotated-cloud", "fake_date_1", "fake_date_2"),
        ("plot-rotated-cloud", "2021-01-01", "fake_date_2"),
        ("plot-rotated-cloud", "2021-01-02", "2021-01-01"),
        ("plot-rotated-and-non-rotated", "fake_date_1", "fake_date_2"),
        ("plot-rotated-and-non-rotated", "2021-01-01", "fake_date_2"),
        ("plot-rotated-and-non-rotated", "2021-01-02", "2021-01-01"),
    ],
)
def test_invalid_dates_exception(command, date_1, date_2):
    result = runner.invoke(app, [command, date_1, date_2])
    assert result.exit_code == 0
    assert result.return_value is None
