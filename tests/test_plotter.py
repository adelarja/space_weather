from datetime import datetime

import matplotlib
import matplotlib.pyplot as plt

import pytest

from solarwindpy.data_manager import MagneticField
from solarwindpy.plotter import plot_mf, plot_rw, plot_rw_and_mf
from solarwindpy.rotation import RotatedWind

X_FAKE_VALUES = [1, 2, 3, 4, 5]
Y_FAKE_VALUES = [0.1, 0.2, 0.3, 0.4, 0.5]
Z_FAKE_VALUES = [0.01, 0.02, 0.03, 0.04, 0.05]

FAKE_MAGNETIC_FIELDS = [
    MagneticField(time=datetime(2021, 1, 1), bgse0=0.0, bgse1=1.0, bgse2=2.0),
    MagneticField(time=datetime(2021, 1, 1), bgse0=0.0, bgse1=0.1, bgse2=0.2),
    MagneticField(
        time=datetime(2021, 1, 1), bgse0=0.0, bgse1=0.01, bgse2=0.02
    ),
    MagneticField(
        time=datetime(2021, 1, 1), bgse0=0.0, bgse1=0.001, bgse2=0.002
    ),
]

FAKE_ROTATED_WIND = RotatedWind(
    bgse0=X_FAKE_VALUES, bgse1=Y_FAKE_VALUES, bgse2=Z_FAKE_VALUES
)


@pytest.mark.parametrize(
    "plot_function, arguments",
    [
        (plot_mf, FAKE_MAGNETIC_FIELDS),
        (plot_rw, FAKE_ROTATED_WIND),
    ],
)
def test_plots(plot_function, arguments):
    ax = plot_function(arguments)
    assert isinstance(ax.get_yaxis(), matplotlib.axis.YAxis)
    assert isinstance(ax.get_xaxis(), matplotlib.axis.XAxis)


def test_plot_rw_and_mf():

    _, axs = plt.subplots(2)
    axs = plot_rw_and_mf(FAKE_MAGNETIC_FIELDS, FAKE_ROTATED_WIND, axs)

    assert isinstance(axs[0].get_yaxis(), matplotlib.axis.YAxis)
    assert isinstance(axs[0].get_xaxis(), matplotlib.axis.XAxis)
    assert isinstance(axs[1].get_yaxis(), matplotlib.axis.YAxis)
    assert isinstance(axs[1].get_xaxis(), matplotlib.axis.XAxis)
