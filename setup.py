# This file is part of the
# Solarwindpy Project(https://github.com/adelarja/space_weather).

# Copyright (c) 2021, Adriana Gulisano, Adel Arja, Ricardo Pafundigit,
# Violeta Bazzano.
# All rights reserved.

# License: BSD 3-Clause License
# Full Text: https://github.com/adelarja/space_weather/blob/main/LICENSE


# =====================================================================
# DOCS
# =====================================================================

"""This file is for distribute and install solarwindpy."""

# ======================================================================
# IMPORTS
# ======================================================================

import os
import pathlib

from setuptools import setup  # noqa

# =============================================================================
# CONSTANTS
# =============================================================================

PATH = pathlib.Path(os.path.abspath(os.path.dirname(__file__)))


REQUIREMENTS = [
    "numpy",
    "pandas",
    "scipy",
    "attrs",
    "matplotlib",
    "heliopy==0.15.4",
    "requests",
    "sunpy",
    "h5netcdf",
    "cdflib",
    "typer",
]

with open(PATH / "solarwindpy" / "__init__.py") as fp:
    for line in fp.readlines():
        if line.startswith("__version__ = "):
            VERSION = line.split("=", 1)[-1].replace('"', "").strip()
            break


with open("README.md") as fp:
    LONG_DESCRIPTION = fp.read()


# =============================================================================
# FUNCTIONS
# =============================================================================

setup(
    name="swindpy",
    version=VERSION,
    description="Package to analyze magnetic storm phenomena.",
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    author="",
    author_email="",
    url="https://github.com/adelarja/space_weather",
    py_modules=[],
    packages=[
        "solarwindpy",
    ],
    license="BSD 3-Clause License",
    install_requires=REQUIREMENTS,
    keywords=[
        "solar",
        "wind",
        "space weather",
        "solarwindpy",
        "magnetic cloud",
    ],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering",
    ],
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "swindpy=solarwindpy.cli:main",
        ],
    },
)
