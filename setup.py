# =====================================================================
# DOCS
# =====================================================================

"""Add a description."""

# ======================================================================
# IMPORTS
# ======================================================================

import os
import pathlib

import ez_setup

ez_setup.use_setuptools()

from setuptools import setup  # noqa

# =============================================================================
# CONSTANTS
# =============================================================================

PATH = pathlib.Path(os.path.abspath(os.path.dirname(__file__)))


REQUIREMENTS = ["numpy", "pandas", "scipy", "attrs", "matplotlib", "heliopy"]

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
    name="Solarwindpy",
    version=VERSION,
    description="Package to make some analysis of solar storm phenomena.",
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    author="",
    author_email="",
    url="",
    py_modules=["ez_setup"],
    packages=[
        "solarwindpy",
    ],
    license="The MIT License",
    install_requires=REQUIREMENTS,
    keywords=["solar", "wind", "space weather", "solarwindpy"],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering",
    ],
    include_package_data=True,
)