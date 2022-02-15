# solarwindpy
<p align="center">
    <img width=200 src="https://raw.githubusercontent.com/adelarja/space_weather/main/res/logo.jpg">
</p>

[![Build status](https://github.com/adelarja/space_weather/actions/workflows/solarwindpy_ci.yml/badge.svg)](https://github.com/adelarja/space_weather/actions)
[![PyPI](https://img.shields.io/pypi/v/swindpy?color=blue)](https://pypi.org/project/swindpy/)

## Motivation
Magnetic clouds (MCs) are highly magnetized plasma structures that have a low proton temperature and a magnetic fieldvector that rotates when seen by a heliospheric observer.

MC are the most geo-effective structures in the Interplanetary medium. Even though MCs have been studied formore than 25 years, there is no agreement about its true magnetic configuration. This is mainly because the magnetic field data retrieved in situ by a spacecraft correspond only to the one-dimensional cut along its trajectory and, thus, it is necessary to make some assumptions to infer the cloud 3D structure from observations.

MCs have been locally considered ascylindrical symmetric structures. From in situ observations and assumptions on the magnetic distribution inside the MC, it ispossible to estimate some global magnetohydrodynamic quantities, such as magnetic fluxes and helicity for instance. In order to obtain good estimations of these quantities, it is necessary tofind the correct MC orientation, and improve the estimation ofits size and components of the cloud frame.

The Minimum Variance method (MV) has been extensively used to find the orientation of structures in the interplanetary medium.   The  Minimum  Variance  method  applied  to  the  observed temporal series of the magnetic field can estimate quitewell the  orientation  of  the  cloud  axis,  when  the  distance  between  the  axis  and  the  spacecraft  trajectory  in  the  MC  (the impact parameter, p) is low with respect to the cloud radius. The MV method has two main advantages with respect to other more sophisticated techniques that are also used to find the orientation  of  an  MC:
- 1 it is relatively easy to apply
- 2 it makes a minimum number of assumptions on the magnetic configuration, only local cylindrical symmetry (so it is model independent).

In this project the MV approach is used to rotate MCs. In the future, another kind of analysis will be available in the project.

## Features
This package aims tu be used to analyze MCs phenomena.

Currently, it uses Minimum Variance method in order to obtain a rotated cloud.

With help of a command line tool it is also possible to export magnetic cloud data into a
csv file.

--------------------------------------------------------------------------------

## Requirements
You need a Python 3.9+ environment to run solarwindpy.

## Installation
The package is available in pypi. You can install it using pip:

        $ pip install swindpy
        
## Test
If you want to run the swindpy tests, you have to clone the repository and use the
pytest module.

        $ pytest tests
