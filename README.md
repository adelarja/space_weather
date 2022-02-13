# solarwindpy
<p align="center">
    <img width=200 src="https://github.com/adelarja/space_weather/blob/main/res/logo.jpg">
</p>

[![Build status](https://github.com/adelarja/space_weather/actions/workflows/solarwindpy_ci.yml/badge.svg)](https://github.com/adelarja/space_weather/actions)
[![PyPI](https://img.shields.io/pypi/v/swindpy?color=blue)](https://pypi.org/project/swindpy/)

## Motivation
Magnetic clouds (MCs) are highly magnetized plasma struc-tures that have a low proton temperature and a magnetic fieldvector that rotates when seen by a heliospheric observer.   Asmotivation MC are the most geo-effective structures in the In-terplanetary medium.  Even though MCs have been studied formore than 25 years, there is no agreement about its true mag-netic configuration.  This is mainly because the magnetic fielddata retrieved in situ by a spacecraft correspond only to the onedimensional cut along its trajectory and, thus, it is necessary tomake some assumptions to infer the cloud 3D structure fromobservations. Magnetic Clouds have been locally considered ascylindrical symmetric structures. From in situ observations andassumptions on the magnetic distribution inside the MC, it ispossible to estimate some global magnetohydrodynamic quan-tities, such as magnetic fluxes and helicity for instance. In orderto obtain good estimations of these quantities, it is necessary tofind the correct MC orientation, and improve the estimation ofits size and components ofBin the cloud frame.The Minimum Variance method (MV) has been extensivelyused to find the orientation of structures in the interplanetarymedium.   The  Minimum  Variance  method  applied  to  the  ob-served temporal series of the magnetic field can estimate quitewell  the  orientation  of  the  cloud  axis,  when  the  distance  be-tween  the  axis  and  the  spacecraft  trajectory  in  the  MC  (theimpact parameter,p) is low with respect to the cloud radius.The MV method has two main advantages with respect to othermore sophisticated techniques that are also used to find the ori-entation  of  an  MC:  (1)  it  is  relatively  easy  to  apply,  and  (2)it makes a minimum number of assumptions on the magneticconfiguration, only local cylindrical symmetry (so it is modelindependent). Whenpis significant, the MV approach providesorientations with more error from the real ones, and these errorscan be estimated as well .In Section 1, we describe briefly theMinimum Variance method applied to magnetic clouds underthe hypothesis of locally cylindrical symmetric configurations



## Features
The input file must be a .cdf file with three columns: Bx, By and Bz (Components of Solar Magnetic Cloud in GSE.)


--------------------------------------------------------------------------------

## Requirements
You need Python 3.9+ to run solarwindpy.

## Installation
Clone this repo and then inside the local directory execute

        $ pip install -e .
        
### Test
Testing solarwindpy.wind

        $ cd tests
        $ python test_wind.py 
        
 If the assert is satisfied then is printed 'This is the espected result'. If the assert is false, then is printed AssertionError.
