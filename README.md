# solarwindpy
<img src="res/logo_SWx.jpg" width="96" height="96" />
A project to rotate an interplanetary magnetic cloud

## Motivation
Magnetic clouds (MCs) are highly magnetized plasma structures that have a low proton
temperature and a magentic field vector that rotates when seen by a heliospheric observer.
More than 25 years of observations of magnetic and plasma properties of MCs at 1 AU
have provided singificant knowledge of their magnetic structure. However, because in situ observations only give information along the trajectory of the spacecraft, their real 3D magnetic configuration remains still partially unknown. Assuming a local cylindrical symmetry it is possible to rotate this structures in the frame of the cloud.

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
