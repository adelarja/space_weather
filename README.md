# solarwindpy
<img src="res/logo_SWx.jpg" width="96" height="96" />
A project to rotate a solar magnetic cloud

## Motivation
COMPLETAR

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
