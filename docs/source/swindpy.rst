**Swindpy**
***********

.. figure:: _static/logo_SWx.png
   :alt: alternate text
   :height: 200
   :width: 200
   :scale: 200
   :align: center
   :figclass: align-center

 
.. sidebar:: Abstract
    :subtitle: swindpy

    In this application we use Minimum Variance 
    Technique to obtain the orientation of magnetic
    clouds as long as the correct borders of the 
    cloud are provided.


**Motivation**
==============

Magnetic clouds (MCs) are highly magnetized
plasma structures that have a low proton
temperature and a magnetic field vector that
rotates when seen by a heliospheric
observer.

As motivation MC are the most geo-effective structures in the 
Interplanetary medium. 
Even though MCs have been studied for more than 25 years,
there is no agreement about its true
magnetic configuration.  This is mainly because
the magnetic field data retrieved in
situ by a spacecraft correspond only to the one dimensional
cut along its trajectory and, thus,
it is necessary to make some assumptions to infer
the cloud 3D structure from observations.
Magnetic Clouds have been locally considered as cylindrical symmetric 
structures.
In order to obtain good estimations of these quantities,
it is necessary to find the correct MC orientation, and
improve the estimation of its size and components
of B in the cloud frame.

The Minimum Variance method (MV) has been extensively
used to find the orientation of structures
in the interplanetary medium.

The Minimum Variance method applied to the observed temporal
series of the magnetic field can estimate quite well
the orientation of the cloud axis, when
the distance between the axis and the spacecraft
trajectory in the MC (the impact parameter, p)
is low with respect to the cloud radius.

The MV method has two main advantages
with respect to other more sophisticated techniques that
are also used to find the orientation of an MC:
(1) it is relatively easy to apply,
and (2) it makes a minimum number of assumptions on
the magnetic configuration, only local cylindrical symmetry  
(so it is model independent). 

When p is significant,the MV approach provides orientations
with more error from the real ones, and these errors can be 
estimated as well

**Requeriments**
================

.. code-block:: bash

    * Python 3.9+


**Dependencies for this project**
=================================

    * `NumPy  <https://numpy.org>`_

    * `Pandas <https://pandas.pydata.org/>`_

    * `SciPy  <https://scipy.org/>`_

    * `matplotlib  <https://matplotlib.org/>`_

    * `typer  <https://typer.tiangolo.com/>`_

    * `heliopy==0.15.4  <https://ui.adsabs.harvard.edu/abs/2019zndo...1009079S/abstract>`_

    * `requests  <https://docs.python-requests.org/en/latest/>`_ 

    * `sunpy  <https://sunpy.org/>`_

    * `h5netcdf  <https://anaconda.org/conda-forge/h5netcdf>`_

    * `cdflib  <https://pypi.org/project/cdflib/>`_

      

**Installation**
================

.. code-block:: bash

    $ pip install swindpy


**Tutorial**
============

In this tutorial we are going to learn how to rotate a Magnetic Cloud using 
the swindpy abstractions.

First of all, we need to set the datetimes we are interesting in:


.. code-block:: python

        from datetime import datetime

        import matplotlib.pyplot as plt

        import swindpy.plotter as plotter
        from swindpy.data_manager import Period, MagneticField, DataManager
        from swindpy.rotation import RotatedWind

For this example, we are going to use this dates: 10-Jan-1997 05:00 
to 11-Jan-1997 02:00

.. code-block:: python

        # We set the datetimes we are interesting in
        date_from = datetime(1997, 1, 10, 5, 0, 0)
        date_to = datetime(1997, 1, 11, 2, 0, 0)

        # We create the Period object
        period = Period(date_from, date_to)

        #Using the DataManager, we retrieve the cdf information
        cdf_data = DataManager.get_gse_magnetic_vector(period)

The obtained data, is a list o MagneticField objects (a swindpy abstraction), that has 
information about the datetime and the gse coordinates measures.

.. code-block:: python

        cdf_data[0]

.. code-block:: bash

     MagneticField(time=Timestamp('1997-01-10 05:00:30'), bgse0=-1.8067849, bgse1=9.860464, bgse2=-8.717464)

Now, using the DataManager, we are able to filter nan and infinite values.

.. code-block:: python

        filtered_data = DataManager.filter_nan_and_inf_values(cdf_data)

Now we can obtain the RotatedWind simply calling a classmethod

.. code-block:: python

        rotated_wind = RotatedWind.get_rotated_wind(filtered_data)

Using the plotter method, we are able to plot non rotated winds and rotated winds (we are also able to add labels, 
change size, etc).   

.. code-block:: python

        # Plotting non rotated
        plotter.plot_mf(filtered_data)

.. figure:: _static/imagen1.png
   :alt: alternate text
   :height: 100
   :width: 200
   :scale: 200
   :align: center
   :figclass: align-center



.. code-block:: python

        # Plotting rotated
        plotter.plot_rw(rotated_wind)

.. figure:: _static/imagen2.png
   :alt: alternate text
   :height: 100
   :width: 200
   :scale: 200
   :align: center
   :figclass: align-center

.. code-block:: python

        # Obtain the rotation angles
        theta, phi = get_rotation_angles(filtered_data)

        # Calculate gamma using calc_gamma
        gamma = Angle(
                "gamma",
                calc_gamma(theta.angle, phi.angle)
        )

        print(theta, phi, gamma)

.. code-block:: bash

        theta
         RAD: -0.31457481937133086
         DEG: -18.02380949106747
        phi
         RAD: 1.6979828788548021
         DEG: 97.28725264385352
        gamma
         RAD: 0.31696887828176834
         DEG: 18.160978962541225


Using swindpy command line interface

We also created a CLI that makes it easier to a user to process magnetic clouds data.

If you want to plot a no rotated cloud, you could do that with the next command:

.. code-block:: bash

        swindpy plot-cloud 2021-01-01 2021-01-02

If you want to plot a rotated cloud, you could do that with the next command:

.. code-block:: bash

        swindpy plot-rotated-cloud 2021-01-01 2021-01-02

To export data about the Magnetic Fields in a period time, use the next command:

.. code-block:: bash

        swindpy to-csv 2021-01-01 2021-01-02 output

The previous command will generate a csv with the period data

You are also able to plot both, no rotated and rotated clouds so you can compare and 
make a quick analysis of the results:

.. code-block:: bash

        swindpy plot-rotated-and-non-rotated 2021-01-01 2021-01-02

**Tools that we use so far**
============================

======================== =========================
**Tools**                **Detail**
------------------------ -------------------------
Python                   Language
Pylint                   Syntax
Pypi                     Publish the library
Sphinx                   Keep record
Github                   Share the code
readthedocs.org          Publish documentation
======================== =========================

**Contact**
===========

You can contact us via email, agulisano@iafe.uba.ar

**Issues**
==========

Please submit bug reports, suggestions for improvements and patches via the issue tracker.

**Links**
=========

    * `Documentation <https://swindpy.readthedocs.io/en/latest/>`_
    * `PyPl Releases <https://pypi.org/project/swindpy/>`_
    * `GitHub <https://github.com/adelarja/space_weather/>`_

**Credits**
===========

We propose to use the open source software Solarwindpy using for the 
calculation the Minimum Variance Technique to obtain the orientation 
of magnetic clouds provided the correct cloud edges are provided.

**License**
===========

 The four essential freedoms 

 A program is free software if users have all four essential freedoms:

The freedom to run the program as desired, for any purpose (freedom 0).
The freedom to study how the program works, and change it to do what you want
(freedom 1). Access to the source code is a necessary condition for this.
The freedom to redistribute copies to help others (freedom 2).
The freedom to distribute copies of your modified versions to third parties
(freedom 3). This allows you to offer the entire community the opportunity
to benefit from the changes. Access to the source code is a necessary condition
for this.

 The two main categories of free software licenses are copyleft and non-copyleft. 
 Copyleft licenses, such as the GNU GPL, insist that modified versions of a free 
 program must also be free software. Non-copyleft licenses do not engage in this.

 -BSD License (Berkeley Software Distribution):

 It is a permissive free software license. In other words, it is in contrast to 
 copyleft licenses, which have share-alike reciprocity requirements. The BSD license 
 allows the use of the source code in non-free software. The original version has 
 already been revised and its variants are called modified BSD licenses.

 Copyright <year> <copyright holder>

 Redistribution and use in source and binary forms, with or without modification, are 
 permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list
of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or other
materials provided with the distribution.
Neither the name of the copyright holder nor the names of its contributors may be
used to endorse or promote products derived from this software without specific
prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT 
SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR 
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
