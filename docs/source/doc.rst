**Swindpy**
***************

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

As motivation MC are the most geo-effective structures in the Interplanetary 
medium. 
Even though MCs have been studied for more than 25 years,
there is no agreement about its true
magnetic configuration.  This is mainly because
the magnetic field data retrieved in
situ by a spacecraft correspond only to the one dimensional
cut along its trajectory and, thus,
it is necessary to make some assumptions to infer
the cloud 3D structure from observations.
Magnetic Clouds have been locally considered as cylindrical symmetric structures.
In order to obtain good estimations of these quantities,
it is necessary to find the correct MC orientation, and
improve the estimation of its size and components
of $\bf B$ in the cloud frame.

The Minimum Variance method (MV) has been extensively
used to find the orientation of structures
in the interplanetary medium.

The Minimum Variance method applied to the observed temporal
series of the magnetic field can estimate quite well
the orientation of the cloud axis, when
the distance between the axis and the spacecraft
trajectory in the MC (the impact parameter, $p$)
is low with respect to the cloud radius.

The MV method has two main advantages
with respect to other more sophisticated techniques that
are also used to find the orientation of an MC:
(1) it is relatively easy to apply,
and (2) it makes a minimum number of assumptions on
the magnetic configuration, only local cylindrical symmetry  
(so it is model independent). 

When $p$ is significant,the MV approach provides orientations
with more error from the real ones, and these errors can be 
estimated as well

**Requeriments**
================

.. code-block:: bash

    * Python 3.10+

**Dependencies for this project**
=================================

    * `NumPy  <https://numpy.org>`_

    * `Pandas <https://pandas.pydata.org/>`_

    * `SciPy  <https://scipy.org/>`_

    * `matplotlib  <https://matplotlib.org/>`_

    * `attrs  <https://attrs.org/en/stable/>`_

    * `heliopy==0.15.4  <https://ui.adsabs.harvard.edu/abs/2019zndo...1009079S/abstract>`_

    * `requests  <https://docs.python-requests.org/en/latest/>`_ 

    * `sunpy  <https://sunpy.org/>`_

    * `h5netcdf  <https://anaconda.org/conda-forge/h5netcdf>`_

    * `cdflib  <https://pypi.org/project/cdflib/>`_

      

**Instalation**
===============

.. code-block:: bash

    $ pip install libreriaswindpy==0.0.1


**Function**
=============

- **multiplicacion()**
    Una función que espera dos números y retorna el resultado

**Algoritmos**
==============




**Que usamos hasta ahora**
==========================

======================== =========================
**Herramientas**         **Detalle**
------------------------ -------------------------
Python                   Lenguaje
Pylint                   Sintaxis
Pypi                     Publicar la libreria
Sphinx                   Documentar
Github                   Compartir el codigo
readthedocs.org          Publicar la documentacion
======================== =========================

**Indexing**
============

    * `Astroalign <https://arxiv.org/abs/1909.02946/>`_ 
    * `Feets <https://arxiv.org/abs/1809.02154/>`_ 
    * `Corral <https://arxiv.org/abs/1701.05566/>`_
    * `Grispy <https://arxiv.org/abs/1912.09585/>`_

**Contact**
===========

you can contact us via email...

**Issues**
==========

Please submit bug reports, suggestions for improevements and patches via the issue tracker.

**Links**
=========

Documentation
Example Application
PyPl Releases
Changelog

**Credits**
===========

We propose using the open source software Solarwindpy for the calculation.........

**License**
===========

This project is licensed under the MIT License (see the LICENSE file for details)


