HiPERCAM/ULTRACAM flux calibration package
===================================


Installation
------------

The software is written as much as possible to make use of core Python
components. The third-party requirements are:

- `astropy <http://astropy.org/>`_, a package for astronomical calculations;

- `hipercam <https://github.com/HiPERCAM/hipercam>`_, the HiPERCAM data reduction pipeline;


Installing with pip will handle all standard python packages required but the HiPERCAM pipeline will need a manual install

 pip install cam-calib

or if you don't have root access::

 pip install --prefix=my_own_installation_directory cam-calib