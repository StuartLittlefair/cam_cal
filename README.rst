HiPERCAM/ULTRACAM flux calibration package
===================================


Installation
------------

The third-party requirements are:

- `numpy <https://numpy.org/>`_,

- `scipy <https://scipy.org/>`_,

- `matplotlib <https://matplotlib.org/>`_,

- `pandas <https://pandas.pydata.org/>`_,

- `astropy <http://astropy.org/>`_, a package for astronomical calculations;

- `mergedeep <https://mergedeep.readthedocs.io/en/latest/>`_, a package for merging dictionaries

- `hipercam <https://github.com/HiPERCAM/hipercam>`_, the HiPERCAM data reduction pipeline;


Installing with pip will handle all standard python packages required but the HiPERCAM pipeline will need a manual install

 pip install cam-calib

or if you don't have root access::

 pip install --prefix=my_own_installation_directory cam-calib

For ease of use, environment variables can be set to tell the package where to look for HiPERCAM or ULTRACAM data::
 
 setenv HCAM_DATA /home/user/hipercam
 setenv UCAM_DATA /home/user/ultracam