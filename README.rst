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


Installing with pip will handle all standard python packages required but the HiPERCAM pipeline will need a manual install.

 cd cam_cal
 pip install .

or if you don't have root access::

 pip install --prefix=my_own_installation_directory cam-calib

For ease of use, environment variables can be set to tell the package where to look for HiPERCAM or ULTRACAM data::
 
 setenv HCAM_DATA /home/user/hipercam
 setenv UCAM_DATA /home/user/ultracam


Use
---

First import the ``Observation`` class and create an instance.

.. code-block:: python

    from cam_cal.observation import Observation

    obs = Observation('ultracam') # or 'hipercam'

Observations or imaging runs are then added with the ``add_observation`` method.
These are given a name for reference, a list of hipercam logfiles, and and obs_type.
An observation can have one of three obs_types: 'atm', 'std', or 'science'.

**atm** - The logfile from a reduction of a selection of stable stars over a decent airmass range.
          More than one logfile can be added for the one observation in case a target is observed more than once in a night (increasing the airmass range covered),
          However, aperture numbering must correspond exactly between both observations so stars can be matched across both runs.
**std** - The logfile from the reduction of a flux standard.
**science** - The logfile of the science data. The target must be aperture 1. It is wise to match the aperture radii with that of the flux standard in order to prevent systematic effects. 

First add a run over a decent airmass range and fit the atmospheric extinction.

.. code-block:: python

    obs.add_observation('SDSSJ0931_atm', logfiles=['2021_01_22/run025_atm.log'], obs_type='atm')
    obs.get_atm_ex(plot=True)

Second, add a run on a flux standard. The calibrated magnitudes in the relevant filter system are required.
These can be given either as a dictionary, or, if the standard is in the list of Gaia spectrophotometric standards,
then either the name or an astropy Skycoord instance can be given.
With a run added, the zeropoint can be calculated.

.. code-block:: python


    obs.add_observation(name='GD108', logfiles=['2021_01_22/run031.log'], obs_type='std', cal_mags='GD108')

    WD1225_006_mags = dict(mean={'us':15.357, 'gs':14.988, 'rs':15.022, 'is':15.135, 'zs':15.310},
                           err={'us':0.02, 'gs':0.02, 'rs':0.02, 'is':0.02, 'zs':0.02})
    obs.add_observation(name='WD1225+006', logfiles=['2021_01_22/run029.log'], obs_type='std', cal_mags=WD1225_006_mags)

    obs.get_zeropoint()

Finally, our science data can be added and calibrated. If the target is a detached eclipsing binary then the data
centred around the eclipse can be extracted and will automatically increase the weighting of the ingress/egress to constitute
an equal portion of the total light curve. This is still experimental though so be careful.

.. code-block:: python

    obs.add_observation(name='SDSSJ1028', logfiles=['2021_01_22/run022.log'], obs_type='science')
    obs.calibrate_science('SDSSJ1028', eclipse=1.5)
    # eclipse=1.5 extracts 1.5x the eclipse width either side of the eclipse midpoint
    # i.e. the eclipse with an eclipse width's worth of out-of-eclipse data either side.
