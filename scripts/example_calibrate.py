from cam_cal.observation import Observation

obs = Observation('ultracam')

obs.add_observation(name='ZTFJ1341_atm', logfiles=['/local/alex/backed_up_on_astro3/Data/photometry/ultracam/2021_01_22/run025_atm.log'], obs_type='atm')
obs.add_observation(name='ZTFJ1404_atm', logfiles=['/local/alex/backed_up_on_astro3/Data/photometry/ultracam/2021_01_22/run027_atm.log'], obs_type='atm')
obs.get_atm_ex()
WD1225_006_mags = dict(mean={'us':15.357, 'gs':14.988, 'rs':15.022, 'is':15.135, 'zs':15.310},
                       err={'us':0.02, 'gs':0.02, 'rs':0.02, 'is':0.02, 'zs':0.02})
obs.add_observation(name='WD1225+006', logfiles=['/local/alex/backed_up_on_astro3/Data/photometry/ultracam/2021_01_22/run029_1.8.log'], obs_type='std', cal_mags=WD1225_006_mags)
obs.get_zeropoint()
obs.add_observation(name='ZTFJ1026', logfiles=['/local/alex/backed_up_on_astro3/Data/photometry/ultracam/2021_01_22/run022_1.8.log'], obs_type='science')
obs.calibrate_science('ZTFJ1026', eclipse=1.5)