import numpy as np
import os
from astropy.io import fits
from astropy.table import Table


class FITS:

    def __init__(self, hdul):
        self.hdul = hdul

    
    @classmethod
    def open(cls, fname):
        hdul = fits.open(fname)
        cls.fname = fname
        return cls(hdul)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.hdul.close()

    
    @classmethod
    def create(cls, data_dict, header_dict):
        header = fits.Header()
        for key, value in header_dict.items():
            header.set(key, value)
        header.comments['TIME_CAL'] = "Time data was calibrated (UTC)."
        header.comments['TIME_ID'] = "Unique ID of calibration."
        primary_hdu = fits.PrimaryHDU(header=header)
        hdul = fits.HDUList([primary_hdu])
        for filt in data_dict['header'].keys():
            hdr = fits.Header()
            for key, value in data_dict['header'][filt].items():
                hdr.set(key, value)
            BMJD = fits.Column(name='BMJD(TDB)', format='D', array=data_dict['data'][filt]['BMJD(TDB)'], unit='d')
            Exp_time = fits.Column(name='Exp_time', format='D', array=data_dict['data'][filt]['exp_time'], unit='d')
            Flux = fits.Column(name='Flux', format='D', array=data_dict['data'][filt]['flux'], unit='Jansky')
            Flux_err = fits.Column(name='Flux_err', format='D', array=data_dict['data'][filt]['flux_err'], unit='Jansky')
            weight = fits.Column(name='Weight', format='D', array=data_dict['data'][filt]['weight'])
            esubd = fits.Column(name='Esubd', format='D', array=data_dict['data'][filt]['esubd'])
            coldefs = fits.ColDefs([BMJD, Exp_time, Flux, Flux_err, weight, esubd])
            hdul.append(fits.BinTableHDU.from_columns(coldefs, header=hdr))
        return cls(hdul)


    def write(self, fname):
        self.hdul.writeto(fname, output_verify='warn', overwrite=True)

    
    def to_lcurve(self, filename=None, filepath=None):
        header = self.hdul[0].header
        for i, hdu in enumerate(self.hdul[1:]):
            hdr = hdu.header
            t = hdu.data['BMJD(TDB)']
            t_exp = hdu.data['Exp_time']
            flux = hdu.data['Flux']
            flux_err = hdu.data['Flux_err']
            weights = hdu.data['Weight']
            esubd = hdu.data['Esubd']

            array = np.column_stack((t, t_exp, flux, flux_err, weights, esubd))

            fname = filename
            fpath = filepath
            if not fname:
                fname = f"{header['TARGET']}_{header['RUN']}_{header['TIME_ID']}_{hdr['FILTER']}_fc.dat"
                if not fpath:
                    try:
                        fpath = os.path.dirname(self.fname)
                    except Exception as e:
                        fpath = ""
                fname = os.path.join(fpath, fname)
            
            np.savetxt(fname, array)