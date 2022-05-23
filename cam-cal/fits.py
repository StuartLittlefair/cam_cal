from astropy.io import fits
from astropy.table import Table


def write_FITS(fname, data_dict, header_dict):

    header = fits.Header()
    for key, value in header_dict.items():
        header.set(key, value)
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
        # hdul.append(fits.BinTableHDU(data_dict['data'][filt], header=hdr))

    hdul.writeto(fname, output_verify='warn', overwrite=True)