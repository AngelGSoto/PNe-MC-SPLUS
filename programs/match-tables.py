import numpy as np
from astropy.table import Table, Column, MaskedColumn
import astropy.coordinates as coord

tab = Table.read('luis-programas/arcs-summary-unique.ecsv', 
                    format='ascii.ecsv')
#                    fill_values=('--', np.nan) ).filled(np.nan)
outcols = ['Object']

arctab['coord'] = coord.SkyCoord(ra=arctab['RA'], dec=arctab['Dec'],
                                 unit=('hourangle', 'deg'))

catalogs = {'RSS': {'file': '../../Halpha-DR4splus-xy/iDR4-SPLUS-psf-STAR-14r16-StN5-err02.csv',
                    'ID': 'ID'},
            'COUP': {'file': '../../Halpha-DR4splus-xy/Getman2005/Getman-2005-Vizier.fits',
                     'ID': 'COUP'},
            'MAX': {'file': '../../Halpha-DR4splus-xy/robberto-mir.fits',
                    'ID': 'MAX'},
            'MIR': {'file': '../../Halpha-DR4splus-xy/',
                    'ID': 'recno'},
}

MAXSEP = 1.0
for cat, metadata in catalogs.items():
    cattab = Table.read(metadata['file'])
    cattab['coord'] = coord.SkyCoord(ra=cattab['RAJ2000'], dec=cattab['DEJ2000'],
                                     unit=('deg', 'deg'))
    sourcenames = []
    sourceseps = []
    for arc in arctab:
        seps = arc['coord'].separation(cattab['coord']).arcsec
        iclosest = seps.argmin()
        sourcenames.append('{} {}'.format(cat, cattab[iclosest][metadata['ID']]))
        sourceseps.append(seps[iclosest])

    mask = np.array(sourceseps) > MAXSEP
    outcols.append(cat)
    arctab.add_column(
        MaskedColumn(name=cat, data=sourcenames, mask=mask, meta=metadata,
                     description='Closest source from the {} catalog'.format(cat)))
    outcols.append(cat + ' sep')
    arctab.add_column(
        MaskedColumn(name=cat + ' sep', data=sourceseps, mask=mask, format='{:.3f}'))



outtab = arctab[outcols]
outtab.write('arc-stellar-sources.ecsv', format='ascii.ecsv')
