'''
This script makes the SDSS spectra overlapped on SPLUS photometry
'''
import astropy.coordinates as coord
import astropy.units as u
from astropy.io import ascii
from astropy.table import Table, hstack
from astropy.coordinates import SkyCoord 
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import seaborn as sn
sn.set_context("poster")
import glob
import argparse
import sys
import os
from astropy.visualization import hist
from astroML.datasets import fetch_imaging_sample, fetch_sdss_S82standards
from astroML.crossmatch import crossmatch_angular

# Read the files
parser = argparse.ArgumentParser(
    description="""Make a spectras""")

parser.add_argument("fileVlt", type=str,
                    default="teste-program",
                    help="Name of file, taken the prefix")

parser.add_argument("TableSplus", type=str,
                    default="teste-program",
                    help="Name of table, taken the prefix")

parser.add_argument("--ymin", required=False, type=float, default=None,
                    help="""Value y-axis min""")

parser.add_argument("--ymax", required=False, type=float, default=None,
                    help="""Value y-axis max""")

cmd_args = parser.parse_args()
file_vlt = cmd_args.fileVlt + ".dat"
file_table = cmd_args.TableSplus + ".ecsv"

#hdu = fits.open(file_spec)
datadir = "../"
datadir1 = "1D/"

table_vlt = Table.read(os.path.join(datadir, file_vlt), format="ascii")
# Coordinates of the VLT spectra
table_vlt['coord'] = coord.SkyCoord(ra=table_vlt['RA'], dec=table_vlt['DEC'],
				    unit=('deg', 'deg'))

# Splus table
table_splus = Table.read(os.path.join(datadir, file_table), format="ascii.ecsv")
table_splus['coord'] = coord.SkyCoord(ra=table_splus['RA_r'], dec=table_splus['DEC_r'],
				    unit=('deg', 'deg'))

# Find the VLT object on the SPLUS list, makes crossmacth using 1 arcsec

MAXSEP = 0.3
# Definition
def match(tab, tab_base):
    sourceseps = []
    for pn in tab:
        seps = pn['coord'].separation(tab_base['coord']).arcsec
        iclosest = seps.argmin()
        sourceseps.append(seps[iclosest])
    return sourceseps  

#mask
mask = np.array(match(table_splus, table_vlt)) <= MAXSEP
tab_splus_match = table_splus[mask]
# print("******************************************************")
# print("Objects VLT source:", tab_splus_match)
# print("******************************************************")
tab_splus_match.remove_column('coord')
tab_splus_match.sort('RA_r')
table_vlt.remove_column('coord')
table_vlt.sort('RA')
table_ = hstack([tab_splus_match, table_vlt])
# Identify rows with missing data and create a mask
missing_data_mask = table_['RA'].mask | table_['VLT'].mask

# Remove rows with missing data from the table
table_mask = table_[~missing_data_mask]
print(table_mask)
table_mask.write('../VLT-XSHOOTER-spectra-splusmag.ecsv', format='ascii.ecsv', overwrite=True)
# sys.exit()
WL, FLUX, fits_ = [], [], []
mag, mag_err = [], []
# Data from the VLT spectra and SPLUS
for table_all in table_mask:
    hdulist= fits.open(os.path.join(datadir1, table_all['VLT']))
    data_spc = hdulist[1].data
    wl = data_spc["WAVE"]
    Flux = data_spc["FLUX"]
    WL.append(wl)
    FLUX.append(Flux)
    fits_.append(table_all['VLT'])
    
    # Now the magnitudes
    mag.append(table_all["u_psf"])
    mag.append(table_all["J0378_psf"])
    mag.append(table_all["J0395_psf"])
    mag.append(table_all["J0410_psf"])
    mag.append(table_all["J0430_psf"])
    mag.append(table_all["g_psf"])
    mag.append(table_all["J0515_psf"])
    mag.append(table_all["r_psf"])
    mag.append(table_all["J0660_psf"])
    mag.append(table_all["i_psf"])
    mag.append(table_all["J0861_psf"])
    mag.append(table_all["z_psf"])
    #ERRO psf                                     
    mag_err.append(table_all["e_u_psf"])
    mag_err.append(table_all["e_J0378_psf"])
    mag_err.append(table_all["e_J0395_psf"])
    mag_err.append(table_all["e_J0410_psf"])
    mag_err.append(table_all["e_J0430_psf"])
    mag_err.append(table_all["e_g_psf"])
    mag_err.append(table_all["e_J0515_psf"])
    mag_err.append(table_all["e_r_psf"])
    mag_err.append(table_all["e_J0660_psf"])
    mag_err.append(table_all["e_i_psf"])
    mag_err.append(table_all["e_J0861_psf"])
    mag_err.append(table_all["e_z_psf"])

# Convert in a array mag and err_mag
mag_arr = np.array(mag).reshape(8, 12)
mag_err_arr = np.array(mag_err).reshape(8, 12)

# Data of the SPLUs list
wl_sp = [3485, 3785, 3950, 4100, 4300, 4803, 5150, 6250, 6600, 7660, 8610, 9110]
# wl_sp = np.array(wl_sp).reshape(1, 12)
# Repeat the row to create a 2D array with multiple identical rows
# Create a 2D array with repeated rows
repeated_wl_sp = np.repeat([np.array(wl_sp)], 8, axis=0)

color = ["#CC00FF", "#9900FF", "#6600FF", "#0000FF", "#009999", "#006600", "#DD8000", "#FF0000", "#CC0066", "#990033", "#660033", "#330034"]
repeated_color =  np.repeat([color], 8, axis=0)
# repeated_color_ = repeated_color.reshape(8, 12)
marker = ["s", "o", "o", "o", "o", "s", "o", "s", "o", "s", "o", "s"] ### tienen todos los filtros
repeated_marker = np.repeat([marker], 8, axis=0)
# repeated_marker_ = repeated_marker.reshape(8, 12)

# sys.exit()
# ff = (10**(-(table["R_psf"][ind] + 2.41) / 2.5)) / 6250.0**2
# print(ff)
# for i, ii in zip(wl, Flux):
#     if i> 6000 and i< 6300:
#          print(i, ii)

# Find scale factor
# m = wl == 6250.289 
# wl_part = wl[m]
# flux_part = Flux[m]
# Fsp = (10**(-(table["r_psf"][ind] + 2.41) / 2.5)) / 6250.0**2
# factor = flux_part / Fsp

# Propagation of error
err_ = []
for wll, magg, magerr in zip(repeated_wl_sp, mag_arr, mag_err_arr):
    c = (10**(-2.41/2.5)) / wll**2
    b = -(1. / 2.5)
    err = np.sqrt(((c*10**(b*magg))**2)*(np.log(10)*b*magerr)**2)
    err_.append(err)

#axis limit
# mask_lim = (wl > 6100.) & (wl < 6900.)
# Flux_lim = Flux[mask_lim]
# if max(Flux_lim) > 5 * np.mean(Flux_lim):
#     max_y_lim = max(Flux_lim) * .9
#     #plt.ylim(ymax=max_y_lim)
#     # min_y_lim = min(Flux_lim) - 0.2
#     # plt.ylim(ymin=min_y_lim,ymax=max_y_lim)

# # set Y-axis range (if applicable)
# if cmd_args.ymin is not None and cmd_args.ymax is not None:
#     plt.ylim(cmd_args.ymin,cmd_args.ymax)
# elif cmd_args.ymin is not None:
#     plt.ylim(ymin=cmd_args.ymin)
# elif cmd_args.ymax is not None:
#     plt.ylim(ymax=cmd_args.ymax)

#plt.ylim(ymin=-50.0,ymax=200)

# if max(Flux) > 9 * np.mean(Flux):
#     max_y_lim = max(Flux) * .9
#     min_y_lim = min(Flux) 
#     plt.ylim(ymin=min_y_lim,ymax=max_y_lim)
# ax.plot(wl, Flux, c = "gray", linewidth=1.3, alpha=0.6, zorder=5)
for fitname, wl_, FLUX_, wl_filter, mag_arr_, mag_err_arr_, repeated_color_, repeated_marker_ in zip(fits_, WL, FLUX, repeated_wl_sp, mag_arr, err_, repeated_color, repeated_marker):
    fig, ax = plt.subplots(figsize=(12, 9))
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)
    ax.set(xlim=[3350,9300])
    ax.set(ylim=[-0.01e-12, 0.07e-12])
    #plt.ylim(ymin=-50.0,ymax=200)
    ax.set(xlabel='Wavelength $(\AA)$')
    ax.set(ylabel=r'F$(\mathrm{erg\ s^{-1} cm^{-2} \AA^{-1}})$')

    ax.plot(wl_, FLUX_, c = "gray", linewidth=2.3, alpha=0.6, zorder=5)
    F = (10**(-(mag_arr_ + 2.41) / 2.5)) / wl_filter**2
    # Crea un nuevo gráfico
        #print(wl_filter, mag_err_arr_)
        # #F *= factor
    for i in range(len(mag_arr_)):
        ax.scatter(wl_filter[i], F[i], c = repeated_color_[i], s=250, marker = repeated_marker_[i], zorder=4)#[]
        ax.errorbar(wl_filter[i], F[i], yerr=mag_err_arr_[i], marker='.', fmt='.', color=repeated_color_[i], ecolor=repeated_color_[i], elinewidth=3.9, markeredgewidth=3.2, capsize=10)
    # ax.legend()
    plt.tight_layout()
    figpdf = fitname.replace(".fits", ".pdf")
    plt.savefig(figpdf)
plt.close()
        

#ax.axvline(4686, color='r', linewidth=0.3, linestyle='-', zorder = 6, label="He II")
# plt.text(0.70, 0.19, table["ID"].split("R3.")[-1]).replace(".", "-"),
#              transform=ax.transAxes, fontsize=25, weight='bold')

# if cmd_args.ymax is not None:
#     ax.annotate(str(table["ID"][ind]).split("R3.")[-1].replace(".", "-"), xy=(9000, cmd_args.ymax),  xycoords='data', size=13,
#             xytext=(-120, -60), textcoords='offset points', 
#             bbox=dict(boxstyle="round4,pad=.5", fc="0.94"),)
#     ax.annotate("r=" + format(float(table["R_psf"][ind]), '.2f'), xy=(9000, 0.86*cmd_args.ymax),  xycoords='data', size=13,
#             xytext=(-120, -60), textcoords='offset points', 
#             bbox=dict(boxstyle="round4,pad=.5", fc="0.94"),)
# else:
#     None
    
# ax.legend()
# plt.tight_layout()
# asciifile = file_spec.replace(".fits", 
#                   "-"+(str(table["ID"][ind]).split("R3.")[-1]).replace(".", "-")+".pdf")
# plt.savefig(asciifile)
