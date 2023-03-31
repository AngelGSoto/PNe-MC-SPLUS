'''
Make photo-spectra from observed SPLUS objects. This program is an updated version of the program: photo-spectra-SPLUSDR2.py.
'''
from __future__ import print_function
import numpy as np
import glob
import json
import matplotlib.pyplot as plt
from astropy.table import Table
#import seaborn as sns
import sys
import argparse
import os
from colour import Color
from pathlib import Path
ROOT_PATH = Path("..")

Number = []

wl = [3485, 3785, 3950, 4100, 4300, 4803, 5150, 6250, 6600, 7660, 8610, 9110]
color = ["#CC00FF", "#9900FF", "#6600FF", "#0000FF", "#009999", "#006600", "#DD8000", "#FF0000", "#CC0066", "#990033", "#660033", "#330034"]
marker = ["s", "o", "o", "o", "o", "s", "o", "s", "o", "s", "o", "s"] ### tienen todos los filtros
facecolor = ['none', '', '', '', '', 'none', '', 'none', '', 'none', '', 'none']

filters = {'uJAVA': ('#CC00FF', 'none', 3485),
           'J0378': ('#9900FF', '', 3785),
           'J0395': ('#6600FF', '', 3940),
           'J0410': ('#0000FF', '', 4100),
           'J0430': ('#009999', '', 4300),
           'gSDSS': ('#006600', 'none', 4803),
           'J0515': ('#DD8000', '', 5150),
           'rSDSS': ('#FF0000', 'none', 6250),
           'J0660': ('#CC0066', '', 6600),
           'iSDSS': ('#990033', 'none', 7660),
           'J0861': ('#330034', '', 9110),
           'zSDSS': ('#660033', 'none', 8610)
           }

wl_b = [3485, 4803,  6250,  7660,  9110]
wl_n = [3785, 3950, 4100, 4300, 5150,  6600, 8610]
color_b = ["#CC00FF", "#006600",  "#FF0000", "#990033",  "#330034"]
color_n = ["#9900FF", "#6600FF", "#0000FF", "#009999",  "#DD8000",  "#CC0066", "#660033"]
marker_b = ["s", "s",  "s", "s", "s"]
marker_n = ["o", "o", "o", "o", "o",  "o", "o"]


parser = argparse.ArgumentParser(
    description="""Write wave and magnitude of a spectrum""")

parser.add_argument("source", type=str,
                    default="known-PN-jplus-idr",
                    help="Name of source, taken the prefix ")

parser.add_argument("--Object", type=str,
                    default=None,
                    help="Id object of the source under interest ")

parser.add_argument("--debug", action="store_true",
                    help="Print out verbose debugging info about each line in region file")

args = parser.parse_args()
file_ = args.source + ".ecsv"

try:
    data = Table.read(ROOT_PATH / file_, format="ascii.ecsv")
except FileNotFoundError:
    file_ = args.source + ".dat"
    data = Table.read(ROOT_PATH / file_, format="ascii")

# If want the spectra of an objetc in particular

idd1 = []
for i in data:
    if i["ID"].endswith(" '"):
        id1 = i["ID"].split("b'")[-1].split(" ")[0]
        idd1.append(id1)
    else:
        id1 = i["ID"].split("b'")[-1].split("'")[0]
        idd1.append(id1)
        
if args.Object is not None:
    Object_ = args.Object
    mask = np.array([source in Object_ for source in idd1])
    data = data[mask]
else:
    data = data
   
# Number of objects
n = len(data)

# Magnitudes list
Number = []
mag_auto  = [[] for _ in range(n)]
mag_petro = [[] for _ in range(n)]
mag_aper_b = [[] for _ in range(n)]
mag_aper_n = [[] for _ in range(n)]
mag_aper = [[] for _ in range(n)]

# Error lists
mag_auto_err  = [[] for _ in range(n)]
mag_petro_err  = [[] for _ in range(n)]
mag_aper_err_b  = [[] for _ in range(n)]
mag_aper_err_n  = [[] for _ in range(n)]
mag_aper_err  = [[] for _ in range(n)]

print("The numbre of objects for plotting is: %d" % n)

for i in range(n):
    mag_aper_b[i].append(data["u_psf"][i]) 
    mag_aper_n[i].append(data["J0378_psf"][i])
    mag_aper_n[i].append(data["J0395_psf"][i])
    mag_aper_n[i].append(data["J0410_psf"][i])
    mag_aper_n[i].append(data["J0430_psf"][i])
    mag_aper_b[i].append(data["g_psf"][i])
    mag_aper_n[i].append(data["J0515_psf"][i]) 
    mag_aper_b[i].append(data["r_psf"][i]) 
    mag_aper_n[i].append(data["J0660_psf"][i])
    mag_aper_b[i].append(data["i_psf"][i]) 
    mag_aper_n[i].append(data["J0861_psf"][i]) 
    mag_aper_b[i].append(data["z_psf"][i])

    #ERRO Aper
    mag_aper_err_b[i].append(float(data["e_u_psf"][i]))
    mag_aper_err_n[i].append(float(data["e_J0378_psf"][i]))
    mag_aper_err_n[i].append(float(data["e_J0395_psf"][i]))
    mag_aper_err_n[i].append(float(data["e_J0410_psf"][i]))
    mag_aper_err_n[i].append(float(data["e_J0430_psf"][i]))
    mag_aper_err_b[i].append(float(data["e_g_psf"][i]))
    mag_aper_err_n[i].append(float(data["e_J0515_psf"][i])) 
    mag_aper_err_b[i].append(float(data["e_r_psf"][i])) 
    mag_aper_err_n[i].append(float(data["e_J0660_psf"][i])) 
    mag_aper_err_b[i].append(float(data["e_i_psf"][i]))
    mag_aper_err_n[i].append(float(data["e_J0861_psf"][i]))
    mag_aper_err_b[i].append(float(data["e_z_psf"][i]))

    # Again for join the point is few smart, I have to think something better
    mag_aper[i].append(data["u_psf"][i]) 
    mag_aper[i].append(data["J0378_psf"][i])
    mag_aper[i].append(data["J0395_psf"][i])
    mag_aper[i].append(data["J0410_psf"][i])
    mag_aper[i].append(data["J0430_psf"][i])
    mag_aper[i].append(data["g_psf"][i])
    mag_aper[i].append(data["J0515_psf"][i]) 
    mag_aper[i].append(data["r_psf"][i]) 
    mag_aper[i].append(data["J0660_psf"][i])
    mag_aper[i].append(data["i_psf"][i]) 
    mag_aper[i].append(data["J0861_psf"][i]) 
    mag_aper[i].append(data["z_psf"][i])

    #ERRO Aper
    mag_aper_err[i].append(float(data["e_u_psf"][i]))
    mag_aper_err[i].append(float(data["e_J0378_psf"][i]))
    mag_aper_err[i].append(float(data["e_J0395_psf"][i]))
    mag_aper_err[i].append(float(data["e_J0410_psf"][i]))
    mag_aper_err[i].append(float(data["e_J0430_psf"][i]))
    mag_aper_err[i].append(float(data["e_g_psf"][i]))
    mag_aper_err[i].append(float(data["e_J0515_psf"][i])) 
    mag_aper_err[i].append(float(data["e_r_psf"][i])) 
    mag_aper_err[i].append(float(data["e_J0660_psf"][i])) 
    mag_aper_err[i].append(float(data["e_i_psf"][i]))
    mag_aper_err[i].append(float(data["e_J0861_psf"][i]))
    mag_aper_err[i].append(float(data["e_z_psf"][i]))
   
    font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }
    ##########################################################################################
    # Plotting -- psf  ######################################################################
    ##########################################################################################
    if data["ID"][i].endswith(" '"):
        if file_.endswith("csv"):
            plotfile = "photopectrum_splus_"+str(data["ID"][i].split("R4.")[-1]).split(" ")[0].replace(".", "-")+"_{}_psf.pdf".format(file_.split('.ecsv')[0])
        else:
            plotfile = "photopectrum_splus_"+str(data["ID"][i].split("R4.")[-1]).split(" ")[0].replace(".", "-")+"_{}_psf.pdf".format(file_.split('.dat')[0])
    else:
        if file_.endswith("csv"):
            plotfile = "photopectrum_splus_"+str(data["ID"][i].split("R4.")[-1]).split("'")[0].replace(".", "-")+"_{}_psf.pdf".format(file_.split('.ecsv')[0])
        else:
            plotfile = "photopectrum_splus_"+str(data["ID"][i].split("R4.")[-1]).split("'")[0].replace(".", "-")+"_{}_psf.pdf".format(file_.split('.dat')[0])
        
    fig = plt.figure(figsize=(15.5, 9.5))
    ax = fig.add_subplot(1,1,1)
    plt.tick_params(axis='x', labelsize=42) 
    plt.tick_params(axis='y', labelsize=42)
    ax.set_xlim(left=3000, right=9700)
    #ax.set_ylim(ymin=17.5,ymax=23)
    #ax1.set_xlabel(r'$\lambda$')
    ax.set_xlabel(r'Wavelength $[\mathrm{\AA]}$', fontsize = 44)
    ax.set_ylabel(r'Magnitude [AB]', fontsize = 44)

    #Masks
    mask = [mag_aper[i][m] != 99.0 for m in range(len(mag_aper[0]))]
    mask_err = [mag_aper_err[i][m] <= 0.9 for m in range(len(mag_aper_err[0]))]
    mask_t = np.array(mask) & np.array(mask_err)

    mask_b = [mag_aper_b[i][m] != 99.0 for m in range(len(mag_aper_b[0]))]
    mask_err_b = [mag_aper_err_b[i][m] <= 0.9 for m in range(len(mag_aper_err_b[0]))]
    mask_t_b = np.array(mask_b) & np.array(mask_err_b)

    mask_n = [mag_aper_n[i][m] != 99.0 for m in range(len(mag_aper_n[0]))]
    mask_err_n = [mag_aper_err_n[i][m] <= 0.9 for m in range(len(mag_aper_err_n[0]))]
    mask_t_n = np.array(mask_n) & np.array(mask_err_n)
    ax.plot(np.array(wl)[mask_err], np.array(mag_aper[i])[mask_err], '-k', alpha=0.2)#, label='Auto')
    for wl1, mag, mag_err, colors, marker_, in zip(np.array(wl_b), np.array(mag_aper_b[i]), np.array(mag_aper_err_b[i]), color_b, marker_b):
        if mag_err <= 0.9:
            ax.scatter(wl1, mag, color = colors, marker=marker_, facecolors="none", s=600, zorder=10)
            ax.errorbar(wl1, mag, yerr=mag_err, fmt='r', marker=None, linestyle=(0, (5, 5)), color=colors, ecolor=colors, elinewidth=5.9, markeredgewidth=5.2,  capsize=20)

    #ax.plot(np.array(wl_n)[mask_err_n], np.array(mag_aper_n[i])[mask_err_n], '-k', alpha=0.2)#, label='Auto')        
    for wl1, mag, mag_err, colors, marker_, in zip(np.array(wl_n), np.array(mag_aper_n[i]), np.array(mag_aper_err_n[i]), color_n, marker_n):
        if mag_err <= 0.9:
            ax.scatter(wl1, mag, color = colors, marker=marker_, s=600, zorder=10)
            ax.errorbar(wl1, mag, yerr=mag_err, fmt='r', marker=None, linestyle=(0, (5, 5)), color=colors, ecolor=colors, elinewidth=5.9, markeredgewidth=5.2,  capsize=20)   
    # plt.text(0.06, 0.1, "Fr 2-21",
    #          transform=ax.transAxes, fontsize=48,  fontdict=font)
    #plt.subplots_adjust(bottom=0.19)
    #plt.legend(fontsize=20.0)
    plt.tight_layout()
    plt.gca().invert_yaxis()
    if args.debug:
        print("Finished plot S-spectra of:", data["ID"][i].split("R3.")[-1].split("\'")[0].replace(".", "-"))
    #save_path = '../../../Dropbox/JPAS/paper-phot/'
    #file_save = os.path.join(save_path, plotfile)
    plt.savefig(plotfile)
    #plt.clf()
    fig.clear()
    plt.close(fig)
