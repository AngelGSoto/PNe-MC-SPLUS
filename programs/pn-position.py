import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from  astropy.table import Table, vstack, hstack
import seaborn as sns
from astropy.coordinates import SkyCoord
from astropy import units as u

#Files
df_lmc = pd.read_csv('TAP_3_J_A+A_456_451_PNe.csv')
df_smc = pd.read_csv('smc-final.csv')
##Tiles
df_tiles = pd.read_csv('Tiles-DR4/iDR4_pointings.csv')

# Selecting only tiles of MC
df_tiles_mc = df_tiles[df_tiles["Field"].str.contains("MC")]

#LMC                
ra_lmc = df_lmc["_RA"]
dec_lmc = df_lmc["_DE"]
icrs_lmc = SkyCoord(ra=ra_lmc*u.degree, dec=dec_lmc*u.degree, frame='icrs')
gal_lmc = icrs_lmc.galactic
l_rad_lmc = gal_lmc.l.radian
l_rad_lmc[l_rad_lmc > np.pi] -= 2. * np.pi
b_rad_lmc = gal_lmc.b.radian
b_deg_lmc = b_rad_lmc * (180/np.pi)#gal.b.degree#b_rad * (180/np.pi)
l_deg_lmc = l_rad_lmc * (180/np.pi)

#SMC
ra_smc = df_smc["RAJ2000"]
dec_smc = df_smc["DEJ2000"]
icrs_smc = SkyCoord(ra=ra_smc*u.degree, dec=dec_smc*u.degree, frame='icrs')
gal_smc = icrs_smc.galactic
l_rad_smc = gal_smc.l.radian
l_rad_smc[l_rad_smc > np.pi] -= 2. * np.pi
b_rad_smc = gal_smc.b.radian
b_deg_smc = b_rad_smc * (180/np.pi)#gal.b.degree#b_rad * (180/np.pi)
l_deg_smc = l_rad_smc * (180/np.pi)

# TILES
ra_tile = df_tiles_mc["RA"]
dec_tile = df_tiles_mc["DEC"]
icrs_tile = SkyCoord(ra=ra_tile*u.degree, dec=dec_tile*u.degree, frame='icrs')
gal_tile = icrs_tile.galactic
l_rad_tile = gal_tile.l.radian
l_rad_tile[l_rad_tile > np.pi] -= 2. * np.pi
b_rad_tile = gal_tile.b.radian
b_deg_tile = b_rad_tile * (180/np.pi)#gal.b.degree#b_rad * (180/np.pi)
l_deg_tile = l_rad_tile * (180/np.pi)


#Plotting
lgd_kws = {'frameon': True, 'fancybox': True, 'shadow': True}
sns.set_style('ticks')
fig, ax = plt.subplots(figsize=(12, 6))
# fig = plt.figure(figsize=(14,7))
# ax = fig.add_subplot(1,1,1, projection='aitoff')
plt.tick_params(axis='x', labelsize=30) 
plt.tick_params(axis='y', labelsize=30)
plt.xlabel(r'$l$(Galactic longitude)[deg]', fontsize=28)
plt.ylabel(r'$b$(Galactic latitude)[deg]', fontsize=28)
# ax.set_xlim(-1.2, 5.7)
# ax.set_ylim(-3.3, 9.2)
#fig = plt.figure(figsize=(10, 8))
#ax = fig.add_subplot(111)
ax.scatter(l_deg_lmc, b_deg_lmc, marker = "o", s = 80, c = sns.xkcd_palette(["forest green"]), edgecolors= "g", zorder = 3, lw=1, alpha = 0.7, label="PNe")
ax.scatter(l_deg_smc, b_deg_smc, marker = "o", s = 80, c = sns.xkcd_palette(["forest green"]), edgecolors= "g", zorder = 4, lw=1, alpha = 0.7)#, label="No confirmed HASH PNe")
ax.scatter(l_deg_tile, b_deg_tile, marker = "s", s = 800, c = sns.xkcd_palette(["grey"]), edgecolors= "none", zorder = 1, lw=1, alpha = 0.3)
#ax.scatter(l_, b_, s = 100, c = sns.xkcd_palette(["purple"]), edgecolors= "purple", zorder = 2, lw=2)
#ax.scatter(l_deg, b_deg, s = 150, c = sns.xkcd_palette(["azure"]), edgecolors= "b", zorder = 3, lw=2, label="New PN")

bbox_props = dict(boxstyle="round", fc="w", ec="0.78", alpha=0.001, pad=0.1)
plt.text(0.1, 0.8, 'LMC',
         transform=ax.transAxes, c="gray", weight='bold', fontsize=32.8, zorder = 1, bbox=bbox_props)

plt.text(0.85, 0.05, 'SMC',
         transform=ax.transAxes, c="gray", weight='bold', fontsize=32.8, zorder = 2, bbox=bbox_props)

ax.grid(True, linestyle='-.', linewidth=0.7)
ax.legend(prop={'family': 'monospace', 'size': 15}, **lgd_kws)
#plt.gca().invert_yaxis()
plt.tight_layout()
fig.savefig("Text/Figs/galactic-coord-pn.pdf")
plt.clf()

