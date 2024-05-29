import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
import seaborn as sns

# Files
df_lmc = pd.read_csv('TAP_3_J_MNRAS_438_2642_table1.csv')
df_smc = pd.read_csv('smc-final.csv')
df_tiles = pd.read_csv('Tiles-DR4/iDR4_pointings.csv')

# Selecting only tiles of MC
df_tiles_mc = df_tiles[df_tiles["Field"].str.contains("MC")]

# LMC
ra_lmc = df_lmc["_RA"]
dec_lmc = df_lmc["_DE"]
icrs_lmc = SkyCoord(ra=ra_lmc.values*u.degree, dec=dec_lmc.values*u.degree, frame='icrs')

# SMC
ra_smc = df_smc["RAJ2000"]
dec_smc = df_smc["DEJ2000"]
icrs_smc = SkyCoord(ra=ra_smc.values*u.degree, dec=dec_smc.values*u.degree, frame='icrs')

# Tiles
ra_tile = df_tiles_mc["RA"]
dec_tile = df_tiles_mc["DEC"]
icrs_tile = SkyCoord(ra=ra_tile.values*u.degree, dec=dec_tile.values*u.degree, frame='icrs')

# Calculate the number of tiles and size of each tile
num_tiles = len(icrs_tile)
total_fov_deg = 2  # Total field of view in degrees
scaling_factor = 12  # Scaling factor to adjust the size of each tile
tile_side_length_deg = np.sqrt(total_fov_deg / num_tiles) * scaling_factor

# Plotting
sns.set_style('ticks')
fig, ax = plt.subplots(figsize=(18, 9))  # Increase figure size
plt.tick_params(axis='x', labelsize=30)
plt.tick_params(axis='y', labelsize=30)
plt.xlabel(r'$l$ (Galactic longitude) [deg]', fontsize=28)
plt.ylabel(r'$b$ (Galactic latitude) [deg]', fontsize=28)

# Scatter plot for LMC and SMC
ax.scatter(icrs_lmc.galactic.l.deg, icrs_lmc.galactic.b.deg, marker="o", s=120, c=sns.xkcd_palette(["forest green"]), edgecolors="g", zorder=3, lw=1, alpha=0.7, label="PNe")
ax.scatter(icrs_smc.galactic.l.deg, icrs_smc.galactic.b.deg, marker="o", s=120, c=sns.xkcd_palette(["forest green"]), edgecolors="g", zorder=4, lw=1, alpha=0.7)

# Plot tiles with larger symbols
for i in range(num_tiles):
    ax.add_patch(plt.Rectangle((icrs_tile.galactic.l.deg[i] - tile_side_length_deg / 2, icrs_tile.galactic.b.deg[i] - tile_side_length_deg / 2), tile_side_length_deg, tile_side_length_deg, color=(0.5, 0.5, 0.5), alpha=0.3))

# Annotations
bbox_props = dict(boxstyle="round", fc="w", ec="0.78", alpha=0.001, pad=0.1)
plt.text(0.1, 0.8, 'LMC', transform=ax.transAxes, c="gray", weight='bold', fontsize=38.8, zorder=1, bbox=bbox_props)
plt.text(0.85, 0.05, 'SMC', transform=ax.transAxes, c="gray", weight='bold', fontsize=38.8, zorder=2, bbox=bbox_props)

# Grid and legend
lgd_kws = {'frameon': True, 'fancybox': True, 'shadow': True}  # Define legend properties
ax.grid(True, linestyle='-.', linewidth=0.7)
ax.legend(prop={'family': 'monospace', 'size': 34}, **lgd_kws)

plt.tight_layout()

# Save figure
fig.savefig("galactic-coord-pn.pdf")

