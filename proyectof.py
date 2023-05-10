# -*- coding: utf-8 -*-
"""
Created on Sat Apr 29 12:17:04 2023

@author: rcorr
"""

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import num2date, date2num, Dataset as NetCDFFile
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import imageio
import webbrowser
import os


input_file= 'MERRA2_400.inst3_3d_asm_Np.20110522.nc4'
nc = NetCDFFile(input_file)

lons= nc.variables['lon'][:]; Nlons=np.size(lons)
lats= nc.variables['lat'][:]; Nlons=np.size(lats)
levs= nc.variables['lev'][:]; Nlons=np.size(levs)
T= nc.variables['T'][:]
ps= nc.variables['PS'][:]
u = nc.variables['U'][0, 0, :, :]
v = nc.variables['V'][0, 0, :, :]
t = np.arange(0, 24, 3)
nt = T.shape[0]
wind_magnitude = np.sqrt(u**2 + v**2)
joplin_lon = -94.5153
joplin_lat = 37.0842







images = []

for i in range(nt):
    ps = nc.variables['PS'][i, :, :]
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent([joplin_lon - 20, joplin_lon + 20, joplin_lat - 20, joplin_lat + 20], crs=ccrs.PlateCarree())

    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.STATES, linestyle=':', linewidth=0.5)

    cmap = plt.get_cmap('viridis')
    norm = plt.Normalize(vmin=np.min(ps), vmax=np.max(ps))
    cax = ax.pcolormesh(lons, lats, ps, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm)
    cbar = fig.colorbar(cax, ax=ax, orientation='horizontal', pad=0.05)

    cbar.set_label('Presión en superficie (Pa)')
    ax.set_title(f'Presión en superficie (MERRA-2) - 22 de mayo de 2011, {t[i]:02d}:00 UTC')
    ax.plot(joplin_lon, joplin_lat, marker='o', markersize=5, color='red', transform=ccrs.PlateCarree())
    ax.text(joplin_lon + 0.1, joplin_lat + 0.1, 'Joplin, Missouri', transform=ccrs.PlateCarree(), fontsize=12, color='red')

    fig.canvas.draw()
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    images.append(image)

    plt.close(fig)
gif_path = os.path.abspath("ps_map_time_series.gif")
imageio.mimsave('ps_map_time_series.gif', images, fps=2)

webbrowser.open("file://" + gif_path)








temperature_surface = T[0, 0, :, :]

fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())

ax.add_feature(cfeature.BORDERS, linewidth=0.5)
ax.add_feature(cfeature.COASTLINE, linewidth=0.5)

cmap = plt.get_cmap('coolwarm')
norm = plt.Normalize(vmin=np.min(temperature_surface), vmax=np.max(temperature_surface))
cax = ax.pcolormesh(lons, lats, temperature_surface, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm)
cbar = fig.colorbar(cax, ax=ax, orientation='horizontal', pad=0.05)


cbar.set_label('Temperatura en superficie (K)')
ax.set_title('Temperatura en superficie (MERRA-2) - 22 de mayo de 2011, 00:00 UTC')


plt.show()




images = []

for i in range(nt):
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent([joplin_lon - 20, joplin_lon + 20, joplin_lat - 20, joplin_lat + 20], crs=ccrs.PlateCarree())

    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.STATES, linestyle=':', linewidth=0.5)

    u = nc.variables['U'][i, 0, :, :]
    v = nc.variables['V'][i, 0, :, :]
    wind_magnitude = np.sqrt(u**2 + v**2)

    cmap = plt.get_cmap('viridis')
    norm = plt.Normalize(vmin=0, vmax=np.max(wind_magnitude))
    cax = ax.pcolormesh(lons, lats, wind_magnitude, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm)
    cbar = fig.colorbar(cax, ax=ax, orientation='horizontal', pad=0.05)

   
    ax.streamplot(lons, lats, u, v, transform=ccrs.PlateCarree(), density=2, color='k', linewidth=0.5)

    cbar.set_label('Magnitud de la velocidad del viento en superficie (m/s)')
    ax.set_title(f'Magnitud y dirección de la velocidad del viento en superficie (MERRA-2) - 22 de mayo de 2011, {t[i]:02d}:00 UTC')
    ax.plot(joplin_lon, joplin_lat, marker='o', markersize=5, color='red', transform=ccrs.PlateCarree())
    ax.text(joplin_lon + 0.1, joplin_lat + 0.1, 'Joplin, Missouri', transform=ccrs.PlateCarree(), fontsize=12, color='red')

    fig.canvas.draw()
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    images.append(image)

    plt.close(fig)



gif_path = os.path.abspath("wind_map_time_series.gif")
imageio.mimsave('wind_map_time_series.gif', images, fps=2)

webbrowser.open("file://" + gif_path)












