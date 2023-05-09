# -*- coding: utf-8 -*-
"""
Created on Tue May  2 09:54:51 2023

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

input_file= 'MERRA2_400.tavg1_2d_flx_Nx.20110522.nc4'
nc1 = NetCDFFile(input_file)
joplin_lon = -94.5153
joplin_lat = 37.0842
time = nc1.variables['time'][:]
nt1 = np.size(time)
lats2 = np.array(nc1.variables['lat'][:])
lons2 = np.array(nc1.variables['lon'][:])
pr = nc1.variables['PREVTOT'][0, :, :]
speedmax = nc1.variables['SPEEDMAX'][:]
cn = nc1.variables['CN'][0, :, :]

fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_extent([joplin_lon - 16, joplin_lon + 16, joplin_lat - 16, joplin_lat + 16], crs=ccrs.PlateCarree())

ax.add_feature(cfeature.BORDERS, linewidth=0.5)
ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax.add_feature(cfeature.STATES, linestyle=':', linewidth=0.5)

cmap = plt.get_cmap('viridis')
norm = plt.Normalize(vmin=np.min(pr), vmax=np.max(pr))
cax = ax.pcolormesh(lons2, lats2, pr, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm)
cbar = fig.colorbar(cax, ax=ax, orientation='horizontal', pad=0.05)

cbar.set_label('(kg/m^2/s)')
ax.set_title(f'Tasa de precipitación total (MERRA-2)')
ax.plot(joplin_lon, joplin_lat, marker='o', markersize=5, color='red', transform=ccrs.PlateCarree())
ax.text(joplin_lon + 0.1, joplin_lat + 0.1, 'Joplin, Missouri', transform=ccrs.PlateCarree(), fontsize=12, color='red')

plt.show()




speedmax_t0 = speedmax[0, :, :]


fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())

# Agregar características geográficas al mapa
ax.add_feature(cfeature.BORDERS, linewidth=0.5)
ax.add_feature(cfeature.COASTLINE, linewidth=0.5)


cmap = plt.get_cmap('viridis')
norm = plt.Normalize(vmin=0, vmax=np.max(speedmax_t0))
cax = ax.pcolormesh(lons2, lats2, speedmax_t0, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm)
cbar = fig.colorbar(cax, ax=ax, orientation='horizontal', pad=0.05)


cbar.set_label('Velocidad máxima del viento en superficie (m/s)')
ax.set_title('Velocidad máxima del viento en superficie (MERRA-2) - 22 de mayo de 2011, 00:00 UTC')


plt.show()


def create_speedmax_image(time_index):
    speedmax_t = speedmax[time_index, :, :]

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent([joplin_lon - 20, joplin_lon + 20, joplin_lat - 20, joplin_lat + 20], crs=ccrs.PlateCarree())

    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.STATES, linestyle=':', linewidth=0.5)

    cmap = plt.get_cmap('viridis')
    norm = plt.Normalize(vmin=0, vmax=np.max(speedmax_t))
    cax = ax.pcolormesh(lons2, lats2, speedmax_t, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm)
    cbar = fig.colorbar(cax, ax=ax, orientation='horizontal', pad=0.05)

    cbar.set_label('Velocidad máxima del viento en superficie (m/s)')
    ax.set_title(f'Velocidad máxima del viento en superficie (MERRA-2) - 22 de mayo de 2011, {time_index:02d}:00 UTC')
    ax.plot(joplin_lon, joplin_lat, marker='o', markersize=5, color='red', transform=ccrs.PlateCarree())
    ax.text(joplin_lon + 0.1, joplin_lat + 0.1, 'Joplin, Missouri', transform=ccrs.PlateCarree(), fontsize=12, color='red')

    fig.canvas.draw()
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    
    plt.close(fig)

    return image


images = [create_speedmax_image(i) for i in range(nt1)]


gif_path = os.path.abspath("speedmax_map_time_series.gif")
imageio.mimsave('speedmax_map_time_series.gif', images, fps=2)


webbrowser.open("file://" + gif_path)


def create_prevtot_image(time_index):
    prevtot_t = nc1.variables['PREVTOT'][time_index, :, :]
    current_time = num2date(time[time_index], nc1.variables['time'].units)

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent([joplin_lon - 20, joplin_lon + 20, joplin_lat - 20, joplin_lat + 20], crs=ccrs.PlateCarree())

    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.STATES, linestyle=':', linewidth=0.5)

    cmap = plt.get_cmap('viridis')
    norm = plt.Normalize(vmin=0, vmax=np.max(prevtot_t))
    cax = ax.pcolormesh(lons2, lats2, prevtot_t, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm)
    cbar = fig.colorbar(cax, ax=ax, orientation='horizontal', pad=0.05)

    cbar.set_label('Tasa de precipitación total (kg/m^2/s)')
    ax.set_title(f'Tasa de precipitación total (MERRA-2) - {current_time.strftime("%Y-%m-%d %H:%M:%S")} UTC')
    ax.plot(joplin_lon, joplin_lat, marker='o', markersize=5, color='red', transform=ccrs.PlateCarree())
    ax.text(joplin_lon + 0.1, joplin_lat + 0.1, 'Joplin, Missouri', transform=ccrs.PlateCarree(), fontsize=12, color='red')

    fig.canvas.draw()
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    plt.close(fig)

    return image

images = [create_prevtot_image(i) for i in range(nt1)]

gif_path = os.path.abspath("prevtot_map_time_series.gif")
imageio.mimsave('prevtot_map_time_series.gif', images, fps=2)

webbrowser.open("file://" + gif_path)


time_index = 0
cn_t = nc1.variables['CN'][time_index, :, :]

fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_extent([joplin_lon - 12, joplin_lon + 12, joplin_lat - 12, joplin_lat + 12], crs=ccrs.PlateCarree())

ax.add_feature(cfeature.BORDERS, linewidth=0.5)
ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax.add_feature(cfeature.STATES, linestyle=':', linewidth=0.5)

cmap = plt.get_cmap('viridis')
norm = plt.Normalize(vmin=0, vmax=np.max(cn_t))
cax = ax.pcolormesh(lons2, lats2, cn_t, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm)
cbar = fig.colorbar(cax, ax=ax, orientation='horizontal', pad=0.05)

cbar.set_label('Concentración de nubes')
ax.set_title(f'Concentración de nubes (MERRA-2) - 22 de mayo de 2011, {time_index:02d}:00 UTC')
ax.plot(joplin_lon, joplin_lat, marker='o', markersize=5, color='red', transform=ccrs.PlateCarree())
ax.text(joplin_lon + 0.1, joplin_lat + 0.1, 'Joplin, Missouri', transform=ccrs.PlateCarree(), fontsize=12, color='red')

plt.show()




def create_tczpbl_image(time_index):
    tczpbl_t = nc1.variables['TCZPBL'][time_index, :, :]

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent([joplin_lon - 20, joplin_lon + 20, joplin_lat - 20, joplin_lat + 20], crs=ccrs.PlateCarree())

    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.STATES, linestyle=':', linewidth=0.5)

    cmap = plt.get_cmap('viridis')
    norm = plt.Normalize(vmin=np.min(tczpbl_t), vmax=np.max(tczpbl_t))
    cax = ax.pcolormesh(lons2, lats2, tczpbl_t, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm)
    cbar = fig.colorbar(cax, ax=ax, orientation='horizontal', pad=0.05)

    cbar.set_label('Altura de la capa límite planetaria (m)')
    ax.set_title(f'Altura de la capa límite planetaria (MERRA-2) - 22 de mayo de 2011, {time_index:02d}:00 UTC')
    ax.plot(joplin_lon, joplin_lat, marker='o', markersize=5, color='red', transform=ccrs.PlateCarree())
    ax.text(joplin_lon + 0.1, joplin_lat + 0.1, 'Joplin, Missouri', transform=ccrs.PlateCarree(), fontsize=12, color='red')

    fig.canvas.draw()
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    plt.close(fig)

    return image

images = [create_tczpbl_image(i) for i in range(nt1)]
gif_path = os.path.abspath("tczpbl_map_time_series.gif")
imageio.mimsave('tczpbl_map_time_series.gif', images, fps=2)

webbrowser.open("file://" + gif_path)



def create_qlml_image(time_index):
    qlml_t = nc1.variables['QLML'][time_index, :, :]

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent([joplin_lon - 20, joplin_lon + 20, joplin_lat - 20, joplin_lat + 20], crs=ccrs.PlateCarree())

    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.STATES, linestyle=':', linewidth=0.5)

    cmap = plt.get_cmap('viridis')
    norm = plt.Normalize(vmin=np.min(qlml_t), vmax=np.max(qlml_t))
    cax = ax.pcolormesh(lons2, lats2, qlml_t, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm)
    cbar = fig.colorbar(cax, ax=ax, orientation='horizontal', pad=0.05)

    cbar.set_label('Contenido de humedad en la capa límite (kg/m^2)')
    ax.set_title(f'Contenido de humedad en la capa límite (MERRA-2) - 22 de mayo de 2011, {time_index:02d}:00 UTC')
    ax.plot(joplin_lon, joplin_lat, marker='o', markersize=5, color='red', transform=ccrs.PlateCarree())
    ax.text(joplin_lon + 0.1, joplin_lat + 0.1, 'Joplin, Missouri', transform=ccrs.PlateCarree(), fontsize=12, color='red')

    fig.canvas.draw()
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    plt.close(fig)

    return image
images = [create_qlml_image(i) for i in range(nt1)]
gif_path = os.path.abspath("qlml_map_time_series.gif")
imageio.mimsave('qlml_map_time_series.gif', images, fps=2)

webbrowser.open("file://" + gif_path)





def create_tlml_image(time_index):
    tlml_t = nc1.variables['TLML'][time_index, :, :]

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent([joplin_lon - 20, joplin_lon + 20, joplin_lat - 20, joplin_lat + 20], crs=ccrs.PlateCarree())

    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.STATES, linestyle=':', linewidth=0.5)

    cmap = plt.get_cmap('coolwarm')
    norm = plt.Normalize(vmin=np.min(tlml_t), vmax=np.max(tlml_t))
    cax = ax.pcolormesh(lons2, lats2, tlml_t, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm)
    cbar = fig.colorbar(cax, ax=ax, orientation='horizontal', pad=0.05)

    cbar.set_label('Temperatura en la capa límite (K)')
    ax.set_title(f'Temperatura en la capa límite (MERRA-2) - 22 de mayo de 2011, {time_index:02d}:00 UTC')
    ax.plot(joplin_lon, joplin_lat, marker='o', markersize=5, color='red', transform=ccrs.PlateCarree())
    ax.text(joplin_lon + 0.1, joplin_lat + 0.1, 'Joplin, Missouri', transform=ccrs.PlateCarree(), fontsize=12, color='red')

    fig.canvas.draw()
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    plt.close(fig)

    return image
images = [create_tlml_image(i) for i in range(nt1)]
gif_path = os.path.abspath("tlml_map_time_series.gif")
imageio.mimsave('tlml_map_time_series.gif', images, fps=2)

webbrowser.open("file://" + gif_path)












nc1.close()












