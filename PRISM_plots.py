# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 12:09:36 2023

@author: EricTien
"""

import os
path=__file__
os.chdir(os.path.dirname(path))
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import xarray as xr
import rioxarray as rxr
import geopandas as gpd
import fiona
import scipy
from datetime import date,timedelta

#%% set PRISM dir
PRISMdir='C:/UCLA/research/PRISM'

#%% Plotting PRISM 2023 water year wet season (Oct to Mar) for western US
# extent 125W to 105W, 30N - 50 N

usstate=gpd.read_file('shp/cb_2019_us_state_20m.shp')
CVbasin=gpd.read_file('shp/Central Valley basin.shp')
Colorado=gpd.read_file('shp/Colorado River Basin in US.shp')
prismlist=['202210','202211','202212','202301','202302','202303']
wusprismtotal,wusprismall=[],[]
for i in range(6):
    datatype='provisional'
    prism=rxr.open_rasterio(PRISMdir+'/PRISM_ppt_'+datatype+'_4kmM3_'+prismlist[i]+'_bil.bil')[0]
    wusprism=prism[prism.y>30,prism.x<-105]
    wusprism=wusprism.where(wusprism>-1)
    wusprism=wusprism.rename({'y':'Lats','x':'Lons'})
    wusprism.rio.write_crs(4326,inplace=True)
    wusprism.rio.set_spatial_dims('Lons','Lats',inplace=True)
    wusprismall.append(wusprism)
    if i==0:
        wusprismtotal=wusprism
    else:
        wusprismtotal=wusprismtotal+wusprism

theLats=wusprism.Lats
theLons=wusprism.Lons
wusprismall=xr.DataArray(wusprismall,dims=['month','Lats','Lons'],coords=[[10,11,12,1,2,3],theLats,theLons])

plt.rcParams['font.size'] = 18
preccolor=plt.colormaps.get_cmap('RdBu')(np.linspace(0,1,12))
purple=np.array([165/256, 55/256, 253/256,1])
preccolor=np.append(preccolor,[purple],axis=0)
preccmap=ListedColormap(preccolor)
fig,ax=plt.subplots(figsize=(12,12))
wusprismtotal.plot(ax=ax,cmap=preccmap,levels=np.arange(0,2200,200))
usstate.boundary.plot(ax=ax,color='gray')
CVbasin.boundary.plot(ax=ax,color='lime')
Colorado.boundary.plot(ax=ax,color='lime')
ax.set_facecolor('silver')
ax.set_title('202210-202303 Precipitation')
ax.set_aspect('equal')
plt.show()

'''plt.rcParams['font.size'] = 12
preccolor=plt.colormaps.get_cmap('RdBu')(np.linspace(0,1,12))
purple=np.array([165/256, 55/256, 253/256,1])
preccolor=np.append(preccolor,[purple],axis=0)
preccmap=ListedColormap(preccolor)
monthlist=[10,11,12,1,2,3]
fig,ax=plt.subplots(1,6)
for i in range(6):
    cl=wusprismall[i].plot(ax=ax[i],cmap=preccmap,levels=np.arange(0,330,30),add_colorbar=False)
    usstate.boundary.plot(ax=ax[i],color='black')
    ax[i].set_title('Month '+str(monthlist[i]))
    ax[i].set_aspect('equal')
fig.colorbar(cl,ax=ax[5],shrink=0.8,ticks=np.arange(0,330,30),label='mm')
plt.show()'''
#%% Plotting PRISM wetseason average
prismall=[[],[],[],[],[],[]]
for year in range(2001,2022):
    for monthno in range(6):
        if monthno>=3:
            ymstr=str(year+1)+'0'+str(monthno-2)
        else:
            ymstr=str(year)+str(monthno+10)
        datatype='stable'
        prism=rxr.open_rasterio(PRISMdir+'/PRISM_ppt_'+datatype+'_4kmM3_'+ymstr+'_bil.bil')[0]
        wusprism=prism[prism.y>30,prism.x<-105]
        wusprism=wusprism.where(wusprism>-1)
        wusprism=wusprism.rename({'y':'Lats','x':'Lons'})
        wusprism.rio.write_crs(4326,inplace=True)
        wusprism.rio.set_spatial_dims('Lons','Lats',inplace=True)
        if year==2001:
            prismall[monthno]=wusprism.values
        else:
            prismall[monthno]=prismall[monthno]+wusprism.values

theLats=wusprism.Lats
theLons=wusprism.Lons
for monthno in range(6):
    prismall[monthno]=xr.DataArray(prismall[monthno],dims=['Lats','Lons'],coords=[theLats,theLons])/21
prismallwet=xr.DataArray(prismall,dims=['month','Lats','Lons'],coords=[[10,11,12,1,2,3],theLats,theLons])

plt.rcParams['font.size'] = 18
preccolor=plt.colormaps.get_cmap('RdBu')(np.linspace(0,1,12))
purple=np.array([165/256, 55/256, 253/256,1])
preccolor=np.append(preccolor,[purple],axis=0)
preccmap=ListedColormap(preccolor)
fig,ax=plt.subplots(figsize=(12,12))
ax.set_facecolor('silver')
CVbasin.boundary.plot(ax=ax,color='lime')
Colorado.boundary.plot(ax=ax,color='lime')
prismallwet.sum(axis=0,skipna=False).plot(ax=ax,cmap=preccmap,levels=np.arange(0,2200,200))
usstate.boundary.plot(ax=ax,color='gray')
ax.set_title('2001~2021 Oct-Mar Precipitation')
ax.set_aspect('equal')
plt.show()

plt.rcParams['font.size'] = 12
preccolor=plt.colormaps.get_cmap('RdBu')(np.linspace(0,1,12))
purple=np.array([165/256, 55/256, 253/256,1])
preccolor=np.append(preccolor,[purple],axis=0)
preccmap=ListedColormap(preccolor)
monthlist=[10,11,12,1,2,3]
fig,ax=plt.subplots(2,6,figsize=(26,8),layout="constrained")
for i in range(6):
    cl=wusprismall[i].plot(ax=ax[0,i],cmap=preccmap,levels=np.arange(0,330,30),add_colorbar=False)
    usstate.boundary.plot(ax=ax[0,i],color='gray')
    ax[0,i].set_facecolor('silver')
    CVbasin.boundary.plot(ax=ax[0,i],color='lime')
    Colorado.boundary.plot(ax=ax[0,i],color='lime')
    ax[0,i].set_title('2022 Month '+str(monthlist[i]))
    ax[0,i].set_aspect('equal')
for i in range(6):
    cl=prismallwet[i].plot(ax=ax[1,i],cmap=preccmap,levels=np.arange(0,330,30),add_colorbar=False)
    usstate.boundary.plot(ax=ax[1,i],color='gray')
    ax[1,i].set_facecolor('silver')
    CVbasin.boundary.plot(ax=ax[1,i],color='lime')
    Colorado.boundary.plot(ax=ax[1,i],color='lime')
    ax[1,i].set_title('2001-2021 Month '+str(monthlist[i]))
    ax[1,i].set_aspect('equal')
fig.colorbar(cl,ax=ax[:,5],shrink=0.8,ticks=np.arange(0,330,30),label='mm')
plt.show()
#%% plotting for 1980 to 2022

prismall=[[],[],[],[],[],[],[],[],[],[],[],[]]
for year in range(1981,2022):
    for monthno in range(12):
        if monthno>=3:
            ymstr=str(year+1)+'0'+str(monthno-2)
            if year>1979:
                M='3'
            else:
                M='2'
        else:
            ymstr=str(year)+str(monthno+10)
            if year>1980:
                M='3'
            else:
                M='2'
        datatype='stable'
        prism=rxr.open_rasterio(PRISMdir+'/PRISM_ppt_'+datatype+'_4kmM'+M+'_'+ymstr+'_bil.bil')[0]
        CAprism=prism[(prism.y>31)&(prism.y<43),(prism.x<-113)&(prism.x>-125)]
        CAprism=CAprism.where(CAprism>-1)
        CAprism=CAprism.rename({'y':'Lats','x':'Lons'})
        CAprism.rio.write_crs(4326,inplace=True)
        CAprism.rio.set_spatial_dims('Lons','Lats',inplace=True)
        if year==1981:
            prismall[monthno]=CAprism.values
        else:
            prismall[monthno]=prismall[monthno]+CAprism.values
        print(ymstr+' done.')
        
theLats=CAprism.Lats
theLons=CAprism.Lons
for monthno in range(12):
    prismall[monthno]=xr.DataArray(prismall[monthno],dims=['Lats','Lons'],coords=[theLats,theLons])
prismallmonths=xr.DataArray(prismall,dims=['month','Lats','Lons'],coords=[[10,11,12,1,2,3,4,5,6,7,8,9],theLats,theLons])/len(list(range(1980,2022)))

plt.rcParams['font.size'] = 18
preccolor=plt.colormaps.get_cmap('RdBu')(np.linspace(0,1,11))
purple=np.array([165/256, 55/256, 253/256,1])
preccolor=np.append(preccolor,[purple],axis=0)
preccmap=ListedColormap(preccolor)
fig,ax=plt.subplots(figsize=(12,12))
ax.set_facecolor('silver')
prismallmonths.sum(axis=0,skipna=False).plot(ax=ax,cmap=preccmap,levels=np.arange(0,2500,250),cbar_kwargs={'label': 'mm/yr'})
usstate.boundary.plot(ax=ax,color='gray')
ax.set_title('1980-2021 Annual Precipitation')
ax.set_xticklabels(['','124W','122W','120W','118W','116W','114W'])
ax.set_yticklabels(['','32N','34N','36N','38N','40N','42N'])
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_aspect('equal')
plt.show()

#%% calculate Central Valley basin and Colorado basin precip of PRISM 1981-present
prismCV,prismC=[],[]
for year in range(1981,2024):
    for monthno in range(12):
        if year==2023 and monthno==3:
            break
        if monthno<9:
            strmonth='0'+str(monthno+1)
        else:
            strmonth=str(monthno+1)
        ymstr=str(year)+strmonth
        datatype='stable'
        if (year == 2022 and monthno>=9) or year==2023:
            datatype='provisional'
        prism=rxr.open_rasterio('C:/UCLA/research/PRISM/PRISM_ppt_'+datatype+'_4kmM3''_'+ymstr+'_bil.bil')[0]
        wusprism=prism[prism.y>30,prism.x<-105]
        wusprism=wusprism.where(wusprism>-1)
        wusprism=wusprism.rename({'y':'Lats','x':'Lons'})
        wusprism.rio.write_crs(4326,inplace=True)
        wusprism.rio.set_spatial_dims('Lons','Lats',inplace=True)
        croppedCV=wusprism.rio.clip(CVbasin['geometry'],crs=4326)
        croppedCV=croppedCV.where(croppedCV>-1).mean()
        croppedC=wusprism.rio.clip(Colorado['geometry'],crs=4326)
        croppedC=croppedC.where(croppedC>-1).mean()
        prismCV.append(croppedCV)
        prismC.append(croppedC)
        print(ymstr+' done.')

pd.DataFrame(np.array(prismCV)).to_csv('result/Central Valley/PRISM_19812023.csv')
pd.DataFrame(np.array(prismC)).to_csv('result/Colorado/PRISM_19812023.csv')

#%% plotting for 1981 to 2022, for multi scale paper figure 1

prismall=[[],[],[],[],[],[],[],[],[],[],[],[]]
for year in range(1981,2022):
    for monthno in range(12):
        if monthno>=3:
            ymstr=str(year+1)+'0'+str(monthno-2)
            if year>1979:
                M='3'
            else:
                M='2'
        else:
            ymstr=str(year)+str(monthno+10)
            if year>1980:
                M='3'
            else:
                M='2'
        datatype='stable'
        prism=rxr.open_rasterio('C:/UCLA/research/PRISM/PRISM_ppt_'+datatype+'_4kmM'+M+'_'+ymstr+'_bil.bil')[0]
        wusprism=prism[prism.y>30,prism.x<-105]
        wusprism=wusprism.where(wusprism>-1)
        wusprism=wusprism.rename({'y':'Lats','x':'Lons'})
        wusprism.rio.write_crs(4326,inplace=True)
        wusprism.rio.set_spatial_dims('Lons','Lats',inplace=True)
        if year==1981:
            prismall[monthno]=wusprism.values
        else:
            prismall[monthno]=prismall[monthno]+wusprism.values
        print(ymstr+' done.')
        
theLats=wusprism.Lats
theLons=wusprism.Lons
for monthno in range(12):
    prismall[monthno]=xr.DataArray(prismall[monthno],dims=['Lats','Lons'],coords=[theLats,theLons])
prismallmonths=xr.DataArray(prismall,dims=['month','Lats','Lons'],coords=[[10,11,12,1,2,3,4,5,6,7,8,9],theLats,theLons])/len(list(range(1980,2022)))

plt.rcParams['font.size'] = 14
preccolor=plt.colormaps.get_cmap('RdBu')(np.linspace(0,1,11))
purple=np.array([165/256, 55/256, 253/256,1])
preccolor=np.append(preccolor,[purple],axis=0)
preccmap=ListedColormap(preccolor)
fig,ax=plt.subplots(figsize=(8,8))
ax.set_facecolor('silver')
prismallmonths[2:6].sum(axis=0,skipna=False).plot(ax=ax,cmap=preccmap,levels=np.arange(0,1600,160),cbar_kwargs={'label': 'mm/yr'})
usstate.boundary.plot(ax=ax,color='gray')
CVbasin.boundary.plot(ax=ax,color='lime',linewidth=2)
ax.set_title('1981-2021 Winter Precipitation')
#ax.set_xticks([])
#ax.set_yticks([])
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_aspect('equal')
plt.show()

plt.rcParams['font.size'] = 14
preccolor=plt.colormaps.get_cmap('RdBu')(np.linspace(0,1,11))
purple=np.array([165/256, 55/256, 253/256,1])
preccolor=np.append(preccolor,[purple],axis=0)
preccmap=ListedColormap(preccolor)
fig,ax=plt.subplots(figsize=(8,8))
ax.set_facecolor('silver')
(prismallmonths[2:6].sum(axis=0,skipna=False)/prismallmonths.sum(axis=0,skipna=False)*100).plot(ax=ax,cmap=preccmap,levels=np.arange(0,110,10),cbar_kwargs={'label': 'winter %'})
usstate.boundary.plot(ax=ax,color='gray')
CVbasin.boundary.plot(ax=ax,color='lime',linewidth=2)
ax.set_title('1981-2021 Winter Precipitation Percentage')
ax.set_xticklabels(['','124W','122W','120W','118W','116W','114W'])
ax.set_yticklabels(['','32N','34N','36N','38N','40N','42N'])
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_aspect('equal')
plt.show()