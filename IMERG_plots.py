# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 18:53:55 2023

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
from datetime import date
#%% set IMERG dir
IMERGdir='C:/UCLA/research/IMERG V7 final'

#%% plotting IMERG
world=gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
Itaipu=gpd.read_file('shp/Itaipu basin.shp')
SF=gpd.read_file('shp/São Francisco basin.shp')

startyear=2001
startmonth=1

precimerg,preccropall,precchirps=[],[],[]
for i in range(270):
    #Loading data
    year=startyear
    monthi = startmonth+i
    while monthi > 12:
        monthi=monthi-12
        year+=1
    strmonth=''
    if monthi < 10:
        strmonth='0'+str(monthi)
    else:
        strmonth=str(monthi)
    f=Dataset(IMERGdir+'/3B-MO.MS.MRG.3IMERG.'+str(year)+strmonth+'01-S000000-E235959.'+strmonth+'.V07B.HDF5')
    #imerg prec unit: mm/hr
    
    #Setting IMERG and forecast data extent
    #extent N5.5,S--34.5;W-74.5,E-35.5
    precip = f['Grid/precipitation'][0][1055:1455,555:955].transpose()
    theLats = f['Grid/lat'][555:955]
    theLons = f['Grid/lon'][1055:1455]
    
    precip = xr.DataArray(precip,dims=['Lats','Lons'],coords=[theLats,theLons])
    precip=precip.where(precip>0,0)
    x, y = np.float32(np.meshgrid(theLons, theLats))
    precip.rio.write_crs(4326,inplace=True)
    precip.rio.set_spatial_dims('Lons','Lats',inplace=True)
    preccropItaipu=precip.rio.clip(Itaipu['geometry'],crs=4326).mean()
    preccropSF=precip.rio.clip(SF['geometry'],crs=4326).mean()
    if i==0:
        precimerg=xr.DataArray(precip)
        preccropallItaipu=xr.DataArray(preccropItaipu)
        preccropallSF=xr.DataArray(preccropSF)
    else:
        precimerg=xr.concat([precimerg,precip],'time')
        preccropallItaipu=xr.concat([preccropallItaipu,preccropItaipu],'time')
        preccropallSF=xr.concat([preccropallSF,preccropSF],'time')
    print(str(year)+' '+strmonth+' done.')
pd.DataFrame(preccropallItaipu).to_csv('result/Itaipu/imerg Itaipu 20012023.csv')
pd.DataFrame(preccropallSF).to_csv('result/Sobradinho/imerg Sobradinho 20012023.csv')

coarseimerg=precimerg.coarsen(Lats=10,Lons=10).mean()
anoimergwetseason=xr.DataArray()
for years in range(1,23):
    if years==1:
        anoimergwetseason=coarseimerg[years*12-2:years*12+3].sum(axis=0)
    else:
        anoimergwetseason=xr.concat([anoimergwetseason,coarseimerg[years*12-2:years*12+3].sum(axis=0)],'time')
anoimergwetseason=anoimergwetseason-anoimergwetseason.mean(axis=0)

plt.rcParams['font.size'] = 14

preccolor=plt.colormaps.get_cmap('RdBu')(np.linspace(0,1,12))
purple=np.array([165/256, 55/256, 253/256,1])
preccolor=np.append(preccolor,[purple],axis=0)
preccmap=ListedColormap(preccolor)
fig,ax=plt.subplots(2,6,figsize=(26,8),layout="constrained")
for month in range(12):
    monthmask=np.full(270,False)
    for year in range(22):
        monthmask[month+year*12]=True
    a=int(month/6)
    b=month%6
    cl=(precimerg[monthmask]*24*30).mean(axis=0).plot(ax=ax[a,b],levels=np.arange(0,390,30),cmap=preccmap,add_colorbar=False)
    world.boundary.plot(ax=ax[a,b],color='black')
    Itaipu.boundary.plot(ax=ax[a,b],color='Lime')
    SF.boundary.plot(ax=ax[a,b],color='Lime')
    ax[a,b].set_aspect('equal')
    ax[a,b].set_title('month '+str(month+1),fontsize=18)
    if a==0:
        ax[a,b].set_xlabel('')
        ax[a,b].set_xticklabels([])
    if b!=0:
        ax[a,b].set_ylabel('')
        ax[a,b].set_yticklabels([])
fig.colorbar(cl,ax=ax[:,5],shrink=0.8,ticks=np.arange(0,390,30),label='mm/month')
#plt.savefig('plots/TE precip spatial IMERG monthly.png')
plt.show()

preccolor=plt.colormaps.get_cmap('RdBu')(np.linspace(0,1,11))
purple=np.array([165/256, 55/256, 253/256,1])
preccolor=np.append(preccolor,[purple],axis=0)
preccmap=ListedColormap(preccolor)
fig,ax=plt.subplots(1)
cl=(precimerg[:240]*24*365).mean(axis=0).plot(ax=ax,levels=np.arange(0,2750,250),cmap=preccmap)
world.boundary.plot(ax=ax,color='black')
Itaipu.boundary.plot(ax=ax,color='Lime')
SF.boundary.plot(ax=ax,color='Lime')
ax.set_aspect('equal')
ax.set_title('Annual',fontsize=18)
#plt.savefig('plots/TE precip spatial IMERG Annual.png')
plt.show()

#%% NMME
#nmmelist=['CanCM4i']
nmmelist=['CanCM4i-IC3','CCSM4','CFSv2','GEM5-NEMO','GEOS5','GFDL-SPEAR']
Novall=[]
Novnmmecorrall=[]
for nmmemodel in nmmelist:
    startyear=2001
    startmonth=1
    precall=[]
    precnmmeavgIL0,precnmmeavgIL1,precnmmeavgIL2,precnmmeavgIL3,precnmmeavgIL4,precnmmeavgIL5,precnmmeavgIL6,precnmmeavgIL7,precnmmeavgIL8=[],[],[],[],[],[],[],[],[]
    precnmmelistI=[precnmmeavgIL0,precnmmeavgIL1,precnmmeavgIL2,precnmmeavgIL3,precnmmeavgIL4,precnmmeavgIL5,precnmmeavgIL6,precnmmeavgIL7,precnmmeavgIL8]
    precnmmeavgSFL0,precnmmeavgSFL1,precnmmeavgSFL2,precnmmeavgSFL3,precnmmeavgSFL4,precnmmeavgSFL5,precnmmeavgSFL6,precnmmeavgSFL7,precnmmeavgSFL8=[],[],[],[],[],[],[],[],[]
    precnmmelistSF=[precnmmeavgSFL0,precnmmeavgSFL1,precnmmeavgSFL2,precnmmeavgSFL3,precnmmeavgSFL4,precnmmeavgSFL5,precnmmeavgSFL6,precnmmeavgSFL7,precnmmeavgSFL8]
    
    for i in range(270):
        #Loading data
        year=startyear
        monthi = startmonth+i
        while monthi > 12:
            monthi=monthi-12
            year+=1
        strmonth=''
        if monthi < 10:
            strmonth='0'+str(monthi)
        else:
            strmonth=str(monthi)
        if nmmemodel=='CFSv2':
            if year >=1981 and year <1985:
                nmme=Dataset('E:/NMME/'+nmmemodel+'/precip M1-14 1981-1984.nc')
                nmme1=Dataset('E:/NMME/'+nmmemodel+'/precip M15-28 1981-1984.nc')
            elif year>=1985 and year<1988:
                nmme=Dataset('E:/NMME/'+nmmemodel+'/precip M1-14 1985-1988.nc')
                nmme1=Dataset('E:/NMME/'+nmmemodel+'/precip M15-28 1985-1988.nc')
            elif year>=1989 and year<1992:
                nmme=Dataset('E:/NMME/'+nmmemodel+'/precip M1-14 1989-1992.nc')
                nmme1=Dataset('E:/NMME/'+nmmemodel+'/precip M15-28 1989-1992.nc')
            elif year>=1993 and year<1996:
                nmme=Dataset('E:/NMME/'+nmmemodel+'/precip M1-14 1993-1996.nc')
                nmme1=Dataset('E:/NMME/'+nmmemodel+'/precip M15-28 1993-1996.nc')
            elif year>=1997 and year<2000:
                nmme=Dataset('E:/NMME/'+nmmemodel+'/precip M1-14 1997-2000.nc')
                nmme1=Dataset('E:/NMME/'+nmmemodel+'/precip M15-28 1997-2000.nc')
            elif year>=2001 and year<2004:
                nmme=Dataset('E:/NMME/'+nmmemodel+'/precip M1-14 2001-2004.nc')
                nmme1=Dataset('E:/NMME/'+nmmemodel+'/precip M15-28 2001-2004.nc')
            elif year>=2005 and year<2008:
                nmme=Dataset('E:/NMME/'+nmmemodel+'/precip M1-14 2005-2008.nc')
                nmme1=Dataset('E:/NMME/'+nmmemodel+'/precip M15-28 2005-2008.nc')
            elif year>=2009 and year<2012:
                nmme=Dataset('E:/NMME/'+nmmemodel+'/precip M1-14 2009-2012.nc')
                nmme1=Dataset('E:/NMME/'+nmmemodel+'/precip M15-28 2009-2012.nc')
            elif year>=2013 and year<2016:
                nmme=Dataset('E:/NMME/'+nmmemodel+'/precip M1-14 2013-2016.nc')
                nmme1=Dataset('E:/NMME/'+nmmemodel+'/precip M15-28 2013-2016.nc')
            elif year>=2017 and year<2020:
                nmme=Dataset('E:/NMME/'+nmmemodel+'/precip M1-14 2017-2020.nc')
                nmme1=Dataset('E:/NMME/'+nmmemodel+'/precip M15-28 2017-2020.nc')
            elif year>=2021:
                nmme=Dataset('E:/NMME/'+nmmemodel+'/precip M1-14 2021-present.nc')
                nmme1=Dataset('E:/NMME/'+nmmemodel+'/precip M15-28 2021-present.nc')
        else:
            if year >=1981 and year <1985:
                nmme=Dataset('E:/NMME/'+nmmemodel+'/precip 1981-1984.nc')
            elif year>=1985 and year<1988:
                nmme=Dataset('E:/NMME/'+nmmemodel+'/precip 1985-1988.nc')
            elif year>=1989 and year<1992:
                if nmmemodel=='GFDL-SPEAR':
                    nmme=Dataset('E:/NMME/'+nmmemodel+'/precip 1991-1992.nc')
                else:
                    nmme=Dataset('E:/NMME/'+nmmemodel+'/precip 1989-1992.nc')
            elif year>=1993 and year<1996:
                nmme=Dataset('E:/NMME/'+nmmemodel+'/precip 1993-1996.nc')
            elif year>=1997 and year<2000:
                nmme=Dataset('E:/NMME/'+nmmemodel+'/precip 1997-2000.nc')
            elif year>=2001 and year<2004:
                if nmmemodel=='GEOS5':
                    nmme=Dataset('E:/NMME/'+nmmemodel+'/precip 2001-2016.nc')
                else:
                    nmme=Dataset('E:/NMME/'+nmmemodel+'/precip 2001-2004.nc')
            elif year>=2005 and year<2008:
                if nmmemodel=='GEOS5':
                    nmme=Dataset('E:/NMME/'+nmmemodel+'/precip 2001-2016.nc')
                else:
                    nmme=Dataset('E:/NMME/'+nmmemodel+'/precip 2005-2008.nc')
            elif year>=2009 and year<2012:
                if nmmemodel=='GEOS5':
                    nmme=Dataset('E:/NMME/'+nmmemodel+'/precip 2001-2016.nc')
                else:
                    nmme=Dataset('E:/NMME/'+nmmemodel+'/precip 2009-2012.nc')
            elif year>=2013 and year<2016:
                if nmmemodel=='GEOS5':
                    nmme=Dataset('E:/NMME/'+nmmemodel+'/precip 2001-2016.nc')
                else:
                    nmme=Dataset('E:/NMME/'+nmmemodel+'/precip 2013-2016.nc')
            elif year>=2017 and year<2020:
                if nmmemodel=='GEOS5':
                    if i==192:#201701 should be in 2001-2016 for GEOS5
                        nmme=Dataset('E:/NMME/'+nmmemodel+'/precip 2001-2016.nc')
                    else:
                        nmme=Dataset('E:/NMME/'+nmmemodel+'/precip 2017-2020.nc')
                else:
                    nmme=Dataset('E:/NMME/'+nmmemodel+'/precip 2017-2020.nc')
            elif year>=2021:
                if nmmemodel=='CanCM4i-IC3' or nmmemodel=='GEM5-NEMO':
                    if i < 251:
                        nmme=Dataset('E:/NMME/'+nmmemodel+'/precip 2017-2020.nc')
                    else:
                        nmme=Dataset('E:/NMME/'+nmmemodel+'/precip 202110-present.nc')
                elif nmmemodel=='GFDL-SPEAR':
                    nmme=Dataset('E:/NMME/'+nmmemodel+'/precip M1-15 2021-present.nc')
                else:
                    nmme=Dataset('E:/NMME/'+nmmemodel+'/precip 2021-present.nc')
        month=np.where(nmme['S'][:]==492+i)[0][0] # month 684=20170101,492=20010101
        Latsouth=int(np.where(nmme['Y'][:]==-34)[0])
        Latnorth=int(np.where(nmme['Y'][:]==(5+1))[0])
        Lonwest=int(np.where(nmme['X'][:]==286)[0])
        Loneast=int(np.where(nmme['X'][:]==(325+1))[0]) #extent N5.5,S--34.5;W-74.5,E-35.5
        Latsnmme = nmme['Y'][Latsouth:Latnorth]
        Lonsnmme = nmme['X'][Lonwest:Loneast]-360
        nmmex,nmmey=np.float32(np.meshgrid(Lonsnmme, Latsnmme))
        if nmme['prec'].dimensions==('S','L','M','Y','X'):
            if nmmemodel == 'CFSv2':
                precnmme=np.ma.concatenate([nmme['prec'][month,0:9,:,Latsouth:Latnorth,Lonwest:Loneast],nmme1['prec'][month,0:9,:,Latsouth:Latnorth,Lonwest:Loneast]],axis=1)
                precnmme=xr.DataArray(precnmme,dims=['L','M','Lats','Lons'],coords=[range(9),range(1,len(nmme['M'])+len(nmme1['M'])+1),Latsnmme,Lonsnmme])
            else:
                precnmme=nmme['prec'][month,0:9,:,Latsouth:Latnorth,Lonwest:Loneast]
                precnmme=xr.DataArray(precnmme,dims=['L','M','Lats','Lons'],coords=[range(9),range(1,len(nmme['M'])+1),Latsnmme,Lonsnmme])
            if i == 0:
                precall=xr.DataArray(precnmme.mean(axis=1))
            else:
                precall=xr.concat([precall,precnmme.mean(axis=1)],'time')
        elif nmme['prec'].dimensions==('S','M','L','Y','X'):
            precnmme=nmme['prec'][month,:,0:9,Latsouth:Latnorth,Lonwest:Loneast]
            precnmme=xr.DataArray(precnmme,dims=['M','L','Lats','Lons'],coords=[range(1,len(nmme['M'])+1),range(9),Latsnmme,Lonsnmme])
            if i == 0:
                precall=xr.DataArray(precnmme.mean(axis=0))
            else:
                precall=xr.concat([precall,precnmme.mean(axis=0)],'time')
        precnmme.rio.write_crs(4326,inplace=True)
        precnmme.rio.set_spatial_dims('Lons','Lats',inplace=True)
        croppedItaipu=precnmme.rio.clip(Itaipu['geometry'],crs=4326)
        croppedSF=precnmme.rio.clip(SF['geometry'],crs=4326)
        
        if nmme['prec'].dimensions==('S','L','M','Y','X'):
            if nmmemodel=='GEOS5':
                for lead in range(len(precnmmelistI)):
                    precnmmelistI[lead].append(croppedItaipu[lead].mean(dim=['Lons','Lats'])[0:4])
                    precnmmelistSF[lead].append(croppedSF[lead].mean(dim=['Lons','Lats'])[0:4])
            elif nmmemodel=='GFDL-SPEAR':
                for lead in range(len(precnmmelistI)):
                    precnmmelistI[lead].append(croppedItaipu[lead].mean(dim=['Lons','Lats'])[0:15])
                    precnmmelistSF[lead].append(croppedSF[lead].mean(dim=['Lons','Lats'])[0:15])
            else:
                for lead in range(len(precnmmelistI)):
                    precnmmelistI[lead].append(croppedItaipu[lead].mean(dim=['Lons','Lats']))
                    precnmmelistSF[lead].append(croppedSF[lead].mean(dim=['Lons','Lats']))
        elif nmme['prec'].dimensions==('S','M','L','Y','X'):
            for lead in range(len(precnmmelistI)):
                precnmmelistI[lead].append(croppedItaipu[:,lead].mean(dim=['Lons','Lats']))
                precnmmelistSF[lead].append(croppedSF[:,lead].mean(dim=['Lons','Lats']))
        
        print(nmmemodel+' '+str(year)+' '+strmonth+' done.')
    precnmmeIalllead=xr.DataArray(precnmmelistI,dims=['lead','time','M'])
    precnmmeSFalllead=xr.DataArray(precnmmelistSF,dims=['lead','time','M'])
    '''for lead in range(len(precnmmelistI)):
        precnmmeIalllead[lead].to_pandas().to_csv('result/Itaipu/'+nmmemodel+'_lead'+str(lead)+'_20012023.csv')
        precnmmeSFalllead[lead].to_pandas().to_csv('result/Sobradinho/'+nmmemodel+'_lead'+str(lead)+'_20012023.csv')
'''
    
    
    Novmask=np.full(270,False)
    Novmaskimerg=np.full(270,False)
    for year in range(1,23):
        Novmask[year*12-2]=True
        Novmaskimerg[(year*12-2):(year*12+3)]=True
    Novall.append(precall[Novmask,0:5].mean(axis=(0,1))*151)
    
    release=11
    #for release in range(6,3,-1):
    forecastprec=xr.DataArray()
    for years in range(22):
        if years==0:
            forecastprec=precall[:][release-1,11-release:11+5-release].sum(axis=0)
        else:
            forecastprec=xr.concat([forecastprec,precall[:][release-1+12*years,11-release:11+5-release].sum(axis=0)],'time')
    forecastprec=forecastprec-forecastprec.mean(axis=0)
    anoimergwetseason['Lats']=forecastprec.Lats
    anoimergwetseason['Lons']=forecastprec.Lons
    forecastcorrimerg=xr.corr(anoimergwetseason[:],forecastprec,dim='time')
    Novnmmecorrall.append(forecastcorrimerg)

preccolor=plt.colormaps.get_cmap('RdBu')(np.linspace(0,1,12))
purple=np.array([165/256, 55/256, 253/256,1])
preccolor=np.append(preccolor,[purple],axis=0)
preccmap=ListedColormap(preccolor)

fig,ax=plt.subplots(3,3,figsize=(11,11),layout="constrained")
counter=0
cl=(coarseimerg[Novmaskimerg]*24*30*5).mean(axis=0).plot(ax=ax[1,0],levels=np.arange(0,1650,150),cmap=preccmap,add_colorbar=False)
world.boundary.plot(ax=ax[1,0],color='black')
Itaipu.boundary.plot(ax=ax[1,0],color='deeppink')
SF.boundary.plot(ax=ax[1,0],color='Lime')
ax[1,0].set_aspect('equal')
ax[1,0].set_title('IMERG',fontsize=24)
ax[1,0].set_ylabel('Latitude',fontsize=18)
ax[1,0].set_xlabel('')
ax[0,0].axis('off')
ax[2,0].axis('off')
for a in range(3):
    for b in range(2):
        cl=Novall[counter].plot(ax=ax[a,b+1],levels=np.arange(0,1650,150),cmap=preccmap,add_colorbar=False)
        world.boundary.plot(ax=ax[a,b+1],color='black')
        Itaipu.boundary.plot(ax=ax[a,b+1],color='deeppink')
        SF.boundary.plot(ax=ax[a,b+1],color='Lime')
        ax[a,b+1].set_title(nmmelist[counter],fontsize=24)
        ax[a,b+1].set_aspect('equal')
        ax[a,b+1].set_xlabel('')
        ax[a,b+1].set_ylabel('')
        if (a!=0 and b!=0) or (a!=2 and b!=0):
            ax[a,b+1].set_yticklabels([])
        if (a!=2 and b!=0) or (a!=2 and b!=1):
            ax[a,b+1].set_xticklabels([])
        counter+=1
ax[1,1].set_yticklabels([])
ax[2,1].set_xlabel('Longitude',fontsize=18)
fig.colorbar(cl,ax=ax[:,2],shrink=0.8,ticks=np.arange(0,1650,150),label='mm/season')
plt.show()

corrcolor=plt.colormaps.get_cmap('RdBu')(np.linspace(0,1,9))
corrcolor[3]=np.array([1,1,1,1])
corrcolor[4]=np.array([1,1,1,1])
corrcmap=ListedColormap(corrcolor)

#new correlation spatial plot, all in one
fig,ax=plt.subplots(2,3,figsize=(12,8),layout="constrained")
counter=0
for a in range(2):
    for b in range(3):
        cl=Novnmmecorrall[counter].plot(ax=ax[a,b],levels=[-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1],cmap=corrcmap,add_colorbar=False)
        world.boundary.plot(ax=ax[a,b],color='black')
        Itaipu.boundary.plot(ax=ax[a,b],color='deeppink')
        SF.boundary.plot(ax=ax[a,b],color='Lime')
        ax[a,b].set_title(nmmelist[counter],fontsize=24)
        ax[a,b].set_aspect('equal')
        ax[a,b].set_xlabel('')
        ax[a,b].set_ylabel('')
        if a==0:
            ax[a,b].set_xticklabels([])
        if b!=0:
            ax[a,b].set_yticklabels([])
        counter+=1
ax[1,1].set_xlabel('Longitude',fontsize=18)
fig.supylabel('Latitude',fontsize=18)
fig.colorbar(cl,ax=ax[:,2],shrink=0.8,ticks=[-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1],label='correlation')
plt.show()