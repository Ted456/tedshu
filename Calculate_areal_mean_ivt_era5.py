# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 11:25:42 2020

@author: ted
"""

import numpy as np
import pandas as pd
import os
import xarray as xr
dirname = os.path.dirname(__file__)
ARfilepath = os.path.join(dirname, 'ars_6h_daily_0to70S_100Eto120W_025025_197909to201808')
Rainfall_path = os.path.join(dirname, 'ClimateDriverRainfallAndARdateData/Rainfall')
os.chdir(Rainfall_path)
f = pd.read_csv('644_station_with_group.csv')
Station_List = f['Agent']
x1 = f['x1Index']
x2 = f['x2Index']
y1 = f['y1Index']
y2 = f['y2Index']
os.chdir(ARfilepath)
IVTfile = ["era5_6h_daily_0to70S_100Eto120W_025025_ivt_%04d%02d%04d%02d.nc" % (Year, 9, Year+1, 8) for Year in range(1950,2020) ]
ARfile = ["ar_%04d%02d%04d%02d_out.nc" % (Year, 9, Year+1, 8)for Year in range(1950,2020)] 
fi = xr.open_mfdataset(IVTfile,concat_dim='time')
fa = xr.open_mfdataset(ARfile,concat_dim='time')
#%%
### Extract IVT,IVTX,IVTY, compute direction, and AR presence (represented by Shape) for grid point 
def GetARivtForPoint(fi,fa,x,y):
    Shape = fa.shape[:,y,x]
    IVT = fi.ivt[:,y,x]
    IVTX = fi.ivtx[:,y,x]
    IVTY = fi.ivty[:,y,x]
    dire = np.arctan2(IVTX,IVTY)*180/np.pi
    return IVT, Shape,dire   

def GetAllARdata(fa,Point1,Point2,Point3,Point4):
    Datetime = fa.time 
    ### Try to convert xarray Dataarray to make up a dataframe, but it is too slow, size of each array is 10,2272 
    ivt1 = Point1[0].to_masked_array()
    ivt2 = Point2[0].to_masked_array()
    ivt3 = Point3[0].to_masked_array()
    ivt4 = Point4[0].to_masked_array()
    dire1 = Point1[2].to_masked_array()
    dire2 = Point2[2].to_masked_array()
    dire3 = Point3[2].to_masked_array()
    dire4 = Point4[2].to_masked_array()
    shape1 = Point1[1].to_masked_array()
    shape2 = Point2[1].to_masked_array()
    shape3 = Point3[1].to_masked_array()
    shape4 = Point4[1].to_masked_array()
    d = {'Date':Datetime,'IVT1':ivt1,'IVT2':ivt2,'IVT3': ivt3, 'IVT4':ivt4,
         'Shape1':shape1,'Shape2': shape2, 'Shape3':shape3,'Shape4': shape4,
        'dire1': dire1,'dire2': dire2,'dire3': dire3, 'dire4':dire4}
    Alldata = pd.DataFrame(data=d)
    Alldata.dropna(inplace=True)
    Alldata['MeanIVT']=Alldata[['IVT1','IVT2','IVT3','IVT4']].mean(axis=1)
    Alldata['MeanDir']=Alldata[['dire1','dire2','dire3','dire4']].mean(axis=1)
    Alldata['MeanDirM']=np.where(Alldata.MeanDir<0,Alldata.MeanDir+360,Alldata.MeanDir)
    AllARdata_path = os.path.join(dirname, 'AR_date_arealMean')
    os.chdir(AllARdata_path)
    WriteFile = str(Station_List[i]) + '_AR_date_arealMean.csv'
    Alldata.to_csv(WriteFile,index=False)
    
for i in range(len(Station_List)): #len(Station_List)14,65
    os.chdir(ARfilepath)
Point4[2]))
    Point1 = GetARivtForPoint(fi,fa,x1[i],y1[i])
    Point2 = GetARivtForPoint(fi,fa,x1[i],y2[i])
    Point3 = GetARivtForPoint(fi,fa,x2[i],y1[i])
    Point4 = GetARivtForPoint(fi,fa,x2[i],y2[i])
    All = GetAllARdata(fa,Point1,Point2,Point3,Point4)

    
