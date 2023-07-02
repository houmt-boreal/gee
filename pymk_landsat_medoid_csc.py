from statsmodels.formula.api import ols
from scipy import stats
import datetime
from osgeo import gdal
import numpy as np
import pymannkendall as mk
import pandas as pd
import glob
import os
import fnmatch
gdal.AllRegister()

inDs = gdal.Open(r'D:\work_uoh2021\onedrive_hel\area_medoid\slcoff_use\medoid_comp_area2\EPSG_32635_medoid_comp_2002.tif.tif')
driver = inDs.GetDriver()
rows = inDs.RasterYSize
cols = inDs.RasterXSize

count_year1= np.empty([rows, cols],dtype=np.int8)


p_value2 = np.empty([rows, cols],dtype=np.float16)
slope2 = np.empty([rows, cols],dtype=np.float16)
# rsquared_adj2= np.empty([rows, cols])

p_value3 = np.empty([rows, cols],dtype=np.float16)
slope3 = np.empty([rows, cols],dtype=np.float16)
# rsquared_adj3= np.empty([rows, cols])
#
# p_value4 = np.empty([rows, cols])
# slope4 = np.empty([rows, cols])
# rsquared_adj4= np.empty([rows, cols])

in_directory = r'D:\\work_uoh2021\\onedrive_hel\\area_medoid\\slcoff_use\\medoid_comp_area2'
files = glob.glob(os.path.join(in_directory, '*.tif'))
def last_17chars(x):
    return (x[-12:])
files = sorted(files, key=last_17chars)
#
in_directory2 = r'D:\\work_uoh2021\\onedrive_hel\\area_medoid\\slcoff_use\\medoid_comp_clearcount_area2'
files2 = glob.glob(os.path.join(in_directory2, '*.tif'))
files2 = sorted(files2, key=last_17chars)

ndvi_list=[]
year_list = []
validcount_list=[]


for filea,fileb in zip(files,files2):
    ds = gdal.Open(filea, gdal.GA_ReadOnly)
    dsb = gdal.Open(fileb, gdal.GA_ReadOnly)
    if ds.RasterYSize==rows and ds.RasterXSize==cols and dsb.RasterYSize==rows and dsb.RasterXSize==cols and \
        filea[-12:-8]==fileb[-12:-8]:
        print(filea[-12:-8],fileb[-12:-8])
        ds = gdal.Open(filea, gdal.GA_ReadOnly)
        ary = ds.GetRasterBand(1).ReadAsArray()
        ss = ary[:, :]
        ss= np.float16(ss)
        year=filea[-12:-8]
        year_new=int(year)-1984
        ndvi_list.append(ss)
        year_list.append(year_new)

        dsb = gdal.Open(fileb, gdal.GA_ReadOnly)
        aryb = dsb.GetRasterBand(1).ReadAsArray()
        ssb = aryb[:, :]
        ssb = np.int8(ssb)
        validcount_list.append(ssb)



ndvi_array = np.asarray(ndvi_list)
validcount_array = np.asarray(validcount_list)

for col in range(cols):
    print (col)
    for row in range(rows):
        cndvi=ndvi_array[:, row, col]
        validcount= validcount_array[:, row, col]


        df = pd.DataFrame()
        df['ndvi'] = cndvi
        df['year'] = year_list
        df['valid_amount'] = validcount


        df = df[df['ndvi'].notna()]
        # df = df.loc[df['ndvi'] > 0]
        df = df[df['ndvi'] > 0]
        df = df[df['valid_amount'] > 2]
        # df = df[df['valid_amount'] > 3] #  for different UTM
        # df = df[df['year'] > 15]  # means  since  2000

        df = df.sort_values('year', ascending=True)
        df = df.reset_index(drop=True)



        if df.empty:
            # print('DataFrame is empty!')
            count_year1[row, col] = 0
            
            # p_value1[row, col] = -999
            # slope1[row, col] = -999
            # rsquared_adj1[row, col] = -999

            p_value2[row, col] = -999
            slope2[row, col] = -999
##            rsquared_adj2[row, col] = -999
##
            p_value3[row, col] = -999
            slope3[row, col] = -999
##            rsquared_adj3[row, col] = -999
##
##            p_value4[row, col] = -999
##            slope4[row, col] = -999
##            rsquared_adj4[row, col] = -999
            
        elif (df.groupby(["year"])['ndvi'].max()).mean()<0.1 or len(df['year'].unique()) < 10:
            # https://www.nature.com/articles/s41467-020-18479-5#Sec7
            # p_value1[row, col] = -999
            # slope1[row, col] = -999
            # rsquared_adj1[row, col] = -999
            count_year1[row, col] = len(df['year'].unique())

            p_value2[row, col] = -999
            slope2[row, col] = -999
##            rsquared_adj2[row, col] = -999
##
            p_value3[row, col] = -999
            slope3[row, col] = -999
##            rsquared_adj3[row, col] = -999
##
##            p_value4[row, col] = -999
##            slope4[row, col] = -999
##            rsquared_adj4[row, col] = -999
            
        else:
            df['kndvi'] = np.tanh(df['ndvi'] * df['ndvi'])


            try:
                count_year1[row, col] = len(df['year'].unique())
            except:
                count_year1[row, col] = 0


            # fitted_model1.params['doy']
            # try:
            #   cc=fitted_model1.params['I(doy ** 2)']
            #   print(cc)
            # except:
            #    cc = 4
            #    print(cc)

            # try:
            #   cc=fitted_model1.params['C(landsat_series)[T.LE07]']
            #   print(cc)
            # except:
            #     cc = 5
            #     print(cc)



            # theil sen, mk

            ysen = np.array(df['ndvi'].values.tolist())
            xsen = np.array(df['year'].values.tolist())
            try:
                slope2[row, col]= stats.theilslopes(ysen, xsen)[0]
            except:
                slope2[row, col] = -999
            try:
                p_value2[row, col] = stats.kendalltau(xsen, ysen)[1]
            except:
                p_value2[row, col] = -999

            if np.isnan(p_value2[row, col]):
                slope2[row, col] = -999
                
            if (slope2[row, col] > 1) or (slope2[row, col] < -1):  # for example, row=171,col=104
                slope2[row, col] = -999

            if np.isnan(slope2[row, col]):
                slope2[row, col] = -999


            ysen = np.array(df['kndvi'].values.tolist())
            xsen = np.array(df['year'].values.tolist())
            try:
                slope3[row, col] = stats.theilslopes(ysen, xsen)[0]
            except:
                slope3[row, col] = -999
            try:
                p_value3[row, col] = stats.kendalltau(xsen, ysen)[1]
            except:
                p_value3[row, col] = -999

            if np.isnan(p_value3[row, col]):
                slope3[row, col] = -999

            if (slope3[row, col] > 1) or (slope3[row, col] < -1):  # for example, row=171,col=104
                slope3[row, col] = -999

            if np.isnan(slope3[row, col]):
                slope3[row, col] = -999

        del df # https://stackoverflow.com/questions/39100971/how-do-i-release-memory-used-by-a-pandas-dataframe
        del cndvi
        del validcount

outDs1 = driver.Create(
    'D:/work_uoh2021/onedrive_hel/area_medoid/medoid_mk_result2/medoid_ndvi_ls_year_count.tif', cols,
    rows, 1, gdal.GDT_Int16)
outDs1.SetGeoTransform(inDs.GetGeoTransform())
outDs1.SetProjection(inDs.GetProjection())
# print outDatas[90,270,i]
outBand1 = outDs1.GetRasterBand(1)
outBand1.WriteArray(count_year1[:, :])  # get count
outBand1 = None
outDs1 = None


outDs2 = driver.Create(
    'D:/work_uoh2021/onedrive_hel/area_medoid/medoid_mk_result2/medoid_ndvi_year_mk_pvaluetest.tif', cols,
    rows, 1, gdal.GDT_Float32)
outDs2.SetGeoTransform(inDs.GetGeoTransform())
outDs2.SetProjection(inDs.GetProjection())
# print outDatas[90,270,i]
outBand2 = outDs2.GetRasterBand(1)
outBand2.WriteArray(p_value2[:, :])  # get count
outBand2.SetNoDataValue(-999)
outBand2 = None
outDs2 = None

outDs3 = driver.Create(
    'D:/work_uoh2021/onedrive_hel/area_medoid/medoid_mk_result2/medoid_ndvi_year_sen_trendtest.tif', cols,
    rows, 1, gdal.GDT_Float32)
outDs3.SetGeoTransform(inDs.GetGeoTransform())
outDs3.SetProjection(inDs.GetProjection())
# print outDatas[90,270,i]
outBand3 = outDs3.GetRasterBand(1)
outBand3.WriteArray(slope2[:, :])
# outBand1.SetNoDataValue(-9999)
outBand3.SetNoDataValue(-999)
outBand3 = None
outDs3 = None




outDs2 = driver.Create(
    'D:/work_uoh2021/onedrive_hel/area_medoid/medoid_mk_result2/medoid_kndvi_year_mk_pvaluetest.tif', cols,
    rows, 1, gdal.GDT_Float32)
outDs2.SetGeoTransform(inDs.GetGeoTransform())
outDs2.SetProjection(inDs.GetProjection())
# print outDatas[90,270,i]
outBand2 = outDs2.GetRasterBand(1)
outBand2.WriteArray(p_value3[:, :])  # get count
outBand2.SetNoDataValue(-999)
outBand2 = None
outDs2 = None

outDs3 = driver.Create(
    'D:/work_uoh2021/onedrive_hel/area_medoid/medoid_mk_result2/medoid_kndvi_year_sen_trendtest.tif', cols,
    rows, 1, gdal.GDT_Float32)
outDs3.SetGeoTransform(inDs.GetGeoTransform())
outDs3.SetProjection(inDs.GetProjection())
# print outDatas[90,270,i]
outBand3 = outDs3.GetRasterBand(1)
outBand3.WriteArray(slope3[:, :])
# outBand1.SetNoDataValue(-9999)
outBand3.SetNoDataValue(-999)
outBand3 = None
outDs3 = None



print('ok')
