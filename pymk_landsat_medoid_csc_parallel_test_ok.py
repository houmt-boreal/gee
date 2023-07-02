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
from concurrent.futures import ProcessPoolExecutor
from joblib import Parallel, delayed
import itertools
import multiprocessing
gdal.AllRegister()

inDs = gdal.Open(r'D:\work_uoh2021\onedrive_hel\area_medoid\slcoff_use\medoid_comp_area2\EPSG_32635_medoid_comp_2002.tif.tif')
driver = inDs.GetDriver()
rows = inDs.RasterYSize
cols = inDs.RasterXSize

count_year1= np.zeros((rows, cols), dtype=np.int8)


p_value2 = np.zeros((rows, cols), dtype=np.float16)
slope2 = np.zeros((rows, cols), dtype=np.float16)
# rsquared_adj2= np.empty([rows, cols])

p_value3 = np.zeros((rows, cols), dtype=np.float16)
slope3 = np.zeros((rows, cols), dtype=np.float16)
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


def process(col, row):


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
                r11 = 0

                r12 = -999
                r13 = -999

                r14 = -999
                r15 = -999


            elif (df.groupby(["year"])['ndvi'].max()).mean()<0.1 or len(df['year'].unique()) < 10:

                r11 = len(df['year'].unique())

                r12 = -999
                r13 = -999

                r14 = -999
                r15 = -999


            else:
                df['kndvi'] = np.tanh(df['ndvi'] * df['ndvi'])
                try:
                    r11= len(df['year'].unique())
                except:
                    r11 = 0



                # theil sen, mk

                ysen = np.array(df['ndvi'].values.tolist())
                xsen = np.array(df['year'].values.tolist())
                try:
                    r13= stats.theilslopes(ysen, xsen)[0]
                except:
                    r13 = -999
                try:
                    r12 = stats.kendalltau(xsen, ysen)[1]
                except:
                    r12 = -999

                if np.isnan(p_value2[row, col]):
                    r13 = -999

                if (slope2[row, col] > 1) or (slope2[row, col] < -1):  # for example, row=171,col=104
                    r13 = -999

                if np.isnan(slope2[row, col]):
                    r13 = -999


                ysen = np.array(df['kndvi'].values.tolist())
                xsen = np.array(df['year'].values.tolist())
                try:
                    r15= stats.theilslopes(ysen, xsen)[0]
                except:
                    r15 = -999
                try:
                    r14 = stats.kendalltau(xsen, ysen)[1]
                except:
                    r14= -999

                if np.isnan(p_value3[row, col]):
                    r14 = -999

                if (slope3[row, col] > 1) or (slope3[row, col] < -1):  # for example, row=171,col=104
                    r14 = -999

                if np.isnan(slope3[row, col]):
                    r14 = -999

            # del df # https://stackoverflow.com/questions/39100971/how-do-i-release-memory-used-by-a-pandas-dataframe
            # del cndvi
            # del validcount
            # return count_year1,p_value2,slope2,p_value3,slope3

            # a.append(r11)
            # b.append(r12)
            # c.append(r13)
            # d.append(r14)
            # e.append(r15)

            return r11,r12,r13,r14,r15

num_cores = multiprocessing.cpu_count()
results = Parallel(n_jobs=num_cores)(delayed(process)(col, row) for col, row in itertools.product(range(cols), range(rows)))
#print (results[3][2])
print (len(results))

a=[]
b=[]
c=[]
d=[]
e=[]
for i in range(len(results)):
   s1=results[i][0]
   s2 = results[i][1]
   s3 = results[i][2]
   s4 = results[i][3]
   s5 = results[i][4]
   a.append(s1)
   b.append(s2)
   c.append(s3)
   d.append(s4)
   e.append(s5)

count_year1 = np.reshape(a, (-1, cols), order='F')   #  new = np.reshape(a, (-1, ncols))   # convert-a-1d-array-to-a-2d-array-in-numpy
p_value2 = np.reshape(b, (-1, cols), order='F')  # https://numpy.org/doc/stable/reference/generated/numpy.reshape.html
slope2 = np.reshape(c, (-1, cols), order='F')
p_value3 = np.reshape(d, (-1, cols), order='F')
slope3 = np.reshape(e, (-1, cols), order='F')

outDs1 = driver.Create(
    'D:/work_uoh2021/onedrive_hel/area_medoid/medoid_mk_result/medoid_ndvi_ls_year_count2.tif', cols,
    rows, 1, gdal.GDT_Int16)
outDs1.SetGeoTransform(inDs.GetGeoTransform())
outDs1.SetProjection(inDs.GetProjection())
# print outDatas[90,270,i]
outBand1 = outDs1.GetRasterBand(1)
outBand1.WriteArray(count_year1[:, :])  # get count
outBand1 = None
outDs1 = None

# outDs2 = driver.Create(
#     'D:/work_uoh2021/onedrive_hel/area_medoid/medoid_mk_result/medoid_ndvi_year_lsfit_pvalue.tif', cols,
#     rows, 1, gdal.GDT_Float32)
# outDs2.SetGeoTransform(inDs.GetGeoTransform())
# outDs2.SetProjection(inDs.GetProjection())
# # print outDatas[90,270,i]
# outBand2 = outDs2.GetRasterBand(1)
# outBand2.WriteArray(p_value1[:, :])  # get count
# outBand2.SetNoDataValue(-999)
# outBand2 = None
# outDs2 = None
#
# outDs3 = driver.Create(
#     'D:/work_uoh2021/onedrive_hel/area_medoid/medoid_mk_result/medoid_ndvi_year_lsfit_trend.tif', cols,
#     rows, 1, gdal.GDT_Float32)
# outDs3.SetGeoTransform(inDs.GetGeoTransform())
# outDs3.SetProjection(inDs.GetProjection())
# # print outDatas[90,270,i]
# outBand3 = outDs3.GetRasterBand(1)
# outBand3.WriteArray(slope1[:, :])
# # outBand1.SetNoDataValue(-9999)
# outBand3.SetNoDataValue(-999)
# outBand3 = None
# outDs3 = None
#
# outDs3 = driver.Create(
#      'D:/work_uoh2021/onedrive_hel/area_medoid/medoid_mk_result/medoid_ndvi_year_lsfit_rsquared_adj.tif', cols,
#     rows, 1, gdal.GDT_Float32)
# outDs3.SetGeoTransform(inDs.GetGeoTransform())
# outDs3.SetProjection(inDs.GetProjection())
# # print outDatas[90,270,i]
# outBand3 = outDs3.GetRasterBand(1)
# outBand3.WriteArray(rsquared_adj1[:, :])
# # outBand1.SetNoDataValue(-9999)
# outBand3.SetNoDataValue(-999)
# outBand3 = None
# outDs3 = None


outDs2 = driver.Create(
    'D:/work_uoh2021/onedrive_hel/area_medoid/medoid_mk_result/medoid_ndvi_year_mk_pvaluetest2.tif', cols,
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
    'D:/work_uoh2021/onedrive_hel/area_medoid/medoid_mk_result/medoid_ndvi_year_sen_trendtest2.tif', cols,
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
    'D:/work_uoh2021/onedrive_hel/area_medoid/medoid_mk_result/medoid_kndvi_year_mk_pvaluetest2.tif', cols,
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
    'D:/work_uoh2021/onedrive_hel/area_medoid/medoid_mk_result/medoid_kndvi_year_sen_trendtest2.tif', cols,
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



##outDs3 = driver.Create(
##    '/scratch/project_2005257/buffer_output/buffer_output_ndsi04_me_corr_median78/ndvi_ls_corrected_year_rsquared_adj_rastigaisa_ndsi04.tif', cols,
##    rows, 1, gdal.GDT_Float32)
##outDs3.SetGeoTransform(inDs.GetGeoTransform())
##outDs3.SetProjection(inDs.GetProjection())
### print outDatas[90,270,i]
##outBand3 = outDs3.GetRasterBand(1)
##outBand3.WriteArray(rsquared_adj2[:, :])
### outBand1.SetNoDataValue(-9999)
##outBand3.SetNoDataValue(-999)
##outBand3 = None
##outDs3 = None

##
##outDs2 = driver.Create(
##    '/scratch/project_2005257/buffer_output/buffer_output_ndsi04_me_corr_median78/kndvi_ls_year_pvalue_rastigaisa_ndsi04.tif', cols,
##    rows, 1, gdal.GDT_Float32)
##outDs2.SetGeoTransform(inDs.GetGeoTransform())
##outDs2.SetProjection(inDs.GetProjection())
### print outDatas[90,270,i]
##outBand2 = outDs2.GetRasterBand(1)
##outBand2.WriteArray(p_value3[:, :])  # get count
##outBand2.SetNoDataValue(-999)
##outBand2 = None
##outDs2 = None
##
##outDs3 = driver.Create(
##    '/scratch/project_2005257/buffer_output/buffer_output_ndsi04_me_corr_median78/kndvi_ls_year_trend_rastigaisa_ndsi04.tif', cols,
##    rows, 1, gdal.GDT_Float32)
##outDs3.SetGeoTransform(inDs.GetGeoTransform())
##outDs3.SetProjection(inDs.GetProjection())
### print outDatas[90,270,i]
##outBand3 = outDs3.GetRasterBand(1)
##outBand3.WriteArray(slope3[:, :])
### outBand1.SetNoDataValue(-9999)
##outBand3.SetNoDataValue(-999)
##outBand3 = None
##outDs3 = None
##
##outDs3 = driver.Create(
##    '/scratch/project_2005257/buffer_output/buffer_output_ndsi04_me_corr_median78/kndvi_ls_year_rsquared_adj_rastigaisa_ndsi04.tif', cols,
##    rows, 1, gdal.GDT_Float32)
##outDs3.SetGeoTransform(inDs.GetGeoTransform())
##outDs3.SetProjection(inDs.GetProjection())
### print outDatas[90,270,i]
##outBand3 = outDs3.GetRasterBand(1)
##outBand3.WriteArray(rsquared_adj3[:, :])
### outBand1.SetNoDataValue(-9999)
##outBand3.SetNoDataValue(-999)
##outBand3 = None
##outDs3 = None
##
##outDs2 = driver.Create(
##    '/scratch/project_2005257/buffer_output/buffer_output_ndsi04_me_corr_median78/kndvi_ls_corrected_year_pvalue_rastigaisa_ndsi04.tif', cols,
##    rows, 1, gdal.GDT_Float32)
##outDs2.SetGeoTransform(inDs.GetGeoTransform())
##outDs2.SetProjection(inDs.GetProjection())
### print outDatas[90,270,i]
##outBand2 = outDs2.GetRasterBand(1)
##outBand2.WriteArray(p_value4[:, :])  # get count
##outBand2.SetNoDataValue(-999)
##outBand2 = None
##outDs2 = None
##
##outDs3 = driver.Create(
##    '/scratch/project_2005257/buffer_output/buffer_output_ndsi04_me_corr_median78/kndvi_ls_corrected_year_trend_rastigaisa_ndsi04.tif', cols,
##    rows, 1, gdal.GDT_Float32)
##outDs3.SetGeoTransform(inDs.GetGeoTransform())
##outDs3.SetProjection(inDs.GetProjection())
### print outDatas[90,270,i]
##outBand3 = outDs3.GetRasterBand(1)
##outBand3.WriteArray(slope4[:, :])
### outBand1.SetNoDataValue(-9999)
##outBand3.SetNoDataValue(-999)
##outBand3 = None
##outDs3 = None
##
##outDs3 = driver.Create(
##    '/scratch/project_2005257/buffer_output/buffer_output_ndsi04_me_corr_median78/kndvi_ls_corrected_year_rsquared_adj_rastigaisa_ndsi04.tif', cols,
##    rows, 1, gdal.GDT_Float32)
##outDs3.SetGeoTransform(inDs.GetGeoTransform())
##outDs3.SetProjection(inDs.GetProjection())
### print outDatas[90,270,i]
##outBand3 = outDs3.GetRasterBand(1)
##outBand3.WriteArray(rsquared_adj4[:, :])
### outBand1.SetNoDataValue(-9999)
##outBand3.SetNoDataValue(-999)
##outBand3 = None
##outDs3 = None

print('ok')
