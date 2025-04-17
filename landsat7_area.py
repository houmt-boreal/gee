
import ee
import numpy as np
import folium
from geetools import batch
import datetime
import geopandas as gpd
import os
import math
from pyproj import CRS
import geemap
from shapely.geometry import shape
ee.Initialize()

def getQABits(image, start, end, mascara):
    # Compute the bits we need to extract.
    pattern = 0
    for i in range(start,end+1):
        pattern += 2**i
    # Return a single band image of the extracted QA bits, giving the     band a new name.
    return image.select([0], [mascara]).bitwiseAnd(pattern).rightShift(start)


# cloudmask option1
def maskClouds(image):
    dilatedCloud = (1 << 1)
    cloud = (1 << 3)
    cloudShadow = (1 << 4)
    qa = image.select('QA_PIXEL')
    fill = (1 << 0)
    maskfill = qa.bitwiseAnd(fill).eq(0).And(qa.bitwiseAnd(dilatedCloud).eq(0))
    # Both flags should be set to zero, indicating clear conditions.
    mask = qa.bitwiseAnd(cloud).eq(0)\
        .And(qa.bitwiseAnd(cloudShadow).eq(0))
    mask = mask.Not().focal_max(radius=150, units='meters').Not()

    saturationMask = image.select('QA_RADSAT').eq(0)

    Cloud_Confidence = getQABits(qa, 8, 9, 'Cloud_Confid')
    Cloud_Shadow_Confidence = getQABits(qa, 10, 11, 'Cloud_Shadow_Confid')
    Cirrus_Confidence = getQABits(qa, 14, 15, 'Cirrus_Confid')
    maskconfid1 = (Cloud_Confidence.gte(2)).focal_max(radius=0, units='meters')
    # maskconfid2 = (Cloud_Confidence.gte(2).And(Cloud_Shadow_Confidence.gte(2))).focal_max(radius=500, units='meters')
    maskconfid2 = (Cloud_Shadow_Confidence.gte(2)).focal_max(radius=0, units='meters')
    maskconfid3 = (Cirrus_Confidence.gte(2)).focal_max(radius=0, units='meters')
    # return image.updateMask(mask).updateMask(anySaturated.eq(0))
    return image.updateMask(mask).updateMask(maskfill).updateMask(saturationMask).updateMask(maskconfid1.Not()).updateMask(
       maskconfid2.Not()).updateMask(maskconfid3.Not())

def snowflag(image):
    snow = (1 << 5)
    qa = image.select('QA_PIXEL')
    masksnow=qa.bitwiseAnd(snow).eq(0)
    ice_snowConfidence = getQABits(qa, 12, 13, 'ice_snow_Confid')
    return image.updateMask(masksnow.Not().focal_max(radius=30, units='meters').Not())\
        .updateMask((ice_snowConfidence.gte(2)).focal_max(radius=30, units='meters').Not())

def scaleBands(img):
    scaling = img.select('SR_B[1-7]')
    x = scaling.multiply(0.0000275).add(-0.2)
    scaling = img.select('ST_B6')
    scaling = scaling.multiply(0.00341802).add(149.0)
    x = x.addBands(scaling)
    notScaling = img.select(['QA_PIXEL','QA_RADSAT'])
    return x.addBands(notScaling)

def valuefilter(img):
    mask_nir=img.select('SR_B4').gt(0.005).And(img.select('SR_B4').lt(1))
    mask_red = img.select('SR_B3').gt(0.005).And(img.select('SR_B3').lt(1))
    mask_b = img.select('SR_B1').gt(0.005).And(img.select('SR_B1').lt(1))
    mask_g = img.select('SR_B2').gt(0.005).And(img.select('SR_B2').lt(1))
    mask_swir = img.select('SR_B5').gt(0).And(img.select('SR_B5').lt(1))
    return img.updateMask(mask_nir).updateMask(mask_red).updateMask(mask_g).updateMask(mask_b).updateMask(mask_swir)


def addIndices(img):
    NDVI = img.normalizedDifference(['SR_B4','SR_B3']).rename('NDVI')
    imgDict = {
    'N': img.select('SR_B4'),
    'R': img.select('SR_B3'),
    'B': img.select('SR_B1'),
    'G': img.select('SR_B2'),
    'SWIR1': img.select('SR_B5')
    }
    formula = '2.5 * (N - R) / (N + 6.0 * R - 7.5 * B + 1.0)'
    formula3 = '(G - SWIR1) / (G + SWIR1 )'

    NDSI = img.expression(formula3, imgDict).rename('NDSI')
    #formula2 = '((N - R) / (N + R ))* ((N - R) / (N + R))'
    EVI = img.expression(formula,imgDict).rename('EVI')
    # GNDVI = img.normalizedDifference(['SR_B','SR_B']).rename('GNDVI')
    #NDVI_square = img.expression(formula2,imgDict).rename('NDVI_square')
    return img.addBands([NDVI,EVI,NDSI])
   # return img.addBands([NDVI,EVI,GNDVI]).unmask(-2)


def ndvivalue(img):
    mask_ndvivalue = img.select('NDVI').gt(0).And(img.select('NDVI').lt(1))
    return img.updateMask(mask_ndvivalue)

# area = ee.Geometry.Rectangle([[20, 68.0],[27, 71]])
# area = ee.Geometry.Point([256360.8,7661089.7],'EPSG:3067')
area = ee.Geometry.Rectangle([[26.07, 69.94],[26.35,70.06]]) # rastigasa   70°N, 26.2°E
# area = ee.Geometry.Rectangle([[20.85, 68.9],[21.15,69.1]]) # bm1 kilpisyavi 69°N, 21°E
#area = ee.Geometry.Polygon([[464463.7, 7761073.5],[474463.7, 7761073.5],[474463.7, 7771073.5],[464463.7, 7771073.5],[464463.7, 7761073.5]], proj = "EPSG:4088")

Lt_collection = (ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
               .filterDate('1999-01-01', '2021-12-31')
               #.filterDate('2011-08-19', '2011-08-20')
               # .filter(ee.Filter.calendarRange(j, field = "month"))
               .filterBounds(area)
               .filterMetadata('IMAGE_QUALITY', 'equals', 9)
               .filterMetadata('SUN_ELEVATION', 'greater_than', 30)  # https://developers.google.com/earth-engine/apidocs/ee-imagecollection-filtermetadata
               #.filterMetadata('CLOUD_COVER', 'less_than', 80)
               # .filter(ee.Filter.neq(property name, None))
               #.filterMetadata('GEOMETRIC_RMSE_MODEL', 'less_than', 10)
               #.filterMetadata('GEOMETRIC_RMSE_MODEL', 'not_less_than', 10)
               .filterMetadata('GEOMETRIC_RMSE_MODEL', 'not_greater_than', 30)
               .map(maskClouds)
               #.map(snowflag)
               .map(scaleBands)
               .map(valuefilter)
               .map(addIndices)
               #.map(ndvivalue)  # carry on later
               #.select('NDVI'))
                 .select(['SR_B4','SR_B3','NDSI','NDVI']))


# Get the number of images of collection
count = Lt_collection.size().getInfo()
print(count)
imageList = Lt_collection.toList(count)

for i in range(count):
  img = ee.Image(imageList.get(i))
  crs_raw = img.projection().getInfo()['crs']
  if crs_raw == 'EPSG:32635':
      theMax = img.select('NDVI').reduceRegion(ee.Reducer.max(), area)
      if theMax.getInfo().get('NDVI') is not None:  # check if theMax.getInfo().get('NDVI') is none
          crs_raw = crs_raw.replace(":", "_")

          imageraw = ee.Image('LANDSAT/LE07/C02/T1_L2/' + img.getInfo()['properties']['system:index'])
          gsw = ee.Image('JRC/GSW1_3/GlobalSurfaceWater')
          # print (imageraw.geometry().getInfo())
          extentGSW = gsw.select('occurrence').gte(50).clip(imageraw.geometry()).unmask(0)
          #extentGSW = gsw.select('max_extent').clip(area).unmask(0)
          imagewater = extentGSW.reproject(imageraw.projection())
          #print (imagewater.getInfo())

          watermask=imagewater.select('occurrence').neq(1)
          id_snow = img.select('NDSI').gt(0.4).And(imagewater.select('occurrence').neq(1))
          snowmask=id_snow.focal_max(radius=60, units='meters').Not()
          img = img.select(['SR_B4','SR_B3','NDVI']).updateMask(watermask).updateMask(snowmask)


          mask_ndvivalue1 = img.select('NDVI').gt(0).And(img.select('NDVI').lt(1))
          img = img.updateMask(mask_ndvivalue1)
          theMax1 = img.select('NDVI').reduceRegion(ee.Reducer.max(), area)
          if theMax1.getInfo().get('NDVI') is not None:

      # export  cloud mask
          #output = ee.Image(0).where(maskClouds(imageraw).select('SR_B3'), ee.Image(1)).reproject(imageraw.projection())

          # rgb_img = geemap.ee_to_numpy(img.select('SR_B3'), region=area)
          # print(rgb_img.shape)

              # countpixels = img.select('NDVI').reduceRegion(ee.Reducer.count(), area, scale=30).getInfo()
              # print(countpixels['NDVI'])
              # if countpixels['NDVI'] > 32743:  # 32743= 363 col * 451 row * 20%
                      filename = crs_raw + '_ndvi_buffer_ndsi04_60m_cloud150_' + img.getInfo()['properties']['system:index'] + '.tif'
                      geemap.ee_export_image_to_drive(img.select('NDVI'), description=filename, folder='north_l7_ndvi_cloudbuffer150m_ndsi04_60m_2022',region=area,scale=30)


