
import ee
import datetime
import pandas as pd
import time
import geemap

ee.Initialize()

aoi = ee.Geometry.Rectangle([[89.6, 47.2], [90.6, 47.5]])  #

startYear = 2020
endYear = 2021
startDay = '07-01'
endDay = '08-31'


def addNDVI(img):
    ndvi = img.expression(
        '(NIR - RED) / (NIR + RED)', {
            'NIR': img.select('SR_B4'),
            'RED': img.select('SR_B3')
        }).select([0], ['NDVI']).set('system:time_start', img.get('system:time_start'))
    return ndvi


runParams = {
    'maxSegments': 6,
    'spikeThreshold': 0.9,
    'vertexCountOvershoot': 3,
    'preventOneYearRecovery': True,
    'recoveryThreshold': 0.25,
    'pvalThreshold': 0.05,
    'bestModelProportion': 0.75,
    'minObservationsNeeded': 6
}


def getQABits(image, start, end, mascara):
    # Compute the bits we need to extract.
    pattern = 0
    for i in range(start, end + 1):
        pattern += 2 ** i
    # Return a single band image of the extracted QA bits, giving the     band a new name.
    return image.select([0], [mascara]).bitwiseAnd(pattern).rightShift(start)


def scaleBands(img):
    scaling = img.select('SR_B[1-7]')
    x = scaling.multiply(0.0000275).add(
        -0.2).float()  # need to cast raw bands to float to make sure that we don't get errors regarding incompatible bands
    # scaling = img.select('ST_B6')
    # scaling = scaling.multiply(0.00341802).add(149.0)
    # x = x.addBands(scaling)
    notScaling = img.select(['QA_PIXEL', 'QA_RADSAT'])
    return x.addBands(notScaling)


def valuefilter(img):
    mask_nir = img.select('SR_B4').gt(0.005).And(img.select('SR_B4').lt(1))
    mask_red = img.select('SR_B3').gt(0.005).And(img.select('SR_B3').lt(1))
    mask_b = img.select('SR_B1').gt(0.005).And(img.select('SR_B1').lt(1))
    mask_g = img.select('SR_B2').gt(0.005).And(img.select('SR_B2').lt(1))
    mask_swir = img.select('SR_B5').gt(0).And(img.select('SR_B5').lt(1))
    return img.updateMask(mask_nir).updateMask(mask_red).updateMask(mask_g).updateMask(mask_b).updateMask(mask_swir)


def valuefilter8(img):
    mask_nir = img.select('SR_B5').gt(0.005).And(img.select('SR_B5').lt(1))
    mask_red = img.select('SR_B4').gt(0.005).And(img.select('SR_B4').lt(1))
    mask_b = img.select('SR_B2').gt(0.005).And(img.select('SR_B2').lt(1))
    mask_g = img.select('SR_B3').gt(0.005).And(img.select('SR_B3').lt(1))
    # above reference  https://github.com/logan-berner/arctic_greening/blob/master/code/0.1_fun_lsat_tools.R
    mask_swir = img.select('SR_B6').gt(0).And(img.select('SR_B6').lt(1))
    return img.updateMask(mask_nir).updateMask(mask_red).updateMask(mask_g).updateMask(mask_b).updateMask(mask_swir)


#  make an image collection from an image with 6 bands all set to 0 and then make them masked values
dummyCollection = ee.ImageCollection([ee.Image([0, 0, 0, 0, 0, 0]).mask(ee.Image(0))])


## slope and intercept citation: Roy, D.P., Kovalskyy, V., Zhang, H.K., Vermote, E.F., Yan, L., Kumar, S.S, Egorov, A., 2016, Characterization of Landsat-7 to Landsat-8 reflective wavelength and normalized difference vegetation index continuity, Remote Sensing of Environment, 185, 57-70.(http:##dx.doi.org/10.1016/j.rse.2015.12.024) Table 2 - reduced major axis (RMA) regression coefficients
def harmonizationRoy(oli):
    slopes = ee.Image.constant([0.9785, 0.9542, 0.9825, 1.0073, 1.0171, 0.9949])
    itcp = ee.Image.constant([-0.0095, -0.0016, -0.0022, -0.0021, -0.0030, 0.0029])
    # for TOA refl https://www.sciencedirect.com/science/article/pii/S0034425718305212#bb0225
    return oli.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'],
                      ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']).subtract(
        itcp).divide(slopes).set('system:time_start', oli.get('system:time_start'))


def getSRcollection(year, startDay, endDay, sensor, aoi):
    # print (sensor)
    if sensor == 'LC08':
        srCollection = (ee.ImageCollection('LANDSAT/' + sensor + '/C02/T1_L2')
                        .filterBounds(aoi)
                        .filterDate(str(year) + '-' + startDay, str(year) + '-' + endDay)
                        .filterMetadata('IMAGE_QUALITY_OLI', 'equals', 9)
                        .filterMetadata('CLOUD_COVER', 'less_than', 80)
                        .filterMetadata('SUN_ELEVATION', 'greater_than', 30)
                        .filterMetadata('GEOMETRIC_RMSE_MODEL', 'not_greater_than', 30)
                        .map(scaleBands).map(valuefilter8))
    # elif sensor == 'LE07':
    #     srCollection = (ee.ImageCollection('LANDSAT/' + sensor + '/C02/T1_L2')
    #                     .filterBounds(aoi)
    #                     .filterDate(str(year) + '-' + startDay, str(year) + '-' + endDay)
    #                     .filter(ee.Filter.date('2003-05-31', '2021-12-31').Not())
    #                     .filterMetadata('IMAGE_QUALITY', 'equals', 9)
    #                     .filterMetadata('CLOUD_COVER', 'less_than', 80)
    #                     .filterMetadata('SUN_ELEVATION', 'greater_than', 30)
    #                     .filterMetadata('GEOMETRIC_RMSE_MODEL', 'not_greater_than', 30)
    #                     .map(scaleBands).map(valuefilter))
    else:
        # print (sensor)
        # print('ok')
        srCollection = (ee.ImageCollection('LANDSAT/' + sensor + '/C02/T1_L2')
                        .filterBounds(aoi)
                        .filterDate(str(year) + '-' + startDay, str(year) + '-' + endDay)
                        .filterMetadata('IMAGE_QUALITY', 'equals', 9)
                        .filterMetadata('CLOUD_COVER', 'less_than', 80)
                        .filterMetadata('SUN_ELEVATION', 'greater_than', 30)
                        .filterMetadata('GEOMETRIC_RMSE_MODEL', 'not_greater_than', 30)
                        # https://courses.spatialthoughts.com/end-to-end-gee.html
                        .filter(ee.Filter.eq('system:index', 'LT05_194011_20040723').Not())  # error geometric
                        .map(scaleBands).map(valuefilter))

    crslist = []
    for ik in range(srCollection.size().getInfo()):
        aa = (ee.Image(srCollection.toList(srCollection.size().getInfo()).get(ik))).getInfo()['properties'][
            'system:index']
        aa2 = (ee.Image(srCollection.toList(srCollection.size().getInfo()).get(ik))).projection().getInfo()['crs']
        # print(aa,aa2)
        if aa2 == 'EPSG:32645':
            crslist.append(aa)
    srCollection = (srCollection.filter(ee.Filter.inList('system:index', crslist)))
    print(crslist, srCollection.size().getInfo())


    def prepImages(img):
        dat = ee.Image(
            ee.Algorithms.If(
                sensor == 'LC08',
                harmonizationRoy(img.unmask()),
                img.select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']) \
                .unmask() \
                # .resample('bicubic') \
                .set('system:time_start', img.get('system:time_start'))
            )
        )

        qa = img.select('QA_PIXEL')
        # print (sensor)

        if sensor == 'LC08':

            dilatedCloud = (1 << 1)
            cirrus = (1 << 2)
            cloud = (1 << 3)
            cloudShadow = (1 << 4)
            fill = (1 << 0)
            snow = (1 << 5)
            water = (1 << 7)

            mask = qa.bitwiseAnd(fill).eq(0) \
                .And(qa.bitwiseAnd(dilatedCloud).eq(0)) \
                .And(qa.bitwiseAnd(cirrus).eq(0)) \
                .And(qa.bitwiseAnd(cloud).eq(0)) \
                .And(qa.bitwiseAnd(cloudShadow).eq(0)) \
                .And(qa.bitwiseAnd(snow).eq(0)) \
                .And(qa.bitwiseAnd(water).eq(0))

        else:
            dilatedCloud = (1 << 1)
            cloud = (1 << 3)
            cloudShadow = (1 << 4)
            fill = (1 << 0)
            snow = (1 << 5)
            water = (1 << 7)

            mask = qa.bitwiseAnd(fill).eq(0) \
                .And(qa.bitwiseAnd(dilatedCloud).eq(0)) \
                .And(qa.bitwiseAnd(cloud).eq(0)) \
                .And(qa.bitwiseAnd(cloudShadow).eq(0)) \
                .And(qa.bitwiseAnd(snow).eq(0)) \
                .And(qa.bitwiseAnd(water).eq(0))

        saturationMask = img.select('QA_RADSAT').eq(0)

        Cloud_Confidence = getQABits(qa, 8, 9, 'Cloud_Confid').gte(2).Not()
        Cloud_Shadow_Confidence = getQABits(qa, 10, 11, 'Cloud_Shadow_Confid').gte(2).Not()
        Cirrus_Confidence = getQABits(qa, 14, 15, 'Cirrus_Confid').gte(2).Not()
        Snow_Ice_Confidence = getQABits(qa, 12, 13, 'snow_ice_Confid').gte(2).Not()

        MappedWater = ee.Image("JRC/GSW1_3/GlobalSurfaceWater")
        MappedWaterBinary = MappedWater.expression("(b('occurrence') > 50) ? 0"
                                                   ": 1").clip(aoi).reproject(img.projection())

        # https://developers.google.com/earth-engine/guides/image_relational#colab-python_3
        # zones_exp = nl_2012.expression("(b('stable_lights') > 62) ? 3 "
        #                                ": (b('stable_lights') > 55) ? 2 "
        #                                ": (b('stable_lights') > 30) ? 1 "
        #                                ": 0")
        # https://gis.stackexchange.com/questions/310564/what-is-the-difference-between-the-mask-and-updatemask
        return dat.updateMask(mask).updateMask(saturationMask).updateMask(MappedWaterBinary).updateMask(
            Cloud_Confidence).updateMask(
            Cloud_Shadow_Confidence).updateMask(Cirrus_Confidence).updateMask(Snow_Ice_Confidence)

    return srCollection.map(prepImages)


def getCombinedSRcollection(year, startDay, endDay, aoi):
    lt5 = getSRcollection(year, startDay, endDay, 'LT05', aoi)
    le7 = getSRcollection(year, startDay, endDay, 'LE07', aoi)
    lc8 = getSRcollection(year, startDay, endDay, 'LC08', aoi)
    return ee.ImageCollection(lt5.merge(le7).merge(lc8))


def medoidMosaic(inCollection, dummyCollection):
    imageCount = inCollection.toList(1).length()  # get the number of images
    # if the number of images in this year is 0, then use the dummy collection, otherwise use the SR collection
    finalCollection = ee.ImageCollection(ee.Algorithms.If(imageCount.gt(0), inCollection, dummyCollection))
    # print((ee.Image(finalCollection.toList(finalCollection.size().getInfo()).get(1))).projection().getInfo()['crs'])
    # output: EPSG:32645

    if finalCollection.size().getInfo() > 0:
        img_proj_based = finalCollection.first()
        # print(img_proj_based.getInfo())  # reproject(img_proj_based.projection())

    #  calculate the median of the annual image collection - returns a single 6 band image - the collection median per band
    median = finalCollection.median().reproject(img_proj_based.projection())

    # print(median.projection().getInfo()['crs'])    # output: EPSG:4326 without .reproject

    # calculate the different between the median and the observation per image per band
    def calcDifFromMed(img):
        # get the difference between each image/band and the corresponding band median and take to power of 2 to make negatives positive and make greater differences weight more
        diff = ee.Image(img).subtract(median).pow(ee.Image.constant(2))
        # per image in collection, sum the powered difference across the bands - set this as the first band add the SR bands to it - now a 7 band image collection
        return diff.reduce('sum').addBands(img)

    difFromMedian = finalCollection.map(calcDifFromMed)
    # print((ee.Image(difFromMedian.toList(difFromMedian.size().getInfo()).get(1))).projection().getInfo()['crs'])
    # output: EPSG:32645
    # print((difFromMedian.size().getInfo()))

    # get the medoid by selecting the image pixel with the smallest difference between median and observation per band
    return ee.ImageCollection(difFromMedian).reduce(ee.Reducer.min(7)).select([1, 2, 3, 4, 5, 6],
                                                                              ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4',
                                                                               'SR_B5',
                                                                               'SR_B7'])#.reproject(img_proj_based.projection())


# if   mosaicType == "medoid"
# FUNCTION TO APPLY MEDOID COMPOSITING FUNCTION TO A COLLECTION
def buildMosaic(year, startDay, endDay, aoi, dummyCollection):
    collection = getCombinedSRcollection(year, startDay, endDay, aoi)
    # apply the medoidMosaic function to reduce the collection to single image per year by medoid
    # add the year to each medoid image - the data is hard-coded Aug 1st
    img = medoidMosaic(collection, dummyCollection).set('system:time_start',
                                                        int(time.mktime(datetime.date(year, 8, 1).timetuple())))
    return ee.Image(img)


# FUNCTION TO BUILD ANNUAL MOSAIC COLLECTION
def buildMosaicCollection(startYear, endYear, startDay, endDay, aoi, dummyCollection):
    imgs = []
    for year in range(startYear, endYear + 1):
        # build the medoid mosaic for a given year
        tmp = buildMosaic(year, startDay, endDay, aoi, dummyCollection)
        # concatenate the annual image medoid to the collection (img) and set the date of the image - hard coded to the year that is being worked on for Aug 1s
        imgs.append(tmp.set('system:time_start', int(time.mktime(datetime.date(year, 8, 1).timetuple()))))
    return ee.ImageCollection(imgs)


# FUNCTION TO COUNT NUMBER OF UNMASKED PIXELS IN AN INTRA ANNUAL COLLECTION
def countClearViewPixels(intraAnnualSRcollection):
    def calbinary(img):
        return img.select(0).multiply(0).add(1).unmask(0)

    binary = intraAnnualSRcollection.map(calbinary)
    return binary.sum()


# FUNCTION TO BUILD ANNAUL COLLECTION OF NUMBER OF UNMASKED PIXELS AVAILABLE TO BUILD COMPOSITE
def buildClearPixelCountCollection(startYear, endYear, startDay, endDay, aoi, dummyCollection):
    imgs = []
    for year in range(startYear, endYear + 1):
        collection = getCombinedSRcollection(year, startDay, endDay, aoi)
        imageCount = collection.toList(1).length()  # get the number of images
        # if the number of images in this year is 0, then use the dummy collection, otherwise use the SR collection
        finalCollection = ee.ImageCollection(ee.Algorithms.If(imageCount.gt(0), collection, dummyCollection))
        notMaskCount = countClearViewPixels(finalCollection)
        imgs.append(notMaskCount.set('system:time_start', int(time.mktime(datetime.date(year, 8, 1).timetuple()))))
    return ee.ImageCollection(imgs)


# # if   mosaicType == "medoid"
# build annual image collection and run LandTrendr
annualSRcollection = buildMosaicCollection(startYear, endYear, startDay, endDay, aoi, dummyCollection)
annualNDVIcollection = (annualSRcollection.map(addNDVI).select('NDVI'))
count = annualNDVIcollection.size().getInfo()
print(count)
imageList = annualNDVIcollection.toList(count)

annualclerapixelcollection = buildClearPixelCountCollection(startYear, endYear, startDay, endDay, aoi, dummyCollection)
count2 = annualclerapixelcollection.size().getInfo()
print(count2)
imageList2 = annualclerapixelcollection.toList(count2)

if count == count2:
    for i in range(count):
        img = ee.Image(imageList.get(i))#.reproject('EPSG:4326', None, 30)
        crs_raw = img.projection().getInfo()['crs']
        aa = img.getInfo()['properties']['system:index']

        img_clear = ee.Image(imageList2.get(i))  # default epsg4326
        aa2 = img_clear.getInfo()['properties']['system:index']

        mask_ndvivalue1 = img.select('NDVI').gt(0).And(img.select('NDVI').lt(1))
        img = img.updateMask(mask_ndvivalue1)

        imgtest = ee.Image(imageList.get(i)).reproject('EPSG:32645', None, 30)
        theMax = imgtest.select('NDVI').reduceRegion(ee.Reducer.max(), aoi, maxPixels=1e13)
        # print (theMax.getInfo().get('NDVI'))
        if theMax.getInfo().get('NDVI') is not None:  # check if theMax.getInfo().get('NDVI') is none
            crs_raw = crs_raw.replace(":", "_")

            filename = crs_raw + '_' + 'medoid_comp_' + str(int(aa) + startYear) + '.tif'
            geemap.ee_export_image_to_drive(img.select('NDVI'), description=filename, folder='north_l7_cloudnnew1te',
                                            region=aoi, scale=30)

            filename2 = crs_raw + '_' + 'medoid_comp_clearcount_' + str(int(aa2) + startYear) + '.tif'
            geemap.ee_export_image_to_drive(img_clear, description=filename2, folder='north_l7_cloudnnew1te', region=aoi,
                                            scale=30)
