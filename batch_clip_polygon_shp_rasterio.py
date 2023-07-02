import os, glob
import fiona
import rasterio
import rasterio.mask
import numpy as np
in_directory = r'D:\xinjiang_chapter3\l578_medoid\result_2003_2020\medoid_mk_result_summer'
files1 = glob.glob(os.path.join(in_directory, '*.tif'))

with fiona.open(r"D:\d\bou4_4m\altai_3zone.shp", "r") as shapefile:
# with fiona.open("D:/jag_modified/fi_rect_small_fin.shp", "r") as shapefile:
    shapes = [feature["geometry"] for feature in shapefile]


for filea in files1:
  #print filea[len(filea) - 46:len(filea) - 29]
  with rasterio.open(filea) as src:
    out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True, nodata=-99,all_touched=True)
    out_meta = src.meta

  out_meta.update({"driver": "GTiff",
                 "height": out_image.shape[1],
                 "width": out_image.shape[2],
                 "transform": out_transform})

  out_raster = r'D:\xinjiang_chapter3\l578_medoid\result_2003_2020\medoid_mk_result_summer_altai\altai_summer_' + filea[len(filea) - 34:len(filea) - 8] +'.tif'
  with rasterio.open(out_raster, "w", **out_meta) as dest:

      dest.nodata = -99
      # out_image = np.where(out_image == -99, -99, out_image)
      # print (np.argwhere(np.isnan(out_image)))
      out_image[np.isnan(out_image)] = -99
      out_image = np.where(out_image > 1, -99, out_image)
      out_image = np.where(out_image < 0, -99, out_image)
      # print (out_image)
      dest.write(out_image)

print ('ok')


