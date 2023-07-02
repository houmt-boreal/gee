import os, glob
import fiona
import rasterio
import rasterio.mask
import numpy as np
in_directory = r'D:\xinjiang_chapter3\l578_medoid\result_2003_2020\medoid_mk_result_summer'
list=os.listdir(in_directory)#read the fold and get content

with fiona.open(r"D:\d\bou4_4m\altai_3zone.shp", "r") as shapefile:
# with fiona.open("D:/jag_modified/fi_rect_small_fin.shp", "r") as shapefile:
    shapes = [feature["geometry"] for feature in shapefile]

# https://rasterio.readthedocs.io/en/latest/api/rasterio.mask.html
for i in range(0,len(list),1):
  filea=str(in_directory+'\\'+list[i])
  #print filea[len(filea) - 46:len(filea) - 29]
  with rasterio.open(filea) as src:
    out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True, nodata=-99,all_touched=True)
    out_meta = src.meta

  out_meta.update({"driver": "GTiff",
                 "height": out_image.shape[1],
                 "width": out_image.shape[2],
                 "transform": out_transform})

  out_raster = r'D:\xinjiang_chapter3\l578_medoid\result_2003_2020\medoid_mk_result_summer_altai\altai_summer_' + list[i]
  with rasterio.open(out_raster, "w", **out_meta) as dest:

      dest.nodata = -99
      # out_image = np.where(out_image == -99, -99, out_image)
      # print (np.argwhere(np.isnan(out_image)))
      out_image[np.isnan(out_image)] = -99
      out_image = np.where(out_image > 50, -99, out_image)
      out_image = np.where(out_image < -10, -99, out_image)
      # print (out_image)
      dest.write(out_image)

print ('ok')

