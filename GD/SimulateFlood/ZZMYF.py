# coding:utf-8
import os.path

from osgeo import gdal
from osgeo.gdalconst import *  # GDAL中常用的一些常量
import geopandas
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from osgeo import ogr
from matplotlib import colors

import SimulateFlood.LongLine, SimulateFlood.CreateDivLine
import CalHydrological.Common as HC
import setting


def getSeedXY(seed, data, dmx, dem_array, field_z, field_h):
    """
    Obtain the longitude and latitude coordinates and other attributes of all seeds, convert them into row and column numbers, and store them in a seed dictionary
    """
    seed_dict = []
    layer = seed.GetLayer(0)
    layer_dmx = dmx.GetLayer(0)
    transform = data.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = transform[5]

    for i in range(layer.GetFeatureCount()):
        feat = layer.GetFeature(i)
        feat_dmx = layer_dmx.GetFeature(i)
        x = feat.GetField('POINT_X')
        y = feat.GetField('POINT_Y')
        x_col = int((x - xOrigin) / pixelWidth)
        y_row = int((y - yOrigin) / pixelHeight)
        high = dem_array[y_row][x_col]

        seed_h = feat_dmx.GetField(field_z) + high
        HC.UpdateField(layer_dmx, feat_dmx, field_h, seed_h)  # Water level + elevation value

        seed_dict.append({
            "FID": feat.GetField('FID_middle'),
            "X": x_col,
            "Y": y_row,
            "done": 0,
            "z": seed_h
        })

    # print(seed_dict)
    return seed_dict


def saveResult(data, floodPart, savePath):
    """
    Generate inundation zones and convert numpy arrays to tif
    """
    fleed = np.zeros((data.RasterYSize, data.RasterXSize))
    for i in floodPart:
        # print(i)
        fleed[i[1], i[0]] = 1

    output_file = savePath

    rows, cols = fleed.shape

    driver = gdal.GetDriverByName("GTiff")
    out_data_set = driver.Create(output_file, cols, rows, 1, gdal.GDT_Float32)

    out_data_set.GetRasterBand(1).WriteArray(fleed)

    out_data_set.SetGeoTransform(data.GetGeoTransform())
    out_data_set.SetProjection(data.GetProjection())

    out_data_set = None


def main(data, dmx, seed, inZField, inHField):
    """
    Input data: DEM, section line, seed point (with flood peak water level)
    data = gdal.Open(r"E:\College\project\geoData\dsz_dem")  # DEM
    dmx = ogr.Open(r"E:\College\project\geoData\dxm.shp") # section line
    seed = ogr.Open(r"E:\College\project\geoData\seed.shp") # seed point
    """
    print("-----------------“计算淹没区”开始执行-----------------")
    # print(inZField)
    temp_dmxLong = os.path.join(setting.output_dir,'temp_dmxLong.shp')
    temp_divLine = os.path.join(setting.output_dir,'temp_divLine.shp')
    temp_divLine_tif =  os.path.join(setting.output_dir,'temp_divLine.tif')
    result_path = os.path.join(setting.output_dir, 'result_' + inZField + '.tif')
    dem_array = data.ReadAsArray()
    print("Extended section line...")
    SimulateFlood.LongLine.main(dmx, temp_dmxLong)
    print("Cut the river...")
    SimulateFlood.CreateDivLine.main(seed, temp_divLine)
    HC.vector2raster(inputfilePath=temp_divLine, outputfile=temp_divLine_tif)
    """
    Read all the cells occupied by the extended section lines and store them in the section line array
    """
    print("Read the extended cross-sectional line...")
    newDmx = gdal.Open(temp_divLine_tif)
    newDmx_raster_array = newDmx.ReadAsArray()
    newDmx_yx = np.where(newDmx_raster_array == 0)  # [[y1,y2,y3,...],[x1,x2,x3....]]
    newDmx_xy = []
    for index, x in enumerate(newDmx_yx[1]):
        temp_xy = [x, newDmx_yx[0][index]]
        newDmx_xy.append(temp_xy)
    # a = pd.DataFrame(dmx_yx)

    # seed_dict = getSeedXY(seed, data)
    seed_dict = getSeedXY(seed, data, dmx, dem_array, inZField, inHField)

    """
    Simulated inundation area
    """
    print("Simulated inundation area...")

    p = inZField  # Frequency field
    floodPart = []  # Flooded array
    newPoint = []  # New origin array
    haveDonePoint = []  # has traversed the origin array
    for i in seed_dict:
        seed_x = i["X"]
        seed_y = i["Y"]

        if i["done"] == 0:

            newPoint.append([seed_x, seed_y])
            for point in newPoint:

                if point not in haveDonePoint:
                    for x in range(point[0] - 1, point[0] + 2):
                        for y in range(point[1] - 1, point[1] + 2):
                            # 1. Determine the grid line array [[y1, y2, y3,...], [x1, x2, x3...]] : if in, don't do b
                            # 2. The origin elevation is ≤ the current seed flood peak water level
                            if [x, y] not in newDmx_xy:
                                if x in range(data.RasterXSize) and y in range(data.RasterYSize):
                                    thisDem = dem_array[y][x]
                                else:
                                    continue
                                if thisDem != -32768 and 0 < thisDem <= i["z"]:
                                    # The flooded cells are stored in a flooded array and a new origin array
                                    floodPart.append([x, y])
                                    newPoint.append([x, y])
                                else:
                                    pass
                            else:
                                floodPart.append([x, y])  # By default, all points on the boundary are submerged
                                break

                        else:
                            i["done"] = 1  # After the eight-neighborhood traversal is completed, the current seed is marked as having been traversed
                            haveDonePoint.append(point)
                            continue
                        break

        newPoint = []

    print("save...")
    saveResult(data, floodPart, result_path)
    print("----------------- running successfully -----------------")
