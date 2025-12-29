# coding:utf-8
import os.path

import numpy as np
from osgeo import ogr, gdal, osr

import setting
import CalHydrological.Common as HC
from setting import *
"""
Calculate the peak flood flow in the basin
Qm_0: The default assumed value is 900 (m3/s)
Qm_1: The initial calculated value is 0 (m3/s)
"""


def hongfeng(F, L, J, n, Sp, u, m):
    """
    Calculate the peak flood flow
    param F: Drainage area
    param L: The length of the longest river within the basin
    param J: Average slope drop along L
    param n: Rainstorm Attenuation Index
    :param Sp: Design Rainstorm
    :param u: Flow generation parameter
    :param m: Convergence parameter
    return: Peak flood flow rate Qm
    """

    Qm_0 = setting.Qm_0
    Qm_1 = 0

    def __getQmWhenYes__():  # The calculation formula for the total confluence flood peak
        # print(pow(t, n) - u)
        out_Qm = 0.278 * (Sp / pow(t, n) - u) * F
        return out_Qm

    def __getQmWhenNo__():  # Calculation formula for partial confluence flood peaks
        out_Qm = 0.278 * (Sp * pow(tc, 1 - n) - u * tc) * F / t
        return out_Qm

    while True:
        t = 0.278 * L / (m * pow(J, 1 / 3.0) * pow(Qm_0, 1 / 4.0))
        tc = pow((1 - n) * Sp / u, 1 / n)
        # print(Qm_0,Qm_1,abs(Qm_0 - Qm_1) / Qm_0)
        if tc >= t:
            Qm_1 = __getQmWhenYes__()
        else:
            Qm_1 = __getQmWhenNo__()

        if abs(Qm_0 - Qm_1) / Qm_0 <= 0.01:
            break
        else:
            Qm_0 = Qm_1

    m1 = 0.278 * L / (t * pow(J, 1 / 3.0) * pow(Qm_0, 1 / 4.0))
    if abs(m1 - m) <= 0.01:
        qm_result = int(Qm_0)
        return qm_result


def main(inFc, inRiver, inSp, inN, inField):
    """
    Calculate the peak flood flow of each sub-basin
    :param inSeed: Seed point
    :param inFlowDir: Flow to the raster directory
    param inFc: Basin surface
    param inRiver: A river line
    :param inSlope: Slope grid
    """
    print("-----------------The Computational Design Peak has begun to be implemented-----------------")
    temp_riverInUnit = os.path.join(setting.output_dir, "temp_riverinunt.shp")  # rasterized river
    temp_river = os.path.join(setting.output_dir, 'temp_river.tif')  # rasterized river

    layer_ds = ogr.Open(setting.temp_units, 1)
    layer_unit = layer_ds.GetLayer()
    # Add fields to the sub-basin
    HC.CreateNewField(layer_unit, setting.unit_fields['Sp'], ogr.OFTReal)
    HC.CreateNewField(layer_unit, setting.unit_fields['n'], ogr.OFTReal)
    print(" Get Basin information...")
    listdata = []
    river_lyr = inRiver.GetLayer()
    HC.CreateNewField(layer_unit, "J", ogr.OFTReal)
    HC.CreateNewField(layer_unit, "L", ogr.OFTReal)
    for i in range(layer_unit.GetFeatureCount()):
        fc_unit = layer_unit.GetFeature(i)
        # fc_unit = layer_unit.GetFeature(layer_unit.GetFeatureCount()-1)
        geom_unit = fc_unit.GetGeometryRef()
        riverInUnit = None  # Rivers within the current basin
        for river_i in range(river_lyr.GetFeatureCount()):
            river_fc = river_lyr.GetFeature(river_i)
            river_geom = river_fc.GetGeometryRef()
            riverInUnit = geom_unit.Intersection(river_geom)
            if riverInUnit.GetPointCount() == 0:
                # print("The two do not intersect.")
                pass
            else:
                # print("The two lines intersect.")
                break
        # print("riverInUnit", type(riverInUnit), riverInUnit)
        # Calculate the river gradient
        HC.createLine(riverInUnit, temp_riverInUnit, ogr.wkbLineString)
        HC.vector2raster(temp_riverInUnit, temp_river)
        river_arr = gdal.Open(temp_river).ReadAsArray()

        temp_raster_yx = np.where(river_arr == 0)
        temp_raster_xy = []  # Rows occupied by rivers
        for index, x in enumerate(temp_raster_yx[1]):
            temp_xy = [x, temp_raster_yx[0][index]]
            temp_raster_xy.append(temp_xy)

        slope_array = (gdal.Open(temp_slope)).ReadAsArray()
        sum = 0
        for n in temp_raster_xy:
            slope_value = slope_array[n[1]][n[0]]
            if slope_value == -3.402823e+38 or slope_value < 0:
                slope_value = 0
            sum += slope_value
        if sum == 0:
            mean = 22.121  # Use the average of all rivers as a substitute
        else:
            mean = sum / len(temp_raster_xy)
        HC.UpdateField(layer_unit, fc_unit, "J", mean)

        # print('mean:', sum, len(temp_raster_xy), mean)
        area = geom_unit.GetArea() / 1000000  # The drainage area is km²

        length = riverInUnit.Length()
        if length <=0:
            length = 100
        HC.UpdateField(layer_unit, fc_unit, "L", length)
        # print("length",length,mean,i)
        if area <= 50:
            m = 1
        elif area <= 100:
            m = 1.5
        else:
            m = 2
        listdata.append([area,  # drainage area km²
                         length / 1000,  # river length, modify the unit to kilometers
                         mean / 100,  # river gradient, for example :90%-->0.9
                         inN,
                         inSp,
                         2.5,  # Flow parameter u
                         m])
    # print("listdata", listdata)
    print(" Calculate and save peak flood flow...")
    HC.CreateNewField(layer_unit, inField, ogr.OFTReal)
    for i in range(layer_unit.GetFeatureCount()):
        fc_unit = layer_unit.GetFeature(i)
        temp = hongfeng(listdata[i][0], listdata[i][1], listdata[i][2], listdata[i][3], listdata[i][4], listdata[i][5],listdata[i][6])
        HC.UpdateField(layer_unit, fc_unit, inField, temp)
    layer_ds = None
    print("-----------------The Computational Design Peak operation was successful-----------------")
