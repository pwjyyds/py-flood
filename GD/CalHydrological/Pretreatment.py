# coding:utf-8
"""
Data preprocessing, including breaking section lines and calculating specific drops, etc
"""
import os
import CalHydrological.Common as HC
import numpy as np
from osgeo import gdal, ogr, osr

import setting
import SimulateFlood.LongLine as ZL


def getJ(inLine, inSlpp, inRiverDiv):
    """
    Calculate the average ratio drop of each cross-sectional line
    :param inLine: Cross-section line
    param inSlpp: Slope
    :param inRiverDiv: A broken river line
    return: Cross-sectional line with average ratio drop
    """
    # 1. Using cross-sectional lines to break the river --> Just input the processed data directly. GDAL is so hard to achieve
    # 2. Calculate the average slope J of each section of the river
    layer = inRiverDiv.GetLayer()
    layer_dmx = inLine.GetLayer()
    # temp_tif = r"E:\College\project\GD\geodata\keshan\temp_vector2raster.tif"
    temp_tif = os.path.join(setting.output_dir,"temp_vector2raster.tif")
    # temp_shp = r'E:\College\project\GD\geodata\keshan\temp_createLine.shp'
    temp_shp =os.path.join(setting.output_dir,"temp_createLine.shp")

    addField = setting.dmx_field['J']
    HC.CreateNewField(layer_dmx, addField, ogr.OFTReal)

    for i in range(layer.GetFeatureCount()):
        # for i in range(1):
        fc = layer.GetFeature(i)
        geom = fc.GetGeometryRef()
        HC.createLine(geom, temp_shp, ogr.wkbLineString)
        HC.vector2raster(inputfilePath=temp_shp,
                         outputfile=temp_tif)
        temp_raster = gdal.Open(temp_tif)
        temp_raster_array = temp_raster.ReadAsArray()
        temp_raster_yx = np.where(temp_raster_array == 0)
        temp_raster_xy = []
        for index, x in enumerate(temp_raster_yx[1]):
            temp_xy = [x, temp_raster_yx[0][index]]
            temp_raster_xy.append(temp_xy)

        slope_array = inSlpp.ReadAsArray()
        sum = 0
        for n in temp_raster_xy:
            slope_value = slope_array[n[1]][n[0]]
            # print("slope_array[n[1]][n[0]]:",slope_value)
            if slope_value == -3.402823e+38 or slope_value < 0:
                slope_value = 0
            sum += slope_value
        mean = sum / len(temp_raster_xy)
        # print(sum, len(temp_raster_xy), mean)

        geom_list = []
        for count in range(geom.GetPointCount()):
            geom_list.append([round(geom.GetPoint(count)[0], 4), round(geom.GetPoint(count)[1], 4)])
        for j in range(layer_dmx.GetFeatureCount()):
            fc_dmx = layer_dmx.GetFeature(j)
            geom_dmx = fc_dmx.GetGeometryRef()

            # Obtain the intersection points of cross-sectional line data and river line data
            for count in range(geom_dmx.GetPointCount()):
                if [geom_dmx.GetPoint(count)[0], geom_dmx.GetPoint(count)[1]] in geom_list or [
                    round(geom_dmx.GetPoint(count)[0], 4), round(geom_dmx.GetPoint(count)[1], 4)] in geom_list:
                    HC.UpdateField(layer_dmx, fc_dmx, addField, mean)

        del temp_raster

        os.remove(temp_tif)
        os.remove(temp_shp)

    for j in range(layer_dmx.GetFeatureCount()):
        fc_dmx = layer_dmx.GetFeature(j)
        if fc_dmx.GetField(addField) is None or fc_dmx.GetField(addField) <= 0:
            HC.UpdateField(layer_dmx, fc_dmx, addField, layer_dmx.GetFeature(j + 1).GetField(addField))


# def main(inLine, inDEM, inSlpp, inRiver, workSpace):
def main(dmx, slope, river, riverDiv, dem, inFlowDir, inSeed):
    """
    Break the cross-section line and obtain the value of each breakpoint
    param inRiver: A river line
    param inSlpp: Slope
    :param inDEM: Digital Elevation model
    :param inLine: Cross-section line
    :return: Breaking point
    """
    print("-----------------The execution of breaking the cross-section line begins-----------------")
    dmxLong = setting.temp_dmxLong
    ZL.main(dmx,dmxLong)
    dmxLongDs = ogr.Open(dmxLong)
    print(" Calculate the average specific drop of the cross-section line...")
    getJ(dmx, slope, riverDiv)
    print(" Calculate roughness...")
    layer = dmx.GetLayer()
    layerLong = dmxLongDs.GetLayer()
    HC.CreateNewField(layer, setting.dmx_field['n0'], ogr.OFTReal)
    # tempDmx = r'E:\College\project\GD\geodata\keshan\temp_dmx.shp'
    tempDmx = os.path.join(setting.output_dir,'temp_dmx.shp')
    # tempDmxRaster = r'E:\College\project\GD\geodata\keshan\temp_dmxRaster.tif'
    tempDmxRaster = os.path.join(setting.output_dir,'temp_dmxRaster.tif')
    print(" Break the cross-section line...")
    dem_arr = dem.ReadAsArray()
    transform = slope.GetGeoTransform()
    cell_size_x = transform[1]
    # print(cell_size_x)

    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds = driver.CreateDataSource(setting.dmx_points)

    output_srs = osr.SpatialReference()
    output_srs.ImportFromEPSG(32650)  # WGS_1984_UTM_Zone_50N的EPSG代码是32650

    layer_points = ds.CreateLayer('temp_dmx_points', geom_type=ogr.wkbPoint, srs=output_srs)

    HC.CreateNewField(layer_points, setting.dmxPoints_field['ObjectID'], ogr.OFTInteger)
    HC.CreateNewField(layer_points, setting.dmxPoints_field['DmxID'], ogr.OFTInteger)
    HC.CreateNewField(layer_points, setting.dmxPoints_field['J'], ogr.OFTReal)
    HC.CreateNewField(layer_points, setting.dmxPoints_field['n0'], ogr.OFTReal)
    HC.CreateNewField(layer_points, setting.dmxPoints_field['DemValue'], ogr.OFTReal)

    for field in setting.fields_z:
        HC.CreateNewField(layer, field, ogr.OFTReal)
    for field in setting.fields_h:
        HC.CreateNewField(layer, field, ogr.OFTReal)



    for j in range(layerLong.GetFeatureCount()):
        fc_ = layerLong.GetFeature(j)
        geom_ = fc_.GetGeometryRef()
        HC.UpdateField(layerLong, fc_, setting.dmx_field['n0'], 0.025)
        HC.createLine(geom_, tempDmx, ogr.wkbLineString)

        HC.vector2raster(tempDmx, tempDmxRaster)
        dmx_raster = gdal.Open(tempDmxRaster)

        raster_arr = dmx_raster.ReadAsArray()
        raster_yx = np.where(raster_arr == 0)
        for index, x in enumerate(raster_yx[1]):
            px, py = HC.imagexy2geo(dmx_raster, raster_yx[0][index], x)
            if px > 0 and py > 0:

                feature = ogr.Feature(layer_points.GetLayerDefn())
                point = ogr.Geometry(ogr.wkbPoint)
                point.AddPoint(px + cell_size_x / 2, py - cell_size_x / 2)
                feature.SetGeometry(point)
                layer_points.CreateFeature(feature)
                # print('ObjectID', index, "DmxID", layer.GetFeature(j).GetField(setting.dmx_field['ObjectID']))
                HC.UpdateField(layer_points, feature, setting.dmxPoints_field['ObjectID'], index)
                HC.UpdateField(layer_points, feature, setting.dmxPoints_field['DmxID'],
                               layer.GetFeature(j).GetField(setting.dmx_field['ObjectID']))
                HC.UpdateField(layer_points, feature, setting.dmxPoints_field['J'],
                               layer.GetFeature(j).GetField(setting.dmx_field['J']))
                HC.UpdateField(layer_points, feature, setting.dmxPoints_field['n0'],
                               layer.GetFeature(j).GetField(setting.dmx_field['n0']))
                if float(dem_arr[raster_yx[0][index]][x]) < 0:
                    HC.UpdateField(layer_points, feature, setting.dmxPoints_field['DemValue'],
                                   500)
                else:
                    HC.UpdateField(layer_points, feature, setting.dmxPoints_field['DemValue'],
                                   float(dem_arr[raster_yx[0][index]][x]))
        del dmx_raster
        os.remove(tempDmx)
        os.remove(tempDmxRaster)

    ds = None
    dmxLongDs = None
    dmx = None
    print(" Generate sub-watershed...")
    HC.CreateBasins(inFlowDir, inSeed, setting.temp_units)
    print("-----------------The breaking section line operation was successful-----------------")
    return setting.dmx_points
