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
import geopandas as gpd
from shapely.ops import unary_union
from setting import *


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
    inRiverDiv = ogr.Open(inRiverDiv, 1)
    inSlpp = gdal.Open(inSlpp)
    layer = inRiverDiv.GetLayer()
    layer_dmx = inLine.GetLayer()
    temp_tif = os.path.join(setting.output_dir, "temp_vector2raster.tif")
    temp_shp = os.path.join(setting.output_dir, "temp_createLine.shp")

    addField = setting.dmx_field['J']
    HC.CreateNewField(layer_dmx, addField, ogr.OFTReal)

    for i in range(layer.GetFeatureCount()):
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
            geom_list.append([round(geom.GetPoint(count)[0], 3), round(geom.GetPoint(count)[1], 3)])
        for j in range(layer_dmx.GetFeatureCount()):
            fc_dmx = layer_dmx.GetFeature(j)
            geom_dmx = fc_dmx.GetGeometryRef()

            # Obtain the intersection points of cross-sectional line data and river line data
            for count in range(geom_dmx.GetPointCount()):
                # print([geom_dmx.GetPoint(count)[0], geom_dmx.GetPoint(count)[1]] )
                # print([round(geom_dmx.GetPoint(count)[0], 2), round(geom_dmx.GetPoint(count)[1], 2)])
                # print(geom_list)
                if [geom_dmx.GetPoint(count)[0], geom_dmx.GetPoint(count)[1]] in geom_list or [
                    round(geom_dmx.GetPoint(count)[0], 3), round(geom_dmx.GetPoint(count)[1], 3)] in geom_list:
                    print(mean)
                    HC.UpdateField(layer_dmx, fc_dmx, addField, mean)

        del temp_raster

        os.remove(temp_tif)
        os.remove(temp_shp)

    for j in range(layer_dmx.GetFeatureCount()):
        fc_dmx = layer_dmx.GetFeature(j)
        if fc_dmx.GetField(addField) is None or fc_dmx.GetField(addField) <= 0:
            # print(addField)
            print(layer_dmx.GetFeature(j + 1).GetField(addField))
            HC.UpdateField(layer_dmx, fc_dmx, addField, layer_dmx.GetFeature(j + 1).GetField(addField))


def split_river_by_cross_sections(river_path, cross_section_path, output_path, buffer_size=0.0001):
    """
    使用断面线（多条线）切割一条河流线，并在每个交点处创建断裂。

    Args:
        river_path (str): 河流线 Shapefile 文件路径 (应只包含一条线)。
        cross_section_path (str): 断面线 Shapefile 文件路径 (可包含多条线)。
        buffer_size (float): 用于切割的缓冲区半径大小。确保此值小于任何两个交点之间的距离。
        output_path (str, optional): 可选的输出文件路径 (例如 'result.shp')。如果为 None，则不保存文件。

    Returns:
        gpd.GeoDataFrame: 包含分割后的河流线段的 GeoDataFrame。

    Raises:
        ValueError: 如果河流线文件不包含恰好一条线段时抛出。
    """

    # 1. 读取数据
    print("split rivers...")
    # river_gdf = gpd.read_file(r"L:\College\project\GD\ningguo_geodata\Input\river.shp")
    river_gdf = gpd.read_file(river_path)
    cross_section_gdf = gpd.read_file(cross_section_path)

    # 验证河流数据
    if len(river_gdf) == 0:
        raise ValueError("河流线数据（river_path）中不包含任何几何体。")
    if len(river_gdf) > 1:
        print("警告: 河流线文件中包含多条线段。函数将只处理第一条线。")

    # 获取唯一的河流线几何体
    river_line = river_gdf.iloc[0].geometry
    # 合并所有断面线
    cross_sections_union = cross_section_gdf.geometry.unary_union

    # 2. 查找交点
    intersection_geoms = river_line.intersection(cross_sections_union)

    # 提取交点列表
    points = []
    if intersection_geoms.geom_type == 'Point':
        points = [intersection_geoms]
    elif intersection_geoms.geom_type == 'MultiPoint':
        points = list(intersection_geoms.geoms)
    # 忽略 LineString 或其他非点状的交集类型
    if points:
        intersection_data = {
            'point_id': range(1, len(points) + 1),
            'POINT_X': [point.x for point in points],  # 新增X坐标字段
            'POINT_Y': [point.y for point in points],  # 新增Y坐标字段
            'geometry': points
        }
        # 使用河流线的坐标系
        intersection_gdf = gpd.GeoDataFrame(
            intersection_data,
            crs=river_gdf.crs,
            geometry='geometry'
        )
    else:
        # 如果没有交点，则创建空的 GeoDataFrame
        intersection_gdf = gpd.GeoDataFrame(
            {'point_id': [],
             'POINT_X': [],
             'POINT_Y': []},
            crs=river_gdf.crs,
            geometry=gpd.points_from_xy([], []),
        )
        print("💡 没有找到交点。河流线不会被分割。")
        # 直接返回原始河流线和空的交点GeoDataFrame
        return river_gdf.copy(), intersection_gdf

    # 3. 创建切割工具（点缓冲区）
    #
    buffers = [point.buffer(buffer_size) for point in points]
    cutting_tool = unary_union(buffers)  # 得到一个 MultiPolygon

    # 4. 执行切割操作
    split_lines = river_line.difference(cutting_tool)

    # 5. 结果整理
    final_geometries = []
    if split_lines.geom_type == 'LineString':
        final_geometries = [split_lines]
    elif split_lines.geom_type == 'MultiLineString':
        # 拆解 MultiLineString 为独立的 LineString
        final_geometries = list(split_lines.geoms)
    else:
        print(f"⚠️ 分割结果类型为 {split_lines.geom_type}，未得到预期的 LineString 或 MultiLineString。")

    # 创建最终的 GeoDataFrame
    result_data = {
        'id': range(1, len(final_geometries) + 1),
        'original_river_name': [river_gdf.iloc[0].get(river_gdf.columns[0], 'Original River')] * len(final_geometries),
        # 尝试使用第一个属性列作为名称
        'geometry': final_geometries
    }
    result_gdf = gpd.GeoDataFrame(
        result_data,
        crs=river_gdf.crs
    )

    # 6. 结果输出
    if output_path:
        result_gdf.to_file(output_path, encoding='utf-8')
        intersection_gdf.to_file(temp_seed, encoding='utf-8')


def main(dmx, dem, createBasin=False):
    """
    Break the cross-section line and obtain the value of each breakpoint
    param inRiver: A river line
    param inSlpp: Slope
    :param createBasin:
    :param inDEM: Digital Elevation model
    :param inLine: Cross-section line
    :return: Breaking point
    """
    print("-----------------The execution of breaking the cross-section line begins-----------------")
    dmxLong = setting.temp_dmxLong
    ZL.main(dmx, dmxLong)
    dmxLongDs = ogr.Open(dmxLong)
    print(" Calculate the average specific drop of the cross-section line...")
    split_river_by_cross_sections(river_dir, dmx_dir, temp_river_div)

    HC.calculate_slope(dem_dir, temp_slope)
    getJ(dmx, temp_slope, temp_river_div)
    print(" Calculate roughness...")
    layer = dmx.GetLayer()
    layerLong = dmxLongDs.GetLayer()
    HC.CreateNewField(layer, setting.dmx_field['n0'], ogr.OFTReal)
    # tempDmx = r'E:\College\project\GD\geodata\keshan\temp_dmx.shp'
    tempDmx = os.path.join(setting.output_dir, 'temp_dmx.shp')
    # tempDmxRaster = r'E:\College\project\GD\geodata\keshan\temp_dmxRaster.tif'
    tempDmxRaster = os.path.join(setting.output_dir, 'temp_dmxRaster.tif')
    print(" Break the cross-section line...")
    dem_arr = dem.ReadAsArray()
    transform = dem.GetGeoTransform()
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
                # layer.GetFeature(j).GetField(setting.dmx_field['ObjectID']))
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

    if createBasin:
        print(" Generate sub-watershed...")
        HC.CreateBasins(dem, temp_seed, setting.temp_units)

    print("-----------------The breaking section line operation was successful-----------------")
    return setting.dmx_points
