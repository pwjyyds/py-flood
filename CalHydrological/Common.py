"""
General method
"""
import os
import pygeodesy
import pyproj
from osgeo import ogr, gdal, osr
import richdem as rd
import geopandas as gpd
import numpy as np
import pyflwdir
import rasterio
from pygeodesy.sphericalNvector import LatLon

import setting


# convenience method for vectorizing a raster
def vectorize(data, nodata, transform, name="value"):
    from rasterio import features
    # read example elevation data and derive background hillslope
    fn = os.path.join(os.path.dirname(__file__), setting.dem_dir)
    with rasterio.open(fn, "r") as src:
        elevtn = src.read(1)
        extent = np.array(src.bounds)[[0, 2, 1, 3]]
        crs = src.crs
    feats_gen = features.shapes(
        data,
        mask=data != nodata,
        transform=transform,
        connectivity=8,
    )
    feats = [
        {"geometry": geom, "properties": {name: val}} for geom, val in list(feats_gen)
    ]

    # parse to geopandas for plotting / writing to file
    gdf = gpd.GeoDataFrame.from_features(feats, crs=crs)
    gdf[name] = gdf[name].astype(data.dtype)
    return gdf


def CreateNewField(layer, fieldName, fieldType):
    """
    New field
    :param layer: Layer
    param fieldName: The name of the field to be added
    :return:
    """
    layerDefinition = layer.GetLayerDefn()
    ifExistField = False
    for i in range(layerDefinition.GetFieldCount()):
        # print(layerDefinition.GetFieldDefn(i).GetName())
        if layerDefinition.GetFieldDefn(i).GetName() == fieldName:
            ifExistField = True
            break
    if not ifExistField:
        fieldDefn = ogr.FieldDefn(fieldName, fieldType)
        layer.CreateField(fieldDefn)


def UpdateField(layer, feature, fieldName, fieldValue):
    """
    Update field
    :param layer: Layer
    :param feature: Feature object
    :param fieldName: Field name
    :param fieldValue: The value to be updated
    :return:
    """
    feature.SetField(fieldName, fieldValue)
    layer.SetFeature(feature)


def vector2raster(inputfilePath, outputfile, bands=[1], burn_values=[0], field="", all_touch="False"):
    """
    inputfilePath input vector files
    outputfile outputs raster files
    """
    import warnings
    warnings.filterwarnings('ignore')
    data = gdal.Open(
        setting.dem_dir)  # Raster template file, determine the metadata of the output raster (coordinate system, etc., raster size, range, etc.)
    # Determine the grid size
    x_res = data.RasterXSize
    y_res = data.RasterYSize
    vector = inputfilePath
    if type(inputfilePath) is str:
        vector = ogr.Open(inputfilePath)

    layer = vector.GetLayer()

    featureCount = layer.GetFeatureCount()

    targetDataset = gdal.GetDriverByName('GTiff').Create(outputfile, x_res, y_res, 1, gdal.GDT_Byte)

    targetDataset.SetGeoTransform(data.GetGeoTransform())
    targetDataset.SetProjection(data.GetProjection())

    band = targetDataset.GetRasterBand(1)

    NoData_value = 255
    band.SetNoDataValue(NoData_value)
    band.FlushCache()
    if field:
        # Call the rasterization function. The RasterizeLayer function has four parameters, namely the raster object, band, vector object, and options
        # options can have multiple attributes, among which the ATTRIBUTE attribute takes the value of a certain field attribute of the vector layer as the converted raster value
        gdal.RasterizeLayer(targetDataset, bands, layer, burn_values=burn_values,
                            options=["ALL_TOUCHED=" + all_touch, "ATTRIBUTE=" + field])
    else:
        gdal.RasterizeLayer(targetDataset, bands, layer, burn_values=burn_values, options=["ALL_TOUCHED=" + all_touch])


def createLine(inGeom, savePath, type):
    # Create a new shapefile based on geom
    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds = driver.CreateDataSource(savePath)
    output_srs = osr.SpatialReference()
    output_srs.ImportFromEPSG(32650)  # WGS_1984_UTM_Zone_50N

    layer = ds.CreateLayer('temp_createLine', geom_type=type, srs=output_srs)
    # print(layer)
    feature = ogr.Feature(layer.GetLayerDefn())
    feature.SetGeometry(inGeom)
    layer.CreateFeature(feature)

    ds = None


def imagexy2geo(dataset, row, col):
    """
    Convert the coordinates (row and column numbers) on the image map to projection coordinates or geographic coordinates (based on the coordinate system conversion of specific data) according to the six-parameter model of GDAL.
    param dataset: GDAL geographic data
    :param row: The row number of the pixel
    :param col: The column number of the pixel
    return: The projected coordinates or geographic coordinates (x, y) corresponding to the row number (row, col)
    """
    trans = dataset.GetGeoTransform()
    px = trans[0] + col * trans[1] + row * trans[2]
    py = trans[3] + col * trans[4] + row * trans[5]
    return px, py


def RasterToPoint(inputRaster, savePath):
    """
    Raster to vector point conversion
    :return:
    """
    transform = inputRaster.GetGeoTransform()

    cell_size_x = transform[1]
    # print(cell_size_x)

    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds = driver.CreateDataSource(savePath)

    output_srs = osr.SpatialReference()
    output_srs.ImportFromEPSG(32650)  # WGS_1984_UTM_Zone_50N

    layer = ds.CreateLayer('temp_createLine', geom_type=ogr.wkbPoint, srs=output_srs)

    raster_arr = inputRaster.ReadAsArray()
    raster_yx = np.where(raster_arr == 1)
    for index, x in enumerate(raster_yx[1]):
        px, py = imagexy2geo(inputRaster, raster_yx[0][index], x)
        if px > 0 and py > 0:
            feature = ogr.Feature(layer.GetLayerDefn())
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(px + cell_size_x / 2,
                           py - cell_size_x / 2)  # Let the dot be generated in the middle of the grid, or in the upper left corner of the grid
            feature.SetGeometry(point)
            layer.CreateFeature(feature)

    ds = None


# def flow(dem_dir):
#     """
#     :param dem_dir:DEM
#     """
#     dem = rd.LoadGDAL(dem_dir)
#     accum_d8 = rd.FlowAccumulation(dem, method='D8')
#     # print(accum_d8, type(accum_d8))
#     rd.SaveGDAL("/geodata/keshan\hydrink\liuxiangFromRd.tif", accum_d8)
#     # d8_fig = rd.rdShow(accum_d8, figsize=(12, 8), axes=False, cmap='jet')


# flow(r"E:\College\project\GD\geodata\Input\dem.tif")


def calculate_d8_flow_direction(input_dem_path, output_flow_path, fill_depressions=True):
    """
    使用 RichDEM 库计算 DEM 的 D8 流向栅格。

    Args:
        input_dem_path (str): 输入 DEM 文件的路径。
        output_flow_path (str): 输出 D8 流向栅格文件的路径。
        fill_depressions (bool): 是否执行洼地填充预处理。
    """

    # 1. 读取 DEM 数据
    # RichDEM 的核心数据结构是 rd.rdarray，可以从文件加载
    dem = rd.LoadGDAL(input_dem_path)

    # 2. (可选) 洼地填充预处理
    # 洼地（Depressions）会阻止水流，填充是计算流向前的标准步骤
    if fill_depressions:
        print("正在执行洼地填充...")
        # rd.FillDepressions 返回一个地形修正后的 rd.rdarray
        dem_filled = rd.FillDepressions(dem, epsilon=True, in_place=False)
    else:
        dem_filled = dem

    # 3. 计算 D8 流向 (Flow Direction)
    # rd.FlowDirection(DEM, method='D8') 是核心函数
    print("正在计算 D8 流向...")
    # D8 算法将流向编码为 1, 2, 4, 8, 16, 32, 64, 128
    # 对应东、东南、南、西南、西、西北、北、东北。
    #
    flow_direction_d8 = rd.FlowDirection(dem_filled, method='D8')

    # 4. 写入输出文件
    # 使用 Rasterio 的 GDAL 驱动将 rd.rdarray 写入文件
    # RichDEM 对象可以直接调用 SaveGDAL 方法
    rd.SaveGDAL(output_flow_path, flow_direction_d8)

    print(f"✅ D8 流向栅格已成功计算并保存到: {output_flow_path}")


def FlowDir(dem_dir, save_path):
    """
    Flow direction
    https://deltares.github.io/pyflwdir/latest/_examples/from_dem.html#Derive-flow-direction
    """
    # read elevation data of the rhine basin using rasterio
    with rasterio.open(dem_dir, "r") as src:
        elevtn = src.read(1, out_dtype='float32')
        nodata = src.nodata
        transform = src.transform
        crs = src.crs
        extent = np.array(src.bounds)[[0, 2, 1, 3]]
        latlon = src.crs.is_geographic
        prof = src.profile
    flw = pyflwdir.from_dem(
        data=elevtn,
        nodata=src.nodata,
        transform=transform,
        latlon=latlon,
        outlets="min",
    )
    d8_data = flw.to_array(ftype="d8")
    prof.update(dtype=d8_data.dtype, nodata=247)
    with rasterio.open(save_path, "w", **prof) as src:
        src.write(d8_data, 1)
    # print("FlowDir ok!")


# FlowDir(r"E:\College\project\GD\geodata\Input\dem.tif")


# RasterToPoint()
def CreateBasins(dem, outlet, outPath):
    """
    Generating basin surface
    :param outPath: Result saving path
    :param fldir: Flow to raster
    param outlet: Data of the water outlet point
    :return:
    """
    FlowDir(setting.dem_dir, setting.temp_flowDir)
    # read and parse data
    with rasterio.open(setting.temp_flowDir, "r") as src:
        flwdir = src.read(1)
        crs = src.crs
        flw = pyflwdir.from_array(
            flwdir,
            ftype="d8",
            transform=src.transform,
            latlon=crs.is_geographic,
            cache=True,
        )

    # define output locations
    outlet = ogr.Open(outlet, 1)
    layer = outlet.GetLayer(0)
    lyr_count = layer.GetFeatureCount()
    x_list = []
    y_list = []
    for i in range(lyr_count):
        # for i in range(1):
        feat = layer.GetFeature(i)
        x_list.append(feat.GetField('POINT_X'))
        y_list.append(feat.GetField('POINT_Y'))
    x = np.array(x_list)
    y = np.array(y_list)

    # gdf_out = gpd.GeoSeries(gpd.points_from_xy(x, y, crs=4326))
    # delineate subbasins
    subbasins = flw.basins(xy=(x, y), streams=flw.stream_order() >= 4)
    # vectorize subbasins using the vectorize convenience method from utils_keshan.py
    gdf_bas = vectorize(subbasins.astype(np.int32), 0, flw.transform, name="basin")
    gdf_bas.crs = pyproj.CRS.from_user_input('EPSG:32650')
    # gdf.rename(columns={'name':'ave_price'},inplace=True)
    # gdf.rename(columns={'addrees':'area_ave_price'},inplace=True)
    gdf_bas.to_file(outPath, driver='ESRI Shapefile', encoding='utf-8')
    # print("OK")


def CreateShapefile(inXY_arr, outPath, epsg):
    """
    Create vectors based on longitude and latitude
    """
    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds = driver.CreateDataSource(outPath)
    output_srs = osr.SpatialReference()
    output_srs.ImportFromEPSG(epsg)

    layer = ds.CreateLayer('newLine', geom_type=ogr.wkbLineString, srs=output_srs)

    feature = ogr.Feature(layer.GetLayerDefn())
    newline = ogr.Geometry(ogr.wkbLineString)
    for pointXY in inXY_arr:
        newline.AddPoint(pointXY[0], pointXY[1])
    feature.SetGeometry(newline)

    # feature.SetField('field_name', 'field_value')
    layer.CreateFeature(feature)

    ds = None


def SimplifyLine(XYlist):
    """
    Simplified line
    """
    points = []
    src_srs = osr.SpatialReference()  # defines the source coordinate system
    src_srs.ImportFromEPSG(32650)  # EPSG code 3857 represents Web Mercator projection
    tgt_srs = osr.SpatialReference()  # defines the target coordinate system
    tgt_srs.ImportFromEPSG(4326)  # EPSG code 4326 represents the WGS84 latitude and longitude coordinate system
    transform = osr.CoordinateTransformation(src_srs, tgt_srs)
    for p in XYlist:
        x, y = p[0], p[1]  #
        lon, lat, z = transform.TransformPoint(x, y)
        p1 = LatLon(lon, lat)
        points.append(p1)
    simplified_points = pygeodesy.simplifyRW(points, radius=1000)
    # print(simplified_points)

    return simplified_points


def calculate_slope(input_dem_path, output_slope_path, scale_factor=1.0, as_percent=True):
    """
    利用 Rasterio 计算 DEM 的坡度栅格。

    Args:
        input_dem_path (str): 输入 DEM 文件的路径。
        output_slope_path (str): 输出坡度栅格文件的路径。
        scale_factor (float): 比例因子。如果DEM是经纬度，使用 ~111120.0。
                              如果DEM是米 (UTM)，使用 1.0 (默认值)。
        as_percent (bool): True 表示输出坡度百分比 (tan(theta)*100)，False 表示角度 (度)。
    """
    if os.path.exists(output_slope_path):
        os.remove(output_slope_path)
    # 打开输入 DEM 文件
    with rasterio.open(input_dem_path) as src:
        dem_data = src.read(1)  # 读取第一个波段的数据
        profile = src.profile  # 获取栅格元数据

        # 获取分辨率
        # 在投影坐标系中，dx 和 dy 是米或 feet。
        # 在地理坐标系中，dx 和 dy 是度。
        dx = profile['transform'][0]
        dy = -profile['transform'][4]  # dy 通常为负，取绝对值

        # 使用 GDAL 的算法计算坡度
        # 计算水平 (dz/dx) 和垂直 (dz/dy) 方向的坡度变化率
        # numpy.gradient 返回数组中元素变化的离散差分
        dz_dy, dz_dx = np.gradient(dem_data, dy * scale_factor, dx * scale_factor)

        # 坡度的数学公式: slope = arctan(sqrt( (dz/dx)^2 + (dz/dy)^2 ))
        # np.sqrt(a**2 + b**2) 是勾股定理计算合向量的长度
        slope_rad = np.arctan(np.sqrt(dz_dx ** 2 + dz_dy ** 2))

        # 将弧度转换为度 (0-90)
        slope_deg = np.degrees(slope_rad)

        # 转换为百分比 (0-1000s)
        if as_percent:
            output_slope = slope_deg / 90.0 * 100.0
        else:
            output_slope = slope_deg  # 默认输出角度 (度)

        # 3. 写入输出文件
        # 更新输出栅格的元数据 (确保数据类型合适，并设置 nodata)
        profile.update(
            dtype=rasterio.float32,
            count=1,
            nodata=profile.get('nodata')  # 保持 NoData 值
        )

        with rasterio.open(output_slope_path, 'w', **profile) as dst:
            dst.write(output_slope.astype(rasterio.float32), 1)

        print(f"✅ 坡度栅格已成功计算并保存到: {output_slope_path}")
