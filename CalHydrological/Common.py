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
            point.AddPoint(px + cell_size_x / 2, py - cell_size_x / 2)  # Let the dot be generated in the middle of the grid, or in the upper left corner of the grid
            feature.SetGeometry(point)
            layer.CreateFeature(feature)

    ds = None


def flow(dem_dir):
    """
    :param dem_dir:DEM
    """
    dem = rd.LoadGDAL(dem_dir)
    accum_d8 = rd.FlowAccumulation(dem, method='D8')
    # print(accum_d8, type(accum_d8))
    rd.SaveGDAL("/geodata/keshan\hydrink\liuxiangFromRd.tif", accum_d8)
    # d8_fig = rd.rdShow(accum_d8, figsize=(12, 8), axes=False, cmap='jet')


# flow(r"E:\College\project\GD\geodata\Input\dem.tif")

def FlowDir(dem_dir):
    """
    Flow direction
    https://deltares.github.io/pyflwdir/latest/_examples/from_dem.html#Derive-flow-direction
    """
    # read elevation data of the rhine basin using rasterio
    with rasterio.open(dem_dir, "r") as src:
        elevtn = src.read(1)
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
    with rasterio.open(r"/geodata/keshan\hydrink\liuxiangFromFlw.tif", "w", **prof) as src:
        src.write(d8_data, 1)
    # print("FlowDir ok!")


# FlowDir(r"E:\College\project\GD\geodata\Input\dem.tif")


# RasterToPoint()
def CreateBasins(fldir, outlet, outPath):
    """
    Generating basin surface
    :param outPath: Result saving path
    :param fldir: Flow to raster
    param outlet: Data of the water outlet point
    :return:
    """

    # read and parse data
    with rasterio.open(fldir, "r") as src:
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

    gdf_out = gpd.GeoSeries(gpd.points_from_xy(x, y, crs=4326))
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
