"""
Extended section line
"""
from osgeo import ogr, osr
import sympy as sp
from scipy.optimize import curve_fit


# 定义线性方程
def linear_eq(x, k, b):
    return k * x + b


def main(oldLine, newLine):
    """
    Extended section line
    param oldLine: Cross-section line
    param newLine: Extended cross-section line
    """
    dataSource = oldLine
    if type(dataSource) is str:
        dataSource = ogr.Open(dataSource)

    sourceLayer = dataSource.GetLayer()

    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds = driver.CreateDataSource(newLine)
    output_srs = osr.SpatialReference()
    output_srs.ImportFromEPSG(32650)  # WGS_1984_UTM_Zone_50N

    layer = ds.CreateLayer('newLine', geom_type=ogr.wkbLineString, srs=output_srs)

    for feature in sourceLayer:
        geom = feature.GetGeometryRef()

        popt, pcov = curve_fit(linear_eq, [geom.GetX(0), geom.GetX(2)], [geom.GetY(0), geom.GetY(2)])  # 打印k和b的值
        k = popt[0]
        b = popt[1]
        # print("y =", k, "* x+", b)

        # Solve for the coordinates of the extended left and right points, that is, the coordinates of the starting point of the new line segment
        x3 = sp.symbols('x3')
        eq = sp.Eq(600, sp.sqrt((x3 - geom.GetX(2)) ** 2 + (k * x3 + b - geom.GetY(2)) ** 2))
        sol = sp.solve(eq, x3)

        feature = ogr.Feature(layer.GetLayerDefn())
        newline = ogr.Geometry(ogr.wkbLineString)
        newline.AddPoint(float(sol[0]), float(k * sol[0] + b))
        newline.AddPoint(float(sol[1]), float(k * sol[1] + b))
        feature.SetGeometry(newline)
        # feature.SetField('field_name', 'field_value')
        layer.CreateFeature(feature)

    ds = None
