# coding:utf-8
"""
Generate a dividing line between two adjacent seeds
"""

from osgeo import ogr


def main(inSeed,outPath):
    s = inSeed
    point_data = s.GetLayer(0)
    output_file = outPath
    driver = ogr.GetDriverByName("ESRI Shapefile")
    out_data_set = driver.CreateDataSource(output_file)
    out_layer = out_data_set.CreateLayer("line_data", geom_type=ogr.wkbLineString)
    # Obtain the quantity of point data
    num_points = point_data.GetFeatureCount()
    # Traverse the point data to generate linear data
    for i in range(num_points - 1):
        # Obtain two adjacent points
        point1 = point_data.GetFeature(i)
        point2 = point_data.GetFeature(i + 1)
        point1_x = point1.GetField('POINT_X')
        point1_y = point1.GetField('POINT_Y')
        point2_x = point2.GetField('POINT_X')
        point2_y = point2.GetField('POINT_Y')

        mid_x = (point1_x + point2_x) / 2
        mid_y = (point1_y + point2_y) / 2
        # print(point1_x, point1_y, point2_x, point2_y, mid_x, mid_y)
        # Fit the cross-sectional line y=k2x+b
        if (point1_x - point2_x) == 0:  # The river channel line section is parallel to the Y-axis, and the cross-section straight line is parallel to the X-axis,y=C
            # print("k=0")
            x1 = mid_x - 300
            x2 = mid_x + 300
            y1 = mid_y
            y2 = mid_y
            k2 = 0
        else:
            k1 = (point2_y - point1_y) / (point2_x - point1_x)  # The product of vertical slopes =-1
            if k1 == 0:  # The slope of the section line does not exist, that is, the section line is parallel to the Y-axis
                x1 = mid_x
                x2 = mid_x
                y1 = mid_y + 300
                y2 = mid_y - 300
                k2 = 0
            else:
                k2 = -1 / k1
                if k2 > 87 or k2 < -87:  # If the slope is greater than 87°, it is determined to be 90°
                    x1 = mid_x
                    x2 = mid_x
                    y1 = mid_y + 300
                    y2 = mid_y - 300
                    k2 = 9999
                else:
                    # print("k=", k2)
                    b = mid_y - k2 * mid_x
                    x1 = mid_x - 300
                    x2 = mid_x + 300
                    y1 = k2 * x1 + b
                    y2 = k2 * x2 + b
        # print(x1, y1, x2, y2)
        line_data = ogr.Geometry(ogr.wkbLineString)
        line_data.AddPoint(x1, y1)
        line_data.AddPoint(x2, y2)
        feature = ogr.Feature(out_layer.GetLayerDefn())
        feature.SetGeometry(line_data)
        out_layer.CreateFeature(feature)
        out_layer.CreateFeature(feature)

    out_data_set = None
