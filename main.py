# coding:utf-8
"""
main file
"""
from osgeo import gdal, ogr

from CalHydrological import Pretreatment, DesignRainstorm, Peak_discharge, FloodPeakWaterLevel
from SimulateFlood import SeedPropagationV2
from setting import *


def main():
    dmx_points = Pretreatment.main(dmx, dem, createBasin=False)
    for index, p in enumerate(p_all):
        Sp, n = DesignRainstorm.main(p, unit)  # Calculate the design rainstorm of the basin
        Peak_discharge.main(unit, river, Sp, n, fields_Qm[index])  # Calculate the design flood peak of the basin
        FloodPeakWaterLevel.main(dmx_points, dmx,fields_z[index], fields_h[index], fields_Qm[index], fields_L[index],fields_S[index])  # Calculate the design water level of the cross-sectional line

        SeedPropagationV2.main(dem_dir, dmx_dir, basin_dir,fields_z[index])


if __name__ == '__main__':
    """
    输入数据：DEM、断面线、种子点（有洪峰水位）
    """
    dem = gdal.Open(dem_dir)  # DEM
    dmx = ogr.Open(dmx_dir, 1)  # Cross-sectional line
    river = ogr.Open(river_dir, 1)  # River line
    unit = ogr.Open(unit_dir, 1)  # The scope of the study area

    geotransform = dem.GetGeoTransform()
    raster_pixel_width = geotransform[1]

    main()

