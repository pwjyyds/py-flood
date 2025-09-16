# coding:utf-8
"""
main file
"""
from osgeo import gdal, ogr

from CalHydrological import Pretreatment, DesignRainstorm, Peak_discharge, FloodPeakWaterLevel
import setting
from setting import p_all, fields_h, fields_z, Qm_0, fields_Qm, fields_L, fields_S
from SimulateFlood import ZZMYF as ZZ


def main():
    dmx_points = Pretreatment.main(dmx, slope, river, riverDiv, data, setting.flowDir_dir, seed)
    for index, p in enumerate(p_all):
        Sp, n = DesignRainstorm.main(p, unit)  # Calculate the design rainstorm of the basin
        Peak_discharge.main(unit, river, slope, Sp, n, fields_Qm[index])  # Calculate the design flood peak of the basin
        FloodPeakWaterLevel.main(dmx_points, dmx,
                                 fields_z[index], fields_h[index], fields_Qm[index], fields_L[index],
                                 fields_S[index])  # Calculate the design water level of the cross-sectional line
        ZZ.main(data, dmx, seed, fields_z[index], fields_h[index])


if __name__ == '__main__':
    """
    输入数据：DEM、断面线、种子点（有洪峰水位）
    """
    data = gdal.Open(setting.dem_dir)  # DEM
    dmx = ogr.Open(setting.dmx_dir, 1)  # Cross-sectional line
    river = ogr.Open(setting.river_dir, 1)  # River line
    riverDiv = ogr.Open(setting.riverDiv_dir, 1)  # The broken river line
    slope = gdal.Open(setting.slope_dir)  # Slope
    unit = ogr.Open(setting.unit_dir, 1)  # The scope of the study area
    seed = ogr.Open(setting.seed_dir, 1)  # Seed point

    geotransform = data.GetGeoTransform()
    setting.raster_pixel_width = geotransform[1]

    main()
