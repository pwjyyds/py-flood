# coding:utf-8
"""
配置文件
"""
import os

p_all = [0.2, 0.1, 0.05, 0.02, 0.01]  # 频率 0.2：5年一遇、0.1：十年、0.05：二十年、0.02：五十年、0.01：百年
fields_z = ['z5', 'z10', 'z20', 'z50', 'z100']  # 每个频率的字段名称，上升水位
fields_h = ['h5', 'h10', 'h20', 'h50', 'h100']  # 上升水位+高程
fields_Qm = ['Qm5', 'Qm10', 'Qm20', 'Qm50', 'Qm100']  # 上升水位+高程
fields_L = ['L5', 'L10', 'L20', 'L50', 'L100']  # 湿周长
fields_S = ['S5', 'S10', 'S20', 'S50', 'S100']  # 断面积
Qm_0 = 1000  # 洪峰流量初试阈值
raster_pixel_width = 30  # 栅格大小，默认30m

# 字段名
dmx_field = {
    "J": "J",
    "n0": "n0",  # 糙度
    "ObjectID": "ObjectID",  # 每个断面线的ID
}
dmxPoints_field = {
    "ObjectID": "ObjectID",  # 每个点的ID
    "DmxID": "DmxID",  # 属于哪一个断面线FID
    "J": "J",
    "n0": "n0",
    "DemValue": "DemValue",
}

unit_fields = {
    'H_6': 'H_6',
    'Cv_6': 'Cv_6',
    'H_24': 'H_24',
    'Cv_24': 'Cv_24',
    'n': 'n',  # 暴雨衰减指数
    'Sp': 'Sp',  # 暴雨雨力
    'Qm': 'Qm',  # 洪峰流量
    'dmxID': 'basin'  # 面内的断面线ID
}

# 当前项目主目录
project_dir = os.getcwd()  # E:\College\project\GD
geodata_dir = os.path.join(project_dir, "ningguo_geodata")
# L:\College\project\GD\ningguo_geodata\Input\river.shp
# 输入数据地址
input_dir = os.path.join(geodata_dir,"Input")
output_dir = os.path.join(geodata_dir,"Output")

dem_dir = os.path.join(input_dir, "dem.tif")  # DEM
dmx_dir = os.path.join(input_dir, "dxm.shp")  # 断面线
river_dir = os.path.join(input_dir, "river.shp")  # 河流线
# riverDiv_dir = os.path.join(input_dir, "river_div.shp")  # 打断的河流线
# slope_dir = os.path.join(input_dir, "slope.tif")  # 坡度（百分比）
unit_dir = os.path.join(input_dir, "ZoneUnit.shp")  # 范围面
basin_dir = os.path.join(input_dir, "basin3.shp")  # 流域面
# flowDir_dir = os.path.join(input_dir, "flowDir.tif")  # 流域面
seed_dir = os.path.join(input_dir, "seed.shp")  # 种子点

# print(river_dir)

# 中间数据地址
dmx_raster_dir = os.path.join(geodata_dir, "dmx_raster1")
dmx_points = os.path.join(geodata_dir, "dmx_points.shp")
temp_flowDir = os.path.join(output_dir, 'temp_flowDir.tif')  # 所有子流域
temp_units = os.path.join(output_dir, 'basin3.shp')  # 所有子流域
# temp_units = os.path.join(output_dir, 'temp_units.shp')  # 所有子流域
temp_seed = os.path.join(output_dir, 'temp_seed.shp')  # 所有子流域
temp_slope = os.path.join(output_dir, 'temp_slope.tif')  # 所有子流域
temp_river_div = os.path.join(output_dir, 'temp_river_div.shp')  # 所有子流域
temp_dmxLong = os.path.join(output_dir, 'temp_dmxLong.shp')  # 延长的断面线
