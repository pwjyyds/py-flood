# coding:utf-8
"""
针对v1的缺点，改进为v2版本
地形连通性：A、B 种子点位于不同支流（不同流域），被分水岭阻隔，即使高程低于阈值，水也无法跨流域流动。
水流方向性：水的流动遵循地形坡度，只能沿下坡 / 平坡方向漫延，而非任意邻域的平面蔓延。
通过水文分析预处理和地形约束的蔓延规则解决上述问题，核心步骤为：
DEM 填洼：消除 DEM 中的洼地（凹陷点），保证水流路径的连续性（水文分析的必要步骤）。
计算水流方向：采用 D8 算法计算每个栅格的水流去向（流向 8 个邻域中的哪一个），反映地形的水流趋势。
划分流域范围：基于水流方向，确定每个种子点的集水流域（分水岭内的地形连通区域），不同流域的种子点互不干扰。
地形约束的种子蔓延：在各流域内，结合高程阈值和水流方向进行 BFS 蔓延，仅允许水沿合理的地形路径流动。
"""

from osgeo import gdal, ogr, osr
import numpy as np
import matplotlib.pyplot as plt
from collections import deque
import os
from setting import *
import geopandas as gpd

# 解决GDAL中文路径和编码问题
gdal.SetConfigOption("GDAL_FILENAME_IS_UTF8", "YES")
gdal.SetConfigOption("SHAPE_ENCODING", "UTF-8")
plt.rcParams["font.family"] = ["SimHei", "Microsoft YaHei", "PingFang SC", "Heiti SC"]  # 支持中文的字体
plt.rcParams["axes.unicode_minus"] = False  # 解决负号显示为方块的问题
plt.rcParams["font.size"] = 10  # 可选：设置默认字体大小

# 4/8邻域定义
FOUR_OFFSETS = [(-1, 0), (1, 0), (0, -1), (0, 1)]
EIGHT_OFFSETS = [(-1, 0), (1, 0), (0, -1), (0, 1),
                 (-1, -1), (-1, 1), (1, -1), (1, 1)]


class DEMFloodAnalysis:
    def __init__(self, dem_path):
        """初始化DEM洪水分析类"""
        self.dem_path = dem_path
        self.ds = gdal.Open(dem_path)
        if self.ds is None:
            raise FileNotFoundError(f"无法打开DEM文件：{dem_path}")

        # 读取DEM基础信息
        self.geotrans = self.ds.GetGeoTransform()
        self.proj = self.ds.GetProjection()
        self.proj_ref = osr.SpatialReference(wkt=self.proj)
        self.rows = self.ds.RasterYSize
        self.cols = self.ds.RasterXSize
        self.band = self.ds.GetRasterBand(1)
        self.no_data = self.band.GetNoDataValue() if self.band.GetNoDataValue() is not None else -9999
        # self.dem_dtype = self.band.DataType  # 记录原DEM的数据类型

        # 读取DEM数组（处理无效值）
        self.dem_array = self.band.ReadAsArray().astype(np.float32)
        self.dem_array[self.dem_array == self.no_data] = np.nan

        # 流域相关变量（外部输入）
        self.watershed_raster = None  # 流域栅格数组（值为流域ID，0为无流域，-1为无效值）
        self.seed_watershed_map = None  # 种子点-流域ID映射：{(x,y,level): watershed_id}


    def geo2pixel(self, lon, lat):
        """地理坐标转栅格行列号"""
        x_origin = self.geotrans[0]
        y_origin = self.geotrans[3]
        pixel_width = self.geotrans[1]
        pixel_height = self.geotrans[5]

        col = int((lon - x_origin) / pixel_width)
        row = int((lat - y_origin) / pixel_height)

        if 0 <= row < self.rows and 0 <= col < self.cols:
            return row, col
        else:
            raise ValueError(f"坐标({lon}, {lat})超出DEM范围")

    def pixel2geo(self, row, col):
        """栅格行列号转地理坐标"""
        x = self.geotrans[0] + col * self.geotrans[1] + self.geotrans[2] * row
        y = self.geotrans[3] + row * self.geotrans[5] + self.geotrans[4] * col
        return x, y

    def rasterize_watershed(self, watershed_shp_path, id_field="ID"):
        """
        将流域面SHP栅格化为与DEM匹配的流域栅格
        :param watershed_shp_path: 流域面SHP文件路径
        :param id_field: 流域ID字段名（SHP面要素的唯一标识）
        :return: 流域栅格数组
        """
        # 1. 打开流域SHP
        shp_ds = ogr.Open(watershed_shp_path)
        if shp_ds is None:
            raise FileNotFoundError(f"无法打开流域SHP文件：{watershed_shp_path}")
        layer = shp_ds.GetLayer(0)

        # 2. 检查ID字段是否存在
        layer_defn = layer.GetLayerDefn()
        field_names = [layer_defn.GetFieldDefn(i).GetName() for i in range(layer_defn.GetFieldCount())]
        if id_field not in field_names:
            raise ValueError(f"流域SHP中不存在{id_field}字段，可用字段：{field_names}")

        # 3. 创建临时栅格（内存中）用于栅格化
        driver = gdal.GetDriverByName("MEM")
        tmp_ds = driver.Create("", self.cols, self.rows, 1, gdal.GDT_Int32)
        tmp_ds.SetGeoTransform(self.geotrans)
        tmp_ds.SetProjection(self.proj)
        tmp_band = tmp_ds.GetRasterBand(1)
        tmp_band.Fill(0)  # 初始值为0（无流域）
        tmp_band.SetNoDataValue(-1)

        # 4. 栅格化流域面（按ID字段赋值）
        gdal.RasterizeLayer(tmp_ds, [1], layer,
                            options=[f"ATTRIBUTE={id_field}",
                                     "ALL_TOUCHED=TRUE"])  # ALL_TOUCHED保证面边缘栅格被包含

        # 5. 读取栅格化结果
        watershed_raster = tmp_band.ReadAsArray()
        # 处理DEM无效值区域
        watershed_raster[np.isnan(self.dem_array)] = -1
        self.watershed_raster = watershed_raster

        print(f"成功栅格化流域面，共识别到{len(np.unique(watershed_raster[watershed_raster > 0]))}个流域")
        return watershed_raster

    def match_seed_to_watershed(self, seed_points):
        """
        匹配种子点到所属流域（空间叠加：点在面内）
        :param seed_points: 种子点列表 [(x1,y1,level1), (x2,y2,level2), ...]
        :return: 种子点-流域ID映射
        """
        if self.watershed_raster is None:
            raise RuntimeError("请先调用rasterize_watershed加载流域栅格")

        seed_watershed_map = {}
        for seed in seed_points:
            x, y, level = seed
            # print(x,y,level)
            try:
                row, col = self.geo2pixel(x, y)
            except ValueError as e:
                print(f"警告：种子点({x}, {y})超出DEM范围，跳过")
                continue

            # 获取种子点所在流域ID
            # print(self.watershed_raster[row, col])
            watershed_id = self.watershed_raster[row, col]
            seed_watershed_map[seed] = watershed_id

        if not seed_watershed_map:
            raise ValueError("无种子点匹配到有效流域")

        self.seed_watershed_map = seed_watershed_map
        return seed_watershed_map

    def read_seed_from_shp(self, shp_path, x_field="POINT_X", y_field="POINT_Y", water_level_field="z100"):
        """从SHP读取种子点（含X/Y坐标和水位字段）"""
        shp_ds = ogr.Open(shp_path)
        if shp_ds is None:
            raise FileNotFoundError(f"无法打开种子点SHP文件：{shp_path}")

        layer = shp_ds.GetLayer(0)

        # 检查字段
        layer_defn = layer.GetLayerDefn()
        field_names = [layer_defn.GetFieldDefn(i).GetName() for i in range(layer_defn.GetFieldCount())]
        required_fields = [x_field, y_field, water_level_field]
        for field in required_fields:
            if field not in field_names:
                raise ValueError(f"种子点SHP中不存在{field}字段，可用字段：{field_names}")

        # 读取种子点
        seed_points = []
        shp_proj = layer.GetSpatialRef()
        for feature in layer:
            x = feature.GetField(x_field)
            y = feature.GetField(y_field)
            level = feature.GetField(water_level_field)

            print("读取种子点",x,y,level)

            # 过滤无效值
            if None in [x, y, level] or np.isnan(x) or np.isnan(y) or np.isnan(level):
                print(f"警告：跳过无效值要素（{x_field}={x}, {y_field}={y}, {water_level_field}={level}）")
                continue

            # 投影转换
            if shp_proj and not shp_proj.IsSame(self.proj_ref):
                try:
                    transform = osr.CoordinateTransformation(shp_proj, self.proj_ref)
                    x, y, _ = transform.TransformPoint(x, y)
                except Exception as e:
                    print(f"警告：种子点({x}, {y})投影转换失败，跳过：{e}")
                    continue

            seed_points.append((x, y, level))

        if not seed_points:
            raise ValueError("种子点SHP中未读取到有效种子点")

        return seed_points

    def flood_by_seed(self, neighbor_type=4, level_type="absolute"):
        """
        基于外部流域面的种子蔓延法
        :param neighbor_type: 邻域类型（4/8）
        :param level_type: 水位类型（absolute/relative）
        :return: 全局淹没范围数组
        """
        if self.watershed_raster is None or self.seed_watershed_map is None:
            raise RuntimeError("请先加载流域栅格并匹配种子点到流域")

        # 1. 初始化全局淹没数组
        flood_array = np.zeros_like(self.dem_array, dtype=np.int8)
        flood_array[np.isnan(self.dem_array)] = -1  # -1：无效值，0：未淹没，1：淹没

        # 2. 定义邻域偏移量
        offsets = FOUR_OFFSETS if neighbor_type == 4 else EIGHT_OFFSETS

        # 3. 按流域分组处理种子点
        watershed_seeds = {}
        for seed, ws_id in self.seed_watershed_map.items():
            if ws_id not in watershed_seeds:
                watershed_seeds[ws_id] = []
            watershed_seeds[ws_id].append(seed)

        # 4. 对每个流域执行淹没计算
        for ws_id, seeds in watershed_seeds.items():
            print(f"处理流域{ws_id}，包含{len(seeds)}个种子点")

            # 生成当前流域的掩码（仅该流域的栅格为True）
            ws_mask = (self.watershed_raster == ws_id)

            # 初始化队列（BFS）
            queue = deque()

            # 处理当前流域的种子点
            for seed in seeds:
                x, y, level = seed
                row, col = self.geo2pixel(x, y)

                # 计算淹没阈值
                if level_type == "absolute":
                    threshold = level
                elif level_type == "relative":
                    # print("相对高程")
                    threshold = self.dem_array[row, col] + level
                else:
                    raise ValueError("水位类型仅支持absolute/relative")

                # 阈值合理性检查
                if threshold < self.dem_array[row, col] - 1e-6:
                    print(f"警告：种子点({x}, {y})阈值低于自身高程，跳过")
                    continue


                # 标记种子点为淹没并加入队列
                if flood_array[row, col] == 0 and ws_mask[row, col]:
                    flood_array[row, col] = 1
                    queue.append((row, col, threshold))

            # 5. 流域内BFS蔓延（仅在流域掩码内）
            while queue:
                row, col, threshold = queue.popleft()

                for dr, dc in offsets:
                    new_row = row + dr
                    new_col = col + dc

                    # 范围检查
                    if 0 <= new_row < self.rows and 0 <= new_col < self.cols:
                        # 仅在当前流域内、未淹没、有效栅格
                        if ws_mask[new_row, new_col] and flood_array[new_row, new_col] == 0 and not np.isnan(
                                self.dem_array[new_row, new_col]):
                            # 高程低于阈值则淹没
                            if self.dem_array[new_row, new_col] <= threshold:
                                flood_array[new_row, new_col] = 1
                                queue.append((new_row, new_col, threshold))

        return flood_array

    def save_flood_result(self, flood_array, output_path):
        """保存淹没结果为TIFF"""

        # if os.path.exists(output_path):
        #     os.remove(output_path)
        #     if os.path.exists(output_path + '.ovr'):
        #         os.remove(output_path + '.ovr')

        driver = gdal.GetDriverByName("GTiff")
        out_ds = driver.Create(
            output_path,
            self.cols,
            self.rows,
            1,
            gdal.GDT_Byte
        )
        out_ds.SetGeoTransform(self.geotrans)
        out_ds.SetProjection(self.proj)

        # 处理无效值（-1→255）
        out_array = flood_array.astype(np.uint8)
        out_array[out_array == -1] = 255
        out_band = out_ds.GetRasterBand(1)
        out_band.WriteArray(out_array)
        out_band.SetNoDataValue(255)

        # 刷新并关闭
        out_band.FlushCache()
        out_ds.FlushCache()
        del out_ds
        print(f"淹没结果已保存至：{output_path}")

    def visualize_result(self, flood_array):
        """可视化结果（DEM+流域+淹没范围）"""
        plt.figure(figsize=(18, 6))

        # 1. DEM高程
        plt.subplot(1, 3, 1)
        im1 = plt.imshow(self.dem_array, cmap="terrain")
        plt.colorbar(im1, label="高程（米）")
        plt.title("DEM高程图")
        plt.axis("off")

        # 2. 流域划分
        plt.subplot(1, 3, 2)
        im2 = plt.imshow(self.watershed_raster, cmap="tab20")
        plt.colorbar(im2, label="流域ID")
        plt.title("外部流域面栅格化")
        plt.axis("off")

        # 3. 淹没范围
        plt.subplot(1, 3, 3)
        im3 = plt.imshow(flood_array, cmap="Blues", vmin=-1, vmax=1)
        plt.colorbar(im3, ticks=[-1, 0, 1], label="淹没状态（-1：无效，0：未淹没，1：淹没）")
        plt.title("洪水淹没范围")
        plt.axis("off")

        plt.tight_layout()
        plt.show()


def line2point_center_with_xy(input_line_shp, output_point_shp):
    """
    将线SHP转换为中心点SHP，保留源数据所有字段，自动补充POINT_X/POINT_Y坐标字段
    :param input_line_shp: 输入线SHP文件路径（如"rivers.shp"）
    :param output_point_shp: 输出点SHP文件路径（如"river_centers.shp"）
    """
    line_gdf = gpd.read_file(input_line_shp, encoding='utf-8')
    line_gdf = line_gdf[~line_gdf.is_empty].reset_index(drop=True)
    point_gdf = line_gdf.copy()
    point_gdf['geometry'] = line_gdf['geometry'].centroid
    x_field = 'POINT_X'
    y_field = 'POINT_Y'
    point_x = point_gdf['geometry'].x
    point_y = point_gdf['geometry'].y
    point_gdf[x_field] = point_x
    point_gdf[y_field] = point_y
    point_gdf.to_file(output_point_shp, encoding='utf-8')


# ---------------------- 测试代码 ----------------------
def main(dem_path, dmx_shp_path, watershed_shp_path, level_field):
    # 1. 初始化分析对象
    # dem_path = r"E:\College\project\GD\geoData\Input\dem.tif"  # 替换为你的DEM文件路径
    # seed_shp_path = r"E:\College\project\GD\geoData\Input\seed.shp"  # 替换为你的SHP点文件路径
    # watershed_shp_path = r"E:\College\project\GD\geoData\Input\basin.shp"  # 替换为你的流域面SHP路径

    line2point_center_with_xy(dmx_shp_path, seed_dir)
    flood_analyzer = DEMFloodAnalysis(dem_path)

    # 2. 加载外部流域面（SHP）
    flood_analyzer.rasterize_watershed(watershed_shp_path, id_field="Id")  # ID为流域面的唯一标识字段

    # 3. 读取种子点（SHP）
    seed_points = flood_analyzer.read_seed_from_shp(
        seed_dir,
        x_field="POINT_X",
        y_field="POINT_Y",
        water_level_field=level_field
    )
    print(f"成功读取{len(seed_points)}个有效种子点")

    # 4. 匹配种子点到所属流域
    flood_analyzer.match_seed_to_watershed(seed_points)

    # 5. 执行淹没分析
    flood_array = flood_analyzer.flood_by_seed(
        neighbor_type=4,
        level_type="relative"  # 若h100为相对水位则改为"relative"
    )

    # 6. 保存结果
    output_path = os.path.join(output_dir, 'result_' + level_field + '.tif')
    flood_analyzer.save_flood_result(flood_array, output_path)

    # 7. 可视化
    # flood_analyzer.visualize_result(flood_array)
