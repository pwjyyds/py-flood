# coding:utf-8
"""
传统的种子蔓延法，该算法有一个不合理之处，就是没有考虑水的流动性，例如A种子点位于支流a，
B种子点位于支流b,两个支流位置很远汇入主流，A种子点的水位是20，B水位是5，理论上B被淹区域很少，
但是这个算法会假设A的水也能演到B处，这显然没有考虑水流的方向性，不贴合真是情况
"""


from osgeo import gdal, ogr, osr
import numpy as np
import matplotlib.pyplot as plt
from collections import deque

# 解决GDAL中文路径和编码问题
gdal.SetConfigOption("GDAL_FILENAME_IS_UTF8", "YES")
gdal.SetConfigOption("SHAPE_ENCODING", "UTF-8")
plt.rcParams["font.family"] = ["SimHei", "Microsoft YaHei", "PingFang SC", "Heiti SC"]  # 支持中文的字体
plt.rcParams["axes.unicode_minus"] = False  # 解决负号显示为方块的问题
plt.rcParams["font.size"] = 10  # 可选：设置默认字体大小


class DEMFloodAnalysis:
    def __init__(self, dem_path):
        """
        初始化DEM洪水分析类
        :param dem_path: DEM数据路径（支持tif、img等格式）
        """
        self.dem_path = dem_path
        self.ds = gdal.Open(dem_path)
        if self.ds is None:
            raise FileNotFoundError(f"无法打开DEM文件：{dem_path}")

        # 读取DEM基本信息
        self.geotrans = self.ds.GetGeoTransform()  # 地理变换参数
        self.proj = self.ds.GetProjection()  # 投影信息（WKT格式）
        self.proj_ref = osr.SpatialReference(wkt=self.proj)  # 投影参考对象
        self.rows = self.ds.RasterYSize  # 行数
        self.cols = self.ds.RasterXSize  # 列数
        self.band = self.ds.GetRasterBand(1)  # 第一波段
        self.no_data = self.band.GetNoDataValue()  # 无效值

        # 读取高程数据为numpy数组
        self.dem_array = self.band.ReadAsArray().astype(np.float32)

    def geo2pixel(self, lon, lat):
        """
        地理坐标（经纬度/平面坐标）转换为栅格行列号
        :param lon: 经度/横坐标
        :param lat: 纬度/纵坐标
        :return: (行号, 列号)
        """
        x_origin = self.geotrans[0]
        y_origin = self.geotrans[3]
        pixel_width = self.geotrans[1]
        pixel_height = self.geotrans[5]

        col = int((lon - x_origin) / pixel_width)
        row = int((lat - y_origin) / pixel_height)

        # 检查行列号是否在有效范围内
        if 0 <= row < self.rows and 0 <= col < self.cols:
            return row, col
        else:
            raise ValueError(f"坐标({lon}, {lat})超出DEM范围")

    def pixel2geo(self, row, col):
        """
        栅格行列号转换为地理坐标
        :param row: 行号
        :param col: 列号
        :return: (经度/横坐标, 纬度/纵坐标)
        """
        x = self.geotrans[0] + col * self.geotrans[1] + self.geotrans[2] * row
        y = self.geotrans[3] + row * self.geotrans[5] + self.geotrans[4] * col
        return x, y

    def read_seed_from_shp(self, shp_path, x_field="POINT_X", y_field="POINT_Y", water_level_field="h100"):
        """
        从SHP点文件读取种子点，提取X/Y坐标和对应的水位字段值
        :param shp_path: SHP点文件路径
        :param x_field: 存储X坐标的字段名（默认POINT_X）
        :param y_field: 存储Y坐标的字段名（默认POINT_Y）
        :param water_level_field: 存储水位值的字段名（默认h100）
        :return: 种子点列表 [(x1, y1, level1), (x2, y2, level2), ...]
        """
        # 打开SHP文件
        shp_ds = ogr.Open(shp_path)
        if shp_ds is None:
            raise FileNotFoundError(f"无法打开SHP文件：{shp_path}")

        # 获取第一个图层
        layer = shp_ds.GetLayer(0)

        # 检查字段是否存在
        layer_defn = layer.GetLayerDefn()
        field_names = [layer_defn.GetFieldDefn(i).GetName() for i in range(layer_defn.GetFieldCount())]
        required_fields = [x_field, y_field, water_level_field]
        for field in required_fields:
            if field not in field_names:
                raise ValueError(f"SHP文件中不存在{field}字段，可用字段：{field_names}")

        # 读取种子点坐标和水位值
        seed_points = []
        for feature in layer:
            # 提取X/Y坐标和水位值
            x = feature.GetField(x_field)
            y = feature.GetField(y_field)
            water_level = feature.GetField(water_level_field)

            # 检查值是否有效
            invalid = False
            for val in [x, y, water_level]:
                if val is None or np.isnan(val) or np.isinf(val):
                    invalid = True
                    break
            if invalid:
                print(f"警告：跳过无效值的要素（{x_field}={x}, {y_field}={y}, {water_level_field}={water_level}）")
                continue

            # 可选：坐标投影转换（若SHP与DEM投影不一致）
            shp_proj = layer.GetSpatialRef()
            if shp_proj and not shp_proj.IsSame(self.proj_ref):
                try:
                    # 创建坐标转换对象
                    transform = osr.CoordinateTransformation(shp_proj, self.proj_ref)
                    # 转换坐标（注意：ogr的TransformPoint返回(x, y, z)，取前两位）
                    x, y, _ = transform.TransformPoint(x, y)
                except Exception as e:
                    print(f"警告：种子点({x}, {y})投影转换失败，已跳过：{e}")
                    continue

            seed_points.append((x, y, water_level))

        if not seed_points:
            raise ValueError("SHP文件中未读取到有效种子点")

        # 打印投影检查信息
        shp_proj = layer.GetSpatialRef()
        if shp_proj and not shp_proj.IsSame(self.proj_ref):
            print(f"提示：SHP与DEM投影不一致，但已尝试自动转换坐标")
        else:
            print(f"提示：SHP与DEM投影一致")

        return seed_points

    def flood_by_seed(self, seed_points, neighbor_type=4, level_type="relative"):
        """
        种子蔓延法实现洪水淹没分析（支持每个种子点自定义水位）
        :param seed_points: 种子点列表，格式为[(x1, y1, level1), (x2, y2, level2), ...]
        :param neighbor_type: 邻域类型，4为四邻域，8为八邻域
        :param level_type: 水位类型，"absolute"为绝对水位（直接用level作为阈值），"relative"为相对水位（种子点高程+level）
        :return: 淹没范围数组（1为淹没，0为未淹没，-1为无效值）
        """
        # 初始化淹没数组
        flood_array = np.zeros_like(self.dem_array, dtype=np.int8)
        flood_array[self.dem_array == self.no_data] = -1  # 标记无效值

        # 定义邻域偏移量
        if neighbor_type == 4:
            offsets = [(-1, 0), (1, 0), (0, -1), (0, 1)]  # 四邻域
        elif neighbor_type == 8:
            offsets = [(-1, 0), (1, 0), (0, -1), (0, 1),
                       (-1, -1), (-1, 1), (1, -1), (1, 1)]  # 八邻域
        else:
            raise ValueError("邻域类型仅支持4或8")

        # 初始化队列（BFS）
        queue = deque()

        # 处理种子点
        for (x, y, level) in seed_points:
            # 转换为栅格行列号
            try:
                row, col = self.geo2pixel(x, y)
            except ValueError as e:
                print(f"警告：种子点({x}, {y})超出DEM范围，已跳过")
                continue

            # 检查种子点是否有效
            if flood_array[row, col] == -1:
                print(f"警告：种子点({x}, {y})位于无效值区域，已跳过")
                continue

            # 计算淹没阈值
            if level_type == "absolute":
                # 绝对水位：直接用h100作为阈值
                flood_threshold = level
            elif level_type == "relative":
                # 相对水位：种子点高程 + h100
                seed_elevation = self.dem_array[row, col]
                flood_threshold = seed_elevation + level
            else:
                raise ValueError("水位类型仅支持absolute（绝对）或relative（相对）")

            # 检查阈值是否合理（避免低于种子点高程）
            if flood_threshold < self.dem_array[row, col] - 1e-6:
                print(f"警告：种子点({x}, {y})的淹没阈值({flood_threshold})低于自身高程({self.dem_array[row, col]})，已跳过")
                continue

            # 标记种子点为淹没，并加入队列
            flood_array[row, col] = 1
            queue.append((row, col, flood_threshold))

        if not queue:
            raise RuntimeError("无有效种子点，无法执行淹没分析")

        # BFS实现种子蔓延
        while queue:
            row, col, threshold = queue.popleft()

            # 遍历邻域
            for dr, dc in offsets:
                new_row = row + dr
                new_col = col + dc

                # 检查行列号是否在有效范围内
                if 0 <= new_row < self.rows and 0 <= new_col < self.cols:
                    # 检查是否为未淹没的有效栅格
                    if flood_array[new_row, new_col] == 0 and self.dem_array[new_row, new_col] != self.no_data:
                        # 检查高程是否低于淹没阈值
                        if self.dem_array[new_row, new_col] <= threshold:
                            flood_array[new_row, new_col] = 1
                            queue.append((new_row, new_col, threshold))

        return flood_array

    def save_flood_result(self, flood_array, output_path):
        """
        保存淹没范围为栅格文件
        :param flood_array: 淹没范围数组
        :param output_path: 输出文件路径（如.tif）
        """
        # 创建输出栅格
        driver = gdal.GetDriverByName("GTiff")
        out_ds = driver.Create(
            output_path,
            self.cols,
            self.rows,
            1,
            gdal.GDT_Byte  # 用Byte类型存储（0-255）
        )

        # 设置地理变换和投影
        out_ds.SetGeoTransform(self.geotrans)
        out_ds.SetProjection(self.proj)

        # 写入数据（将-1转换为255作为无效值）
        out_band = out_ds.GetRasterBand(1)
        out_array = flood_array.astype(np.uint8)
        out_array[out_array == -1] = 255
        out_band.WriteArray(out_array)
        out_band.SetNoDataValue(255)  # 设置无效值

        # 刷新并关闭
        out_band.FlushCache()
        out_ds.FlushCache()
        del out_ds

    def visualize_result(self, flood_array):
        """
        可视化淹没结果
        :param flood_array: 淹没范围数组
        """
        plt.figure(figsize=(12, 8))

        # 绘制DEM高程
        plt.subplot(1, 2, 1)
        im1 = plt.imshow(self.dem_array, cmap="terrain")
        plt.colorbar(im1, label="高程（米）")
        plt.title("DEM高程图")
        plt.axis("off")

        # 绘制淹没范围
        plt.subplot(1, 2, 2)
        im2 = plt.imshow(flood_array, cmap="Blues", vmin=-1, vmax=1)
        plt.colorbar(im2, ticks=[-1, 0, 1], label="淹没状态（-1：无效值，0：未淹没，1：淹没）")
        plt.title("洪水淹没范围")
        plt.axis("off")

        plt.tight_layout()
        plt.show()


# ---------------------- 测试代码 ----------------------
if __name__ == "__main__":
    # 1. 初始化DEM洪水分析对象
    dem_path = r"E:\College\project\GD\geoData\Input\dem.tif"  # 替换为你的DEM文件路径
    shp_path = r"E:\College\project\GD\geoData\Input\seed.shp"  # 替换为你的SHP点文件路径
    # dem_path = r"E:\College\project\GD\geoData\testdem.tif"  # 替换为你的DEM文件路径
    # shp_path = r"E:\College\project\GD\geoData\testseed.shp"  # 替换为你的SHP点文件路径
    flood_analyzer = DEMFloodAnalysis(dem_path)

    # 2. 从SHP文件读取种子点（含h100水位字段）
    seed_points = flood_analyzer.read_seed_from_shp(
        shp_path=shp_path,
        x_field="POINT_X",
        y_field="POINT_Y",
        water_level_field="z100"  # 读取h100字段作为水位
    )
    print(f"成功读取{len(seed_points)}个有效种子点")

    # 3. 执行种子蔓延法淹没分析
    flood_array = flood_analyzer.flood_by_seed(
        seed_points=seed_points,
        neighbor_type=4,
         # 若h100是相对水位，改为"relative"
    )

    # 4. 保存淹没结果
    output_path = "flood_result.tif"
    flood_analyzer.save_flood_result(flood_array, output_path)
    print(f"淹没结果已保存至：{output_path}")

    # 5. 可视化结果
    flood_analyzer.visualize_result(flood_array)
