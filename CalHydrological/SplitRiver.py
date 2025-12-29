import geopandas as gpd
from shapely.ops import unary_union


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
            {'point_id': []},
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
        intersection_gdf.to_file(r"L:\College\project\GD\ningguo_geodata\Input\seeee.shp", encoding='utf-8')


# ----------------------------------------------------------------------


RIVER_FILE = r"L:\College\project\GD\ningguo_geodata\Input\river.shp"  # 包含一条河流的LineString
CROSS_SECTION_FILE = r"L:\College\project\GD\ningguo_geodata\Input\dxm.shp"  # 包含三条断面线的MultiLineString/LineString集合
#
OUTPUT_FILE = r"L:\College\project\GD\ningguo_geodata\Input\river22.shp"
#
result_gdf = split_river_by_cross_sections(
    river_path=RIVER_FILE,
    cross_section_path=CROSS_SECTION_FILE,
    output_path=OUTPUT_FILE
)
