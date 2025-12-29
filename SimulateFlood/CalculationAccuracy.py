
import geopandas as gpd
from shapely.ops import unary_union
from shapely.validation import make_valid


# ===================== 1. 配置文件路径 =====================
SIM_SHP = r"L:\College\project\GD\ningguo_geodata\Output\result_z1002.shp" # 算法模拟的shp
TRUE_SHP = r"L:\floods\data\sen1_analyse\ningguo_0810_water.shp"  # 真实观测的shp
STUDY_SHP = r"L:\College\project\GD\ningguo_geodata\Input\basin2.shp"  # 例："研究区边界.shp"
TARGET_CRS = "EPSG:32650"  # 你的UTM投影（已确定的）

# ===================== 2. 读取并预处理研究区 =====================
# 读取研究区，合并多个面为一个整体几何
study_gdf = gpd.read_file(STUDY_SHP)
study_gdf = study_gdf.to_crs(TARGET_CRS)  # 统一到UTM投影
study_union = unary_union(study_gdf['geometry'].apply(make_valid))  # 合并多个面

# ===================== 3. 裁剪真实/模拟数据到研究区 =====================
# 读取真实数据 → 裁剪到研究区 → 修复几何
true_gdf = gpd.read_file(TRUE_SHP).to_crs(TARGET_CRS)
true_clipped = gpd.clip(true_gdf, study_union)  # 只保留研究区内的真实范围
true_clipped['geometry'] = true_clipped['geometry'].apply(make_valid)
true_final = unary_union(true_clipped['geometry'])  # 合并裁剪后的真实范围

# 读取模拟数据 → 裁剪到研究区 → 修复几何
sim_gdf = gpd.read_file(SIM_SHP).to_crs(TARGET_CRS)
sim_clipped = gpd.clip(sim_gdf, study_union)  # 只保留研究区内的模拟范围
sim_clipped['geometry'] = sim_clipped['geometry'].apply(make_valid)
sim_final = unary_union(sim_clipped['geometry'])  # 合并裁剪后的模拟范围

# ===================== 4. 计算研究区内的指标（单位：平方千米） =====================
true_area = true_final.area / 10**6
sim_area = sim_final.area / 10**6
inter_area = true_final.intersection(sim_final).area / 10**6
union_area = true_final.union(sim_final).area / 10**6

# 核心指标
iou = inter_area / union_area if union_area != 0 else 0
pa = inter_area / true_area if true_area != 0 else 0
ua = inter_area / sim_area if sim_area != 0 else 0
f1 = 2*(pa*ua)/(pa+ua) if (pa+ua)!=0 else 0

# Kappa系数（基于研究区）
study_area = study_union.area / 10**6
tn = study_area - union_area  # 研究区内“未淹没且未模拟”的区域
tp = inter_area
fp = sim_area - tp
fn = true_area - tp
total = tp + fp + fn + tn
po = (tp + tn)/total if total !=0 else 0
pe = ((tp+fn)*(tp+fp) + (fn+tn)*(fp+tn))/(total**2) if total !=0 else 0
kappa = (po - pe)/(1 - pe) if (1 - pe)!=0 else 0

# 输出结果
print("===== 研究区内的洪水模拟精度指标 =====")
print(f"研究区面积: {study_area:.2f} 平方千米")
print(f"真实淹没面积（研究区内）: {true_area:.2f} 平方千米")
print(f"模拟淹没面积（研究区内）: {sim_area:.2f} 平方千米")
print(f"IOU: {iou:.4f}")
print(f"生产者精度（PA）: {pa:.4f}")
print(f"用户精度（UA）: {ua:.4f}")
print(f"F1分数: {f1:.4f}")
print(f"Kappa系数: {kappa:.4f}")