"""
Basic 11: 对称性与监视器 — 用更少网格得到同样物理量

本示例演示三个常用技巧：
1. 使用 mp.Mirror 对称性减少计算区域；
2. 用 at_every 回调记录点场时间序列；
3. 用 get_array 导出二维场切片和介电常数分布。

物理模型：2D 真空中的线偏振高斯脉冲。结构和光源关于 y=0 偶对称，
因此只计算上半空间也能得到完整结果。为了让脚本对新手友好，默认仍
使用完整区域；把 USE_SYMMETRY 改为 True 可开启对称性。
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# ------------------------- 可调参数 -------------------------
resolution = 20              # 每个长度单位的网格数
fcen = 0.25                  # 中心频率，对应真空波长 4
fwidth = 0.12                # 脉冲带宽
sx, sy = 12, 8               # 计算区域尺寸
pml_thickness = 1.0
USE_SYMMETRY = False         # 改为 True 可使用 y 镜面对称性

# ------------------------- 建立仿真 -------------------------
cell = mp.Vector3(sx, sy, 0)
boundary_layers = [mp.PML(pml_thickness)]

sources = [mp.Source(
    src=mp.GaussianSource(frequency=fcen, fwidth=fwidth),
    component=mp.Ez,                         # 2D 标量极化，便于画图
    center=mp.Vector3(-4, 0),
    size=mp.Vector3(0, 2.5),                 # y 方向展开，近似平面波束
)]

symmetries = []
if USE_SYMMETRY:
    # Ez 对 y=0 镜像为偶对称；若源或结构破坏对称性，不能使用该设置。
    symmetries = [mp.Mirror(mp.Y, phase=+1)]

sim = mp.Simulation(
    cell_size=cell,
    boundary_layers=boundary_layers,
    sources=sources,
    symmetries=symmetries,
    resolution=resolution,
)

# ------------------------- 时间监视器 -------------------------
probe_point = mp.Vector3(1.0, 0.0)
time_trace = []

def record_probe(sim):
    """每隔一段时间记录探针点 Ez；回调函数由 MEEP 在 run 中调用。"""
    time_trace.append(sim.get_field_point(mp.Ez, probe_point))

# 等源结束后再等场衰减，避免固定 until 过短或过长。
sim.run(
    mp.at_every(0.5, record_probe),
    until_after_sources=mp.stop_when_fields_decayed(40, mp.Ez, probe_point, 1e-7),
)

print(f"记录到 {len(time_trace)} 个时间采样点")

# ------------------------- 导出空间切片 -------------------------
ez = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Ez)
eps = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)

# ------------------------- 可视化 -------------------------
fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))

axes[0].plot(np.arange(len(time_trace)) * 0.5, np.real(time_trace))
axes[0].set_xlabel("时间 / MEEP 单位")
axes[0].set_ylabel("Ez")
axes[0].set_title("探针点时间信号")
axes[0].grid(alpha=0.3)

im1 = axes[1].imshow(eps.T, origin="lower", cmap="gray", extent=[-sx/2, sx/2, -sy/2, sy/2])
axes[1].set_title("介电常数分布")
axes[1].set_xlabel("x")
axes[1].set_ylabel("y")
fig.colorbar(im1, ax=axes[1], fraction=0.046)

im2 = axes[2].imshow(np.real(ez).T, origin="lower", cmap="RdBu", extent=[-sx/2, sx/2, -sy/2, sy/2])
axes[2].set_title("最终 Ez 场切片")
axes[2].set_xlabel("x")
axes[2].set_ylabel("y")
fig.colorbar(im2, ax=axes[2], fraction=0.046)

plt.tight_layout()
plt.savefig("symmetry_and_monitors.png", dpi=160)
print("图像已保存到 symmetry_and_monitors.png")
plt.show()
