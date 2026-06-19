"""
Example 03: 介质平板的传输与反射

计算一个电介质平板（ε=12）的传输和反射系数。
使用两个通量监视器分别测量透射和反射功率。
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# --- 参数 ---
resolution = 40
freq       = 0.4     # 中心频率
fwidth     = 0.3     # 频率范围
dpml       = 1.0     # PML 厚度
pad        = 2.0     # 源/PML 间距

# 几何体
slab_center = mp.Vector3(0, 0, 0)
slab_thickness = 1.0   # 1 μm 厚
slab_epsilon  = 12.0   # ε = 12 (类似硅)

cell_length = 2 * (dpml + pad + 2.0)
cell = mp.Vector3(cell_length, 0, 0)

# 通量监视器位置
refl_pos = -dpml - pad + 0.5
tran_pos =  dpml + pad - 0.5

# --- 几何体 ---
geometry = [mp.Block(
    size=mp.Vector3(slab_thickness, mp.inf, mp.inf),
    center=slab_center,
    material=mp.Medium(epsilon=slab_epsilon)
)]

# --- 光源 (从左侧入射) ---
source_center = mp.Vector3(-dpml - pad + 2.5, 0, 0)
sources = [mp.Source(
    mp.GaussianSource(frequency=freq, fwidth=fwidth),
    component=mp.Ey,
    center=source_center
)]

# --- 监视器 ---
# 反射功率：光源左边的通量区域
refl_fr = mp.FluxRegion(
    center=mp.Vector3(refl_pos, 0, 0),
    size=mp.Vector3(0, mp.inf, mp.inf)
)

# 透射功率：平版右边的通量区域
tran_fr = mp.FluxRegion(
    center=mp.Vector3(tran_pos, 0, 0),
    size=mp.Vector3(0, mp.inf, mp.inf)
)

refl = sim.add_flux(freq, fwidth, 1, refl_fr)
tran = sim.add_flux(freq, fwidth, 1, tran_fr)

# --- 运行 ---
sim = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    sources=sources,
    boundary_layers=[mp.PML(dpml)],
    resolution=resolution,
)

sim.run(until=200)

# --- 计算反射和透射 ---
freqs = np.linspace(freq - fwidth, freq + fwidth, 100)
refl_flux = np.array(mp.get_fluxes(refl))
tran_flux = np.array(mp.get_fluxes(tran))

R = -refl_flux / (refl_flux + tran_flux)  # 反射率
T =  tran_flux / (refl_flux + tran_flux)  # 透射率

# --- 空仿真（无平板）作为参考 ---
# 先跑一次空的仿真来归一化
geometry_empty = []
sim_empty = mp.Simulation(
    cell_size=cell,
    geometry=geometry_empty,
    sources=sources,
    boundary_layers=[mp.PML(dpml)],
    resolution=resolution,
)
refl0 = sim_empty.add_flux(freq, fwidth, 1, refl_fr)
tran0 = sim_empty.add_flux(freq, fwidth, 1, tran_fr)
sim_empty.run(until=200)
tran0_flux = np.array(mp.get_fluxes(tran0))

# 归一化透射系数
T_norm = tran_flux / tran0_flux
R_norm = -refl_flux / tran0_flux

print(f"频率 {freq:.2f}: 反射率 R = {R_norm[0]:.4f}, 透射率 T = {T_norm[0]:.4f}")

# --- 解析解验证 ---
# 法布里-珀罗公式
n1, n2 = 1.0, np.sqrt(slab_epsilon)
k = 2 * np.pi * freq     # 波矢
eta = n1 / n2             # 折射率比
delta = k * slab_thickness * n2

# 反射系数
r = (1 - eta**2) * (1 - np.exp(2j * delta)) / \
    ((1 + eta)**2 - (1 - eta)**2 * np.exp(2j * delta))
R_ana = abs(r)**2

print(f"解析解：R(ε=12) = {R_ana:.4f}")
print(f"   注意：meep 的 2D 仿真包含了倾斜入射，结果会略有偏差")
