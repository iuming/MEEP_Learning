"""
Example 07: 介质波导模式分析

使用 MEEP 计算平板波导的导模。
通过本征模式求解器找到 TE 和 TM 模式的有效折射率。
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# --- 波导参数 ---
resolution = 64
n_core     = 2.0     # 芯层折射率
n_clad     = 1.0     # 包层折射率
w          = 1.0     # 波导宽度 (μm)
wavelength = 1.55    # 波长 (μm) — 通信波段
freq       = 1.0 / wavelength

# --- 计算区域 ---
cell_size = mp.Vector3(6, 6, 0)

geometry = [
    mp.Block(
        size=mp.Vector3(mp.inf, w, mp.inf),
        center=mp.Vector3(0, 0, 0),
        material=mp.Medium(epsilon=n_core ** 2),
    )
]

# --- 使用模式求解器 ---
# TM 模式 (Ex dominant)
mode_tm = mp.solve_freq(
    mp.Simulation(
        cell_size=cell_size,
        geometry=geometry,
        resolution=resolution,
        default_material=mp.Medium(epsilon=n_clad ** 2),
    ),
    freq=freq,
    k=mp.Vector3(1, 0, 0),  # 传播方向
    direction=mp.X,          # 传播方向
    n=5,                     # 找 5 个模式
)

print("=== TM Modes ===")
for i, m in enumerate(mode_tm):
    neff = m.k.x / (2 * np.pi * freq)
    print(f"  Mode {i}: neff = {neff:.4f}, k = {m.k.x:.4f}")

# TE 模式 (Hz dominant)
mode_te = mp.solve_freq(
    mp.Simulation(
        cell_size=cell_size,
        geometry=geometry,
        resolution=resolution,
        default_material=mp.Medium(epsilon=n_clad ** 2),
    ),
    freq=freq,
    k=mp.Vector3(1, 0, 0),
    direction=mp.X,
    n=5,
    hz=True,
)

print("\n=== TE Modes ===")
for i, m in enumerate(mode_te):
    neff = m.k.x / (2 * np.pi * freq)
    print(f"  Mode {i}: neff = {neff:.4f}, k = {m.k.x:.4f}")

# --- 可视化其中一个模式的场分布 ---
mode_to_view = mode_tm[0]  # 基模

sim = mp.Simulation(
    cell_size=cell_size,
    geometry=geometry,
    resolution=resolution,
    default_material=mp.Medium(epsilon=n_clad ** 2),
)

sim.init_sim()

# 提取基模的 Ey 场分布
ey_field = sim.get_eigenmode(
    freq=freq,
    k=mode_to_view.k,
    direction=mp.X,
    where=mp.Volume(center=mp.Vector3(0, 0, 0), size=cell_size),
    component=mp.Ey,
    match_frequency=True,
)

# 获取场数据并重塑为 2D
ey_2d = ey_field.get_array().reshape(6 * resolution, 6 * resolution)

plt.figure(figsize=(6, 5))
plt.imshow(ey_2d, cmap='RdBu', aspect='auto',
           extent=[-3, 3, -3, 3], origin='lower')
plt.colorbar(label='|Ey|')
plt.axhline(y=w/2, color='green', linestyle='--', linewidth=1)
plt.axhline(y=-w/2, color='green', linestyle='--', linewidth=1,
            label='Waveguide boundary')
plt.xlabel('x (μm)')
plt.ylabel('y (μm)')
plt.title(f'Waveguide TM₀ Mode @ λ={wavelength}μm')
plt.legend()
plt.tight_layout()
plt.savefig('waveguide_mode.png', dpi=150)
plt.show()
