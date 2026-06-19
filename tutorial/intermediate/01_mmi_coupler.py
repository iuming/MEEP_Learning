"""
Intermediate 01: 多模干涉耦合器 (MMI)

设计 1×2 多模干涉耦合器，将输入光功率等分到两个输出端口。
这是集成光子学中最基本的功率分配器件。
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# --- 参数 ---
resolution = 40
n_wg       = 2.8     # 芯层折射率
n_clad     = 1.0     # 包层折射率
wavelength = 1.55
freq       = 1.0 / wavelength

# MMI 参数
wg_width    = 0.5    # 输入/输出波导宽度
mmi_width   = 3.0    # MMI 宽度
mmi_length  = 11.0   # MMI 长度 (估算 3Lπ/2)
taper_length = 2.0   # 锥形过渡区长度
sep         = 1.0    # 输出波导间距
dpml        = 1.0
pad         = 2.0

cell_x = mmi_length + 2 * taper_length + 2 * pad + 4
cell_y = mmi_width + 2 * pad + 2
cell = mp.Vector3(cell_x, cell_y, 0)

# --- 几何体 ---
geometry = []

# MMI 主体
geometry.append(mp.Block(
    size=mp.Vector3(mmi_length, mmi_width, mp.inf),
    center=mp.Vector3(0, 0, 0),
    material=mp.Medium(epsilon=n_wg ** 2),
))

# 输入波导 (左侧)
geometry.append(mp.Block(
    size=mp.Vector3(pad + taper_length, wg_width, mp.inf),
    center=mp.Vector3(-mmi_length / 2 - taper_length / 2 - pad / 2, 0, 0),
    material=mp.Medium(epsilon=n_wg ** 2),
))

# 输入锥形
taper_points = [
    mp.Vector3(-mmi_length / 2, wg_width / 2, 0),
    mp.Vector3(-mmi_length / 2 - taper_length, wg_width / 2, 0),
    mp.Vector3(-mmi_length / 2 - taper_length, -wg_width / 2, 0),
    mp.Vector3(-mmi_length / 2, -mmi_width / 2, 0),
]
geometry.append(mp.Prism(taper_points, height=mp.inf,
                          material=mp.Medium(epsilon=n_wg ** 2)))

# 输出波导 (右侧，两个)
for s in [-sep / 2, sep / 2]:
    geometry.append(mp.Block(
        size=mp.Vector3(pad + taper_length, wg_width, mp.inf),
        center=mp.Vector3(mmi_length / 2 + taper_length / 2 + pad / 2, s, 0),
        material=mp.Medium(epsilon=n_wg ** 2),
    ))

# 输出锥形
for sy, sign in [(-sep / 2, -1), (sep / 2, 1)]:
    taper_out = [
        mp.Vector3(mmi_length / 2, sy + sign * mmi_width / 2, 0),
        mp.Vector3(mmi_length / 2 + taper_length, sy + sign * wg_width / 2, 0),
        mp.Vector3(mmi_length / 2 + taper_length, sy - sign * wg_width / 2, 0),
        mp.Vector3(mmi_length / 2, sy - sign * mmi_width / 2, 0),
    ]
    geometry.append(mp.Prism(taper_out, height=mp.inf,
                              material=mp.Medium(epsilon=n_wg ** 2)))

# --- 光源 ---
sources = [mp.Source(
    mp.GaussianSource(frequency=freq, fwidth=0.1 * freq),
    component=mp.Ez,
    center=mp.Vector3(-mmi_length / 2 - taper_length - pad + 0.5, 0, 0),
    size=mp.Vector3(0, wg_width * 2, 0),
)]

# --- 仿真 ---
sim = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    sources=sources,
    boundary_layers=[mp.PML(dpml)],
    resolution=resolution,
    default_material=mp.Medium(epsilon=n_clad ** 2),
)

# 输出端口通量监视
out1_fr = mp.FluxRegion(
    center=mp.Vector3(mmi_length / 2 + taper_length + pad - 0.5, -sep / 2, 0),
    size=mp.Vector3(0, wg_width * 2, 0),
)
out2_fr = mp.FluxRegion(
    center=mp.Vector3(mmi_length / 2 + taper_length + pad - 0.5, sep / 2, 0),
    size=mp.Vector3(0, wg_width * 2, 0),
)

out1 = sim.add_flux(freq, 0.1 * freq, 1, out1_fr)
out2 = sim.add_flux(freq, 0.1 * freq, 1, out2_fr)

sim.run(until=300)

# --- 结果 ---
P1 = np.array(mp.get_fluxes(out1))[0]
P2 = np.array(mp.get_fluxes(out2))[0]
P_total = P1 + P2

imb = -10 * np.log10(max(P1, P2) / min(P1, P2))

print(f"输出端口 1: {P1:.4f}  ({P1 / P_total * 100:.1f}%)")
print(f"输出端口 2: {P2:.4f}  ({P2 / P_total * 100:.1f}%)")
print(f"不平衡度:    {imb:.2f} dB")
print(f"目标: < 0.5 dB")

# --- 可视化 ---
eps = sim.get_array(component=mp.Dielectric)

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# 结构图
axes[0].imshow(eps.transpose(), interpolation='spline36', cmap='binary',
               extent=[-cell_x/2, cell_x/2, -cell_y/2, cell_y/2])
axes[0].set_title(f'1×2 MMI Coupler\nImbalance: {imb:.2f} dB')
axes[0].set_xlabel('x (μm)')
axes[0].set_ylabel('y (μm)')

# 输出比例饼图
axes[1].pie([P1, P2], labels=[f'Port 1\n{P1/P_total*100:.1f}%',
                               f'Port 2\n{P2/P_total*100:.1f}%'],
            colors=['#3498db', '#e74c3c'], autopct='%1.1f%%',
            explode=(0, 0.05))
axes[1].set_title('Power Split Ratio')

plt.tight_layout()
plt.savefig('mmi_coupler.png', dpi=150)
plt.show()
