"""
Intermediate 02: 光子晶体波导与慢光效应

在二维光子晶体中引入线缺陷形成波导。
分析导模的色散关系，观察慢光（低群速度）区域。
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# --- 参数 ---
resolution = 32
a          = 1.0     # 晶格常数
r          = 0.3 * a # 孔半径
n_slab     = 3.4     # 底板折射率 (如硅)
n_hole     = 1.0     # 孔内折射率
wavelength = 1.55 / a

# 光子晶体尺寸
nx = 31  # x 方向晶胞数
ny = 15  # y 方向晶胞数

cell = mp.Vector3(nx * a, ny * a, 0)

# --- 构建光子晶体 ---
geometry = []

for i in range(nx):
    for j in range(ny):
        # 跳过后 5 行形成线缺陷波导
        if j == ny // 2 or j == ny // 2 + 1 or j == ny // 2 - 1:
            continue
        x = (i - nx / 2 + 0.5) * a
        y = (j - ny / 2 + 0.5) * a
        geometry.append(mp.Cylinder(
            radius=r,
            center=mp.Vector3(x, y, 0),
            height=mp.inf,
            material=mp.Medium(epsilon=n_hole ** 2),
        ))

# --- 光源 ---
freq = 0.3 / a
fwidth = freq * 0.5

sources = [mp.Source(
    mp.GaussianSource(frequency=freq, fwidth=fwidth),
    component=mp.Ez,
    center=mp.Vector3(-nx * a / 2 + 2 * a, 0, 0),
    size=mp.Vector3(0, a * 2, 0),
)]

# --- 仿真 ---
sim = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    sources=sources,
    boundary_layers=[mp.PML(1.0)],
    resolution=resolution,
    default_material=mp.Medium(epsilon=n_slab ** 2),
)

# 在波导中放置多个监视器观察传播
monitor_ys = [0]
monitor_xs = np.linspace(-nx * a / 2 + 3 * a, nx * a / 2 - 3 * a, nx - 6)

# 记录每个时间步各个监视点的场值
field_data = {x: [] for x in monitor_xs}

def record_fields(sim):
    for x in monitor_xs:
        field_data[x].append(sim.get_field_point(mp.Ez, mp.Vector3(x, 0, 0)))

sim.run(
    mp.at_every(0.5, record_fields),
    until=200,
)

# --- 分析传播 ---
# 计算群速度：跟踪波包峰值传播
peak_times = []
for x in monitor_xs:
    data = np.array(field_data[x])
    peak_idx = np.argmax(np.abs(data))
    peak_times.append(peak_idx * 0.5)

peak_times = np.array(peak_times)
# 用线性拟合求群速度
coeffs = np.polyfit(monitor_xs, peak_times, 1)
vg = 1.0 / coeffs[0]  # dx/dt

print(f"群速度 vg: {vg:.4f} c")
print(f"禁带外群速度通常接近 c/n_eff")
print(f"慢光条件: vg << c/n")

# --- 能带计算 (沿波导方向) ---
# Bloch 周期边界 + k 扫描简化版
k_points = np.linspace(0, np.pi / a, 15)
band_freqs = []

for kx in k_points:
    sim_k = mp.Simulation(
        cell_size=mp.Vector3(a, ny * a, 0),
        geometry=[mp.Cylinder(r, material=mp.Medium(epsilon=n_hole ** 2))
                  for _ in range(1)],  # 简化
        resolution=resolution,
        default_material=mp.Medium(epsilon=n_slab ** 2),
        k_point=mp.Vector3(kx, 0, 0),
        ensure_periodicity=True,
    )
    # 此处仅示意，完整能带计算见 basic/06
    sim_k.reset_meep()

# --- 可视化 ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# 结构
eps = sim.get_array(component=mp.Dielectric)
axes[0].imshow(np.log(eps.transpose() + 1),
               interpolation='spline36', cmap='hot',
               extent=[-cell.x/2, cell.x/2, -cell.y/2, cell.y/2])
axes[0].axhline(y=0, color='cyan', linestyle='--', linewidth=1,
                label='Waveguide (line defect)')
axes[0].set_title('Photonic Crystal Waveguide (W1)')
axes[0].set_xlabel('x (a)')
axes[0].set_ylabel('y (a)')
axes[0].legend()

# 时间-空间传播图
time_grid = np.array([field_data[x] for x in monitor_xs])
axes[1].imshow(time_grid, aspect='auto', cmap='RdBu',
               extent=[0, 200 * 0.5, monitor_xs[0], monitor_xs[-1]])
axes[1].set_title(f'Pulse Propagation (vg ≈ {vg:.3f}c)')
axes[1].set_xlabel('Time')
axes[1].set_ylabel('x position (a)')

plt.tight_layout()
plt.savefig('photonic_crystal_wg.png', dpi=150)
plt.show()
