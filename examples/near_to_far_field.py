"""
Examples: 近场-远场变换 (Near-to-Far Field Transformation)

MEEP 提供的 near2far 功能：
从近场仿真结果计算远场辐射图样。
常用于天线设计和散射截面计算。
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# --- 参数 ---
resolution = 40
wavelength = 1.0
freq       = 1.0 / wavelength
dpml       = 1.0
pad        = 1.0

# 远场计算参数
far_field_radius = 1000  # 远场观测距离 (>> 波长)
n_angles         = 360   # 角度采样数

# 偶极子天线参数
dipole_length = 0.5  # 半波偶极子

cell_size = mp.Vector3(6, 6, 0)

# --- 偶极子天线 ---
geometry = [
    mp.Block(
        size=mp.Vector3(dipole_length, 0.05, mp.inf),
        center=mp.Vector3(0, 0, 0),
        material=mp.Medium(epsilon=1.0),  # 理想导体
    ),
]

# 使用电流源激励天线中心
sources = [mp.Source(
    mp.GaussianSource(frequency=freq, fwidth=0.1 * freq),
    component=mp.Ez,
    center=mp.Vector3(0, 0.025, 0),
    size=mp.Vector3(dipole_length * 0.5, 0, 0),
)]

# --- 仿真 ---
sim = mp.Simulation(
    cell_size=cell_size,
    geometry=geometry,
    sources=sources,
    boundary_layers=[mp.PML(dpml)],
    resolution=resolution,
)

# 近场-远场变换表面（包围天线的盒子）
near_field_box = sim.add_near2far(
    freq, 0.1 * freq, 1,
    # 四条边组成闭合表面
    mp.Near2FarRegion(
        center=mp.Vector3(0, 0.5, 0),
        size=mp.Vector3(2.0, 0, 0),
    ),
    mp.Near2FarRegion(
        center=mp.Vector3(0, -0.5, 0),
        size=mp.Vector3(2.0, 0, 0),
    ),
    mp.Near2FarRegion(
        center=mp.Vector3(1.0, 0, 0),
        size=mp.Vector3(0, 1.0, 0),
    ),
    mp.Near2FarRegion(
        center=mp.Vector3(-1.0, 0, 0),
        size=mp.Vector3(0, 1.0, 0),
    ),
)

sim.run(until=200)

# --- 计算远场 ---
angles = np.linspace(0, 2 * np.pi, n_angles)

# 远场点列表（在远场半径的圆上）
far_field_points = [
    mp.Vector3(
        far_field_radius * np.cos(theta),
        far_field_radius * np.sin(theta),
        0
    )
    for theta in angles
]

# 计算远场电场
far_field = sim.get_farfields(near_field_box, far_field_points)

# 提取 Ez 分量
Ez_far = [abs(ff[2]) for ff in far_field]  # Ez 是 z 分量
Ez_far = np.array(Ez_far)

# 归一化
Ez_db = 10 * np.log10(Ez_far / Ez_far.max() + 1e-10)

# --- 绘制辐射方向图 ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6),
                           subplot_kw={'projection': 'polar'})

# 线性标度
axes[0].plot(angles, Ez_far, linewidth=2, color='C0')
axes[0].fill(angles, Ez_far, alpha=0.2, color='C0')
axes[0].set_title('Radiation Pattern (Linear)', va='bottom')
axes[0].set_rlabel_position(45)

# dB 标度
axes[1].plot(angles, Ez_db, linewidth=2, color='C1')
axes[1].fill(angles, Ez_db, alpha=0.2, color='C1')
axes[1].set_ylim(-40, 5)
axes[1].set_title('Radiation Pattern (dB)', va='bottom')
axes[1].set_rlabel_position(45)

plt.suptitle(f'Half-Wave Dipole Antenna\nλ = {wavelength} μm',
             fontsize=14)
plt.tight_layout()
plt.savefig('antenna_radiation_pattern.png', dpi=150)
plt.show()

# --- 方向性计算 ---
# 数值积分求总辐射功率
dtheta = 2 * np.pi / n_angles
P_total = np.sum(Ez_far ** 2) * dtheta
P_max = Ez_far.max() ** 2
directivity = 2 * np.pi * P_max / P_total
D_dB = 10 * np.log10(directivity)

# 主瓣方向
main_lobe_idx = np.argmax(Ez_far)
main_lobe_angle = np.degrees(angles[main_lobe_idx])

# 波束宽度
half_power = Ez_far.max() / np.sqrt(2)
# 找 -3dB 点
above_hp = Ez_far >= half_power
transitions = np.diff(above_hp.astype(int))
# HPBW 近似
hpbw_deg = np.sum(above_hp) * np.degrees(dtheta)

print("=== 天线辐射特性 ===")
print(f"方向性: {directivity:.2f} ({D_dB:.1f} dBi)")
print(f"主瓣方向: {main_lobe_angle:.1f}°")
print(f"HPBW 近似: {hpbw_deg:.1f}°")
print("\n参考: 理想半波偶极子方向性约为 2.15 dBi")
