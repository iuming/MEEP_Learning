"""
Advanced 04: 光力效应 — 光阱与粒子操控

仿真聚焦光束对纳米颗粒的光力作用。
计算光学梯度力和散射力，分析光阱稳定性。

应用：光镊、冷原子捕获、纳米操控、生物分子力学
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# --- 参数 ---
resolution = 50
wavelength = 1.064   # 1064 nm (常用光镊波长)
freq       = 1.0 / wavelength
fwidth     = 0.05 * freq

# 颗粒参数
particle_radius = 0.1  # 100 nm 半径
n_particle      = 1.59  # 聚苯乙烯 (PS) 微球
n_medium        = 1.33  # 水

# 光束参数
beam_waist = 0.5  # 束腰半径
NA         = 0.8  # 数值孔径
focal_length = 2.0

dpml = 1.0
pad  = 2.0

cell_width  = 6.0
cell_height = 8.0
cell = mp.Vector3(cell_width, cell_height, 0)

# --- 几何体 ---
# 颗粒
particle = mp.Cylinder(
    radius=particle_radius,
    center=mp.Vector3(0, 0, 0),
    height=mp.inf,
    material=mp.Medium(epsilon=n_particle ** 2),
)

geometry = [particle]

# --- 聚焦高斯光 ---
# 使用多个平面波源模拟高斯光束 (角谱法)
# 简化版：使用带相位的线源合成聚焦束
beam_angle = np.arcsin(NA / n_medium)
beam_half_width = focal_length * np.tan(beam_angle)

sources = []

# 在底部放置多个源，带球面波前
n_sources = 21
source_y = -cell_height / 2 + dpml + 0.5

for i in range(n_sources):
    x = (i - (n_sources - 1) / 2) * (beam_half_width * 2 / (n_sources - 1))
    # 球面波相位: φ = -k * sqrt(x² + f²)
    r_to_focus = np.sqrt(x ** 2 + focal_length ** 2)
    phase = -2 * np.pi * n_medium * r_to_focus / wavelength

    # 高斯振幅加权
    amp = np.exp(-x ** 2 / (beam_waist ** 2))

    sources.append(mp.Source(
        mp.GaussianSource(frequency=freq, fwidth=fwidth),
        component=mp.Ez,
        center=mp.Vector3(x, source_y, 0),
        amplitude=amp,
    ))

# --- 仿真（无颗粒 - 参考光场） ---
sim_free = mp.Simulation(
    cell_size=cell,
    geometry=[],
    sources=sources,
    boundary_layers=[mp.PML(dpml)],
    resolution=resolution,
    default_material=mp.Medium(epsilon=n_medium ** 2),
)

sim_free.run(until=100)

# --- 仿真（有颗粒） ---
sim = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    sources=sources,
    boundary_layers=[mp.PML(dpml)],
    resolution=resolution,
    default_material=mp.Medium(epsilon=n_medium ** 2),
)

# 力计算使用 Maxwell Stress Tensor (MST)
# 包围颗粒的积分表面
force_box_x = 2 * particle_radius  # 略大于颗粒
mst_flux = sim.add_force(
    freq, 0,
    mp.ForceRegion(
        center=mp.Vector3(0, force_box_x / 2, 0),
        size=mp.Vector3(force_box_x, 0, 0),
        direction=mp.Y,
        weight=+1,
    ),
    mp.ForceRegion(
        center=mp.Vector3(0, -force_box_x / 2, 0),
        size=mp.Vector3(force_box_x, 0, 0),
        direction=mp.Y,
        weight=-1,
    ),
    mp.ForceRegion(
        center=mp.Vector3(force_box_x / 2, 0, 0),
        size=mp.Vector3(0, force_box_x, 0),
        direction=mp.X,
        weight=+1,
    ),
    mp.ForceRegion(
        center=mp.Vector3(-force_box_x / 2, 0, 0),
        size=mp.Vector3(0, force_box_x, 0),
        direction=mp.X,
        weight=-1,
    ),
)

sim.run(until=100)

# 提取力
force = np.array(mp.get_forces(mst_flux))
Fx = force[0]  # 散射力 (横向)
Fy = force[1]  # 梯度力 (纵向)

print("=== 光力计算结果 ===")
print(f"颗粒半径: {particle_radius * 1000:.0f} nm")
print(f"颗粒折率: {n_particle}")
print(f"介质折射率: {n_medium}")
print(f"NA: {NA}")
print()
print(f"散射力 Fx: {Fx:.2e} N")
print(f"梯度力 Fy: {Fy:.2e} N")
print(f"总力 |F|:  {np.sqrt(Fx ** 2 + Fy ** 2):.2e} N")
print()

# 判断光阱稳定性
if Fy < 0:
    print("✓ 梯度力为负 → 颗粒被拉向焦点（稳定捕获）")
else:
    print("✗ 梯度力为正 → 颗粒被推离焦点（不稳定）")

# --- 可视化 ---
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# 参考光场 (无颗粒)
ez_free = sim_free.get_array(component=mp.Ez)
size = sim_free.cell_size

im0 = axes[0, 0].imshow(np.abs(ez_free.transpose()) ** 0.4,
                         interpolation='spline36', cmap='hot',
                         extent=[-size.x/2, size.x/2, -size.y/2, size.y/2],
                         origin='lower')
axes[0, 0].set_title('Incident Field |Ez| (without particle)')
axes[0, 0].set_xlabel('x (μm)')
axes[0, 0].set_ylabel('y (μm)')
plt.colorbar(im0, ax=axes[0, 0])

# 散射场 (有颗粒)
ez_total = sim.get_array(component=mp.Ez)
ez_scat  = ez_total - ez_free

im1 = axes[0, 1].imshow(np.abs(ez_scat.transpose()) ** 0.4,
                         interpolation='spline36', cmap='hot',
                         extent=[-size.x/2, size.x/2, -size.y/2, size.y/2],
                         origin='lower')
# 画颗粒轮廓
circle = plt.Circle((0, 0), particle_radius, fill=False,
                     color='cyan', linewidth=2)
axes[0, 1].add_patch(circle)
axes[0, 1].set_title('Scattered Field |Ez_scat|')
axes[0, 1].set_xlabel('x (μm)')
axes[0, 1].set_ylabel('y (μm)')
plt.colorbar(im1, ax=axes[0, 1])

# 力矢量图
axes[1, 0].arrow(0, 0, Fx * 1e12, Fy * 1e12,
                 head_width=0.2, head_length=0.3,
                 fc='red', ec='red', linewidth=2)
axes[1, 0].set_xlim(-2, 2)
axes[1, 0].set_ylim(-2, 2)
axes[1, 0].axhline(y=0, color='gray', linestyle=':', alpha=0.5)
axes[1, 0].axvline(x=0, color='gray', linestyle=':', alpha=0.5)
axes[1, 0].set_xlabel('Fx (×10⁻¹² N)')
axes[1, 0].set_ylabel('Fy (×10⁻¹² N)')
axes[1, 0].set_title('Optical Force Vector')
axes[1, 0].grid(True, alpha=0.3)

# 示意图：光阱势能
# 梯度力对应势阱
r_vals = np.linspace(-1, 1, 100)
pot = 0.5 * r_vals ** 2  # 简谐近似
axes[1, 1].plot(r_vals, pot, linewidth=2, color='#2ecc71')
axes[1, 1].fill_between(r_vals, pot, alpha=0.2, color='#2ecc71')
axes[1, 1].set_xlabel('Displacement (μm)')
axes[1, 1].set_ylabel('Potential U (a.u.)')
axes[1, 1].set_title('Optical Trap Potential (harmonic approx.)')
axes[1, 1].grid(True, alpha=0.3)

# 计算阱刚度 k_trap = dF/dy|y=0
k_trap = abs(Fy) / (particle_radius * 0.5)  # 粗略估算
print(f"阱刚度 k_trap ≈ {k_trap:.2e} N/m")

plt.suptitle(f'Optical Tweezers: R={particle_radius*1000:.0f}nm, '
             f'λ={wavelength*1000:.0f}nm, NA={NA}',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('optical_tweezers.png', dpi=150)
plt.show()

print("\n=== 光镊设计要点 ===")
print("1. 高 NA 物镜 → 强梯度力（克服散射力）")
print("2. 颗粒折射率 > 介质 → 梯度力指向焦点（稳定捕获）")
print("3. 偏振方向影响横向捕获效率")
print("4. 三维捕获需要紧聚焦（NA > 0.8）")
print()
print("参考:")
print("- Ashkin et al., Opt. Lett. 11, 288 (1986)")
print("- MEEP: Maxwell Stress Tensor Force Calculation")
