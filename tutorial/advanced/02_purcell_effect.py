"""
Advanced 02: 量子发射体 — 自发辐射增强 (Purcell 效应)

将一个量子发射体（偶极子）放置在光学微腔中，
计算 Purcell 因子——自发辐射的增强倍数。

应用：单光子源、量子光学、腔量子电动力学 (Cavity QED)
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# --- 参数 ---
resolution = 60
wavelength = 1.55    # 发射波长 (μm)
freq       = 1.0 / wavelength
fwidth     = 0.3 * freq

# 微腔参数
n_cavity  = 3.4     # 腔折射率 (如 Si)
n_clad    = 1.0     # 空气
cavity_Lx = 0.5     # 腔长度
Q_factor  = 1000    # 目标品质因子

dpml = 1.0
pad  = 2.0

cell = mp.Vector3(8, 4, 4)  # 3D 仿真！

# --- 构建腔 ---
cavity = mp.Block(
    size=mp.Vector3(cavity_Lx, 0.5, 0.5),
    center=mp.Vector3(0, 0, 0),
    material=mp.Medium(epsilon=n_cavity ** 2),
)

# --- 量子发射体 (电偶极子) ---
# 放置在腔中心
dipole_pos = mp.Vector3(0, 0, 0)
dipole_moment = mp.Vector3(1, 0, 0)  # x 极化

# 偶极子源 (对于 Purcell 计算)
sources = [mp.Source(
    mp.GaussianSource(frequency=freq, fwidth=fwidth),
    component=mp.Ex,
    center=dipole_pos,
)]

# --- 仿真 ---
sim = mp.Simulation(
    cell_size=cell,
    geometry=[cavity],
    sources=sources,
    boundary_layers=[mp.PML(dpml)],
    resolution=resolution,
    default_material=mp.Medium(epsilon=n_clad ** 2),
)

# 记录发射体位置的近场
field_data = []

def record_field(sim):
    ex = sim.get_field_point(mp.Ex, dipole_pos)
    field_data.append(ex)

# 放置通量监视器包围腔，计算辐射功率
flux_box = sim.add_flux(
    freq, fwidth, 1,
    mp.FluxRegion(center=mp.Vector3(0.5, 0, 0), size=mp.Vector3(0, 1, 1)),
    mp.FluxRegion(center=mp.Vector3(-0.5, 0, 0), size=mp.Vector3(0, 1, 1)),
)

sim.run(
    mp.at_every(0.1, record_field),
    until=150,
)

# 计算辐射功率
P_cavity = abs(np.array(mp.get_fluxes(flux_box))[0])

# --- 参考仿真：自由空间中的偶极子 ---
sim_free = mp.Simulation(
    cell_size=cell,
    geometry=[],
    sources=sources,
    boundary_layers=[mp.PML(dpml)],
    resolution=resolution,
)

flux_free = sim_free.add_flux(
    freq, fwidth, 1,
    mp.FluxRegion(center=mp.Vector3(0.5, 0, 0), size=mp.Vector3(0, 1, 1)),
    mp.FluxRegion(center=mp.Vector3(-0.5, 0, 0), size=mp.Vector3(0, 1, 1)),
)

sim_free.run(until=150)
P_free = abs(np.array(mp.get_fluxes(flux_free))[0])

# --- Purcell 因子 ---
F_p = P_cavity / P_free

# --- Q 因子估算 (harminv) ---
harm = mp.Harminv(mp.Ex, dipole_pos, freq, fwidth)
sim2 = mp.Simulation(
    cell_size=cell,
    geometry=[cavity],
    sources=sources,
    boundary_layers=[mp.PML(dpml)],
    resolution=resolution,
    default_material=mp.Medium(epsilon=n_clad ** 2),
)
sim2.run(mp.after_sources(harm), until_after_sources=300)

modes = harm.modes
if modes:
    mode = modes[0]
    Q_estimated = mode.Q
    resonance_freq = mode.freq
    resonance_wl = 1.0 / resonance_freq
else:
    Q_estimated = np.nan
    resonance_wl = np.nan

# --- 理论对比 ---
# 对于理想微腔: F_p = (3/4π²) * (λ/n)³ * (Q/V)
V_mode = cavity_Lx * 0.5 * 0.5  # 模式体积近似
F_p_theory = (3 / (4 * np.pi ** 2)) * \
             (wavelength / n_cavity) ** 3 * Q_estimated / V_mode

# --- 输出结果 ---
print("=== Purcell Enhancement ===")
print(f"波长: {wavelength:.3f} μm")
print(f"腔材料: n = {n_cavity}")
print(f"腔尺寸: {cavity_Lx}×0.5×0.5 μm³")
print(f"模式体积 V: {V_mode:.4f} μm³")
print(f"品质因子 Q: {Q_estimated:.1f}" if not np.isnan(Q_estimated) else "N/A")
print(f"谐振波长: {resonance_wl:.4f} μm" if not np.isnan(resonance_wl) else "N/A")
print()
print(f"自由空间辐射功率 P_free: {P_free:.6f}")
print(f"腔中辐射功率 P_cavity:   {P_cavity:.6f}")
print(f"Purcell Factor (仿真):     {F_p:.2f}")
print(f"Purcell Factor (理论):     {F_p_theory:.2f}")

# --- 可视化 ---
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# 时域信号
field_arr = np.array(field_data)
axes[0].plot(np.linspace(0, 15, len(field_arr)), field_arr, linewidth=1)
axes[0].set_xlabel('Time')
axes[0].set_ylabel('Ex @ dipole')
axes[0].set_title('Dipole Near-Field Decay')
axes[0].grid(True, alpha=0.3)

# Purcell 因子
axes[1].bar(['Free Space', 'Cavity (sim)', f'Cavity (theory)\nQ≈{Q_estimated:.0f}'],
            [1, F_p, F_p_theory],
            color=['gray', '#3498db', '#e74c3c'])
axes[1].axhline(y=1, color='gray', linestyle=':', alpha=0.5)
axes[1].set_ylabel('Purcell Factor Fp')
axes[1].set_title('Emission Enhancement')

# 腔模式剖面 (示意)
x = np.linspace(-cavity_Lx, cavity_Lx, 200)
mode_profile = np.cos(np.pi * x / cavity_Lx) * np.exp(-abs(x) / 2)
axes[2].plot(x, mode_profile ** 2, linewidth=1.5)
axes[2].axvspan(-cavity_Lx / 2, cavity_Lx / 2, alpha=0.1, color='red',
                label='Cavity region')
axes[2].set_xlabel('x (μm)')
axes[2].set_ylabel('|E|²')
axes[2].set_title('Fundamental Mode Profile')
axes[2].legend()
axes[2].grid(True, alpha=0.3)

plt.suptitle(f'Cavity QED: Purcell Effect (Fp = {F_p:.2f})',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('purcell_effect.png', dpi=150)
plt.show()
