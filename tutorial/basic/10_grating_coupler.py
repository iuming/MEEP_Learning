"""
Example 10: 亚波长光栅耦合器

设计一个亚波长光栅，将自由空间光耦合进入波导。
演示如何优化光栅周期和占空比以达到最大耦合效率。
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# --- 参数 ---
resolution = 40
n_wg       = 2.8     # 波导折射率 (如 GaAs)
n_clad     = 1.0     # 空气
wg_height  = 0.22    # 波导厚度
wavelength = 1.55
freq       = 1.0 / wavelength
fwidth     = 0.2 * freq

# 光栅参数
n_teeth    = 15      # 光栅齿数
period     = 0.58    # 光栅周期 (初始猜测)
duty_cycle = 0.5     # 占空比 (初始猜测)
etch_depth = 0.1     # 刻蚀深度
dpml       = 1.0
pad        = 3.0

# --- 几何体 ---
cell_height = 2 * dpml + 2 * pad + 3
cell_width  = 2 * dpml + 2 * pad + n_teeth * period + 4

cell = mp.Vector3(cell_width, cell_height, 0)

geometry = []

# 底层波导
geometry.append(mp.Block(
    size=mp.Vector3(mp.inf, wg_height, mp.inf),
    center=mp.Vector3(0, -wg_height / 2, 0),
    material=mp.Medium(epsilon=n_wg ** 2),
))

# 光栅齿 (在波导上面)
for i in range(n_teeth):
    x_pos = -n_teeth * period / 2 + i * period + period / 2
    if i % 2 == 0:  # 占空比 50%
        geometry.append(mp.Block(
            size=mp.Vector3(period * duty_cycle, etch_depth, mp.inf),
            center=mp.Vector3(x_pos, wg_height / 2 + etch_depth / 2, 0),
            material=mp.Medium(epsilon=n_wg ** 2),
        ))

# --- 光源 (高斯光束从上方入射) ---
sources = [mp.Source(
    mp.GaussianSource(frequency=freq, fwidth=fwidth),
    component=mp.Ez,
    center=mp.Vector3(0, dpml + pad, 0),
    size=mp.Vector3(period * n_teeth * 0.6, 0, 0),
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

# 透射监视器 (在波导末端)
tran_fr = mp.FluxRegion(
    center=mp.Vector3(4, 0, 0),
    size=mp.Vector3(0, 2, 0),
)
tran = sim.add_flux(
    np.linspace(freq - fwidth, freq + fwidth, 100),
    1, tran_fr,
)

# 反射监视器 (在上方)
refl_fr = mp.FluxRegion(
    center=mp.Vector3(0, dpml + pad - 0.5, 0),
    size=mp.Vector3(period * n_teeth * 0.6, 0, 0),
)
refl = sim.add_flux(
    np.linspace(freq - fwidth, freq + fwidth, 100),
    1, refl_fr,
)

# 空仿真 (无光栅)
sim_empty = mp.Simulation(
    cell_size=cell,
    geometry=[geometry[0]],  # 只有波导层
    sources=sources,
    boundary_layers=[mp.PML(dpml)],
    resolution=resolution,
    default_material=mp.Medium(epsilon=n_clad ** 2),
)
tran0 = sim_empty.add_flux(
    np.linspace(freq - fwidth, freq + fwidth, 100),
    1, tran_fr,
)

sim.run(until=400)
sim_empty.run(until=400)

# --- 结果 ---
freqs = np.linspace(freq - fwidth, freq + fwidth, 100)
tran_flux = np.array(mp.get_fluxes(tran))
tran0_flux = np.array(mp.get_fluxes(tran0))
refl_flux = np.array(mp.get_fluxes(refl))

T = tran_flux / tran0_flux
R = refl_flux / tran0_flux
coupling = T - R  # 近似耦合效率

wavelengths = 1.0 / freqs

idx_center = np.argmin(np.abs(freqs - freq))
print(f"中心波长 {wavelength:.2f} μm:")
print(f"  耦合效率: {coupling[idx_center]:.2%}")
print(f"  透射:     {T[idx_center]:.2%}")
print(f"  反射:     {R[idx_center]:.2%}")

# --- 绘图 ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# 光谱
axes[0].plot(wavelengths, T, label='Transmission', linewidth=1.5)
axes[0].plot(wavelengths, R, label='Reflection', linewidth=1.5)
axes[0].plot(wavelengths, coupling, label='Net Coupling', linewidth=2, color='black')
axes[0].axvline(x=wavelength, color='red', linestyle='--', alpha=0.5)
axes[0].set_xlabel('Wavelength (μm)')
axes[0].set_ylabel('Efficiency')
axes[0].set_title('Grating Coupler Spectrum')
axes[0].legend()
axes[0].grid(True, alpha=0.3)
axes[0].invert_xaxis()

# 场分布
eps = sim.get_array(component=mp.Dielectric)
axes[1].imshow(eps.transpose(), interpolation='spline36', cmap='binary',
               extent=[-cell_width/2, cell_width/2, -cell_height/2, cell_height/2])
axes[1].set_title('Dielectric Structure')
axes[1].set_xlabel('x (μm)')
axes[1].set_ylabel('y (μm)')

plt.tight_layout()
plt.savefig('grating_coupler.png', dpi=150)
plt.show()
