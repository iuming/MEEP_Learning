"""
Example 08: 环形谐振腔

仿真一个与直波导耦合的环形谐振腔。
通过扫描波长观察谐振特性（透射谱中的 dips）。
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# --- 参数 ---
resolution = 30
n_wg       = 2.0     # 波导折射率
n_ring     = 2.5     # 环折射率（略高确保模式约束）
radius     = 2.5     # 环半径
gap        = 0.3     # 耦合间隙
wg_width   = 0.5     # 波导宽度
wavelength = 1.55    # 中心波长

freq_center = 1.0 / wavelength
freq_width  = 0.15 * freq_center  # 扫描范围
dpml        = 1.0
pad         = 2.0

cell_size = mp.Vector3(
    2 * (dpml + pad + radius + 3),
    2 * (radius + 3),
    0
)

# --- 几何体 ---
# 环形谐振腔
ring = mp.Cylinder(
    radius=radius,
    material=mp.Medium(epsilon=n_ring ** 2),
    center=mp.Vector3(0, 0, 0),
    height=mp.inf,
)

# 直波导
bus_wg = mp.Block(
    size=mp.Vector3(mp.inf, wg_width, mp.inf),
    center=mp.Vector3(0, radius + gap + wg_width / 2, 0),
    material=mp.Medium(epsilon=n_wg ** 2),
)

geometry = [ring, bus_wg]

# --- 光源 ---
sources = [mp.Source(
    mp.GaussianSource(frequency=freq_center, fwidth=freq_width),
    component=mp.Ez,
    center=mp.Vector3(-radius - 2, radius + gap + wg_width / 2, 0),
    size=mp.Vector3(0, wg_width * 2, 0),
)]

# --- 仿真 ---
sim = mp.Simulation(
    cell_size=cell_size,
    geometry=geometry,
    sources=sources,
    boundary_layers=[mp.PML(dpml)],
    resolution=resolution,
)

# 透射通量监视器
nfreq = 200
tran_fr = mp.FluxRegion(
    center=mp.Vector3(radius + 2, radius + gap + wg_width / 2, 0),
    size=mp.Vector3(0, wg_width * 3, 0),
)
tran = sim.add_flux(
    np.linspace(freq_center - freq_width, freq_center + freq_width, nfreq),
    1,
    tran_fr,
)

# 场监视器（用于可视化）
sim.run(
    until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(
        radius + 2, radius + gap + wg_width / 2, 0
    ), 1e-6),
)

# --- 计算透射谱 ---
freqs = np.linspace(freq_center - freq_width, freq_center + freq_width, nfreq)
tran_flux = np.array(mp.get_fluxes(tran))

# 归一化：跑空仿真
sim_empty = mp.Simulation(
    cell_size=cell_size,
    geometry=[bus_wg],  # 只有直波导，没有环
    sources=sources,
    boundary_layers=[mp.PML(dpml)],
    resolution=resolution,
)
tran0 = sim_empty.add_flux(
    np.linspace(freq_center - freq_width, freq_center + freq_width, nfreq),
    1, tran_fr,
)
sim_empty.run(
    until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(
        radius + 2, radius + gap + wg_width / 2, 0
    ), 1e-6),
)
tran0_flux = np.array(mp.get_fluxes(tran0))

T = tran_flux / tran0_flux
wavelengths = 1.0 / freqs

# --- 绘图 ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# 透射谱
ax1.plot(wavelengths, T, linewidth=1.5)
ax1.set_xlabel('Wavelength (μm)')
ax1.set_ylabel('Transmission')
ax1.set_title(f'Ring Resonator Transmission (R={radius}μm, gap={gap}μm)')
ax1.grid(True, alpha=0.3)
ax1.invert_xaxis()

# 找到谐振波长（dips）
from scipy.signal import argrelextrema
min_indices = argrelextrema(T, np.less)[0]
resonant_wavelengths = wavelengths[min_indices]
for wl in resonant_wavelengths:
    ax1.axvline(x=wl, color='red', linestyle='--', alpha=0.4)

# 场分布
eps_data = sim.get_array(center=mp.Vector3(0, 0, 0),
                          size=mp.Vector3(cell_size.x, cell_size.y, 0),
                          component=mp.Dielectric)
ez_data = sim.get_array(center=mp.Vector3(0, 0, 0),
                         size=mp.Vector3(cell_size.x, cell_size.y, 0),
                         component=mp.Ez)

ax2.imshow(eps_data.transpose(),
           interpolation='spline36', cmap='binary',
           extent=[-cell_size.x/2, cell_size.x/2, -cell_size.y/2, cell_size.y/2])
ax2.set_title('Dielectric Structure')
ax2.set_xlabel('x (μm)')
ax2.set_ylabel('y (μm)')

plt.tight_layout()
plt.savefig('ring_resonator.png', dpi=150)
plt.show()

print(f"发现 {len(resonant_wavelengths)} 个谐振峰:")
for wl in resonant_wavelengths:
    print(f"  λ = {wl:.4f} μm")
