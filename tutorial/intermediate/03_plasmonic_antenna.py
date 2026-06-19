"""
Intermediate 03: 等离子体（Plasmonic）纳米天线

仿真金属纳米颗粒的局域表面等离子体共振 (LSPR)。
演示如何使用 Drude-Lorentz 色散模型。

应用场景：增强拉曼散射、生物传感、光热治疗
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# --- 参数 ---
resolution = 80
wavelength_min = 0.4    # 400 nm
wavelength_max = 1.0    # 1000 nm (可见-近红外)
freq_max = 1.0 / wavelength_min
freq_min = 1.0 / wavelength_max
freq_center = (freq_min + freq_max) / 2
freq_width  = (freq_max - freq_min) / 2
nfreq = 200

# --- 金的 Drude-Lorentz 模型 ---
# ε(ω) = 1 - ωp²/(ω² + iγω) + Σ Lorentz terms
gold = mp.Medium(
    epsilon=1.0,
    # Drude: plasma freq, damping rate
    E_susceptibilities=[
        mp.DrudeSusceptibility(
            frequency=9.03,      # ωp (eV → meep units)
            gamma=0.053,         # γ (eV)
            sigma=1.0,
        ),
        # Lorentz terms (interband transitions)
        mp.LorentzianSusceptibility(
            frequency=2.5,
            gamma=0.5,
            sigma=1.0,
        ),
    ],
    valid_freq_range=(freq_min, freq_max),
)

# 替代：简化的 Drude 金属模型
metal = mp.Medium(
    epsilon=1.0,
    E_susceptibilities=[
        mp.DrudeSusceptibility(
            frequency=1.0,       # ωp / (2πc) in 1/μm
            gamma=0.01,          # collision rate
            sigma=1.0,
        ),
    ],
)

# --- 几何体 ---
# 金纳米棒 (rod)
rod_length = 0.08  # 80 nm
rod_radius = 0.02  # 20 nm
cell_size = mp.Vector3(0.5, 0.5, 0.5)  # 3D 仿真！

dpml = 0.1

geometry = [
    mp.Block(
        size=mp.Vector3(rod_length, rod_radius * 2, rod_radius * 2),
        center=mp.Vector3(0, 0, 0),
        material=metal,
    ),
]

# --- 光源 (平面波，偏振沿 rod 长轴) ---
sources = [mp.Source(
    mp.GaussianSource(frequency=freq_center, fwidth=freq_width),
    component=mp.Ex,       # x 偏振
    center=mp.Vector3(-0.2, 0, 0),
    size=mp.Vector3(0, 0.3, 0.3),
)]

# --- 仿真 ---
sim = mp.Simulation(
    cell_size=cell_size,
    geometry=geometry,
    sources=sources,
    boundary_layers=[mp.PML(dpml)],
    resolution=resolution,
)

# 散射截面：通过通量监视器
# 总场/散射场 (TFSF) 配置
flux_box_x = mp.FluxRegion(
    center=mp.Vector3(0.15, 0, 0),
    size=mp.Vector3(0, cell_size.y, cell_size.z),
)
flux_box_y = mp.FluxRegion(
    center=mp.Vector3(0, 0.15, 0),
    size=mp.Vector3(cell_size.x, 0, cell_size.z),
)

# 吸收和散射通量
absorb = sim.add_flux(
    np.linspace(freq_min, freq_max, nfreq),
    1,
    flux_box_x,
)

scat = sim.add_flux(
    np.linspace(freq_min, freq_max, nfreq),
    1,
    flux_box_y,
)

# 运行
sim.run(until=500)

# --- 结果 ---
freqs = np.linspace(freq_min, freq_max, nfreq)
absorb_flux = np.array(mp.get_fluxes(absorb))
scat_flux = np.array(mp.get_fluxes(scat))
extinction = absorb_flux + scat_flux

wavelengths = 1.0 / freqs

# 找到共振峰
peak_idx = np.argmax(extinction)
resonance_wl = wavelengths[peak_idx]

print(f"LSPR 共振波长: {resonance_wl:.3f} μm ({resonance_wl * 1000:.0f} nm)")
print(f"纳米棒长: {rod_length * 1000:.0f} nm, 直径: {rod_radius * 2 * 1000:.0f} nm")

# --- 可视化 ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# 消光光谱
axes[0].plot(wavelengths * 1000, extinction / extinction.max(),
             linewidth=1.5, color='darkorange')
axes[0].axvline(x=resonance_wl * 1000, color='red', linestyle='--',
                label=f'LSPR: {resonance_wl * 1000:.0f} nm')
axes[0].set_xlabel('Wavelength (nm)')
axes[0].set_ylabel('Extinction (a.u.)')
axes[0].set_title('Plasmonic Nanorod Extinction Spectrum')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# 场分布 (x-y 平面)
eps = sim.get_array(center=mp.Vector3(0, 0, 0),
                     size=mp.Vector3(0.3, 0.3, 0),
                     component=mp.Dielectric)
ex = sim.get_array(center=mp.Vector3(0, 0, 0),
                    size=mp.Vector3(0.3, 0.3, 0),
                    component=mp.Ex)

axes[1].imshow(np.abs(ex.transpose()), cmap='hot',
               extent=[-0.15, 0.15, -0.15, 0.15])
axes[1].set_title(f'|Ex| Field Enhancement @ {resonance_wl * 1000:.0f} nm')
axes[1].set_xlabel('x (μm)')
axes[1].set_ylabel('y (μm)')

plt.tight_layout()
plt.savefig('plasmonic_nanorod.png', dpi=150)
plt.show()
