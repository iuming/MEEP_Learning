"""
Project 01: 波分复用器 (WDM) — 基于级联环形谐振腔

设计一个 4 通道 DWDM (密集波分复用) 器件。
使用不同半径的微环谐振腔选择不同波长。

原理：
- 每个环谐振腔的谐振条件: 2πR·n_eff = m·λ
- 通过改变半径 R 选择不同波长通道
- ITU 标准：100 GHz 间隔 → Δλ ≈ 0.8 nm @ 1550 nm
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# --- 设计参数 ---
resolution = 40
n_wg       = 2.5
n_ring     = 2.8
wg_width   = 0.5
gap        = 0.3
dpml       = 1.0

# ITU DWDM 波长 (100 GHz 间隔, 1550 nm 附近)
channels = np.array([1550.92, 1551.72, 1552.52, 1553.32])  # nm
freqs_ch = 1.0 / (channels * 1e-3)  # 归一化到 1/μm

freq_center = np.mean(freqs_ch)
freq_width  = (max(freqs_ch) - min(freqs_ch)) * 1.1

# 各环的谐振半径
m = 60  # 模式阶数
n_eff = n_ring  # 近似有效折射率
radii = m * 1.0 / (2 * np.pi * n_eff * freqs_ch)  # 2πR n_eff = mλ

print("DWDM 通道设计:")
for i, (wl, r) in enumerate(zip(channels, radii)):
    print(f"  Ch{i + 1}: λ = {wl:.2f} nm, R = {r:.3f} μm")

# --- 构建几何体 ---
ring_spacing = max(radii) * 2 + gap + wg_width + 2
cell_width  = 2 * dpml + 4 * ring_spacing + 6
cell_height = 8

cell = mp.Vector3(cell_width, cell_height, 0)

geometry = []

# 总线波导
geometry.append(mp.Block(
    size=mp.Vector3(mp.inf, wg_width, mp.inf),
    center=mp.Vector3(0, radii[0] + gap + wg_width / 2 + 2, 0),
    material=mp.Medium(epsilon=n_wg ** 2),
))

# 输出波导（在下方）
for i in range(4):
    x_pos = -cell_width / 2 + dpml + 2 + i * ring_spacing
    out_y = -cell_height / 2 + dpml + 2
    geometry.append(mp.Block(
        size=mp.Vector3(wg_width, 3, mp.inf),
        center=mp.Vector3(x_pos, out_y + 0.5, 0),
        material=mp.Medium(epsilon=n_wg ** 2),
    ))

# 微环谐振腔
for i, r in enumerate(radii):
    x_pos = -cell_width / 2 + dpml + 2 + i * ring_spacing
    ring_cy = radii[0] + gap + wg_width / 2 + 2 - r - gap
    geometry.append(mp.Cylinder(
        radius=r + wg_width / 2,
        center=mp.Vector3(x_pos, ring_cy, 0),
        height=mp.inf,
        material=mp.Medium(epsilon=n_ring ** 2),
    ))
    # 内孔 (空心)
    geometry.append(mp.Cylinder(
        radius=r - wg_width / 2,
        center=mp.Vector3(x_pos, ring_cy, 0),
        height=mp.inf,
        material=mp.Medium(epsilon=1.0),
    ))

# --- 光源 ---
sources = [mp.Source(
    mp.GaussianSource(frequency=freq_center, fwidth=freq_width),
    component=mp.Ez,
    center=mp.Vector3(-cell_width / 2 + dpml + 1,
                        radii[0] + gap + wg_width / 2 + 2, 0),
    size=mp.Vector3(0, wg_width * 2, 0),
)]

# --- 仿真 ---
sim = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    sources=sources,
    boundary_layers=[mp.PML(dpml)],
    resolution=resolution,
)

# 在每个输出端口放置监视器
nfreq = 200
monitors = []
for i in range(4):
    x_pos = -cell_width / 2 + dpml + 2 + i * ring_spacing
    out_y = -cell_height / 2 + dpml + 0.5
    fr = mp.FluxRegion(
        center=mp.Vector3(x_pos, out_y, 0),
        size=mp.Vector3(wg_width * 2, 0, 0),
    )
    monitors.append(sim.add_flux(
        np.linspace(freq_center - freq_width, freq_center + freq_width, nfreq),
        1, fr,
    ))

# 运行
sim.run(until=500)

# --- 透射谱分析 ---
freqs = np.linspace(freq_center - freq_width, freq_center + freq_width, nfreq)
wavelengths = 1.0 / freqs * 1000  # → nm

fig, axes = plt.subplots(2, 1, figsize=(12, 8))

colors = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12']

for i in range(4):
    flux_i = np.array(mp.get_fluxes(monitors[i]))
    T_norm = flux_i / flux_i.max()  # 归一化
    axes[0].plot(wavelengths, T_norm, color=colors[i],
                 label=f'Ch{i + 1}: {channels[i]:.2f} nm', linewidth=1.5)
    axes[0].axvline(x=channels[i], color=colors[i], linestyle='--', alpha=0.4)

axes[0].set_xlabel('Wavelength (nm)')
axes[0].set_ylabel('Normalized Transmission')
axes[0].set_title('DWDM 4-Channel Demultiplexer Spectrum')
axes[0].legend(loc='upper right')
axes[0].grid(True, alpha=0.3)
axes[0].invert_xaxis()

# 串扰矩阵
crosstalk = np.zeros((4, 4))
for i in range(4):
    flux_i = np.array(mp.get_fluxes(monitors[i]))
    for j in range(4):
        # 在通道 j 的波长处评估通道 i 的响应
        idx = np.argmin(np.abs(freqs - freqs_ch[j]))
        crosstalk[i, j] = flux_i[idx]

# 归一化
crosstalk_norm = np.zeros_like(crosstalk)
for j in range(4):
    crosstalk_norm[:, j] = crosstalk[:, j] / np.max(crosstalk[:, j])

im = axes[1].imshow(crosstalk_norm, cmap='YlOrRd', vmin=0, vmax=1)
axes[1].set_xticks(range(4))
axes[1].set_xticklabels([f'Input\n{channels[i]:.1f}nm' for i in range(4)])
axes[1].set_yticks(range(4))
axes[1].set_yticklabels([f'Output\nCh{i + 1}' for i in range(4)])
axes[1].set_title('Crosstalk Matrix')

for i in range(4):
    for j in range(4):
        text_color = 'white' if crosstalk_norm[i, j] > 0.5 else 'black'
        axes[1].text(j, i, f'{crosstalk_norm[i, j]:.2f}',
                     ha='center', va='center', color=text_color)

plt.colorbar(im, ax=axes[1])
plt.tight_layout()
plt.savefig('wdm_demux.png', dpi=150)
plt.show()

# 性能汇总
print("\n=== DWDM Demux Performance ===")
for j in range(4):
    isolation = -10 * np.log10(np.max(crosstalk_norm[:, j]) / \
                                (crosstalk_norm[j, j] + 1e-10) + 1e-10)
    print(f"  Ch{j + 1} ({channels[j]:.2f} nm): isolation ≈ {isolation:.1f} dB")
