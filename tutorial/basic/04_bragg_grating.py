"""
Example 04: 布拉格光栅与光子带隙

周期交替的电介质层（高/低折射率）形成布拉格光栅。
计算传输谱，观察光子带隙。
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# --- 参数 ---
resolution = 50
n_high     = 2.0     # 高折射率
n_low      = 1.5     # 低折射率
d_period   = 1.0     # 周期长度 (a)
n_periods  = 10      # 周期数
d_high     = 0.5     # 高折射率层厚度
d_low      = 0.5     # 低折射率层厚度

# 布拉格频率 (中心频率): c/(2na)
freq_bragg = 1.0 / (2.0 * d_period)  # ≈ 0.5
freq_min   = freq_bragg * 0.1
freq_max   = freq_bragg * 1.9
freq_range = freq_max - freq_min
num_freqs  = 200

dpml = 1.0
pad  = 1.0

cell_length = 2 * (dpml + pad) + n_periods * d_period
cell = mp.Vector3(cell_length, 0, 0)

# --- 构建光栅几何 ---
geometry = []
for i in range(n_periods):
    center = mp.Vector3(
        -n_periods * d_period / 2 + i * d_period + d_high / 2,
        0, 0
    )
    geometry.append(mp.Block(
        size=mp.Vector3(d_high, mp.inf, mp.inf),
        center=center,
        material=mp.Medium(epsilon=n_high**2)
    ))

# --- 光源 ---
source_pos = -dpml - pad - 1.0
sources = [mp.Source(
    mp.GaussianSource(frequency=freq_bragg, fwidth=freq_range),
    component=mp.Ey,
    center=mp.Vector3(source_pos, 0, 0)
)]

# --- 监视器 ---
tran_pos = dpml + pad + 1.0
refl_pos = -dpml - pad + 0.5

sim = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    sources=sources,
    boundary_layers=[mp.PML(dpml)],
    resolution=resolution,
)

refl_fr = mp.FluxRegion(center=mp.Vector3(refl_pos, 0, 0))
tran_fr = mp.FluxRegion(center=mp.Vector3(tran_pos, 0, 0))

refl = sim.add_flux(
    np.linspace(freq_min, freq_max, num_freqs),
    1, refl_fr
)
tran = sim.add_flux(
    np.linspace(freq_min, freq_max, num_freqs),
    1, tran_fr
)

# 空仿真
sim_empty = mp.Simulation(
    cell_size=cell,
    geometry=[],
    sources=sources,
    boundary_layers=[mp.PML(dpml)],
    resolution=resolution,
)
tran0 = sim_empty.add_flux(
    np.linspace(freq_min, freq_max, num_freqs),
    1, tran_fr
)

# 运行
sim.run(until=300)
sim_empty.run(until=300)

# --- 结果 ---
freqs = np.linspace(freq_min, freq_max, num_freqs)
tran_flux = np.array(mp.get_fluxes(tran))
tran0_flux = np.array(mp.get_fluxes(tran0))
T = tran_flux / tran0_flux

# 绘图
plt.figure(figsize=(10, 5))
plt.plot(freqs, T, linewidth=1.5)
plt.axvline(x=freq_bragg, color='red', linestyle='--',
            label=f'Bragg freq: {freq_bragg:.3f}')
plt.axvspan(freq_bragg * 0.8, freq_bragg * 1.2,
            alpha=0.1, color='gray', label='Band Gap Region')
plt.xlabel('Frequency (c/a)')
plt.ylabel('Transmission')
plt.title(f'Bragg Mirror: {n_periods} periods, n={n_high}/{n_low}')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('bragg_grating_result.png', dpi=150)
plt.show()

print("布拉格光栅仿真完成！")
print(f"中心频率: {freq_bragg:.3f} c/a")
print(f"最低透射率: {T.min():.4f}")
