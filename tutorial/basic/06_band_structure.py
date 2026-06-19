"""
Example 06: 光子晶体能带结构

使用 MEEP 计算一维光子晶体的能带结构。
通过 Bloch 周期边界条件扫描 k 向量，找到本征频率。
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# --- 参数 ---
resolution = 64
n_high     = 3.0     # 高折射率
n_low      = 1.0     # 低折射率
a          = 1.0     # 晶格常数
d_high     = 0.5 * a
d_low      = 0.5 * a

# 扫描 k 向量 (0 到 π/a)
num_k = 20
k_points = np.linspace(0, np.pi / a, num_k)
num_bands = 6  # 计算 6 个能带

geometry_lattice = mp.Lattice(
    size=mp.Vector3(a, 0, 0)
)

geometry = [
    mp.Block(
        size=mp.Vector3(d_high, mp.inf, mp.inf),
        center=mp.Vector3(-a / 2 + d_high / 2, 0, 0),
        material=mp.Medium(epsilon=n_high ** 2),
    )
]

bands = []

for kx_val in k_points:
    k_vec = mp.Vector3(kx_val, 0, 0)

    sim = mp.Simulation(
        cell_size=mp.Vector3(a, 0, 0),
        geometry=geometry,
        resolution=resolution,
        default_material=mp.Medium(epsilon=n_low ** 2),
        k_point=k_vec,
        ensure_periodicity=True,
    )

    sim.init_sim()

    # 使用 harminv 找本征频率
    # 在多个位置放置随机源来激发所有模式
    for i in range(-2, 3):
        src_center = mp.Vector3(i * a / 10, 0, 0)
        sim.add_source(mp.Source(
            mp.GaussianSource(frequency=1.0, fwidth=2.0),
            component=mp.Ey,
            center=src_center,
        ))

    harm = mp.Harminv(mp.Ey, mp.Vector3(0, 0, 0), 0.01, 5.0)
    sim.run(mp.after_sources(harm), until_after_sources=200)

    # 收集本征频率
    modes = harm.modes
    freqs = [m.freq for m in modes[:num_bands]]
    freqs.sort()
    bands.append(freqs)

    print(f"  k = {kx_val:.3f} π/a, found {len(freqs)} modes")

# --- 绘图 ---
bands = np.array(bands)

plt.figure(figsize=(8, 6))
for n in range(bands.shape[1]):
    plt.plot(k_points * a / np.pi, bands[:, n], 'o-', markersize=3,
             label=f'Band {n + 1}')

plt.xlabel('Wave Vector k (π/a)')
plt.ylabel('Frequency ω (c/a)')
plt.title('1D Photonic Crystal Band Structure')
plt.legend(loc='upper right', ncol=2)
plt.grid(True, alpha=0.3)

# 标注光子带隙区域
plt.axhline(y=1.0 / (2 * a), color='gray', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig('band_structure.png', dpi=150)
plt.show()

print(f"\nBand structure calculated for {num_k} k-points.")
print(f"Gapless regions appear where bands overlap in frequency.")
