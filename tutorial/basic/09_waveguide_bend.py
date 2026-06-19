"""
Example 09: 弯曲波导与传播损耗

对比直波导和弯曲波导的传输效率。
演示如何在 MEEP 中构建弯曲结构并计算弯曲损耗。
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# --- 参数 ---
resolution = 30
n_wg       = 2.0
n_clad     = 1.0
wg_width   = 0.5
wavelength = 1.55
freq       = 1.0 / wavelength
bend_radius = 3.0
dpml       = 1.0
pad        = 2.0

# 使用 GaussianSource 作为源
fwidth = 0.1 * freq

# --- 直波导（参考） ---
def make_straight_wg():
    cell = mp.Vector3(12, 6, 0)
    wg = mp.Block(
        size=mp.Vector3(mp.inf, wg_width, mp.inf),
        center=mp.Vector3(0, 0, 0),
        material=mp.Medium(epsilon=n_wg ** 2),
    )
    src = [mp.Source(
        mp.GaussianSource(frequency=freq, fwidth=fwidth),
        component=mp.Ez,
        center=mp.Vector3(-4, 0, 0),
        size=mp.Vector3(0, wg_width * 2, 0),
    )]
    sim = mp.Simulation(
        cell_size=cell,
        geometry=[wg],
        sources=src,
        boundary_layers=[mp.PML(dpml)],
        resolution=resolution,
        default_material=mp.Medium(epsilon=n_clad ** 2),
    )
    return sim

# --- 90度弯曲波导 ---
def make_bent_wg():
    cell = mp.Vector3(10, 10, 0)

    # 水平段
    horz = mp.Block(
        size=mp.Vector3(mp.inf, wg_width, mp.inf),
        center=mp.Vector3(-10, -bend_radius - wg_width / 2, 0),
        material=mp.Medium(epsilon=n_wg ** 2),
    )
    # 弯曲段 (使用环形截面)
    bend = mp.Cylinder(
        radius=bend_radius + wg_width,
        center=mp.Vector3(
            -bend_radius * 2, -bend_radius - wg_width / 2, 0
        ),
        height=mp.inf,
        material=mp.Medium(epsilon=n_wg ** 2),
    )
    inner = mp.Cylinder(
        radius=bend_radius,
        center=mp.Vector3(
            -bend_radius * 2, -bend_radius - wg_width / 2, 0
        ),
        height=mp.inf,
        material=mp.Medium(epsilon=n_clad ** 2),
    )
    # 垂直段
    vert = mp.Block(
        size=mp.Vector3(wg_width, mp.inf, mp.inf),
        center=mp.Vector3(-bend_radius * 2 + bend_radius + wg_width / 2, -10, 0),
        material=mp.Medium(epsilon=n_wg ** 2),
    )

    geometry = [horz, bend, inner, vert]

    src = [mp.Source(
        mp.GaussianSource(frequency=freq, fwidth=fwidth),
        component=mp.Ez,
        center=mp.Vector3(
            -bend_radius * 2 - 2, -bend_radius - wg_width / 2, 0
        ),
        size=mp.Vector3(0, wg_width * 2, 0),
    )]

    sim = mp.Simulation(
        cell_size=cell,
        geometry=geometry,
        sources=src,
        boundary_layers=[mp.PML(dpml)],
        resolution=resolution,
        default_material=mp.Medium(epsilon=n_clad ** 2),
    )
    return sim

# --- 运行直波导 ---
sim_straight = make_straight_wg()
flux_str = sim_straight.add_flux(
    freq, fwidth, 1,
    mp.FluxRegion(center=mp.Vector3(4, 0, 0), size=mp.Vector3(0, wg_width * 3, 0))
)
sim_straight.run(until=200)

# --- 运行弯曲波导 ---
sim_bent = make_bent_wg()
# 在弯曲输出端放置监视器
mon_pos_y = -bend_radius - wg_width / 2 - 4
flux_bent = sim_bent.add_flux(
    freq, fwidth, 1,
    mp.FluxRegion(
        center=mp.Vector3(
            -bend_radius * 2 + bend_radius + wg_width / 2,
            mon_pos_y, 0
        ),
        size=mp.Vector3(wg_width * 3, 0, 0),
    )
)
sim_bent.run(until=300)

# --- 结果 ---
P_straight = np.array(mp.get_fluxes(flux_str))[0]
P_bent = np.array(mp.get_fluxes(flux_bent))[0]

# 弯曲损耗 (dB)
bend_loss = -10 * np.log10(abs(P_bent / P_straight))

print(f"直波导输出功率: {P_straight:.6f}")
print(f"弯曲波导输出功率: {P_bent:.6f}")
print(f"弯曲损耗: {bend_loss:.2f} dB")

# --- 可视化 ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, sim, title in zip(
    axes,
    [sim_straight, sim_bent],
    ['Straight Waveguide', f'90° Bend (R={bend_radius}μm)']
):
    eps = sim.get_array(component=mp.Dielectric)
    size = sim.cell_size
    ax.imshow(eps.transpose(), interpolation='spline36', cmap='binary',
              extent=[-size.x/2, size.x/2, -size.y/2, size.y/2])
    ax.set_title(title)
    ax.set_xlabel('x (μm)')
    ax.set_ylabel('y (μm)')

plt.tight_layout()
plt.savefig('bend_waveguide.png', dpi=150)
plt.show()
