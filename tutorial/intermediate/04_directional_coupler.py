"""
Intermediate 04: 定向耦合器与模式拍长

两个平行波导之间的倏逝波耦合。
计算耦合长度 (beat length)，实现完整功率转移。
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# --- 参数 ---
resolution = 40
n_core     = 3.0     # 波导折射率
n_clad     = 1.0     # 包层折射率
wg_width   = 0.5     # 波导宽度
gap        = 0.3     # 波导间距
wavelength = 1.55
freq       = 1.0 / wavelength
fwidth     = 0.05 * freq

# 耦合器尺寸
coupler_length = 20.0  # 总长度
port_length    = 2.0   # 输入/输出直波导长度
tot_length     = coupler_length + 2 * port_length

cell_x = tot_length + 4
cell_y = 8
cell = mp.Vector3(cell_x, cell_y, 0)

# --- 几何体 ---
y_top = gap / 2 + wg_width / 2
y_bot = -gap / 2 - wg_width / 2

geometry = [
    # 上方波导
    mp.Block(
        size=mp.Vector3(mp.inf, wg_width, mp.inf),
        center=mp.Vector3(0, y_top, 0),
        material=mp.Medium(epsilon=n_core ** 2),
    ),
    # 下方波导
    mp.Block(
        size=mp.Vector3(mp.inf, wg_width, mp.inf),
        center=mp.Vector3(0, y_bot, 0),
        material=mp.Medium(epsilon=n_core ** 2),
    ),
]

# --- 光源（只激励下方波导） ---
sources = [mp.Source(
    mp.GaussianSource(frequency=freq, fwidth=fwidth),
    component=mp.Ez,
    center=mp.Vector3(-tot_length / 2 + port_length / 2, y_bot, 0),
    size=mp.Vector3(0, wg_width * 2, 0),
)]

# --- 仿真 ---
sim = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    sources=sources,
    boundary_layers=[mp.PML(1.0)],
    resolution=resolution,
    default_material=mp.Medium(epsilon=n_clad ** 2),
)

# 输出端口监视器
out_bot = sim.add_flux(
    freq, fwidth, 1,
    mp.FluxRegion(
        center=mp.Vector3(tot_length / 2 - port_length / 2, y_bot, 0),
        size=mp.Vector3(0, wg_width * 2, 0),
    )
)
out_top = sim.add_flux(
    freq, fwidth, 1,
    mp.FluxRegion(
        center=mp.Vector3(tot_length / 2 - port_length / 2, y_top, 0),
        size=mp.Vector3(0, wg_width * 2, 0),
    )
)

# 沿耦合区域放置功率监视点
nx_sample = 50
sample_xs = np.linspace(-coupler_length / 2, coupler_length / 2, nx_sample)
P_top = np.zeros(nx_sample)
P_bot = np.zeros(nx_sample)

def record_power(sim):
    for i, x in enumerate(sample_xs):
        top_ey = sim.get_field_point(mp.Ez, mp.Vector3(x, y_top, 0))
        bot_ey = sim.get_field_point(mp.Ez, mp.Vector3(x, y_bot, 0))
        P_top[i] += abs(top_ey) ** 2
        P_bot[i] += abs(bot_ey) ** 2

sim.run(mp.at_every(1.0, record_power), until=300)

# --- 结果 ---
P_total_top = np.array(mp.get_fluxes(out_top))[0]
P_total_bot = np.array(mp.get_fluxes(out_bot))[0]

ratio = P_total_top / (P_total_top + P_total_bot)

print(f"上方波导输出功率: {P_total_top:.4f}  ({ratio * 100:.1f}%)")
print(f"下方波导输出功率: {P_total_bot:.4f}  ({(1 - ratio) * 100:.1f}%)")

# 估算拍长（通过寻找功率振荡周期）
from scipy.optimize import curve_fit

# 归一化
P_top_norm = P_top / (P_top + P_bot + 1e-10)

# 拟合正弦振荡
def sinusoid(x, A, k, phi, offset):
    return A * np.sin(k * x + phi) + offset

try:
    popt, _ = curve_fit(sinusoid, sample_xs, P_top_norm,
                         p0=[0.5, 2 * np.pi / 10, 0, 0.5])
    beat_length = 2 * np.pi / popt[1]
    print(f"拍长 (beat length): {beat_length:.2f} μm")
except:
    beat_length = coupler_length / 2
    print("拍长拟合失败，使用近似值")

# --- 可视化 ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# 功率演化
axes[0].plot(sample_xs, P_top_norm, label='Top waveguide', linewidth=1.5)
axes[0].plot(sample_xs, 1 - P_top_norm, '--', label='Bottom waveguide', linewidth=1.5)
axes[0].axhline(y=0.5, color='gray', linestyle=':', alpha=0.5)
axes[0].set_xlabel('Position x (μm)')
axes[0].set_ylabel('Normalized Power')
axes[0].set_title(f'Directional Coupler: gap={gap}μm, L={coupler_length}μm\n'
                  f'Beat Length ≈ {beat_length:.1f}μm')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# 结构
eps = sim.get_array(component=mp.Dielectric)
size = sim.cell_size
axes[1].imshow(eps.transpose(), interpolation='spline36', cmap='binary',
               extent=[-size.x/2, size.x/2, -size.y/2, size.y/2])
axes[1].set_title('Dielectric Structure')
axes[1].set_xlabel('x (μm)')
axes[1].set_ylabel('y (μm)')

plt.tight_layout()
plt.savefig('directional_coupler.png', dpi=150)
plt.show()
