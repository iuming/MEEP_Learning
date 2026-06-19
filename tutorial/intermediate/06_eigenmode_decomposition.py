"""
Intermediate 06: 本征模源与模式分解 — 波导端口 S 参数的基础

MEEP 的 EigenModeSource 可以直接激发波导本征模式，比普通点源/线源更适合
集成光子器件。add_mode_monitor + get_eigenmode_coefficients 可把端口处的场
分解到前向/后向模式，得到更接近实验端口定义的透射和反射。

模型：一条 2D 介质直波导。先进行空波导仿真，并在输出端口读取基模前向
系数。理想情况下透射接近 1；偏差主要来自有限仿真时间、PML 和离散误差。
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# ------------------------- 参数 -------------------------
resolution = 24
sx, sy = 16, 8
dpml = 1.0
waveguide_width = 0.5
waveguide_index = 3.4
fcen = 0.645                 # 约等于 1.55 μm 波长
fwidth = 0.08
nfreq = 31

cell = mp.Vector3(sx, sy, 0)

# 直波导贯穿计算区域；x 方向很长，避免端面散射。
geometry = [mp.Block(
    size=mp.Vector3(mp.inf, waveguide_width, mp.inf),
    center=mp.Vector3(),
    material=mp.Medium(index=waveguide_index),
)]

sources = [mp.EigenModeSource(
    src=mp.GaussianSource(frequency=fcen, fwidth=fwidth),
    center=mp.Vector3(-0.5 * sx + dpml + 1.0, 0),
    size=mp.Vector3(0, 3.0),                  # 截面要覆盖模式尾场
    direction=mp.X,
    eig_band=1,                               # 基模
    eig_parity=mp.ODD_Z,                      # 2D Ez-like 模式常用设置
)]

sim = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    sources=sources,
    boundary_layers=[mp.PML(dpml)],
    resolution=resolution,
)

# 输出端模式监视器：放在右侧，离 PML 和源都要有距离。
monitor_center = mp.Vector3(0.5 * sx - dpml - 1.5, 0)
mode_region = mp.ModeRegion(center=monitor_center, size=mp.Vector3(0, 3.0))
mode_monitor = sim.add_mode_monitor(fcen, fwidth, nfreq, mode_region)

# 同时放一个普通通量监视器，便于比较“总功率”和“基模功率”。
flux_monitor = sim.add_flux(fcen, fwidth, nfreq, mp.FluxRegion(center=monitor_center, size=mp.Vector3(0, 3.0)))

sim.run(until_after_sources=mp.stop_when_fields_decayed(60, mp.Ez, monitor_center, 1e-8))

# alpha[频率编号, 模式编号, 方向]；方向 0 为 +x，方向 1 为 -x。
coeffs = sim.get_eigenmode_coefficients(mode_monitor, [1], eig_parity=mp.ODD_Z)
alpha = coeffs.alpha[:, 0, 0]
mode_power = np.abs(alpha) ** 2

freqs = np.array(mp.get_flux_freqs(flux_monitor))
flux_power = np.array(mp.get_fluxes(flux_monitor))
wavelengths = 1 / freqs

# 归一化到中心频率附近的最大值，得到直观曲线；真实器件需先做空结构参考。
mode_power_norm = mode_power / np.max(mode_power)
flux_power_norm = flux_power / np.max(np.abs(flux_power))

print("中心频率附近的基模归一化透射：", mode_power_norm[len(mode_power_norm)//2])

plt.figure(figsize=(7, 4))
plt.plot(wavelengths, mode_power_norm, "o-", label="基模功率 |alpha|²")
plt.plot(wavelengths, flux_power_norm, "s--", label="总通量功率")
plt.gca().invert_xaxis()
plt.xlabel("波长 / μm（若长度单位为 μm）")
plt.ylabel("归一化功率")
plt.title("EigenModeSource + 模式分解示例")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("eigenmode_decomposition.png", dpi=160)
print("图像已保存到 eigenmode_decomposition.png")
plt.show()
