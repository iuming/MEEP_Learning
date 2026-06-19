"""
Example: 二端口波导模式分解模板

这个脚本把“端口 S 参数”的核心步骤封装成函数，便于复制到弯曲波导、
MMI、定向耦合器等器件中。示例器件是一段直波导，因此 S21 应接近 1，
S11 应接近 0。

注意：真实多端口器件通常要为每个输入端口分别做一次归一化仿真。
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

resolution = 24
sx, sy = 16, 8
dpml = 1.0
wg_width = 0.5
wg_index = 3.4
fcen = 0.645
fwidth = 0.08
nfreq = 21
port_size = mp.Vector3(0, 3.0)
left_port = mp.Vector3(-5.5, 0)
right_port = mp.Vector3(5.5, 0)


def build_sim(include_source=True):
    """创建直波导仿真；include_source=False 可用于后续扩展。"""
    geometry = [mp.Block(
        size=mp.Vector3(mp.inf, wg_width, mp.inf),
        center=mp.Vector3(),
        material=mp.Medium(index=wg_index),
    )]

    sources = []
    if include_source:
        sources = [mp.EigenModeSource(
            src=mp.GaussianSource(frequency=fcen, fwidth=fwidth),
            center=left_port,
            size=port_size,
            direction=mp.X,
            eig_band=1,
            eig_parity=mp.ODD_Z,
        )]

    return mp.Simulation(
        cell_size=mp.Vector3(sx, sy, 0),
        geometry=geometry,
        sources=sources,
        boundary_layers=[mp.PML(dpml)],
        resolution=resolution,
    )


def add_ports(sim):
    """添加左右两个模式监视器。"""
    left = sim.add_mode_monitor(fcen, fwidth, nfreq, mp.ModeRegion(center=left_port, size=port_size))
    right = sim.add_mode_monitor(fcen, fwidth, nfreq, mp.ModeRegion(center=right_port, size=port_size))
    return left, right


def get_port_coefficients(sim, monitor):
    """返回某个端口基模的前向/后向复振幅。"""
    coeffs = sim.get_eigenmode_coefficients(monitor, [1], eig_parity=mp.ODD_Z)
    forward = coeffs.alpha[:, 0, 0]
    backward = coeffs.alpha[:, 0, 1]
    return forward, backward

sim = build_sim(include_source=True)
left_monitor, right_monitor = add_ports(sim)
sim.run(until_after_sources=mp.stop_when_fields_decayed(60, mp.Ez, right_port, 1e-8))

left_forward, left_backward = get_port_coefficients(sim, left_monitor)
right_forward, right_backward = get_port_coefficients(sim, right_monitor)
freqs = np.array(mp.get_eigenmode_freqs(right_monitor)) if hasattr(mp, "get_eigenmode_freqs") else np.linspace(fcen - fwidth/2, fcen + fwidth/2, nfreq)

# 用输入端口的前向功率归一化。直波导中，左端口 backward 是反射，右端口 forward 是透射。
input_power = np.maximum(np.abs(left_forward) ** 2, 1e-30)
S11_power = np.abs(left_backward) ** 2 / input_power
S21_power = np.abs(right_forward) ** 2 / input_power
wavelengths = 1 / freqs

print(f"S21 平均值：{np.mean(S21_power):.3f}，S11 平均值：{np.mean(S11_power):.3e}")

plt.figure(figsize=(7, 4))
plt.plot(wavelengths, S21_power, "o-", label="|S21|² 透射")
plt.plot(wavelengths, S11_power, "s-", label="|S11|² 反射")
plt.gca().invert_xaxis()
plt.xlabel("波长 / μm（若长度单位为 μm）")
plt.ylabel("模式功率比例")
plt.title("二端口波导模式分解模板")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("mode_decomposition_ports.png", dpi=160)
print("图像已保存到 mode_decomposition_ports.png")
plt.show()
