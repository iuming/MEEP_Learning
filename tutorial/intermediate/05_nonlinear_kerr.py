"""
Intermediate 05: 非线性光学 — Kerr 效应

演示三阶非线性光学效应（Kerr 效应）。
当光强足够大时，介质折射率随光强变化：n = n0 + n2*I

仿真自相位调制 (SPM) 导致的光谱展宽。
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# --- 参数 ---
resolution = 40
n_linear   = 1.5     # n0 (线性折射率)
n2         = 1e-17   # Kerr 系数 (m²/W) — 典型值
# meep 使用 χ(3) 而非 n2，转换关系:
# χ(3) [SI] = (4/3) n0² ε0 c n2
# 简化: 设置 chi3 = n2 * (4/3) * n_linear²
chi3 = n2 * (4 / 3) * n_linear ** 2

wavelength = 1.55
freq       = 1.0 / wavelength
fwidth     = 0.1 * freq

cell = mp.Vector3(15, 2, 0)
dpml = 1.0

# --- 非线性介质 ---
kerr_material = mp.Medium(
    epsilon=n_linear ** 2,
    chi3=chi3,  # 三阶非线性极化率
)

# 非线性区域
geometry = [
    mp.Block(
        size=mp.Vector3(8, mp.inf, mp.inf),
        center=mp.Vector3(0, 0, 0),
        material=kerr_material,
    ),
]

# --- 高功率脉冲 ---
sources = [mp.Source(
    mp.GaussianSource(frequency=freq, fwidth=fwidth),
    component=mp.Ey,
    center=mp.Vector3(-5, 0, 0),
    amplitude=10.0,  # 高振幅 → 高光强 → 非线性显著
)]

# --- 仿真 (非线性) ---
sim_nl = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    sources=sources,
    boundary_layers=[mp.PML(dpml)],
    resolution=resolution,
)

# 记录输出端时域信号
t_nl = []
ey_nl = []

def record_nl(sim):
    t_nl.append(sim.meep_time())
    ey_nl.append(sim.get_field_point(mp.Ey, mp.Vector3(5, 0, 0)))

sim_nl.run(
    mp.at_every(0.1, record_nl),
    until=60
)

# --- 仿真 (线性参考，低功率) ---
sim_lin = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    sources=[mp.Source(
        mp.GaussianSource(frequency=freq, fwidth=fwidth),
        component=mp.Ey,
        center=mp.Vector3(-5, 0, 0),
        amplitude=0.01,  # 低振幅 → 线性区
    )],
    boundary_layers=[mp.PML(dpml)],
    resolution=resolution,
)

t_lin = []
ey_lin = []

def record_lin(sim):
    t_lin.append(sim.meep_time())
    ey_lin.append(sim.get_field_point(mp.Ey, mp.Vector3(5, 0, 0)))

sim_lin.run(
    mp.at_every(0.1, record_lin),
    until=60
)

# --- 频谱分析 ---
ey_nl_arr = np.array(ey_nl)
ey_lin_arr = np.array(ey_lin)

# FFT
nfft = 2048
spec_nl = np.abs(np.fft.fft(ey_nl_arr * np.hanning(len(ey_nl_arr)), nfft))
spec_lin = np.abs(np.fft.fft(ey_lin_arr * np.hanning(len(ey_lin_arr)), nfft))

dt = 0.1
freqs_fft = np.fft.fftfreq(nfft, dt)
# 只显示正频率
pos_idx = freqs_fft > 0

# --- 可视化 ---
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# 时域 - 线性
axes[0, 0].plot(t_lin, ey_lin_arr, linewidth=1)
axes[0, 0].set_title('Linear (Low Power) — No SPM')
axes[0, 0].set_xlabel('Time')
axes[0, 0].set_ylabel('Ey')
axes[0, 0].grid(True, alpha=0.3)

# 时域 - 非线性
axes[0, 1].plot(t_nl, ey_nl_arr, linewidth=1, color='C3')
axes[0, 1].set_title('Nonlinear (High Power) — SPM Broadening')
axes[0, 1].set_xlabel('Time')
axes[0, 1].set_ylabel('Ey')
axes[0, 1].grid(True, alpha=0.3)

# 频谱 - 线性
axes[1, 0].plot(freqs_fft[pos_idx], spec_lin[pos_idx], linewidth=1.5)
axes[1, 0].axvline(x=freq, color='red', linestyle='--', alpha=0.5,
                   label=f'Input: {freq:.3f}')
axes[1, 0].set_title('Linear Spectrum')
axes[1, 0].set_xlabel('Frequency')
axes[1, 0].set_ylabel('|FFT|')
axes[1, 0].legend()
axes[1, 0].grid(True, alpha=0.3)

# 频谱 - 线性
axes[1, 1].plot(freqs_fft[pos_idx], spec_nl[pos_idx], linewidth=1.5, color='C3')
axes[1, 1].axvline(x=freq, color='red', linestyle='--', alpha=0.5,
                   label=f'Input: {freq:.3f}')
axes[1, 1].set_title('Nonlinear Spectrum — Broadened!')
axes[1, 1].set_xlabel('Frequency')
axes[1, 1].set_ylabel('|FFT|')
axes[1, 1].legend()
axes[1, 1].grid(True, alpha=0.3)

# 计算光谱展宽因子 (FWHM)
def fwhm(freqs, spec):
    half_max = np.max(spec) / 2
    above = spec >= half_max
    fwhm_ = np.sum(above) * (freqs[1] - freqs[0])
    return fwhm_

fwhm_lin = fwhm(freqs_fft[pos_idx], spec_lin[pos_idx])
fwhm_nl  = fwhm(freqs_fft[pos_idx], spec_nl[pos_idx])
broadening = fwhm_nl / fwhm_lin

print(f"FWHM 线性: {fwhm_lin:.4f}")
print(f"FWHM 非线性: {fwhm_nl:.4f}")
print(f"展宽因子: {broadening:.2f}×")

plt.suptitle(f'SPM Spectral Broadening: {broadening:.2f}×',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('kerr_nonlinear.png', dpi=150)
plt.show()
