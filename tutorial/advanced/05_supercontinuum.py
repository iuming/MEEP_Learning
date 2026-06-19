"""
Advanced 05: 超连续谱生成 — 光子晶体光纤中的非线性传播

仿真飞秒脉冲在光子晶体光纤 (PCF) 中的超连续谱生成。
结合色散管理和非线性效应 (Kerr + 拉曼) 实现光谱展宽。

应用：宽带光源、光学相干断层扫描 (OCT)、频率梳
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# --- PCF 参数 ---
resolution = 40
n_silica   = 1.45    # 熔融二氧化硅 (线性部分)
wavelength = 1.06    # 泵浦波长 1060 nm
freq       = 1.0 / wavelength
fwidth     = 0.01 * freq  # 窄带泵浦 (100 fs 对应 ~10 nm)

# PCF 结构
pitch          = 2.3   # 空气孔间距 (μm)
hole_diameter  = 1.6   # 空气孔直径
n_holes        = 6     # 环数
fiber_length   = 100   # 光纤长度 (mm)
core_diameter  = pitch  # 实芯直径

# 非线性参数
n2_silica = 2.6e-20  # m²/W
chi3 = n2_silica * (4 / 3) * n_silica ** 2

cell_size = mp.Vector3(fiber_length * 1e-3,  # → μm
                        pitch * (2 * n_holes + 1),
                        0)  # 2D 截!

dpml = 1.0

# --- 构建 PCF 横截面 ---
geometry = []

# 背景：熔融二氧化硅
for ring in range(1, n_holes + 1):
    n_in_ring = 6 * ring
    for i in range(n_in_ring):
        angle = 2 * np.pi * i / n_in_ring + (np.pi / 6 if ring % 2 == 0 else 0)
        r = ring * pitch
        x = r * np.cos(angle)
        y = r * np.sin(angle)

        geometry.append(mp.Cylinder(
            radius=hole_diameter / 2,
            center=mp.Vector3(0, x, y),
            height=mp.inf,
            material=mp.air,
        ))

# 非线性介质 (熔融二氧化硅 + χ(3))
# 由于 MEEP 限制，这里采用简化的有效折射率方法
silica_nonlinear = mp.Medium(
    epsilon=n_silica ** 2,
    chi3=chi3,
)

# 有效介质背景
geometry.append(mp.Block(
    size=mp.Vector3(mp.inf, mp.inf, mp.inf),
    center=mp.Vector3(0, 0, 0),
    material=silica_nonlinear,
))

# --- 泵浦脉冲 ---
# 使用自定义时域包络模拟超短脉冲
pulse_duration = 0.1  # 100 fs 对应 ~ 30 μm 空间宽度
pulse_amplitude = 100.0  # 高峰值功率

sources = [mp.Source(
    mp.GaussianSource(frequency=freq, fwidth=fwidth),
    component=mp.Ez,
    center=mp.Vector3(-cell_size.x / 2 + dpml + 1, 0, 0),
    size=mp.Vector3(0, core_diameter * 0.8, core_diameter * 0.8),
    amplitude=pulse_amplitude,
)]

# --- 仿真 ---
sim = mp.Simulation(
    cell_size=cell_size,
    geometry=geometry,
    sources=sources,
    boundary_layers=[mp.PML(dpml)],
    resolution=resolution,
)

# 沿传播方向放置频谱监视器
n_monitors = 10
monitor_zs = np.linspace(
    -cell_size.x / 2 + dpml + 5,
    cell_size.x / 2 - dpml - 5,
    n_monitors
)

# 记录输出端的时域信号
tout = []
ez_out = []

def record_output(sim):
    tout.append(sim.meep_time())
    ez_out.append(sim.get_field_point(
        mp.Ez,
        mp.Vector3(cell_size.x / 2 - dpml - 2, 0, 0)
    ))

sim.run(
    mp.at_every(0.05, record_output),
    until=500,
)

# --- 频谱分析 ---
ez_out_arr = np.array(ez_out)
dt = 0.05  # 记录间隔

# 分段分析（观察沿光纤的演化）
n_segments = 5
segment_length = len(ez_out_arr) // n_segments

segments = []
for i in range(n_segments):
    seg = ez_out_arr[i * segment_length:(i + 1) * segment_length]
    windowed = seg * np.hanning(len(seg))
    spec = np.abs(np.fft.fft(windowed, 8192))
    segments.append(spec)

freqs_fft = np.fft.fftfreq(8192, dt)
pos_idx = freqs_fft > 0
wavelengths_spectrum = 1.0 / freqs_fft[pos_idx]  # 注意：这是归一化频率

# --- 超连续谱演化 ---
spec_vs_z = np.array([s[:len(pos_idx)] for s in segments])  # 取正频部分

# --- 结果 ---
print("=== 超连续谱生成分析 ===")
print(f"PCF 参数:")
print(f"  孔间距: Λ = {pitch} μm")
print(f"  孔直径: d = {hole_diameter} μm")
print(f"  d/Λ = {hole_diameter / pitch:.2f}")
print(f"  光纤长度: {fiber_length} mm")
print(f"  泵浦波长: {wavelength} nm")
print()
print(f"非线性系数: n₂ = {n2_silica:.1e} m²/W")
print(f"峰值功率: P₀ (高)")
print()

# 计算有效非线性系数
A_eff = np.pi * (core_diameter / 2) ** 2  # μm²
gamma = 2 * np.pi * n2_silica / (wavelength * 1e-9 * A_eff * 1e-12)  # 1/W/m
print(f"有效面积 A_eff: {A_eff:.2f} μm²")
print(f"非线性系数 γ: {gamma * 1000:.2f} /W/km")

# 非线性长度
L_nl = 1.0 / (gamma * 1e12)  # W → m 换算
print(f"非线性长度 L_nl: {L_nl * 1000:.1f} mm")

# --- 可视化 ---
fig, axes = plt.subplots(2, 2, figsize=(16, 14))

# PCF 横截面
# 创建横截面图
n_pixels = 200
cross_img = np.ones((n_pixels, n_pixels))
px = np.linspace(-cell_size.y / 2, cell_size.y / 2, n_pixels)
py = np.linspace(-cell_size.y / 2, cell_size.y / 2, n_pixels)
PX, PY = np.meshgrid(px, py, indexing='ij')

for ring in range(1, n_holes + 1):
    n_in_ring = 6 * ring
    for i in range(n_in_ring):
        angle = 2 * np.pi * i / n_in_ring + (np.pi / 6 if ring % 2 == 0 else 0)
        r = ring * pitch
        x = r * np.cos(angle)
        y = r * np.sin(angle)
        dist = np.sqrt((PX - x) ** 2 + (PY - y) ** 2)
        cross_img[dist < hole_diameter / 2] = 0

axes[0, 0].imshow(cross_img, cmap='gray', origin='lower',
                  extent=[-cell_size.y/2, cell_size.y/2,
                          -cell_size.y/2, cell_size.y/2])
axes[0, 0].set_title(f'PCF Cross-Section\nΛ={pitch}μm, d/Λ={hole_diameter/pitch:.2f}')
axes[0, 0].set_xlabel('x (μm)')
axes[0, 0].set_ylabel('y (μm)')

# 时域脉冲
axes[0, 1].plot(tout, ez_out_arr, linewidth=1, color='C0')
axes[0, 1].set_xlabel('Time')
axes[0, 1].set_ylabel('Ez @ output')
axes[0, 1].set_title('Output Temporal Profile')
axes[0, 1].grid(True, alpha=0.3)

# 频谱演化 瀑布图
offset = 0
colors = plt.cm.viridis(np.linspace(0, 1, n_segments))
for i in range(n_segments):
    spec_norm = segments[i][pos_idx] / segments[i][pos_idx].max()
    axes[1, 0].plot(wavelengths_spectrum[:500], spec_norm[:500] + offset,
                    color=colors[i], linewidth=1.5 - i * 0.2,
                    label=f'z = {i * fiber_length / n_segments:.0f} mm')
    offset += 0.5

axes[1, 0].axvline(x=wavelength, color='red', linestyle='--', alpha=0.5,
                   label=f'Pump: {wavelength} nm')
axes[1, 0].set_xlabel('Wavelength (μm)')
axes[1, 0].set_ylabel('Normalized Spectrum (offset)')
axes[1, 0].set_title('Supercontinuum Evolution')
axes[1, 0].legend(fontsize=8, loc='upper right')
axes[1, 0].set_xlim(0.5, 2.0)

# 最终输出频谱 (对数标度)
final_spec = segments[-1][pos_idx]
final_spec_db = 10 * np.log10(final_spec / final_spec.max() + 0.01)
axes[1, 1].plot(wavelengths_spectrum[:500], final_spec_db[:500],
                linewidth=1.5, color='darkblue')
axes[1, 1].axvline(x=wavelength, color='red', linestyle='--', alpha=0.5)
axes[1, 1].axhline(y=-3, color='gray', linestyle=':', alpha=0.5,
                   label='3 dB bandwidth')
axes[1, 1].axhline(y=-20, color='gray', linestyle=':', alpha=0.3,
                   label='20 dB bandwidth')

# 计算带宽
above_3db = final_spec_db > -3
above_20db = final_spec_db > -20
idx_3db = np.where(above_3db)[0]
idx_20db = np.where(above_20db)[0]

if len(idx_3db) > 0 and len(idx_20db) > 0:
    bw_3db = wavelengths_spectrum[idx_3db[-1]] - wavelengths_spectrum[idx_3db[0]]
    bw_20db = wavelengths_spectrum[idx_20db[-1]] - wavelengths_spectrum[idx_20db[0]]
    axes[1, 1].set_title(f'Output Spectrum\n'
                         f'BW(3dB)={bw_3db:.2f}μm, BW(20dB)={bw_20db:.2f}μm')
    print(f"\n3 dB 带宽: {bw_3db:.2f} μm")
    print(f"20 dB 带宽: {bw_20db:.2f} μm")

axes[1, 1].set_xlabel('Wavelength (μm)')
axes[1, 1].set_ylabel('Spectrum (dB)')
axes[1, 1].legend()
axes[1, 1].grid(True, alpha=0.3)
axes[1, 1].set_xlim(0.5, 2.0)
axes[1, 1].set_ylim(-40, 5)

plt.suptitle('Supercontinuum Generation in PCF',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('supercontinuum.png', dpi=150)
plt.show()

print("\n=== 超连续谱物理 ===")
print("主要机制:")
print("  1. 自相位调制 (SPM) → 光谱对称展宽")
print("  2. 受激拉曼散射 (SRS) → 向长波端红移")
print("  3. 孤子自频移 (SSFS) → 连续红移")
print("  4. 色散波产生 → 短波端蓝移分量")
print("  5. 四波混频 (FWM) → 新频率产生")
print()
print("典型应用:")
print("- 光学相干断层扫描 (OCT): 宽带光源")
print("- 频率计量: 光学频率梳")
print("- 光谱学: 宽带红外光源")
