"""
Project 02: 超表面透镜 (Metalens)

设计一个基于介质纳米柱的透射式超表面透镜。
通过调整纳米柱直径实现 0~2π 的相位覆盖，
构造双曲相位剖面聚焦入射平面波。

应用：超薄平面光学元件、AR/VR 显示、手机摄像头
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# --- 设计参数 ---
resolution = 40
n_pillar   = 3.5    # 纳米柱折射率 (如 TiO2)
n_sub      = 1.45   # 基底折射率 (如 SiO2)
wavelength = 0.55   # 设计波长 550 nm (绿光)
freq       = 1.0 / wavelength
focal_length = 10.0 # 焦距 (μm)
lens_diameter = 8.0 # 透镜直径 (μm)
pillar_height = 0.6 # 纳米柱高度 (μm)

# 超表面单元参数
unit_cell_size = 0.35  # 单元大小 (亚波长)
n_pillars = int(lens_diameter / unit_cell_size)

# 所需相位剖面: φ(r) = -2π/λ (√(r² + f²) - f)
# 通过调整直径实现

# 直径-相位映射 (通过扫描单个纳米柱得到，此处为示意值)
# 实际应用中需要先跑参数扫描
min_diameter = 0.08  # 最小直径
max_diameter = 0.30  # 最大直径

def required_phase(x, y):
    """双曲相位剖面"""
    r2 = x**2 + y**2
    return -2 * np.pi / wavelength * (np.sqrt(r2 + focal_length**2) - focal_length)

def phase_to_diameter(phase):
    """
    将所需相位映射到纳米柱直径。
    简化模型：线性映射 0→2π → 直径范围 [min, max]
    """
    phase_normalized = (phase - phase.min()) % (2 * np.pi) / (2 * np.pi)
    return min_diameter + phase_normalized * (max_diameter - min_diameter)

# --- 构建超表面 ---
dpml = 1.0
pad  = 2.0
cell_x = lens_diameter + 2 * pad
cell_y = focal_length + 2 * pad + 4
cell = mp.Vector3(cell_x, cell_y, 0)

geometry = []

# SiO2 基底
geometry.append(mp.Block(
    size=mp.Vector3(cell_x, 1, mp.inf),
    center=mp.Vector3(0, -cell_y / 2 + 0.5, 0),
    material=mp.Medium(epsilon=n_sub ** 2),
))

# TiO2 纳米柱阵列
xs = np.linspace(-lens_diameter / 2, lens_diameter / 2, n_pillars)
phases = np.array([[required_phase(x, 0) for x in xs]])

for i, x in enumerate(xs):
    r = abs(x)
    if r > lens_diameter / 2:
        continue

    phase_i = required_phase(x, 0)
    d = phase_to_diameter(np.array([phase_i]))[0]

    geometry.append(mp.Block(
        size=mp.Vector3(d, pillar_height, mp.inf),
        center=mp.Vector3(x, -cell_y / 2 + 1 + pillar_height / 2, 0),
        material=mp.Medium(epsilon=n_pillar ** 2),
    ))

# --- 光源 (平面波从下方入射) ---
source_y = -cell_y / 2 + dpml + 0.5
sources = [mp.Source(
    mp.GaussianSource(frequency=freq, fwidth=0.1 * freq),
    component=mp.Ez,
    center=mp.Vector3(0, source_y, 0),
    size=mp.Vector3(lens_diameter + 2, 0, 0),
)]

# --- 仿真 ---
sim = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    sources=sources,
    boundary_layers=[mp.PML(dpml)],
    resolution=resolution,
)

# 在焦平面放置监视器
focal_plane_y = -cell_y / 2 + 1 + pillar_height + focal_length

sim.run(until=300)

# --- 焦斑分析 ---
# 提取焦平面的场
ez_focal = []
nx_sample = 200
for x in np.linspace(-lens_diameter / 2, lens_diameter / 2, nx_sample):
    ez_focal.append(sim.get_field_point(mp.Ez, mp.Vector3(x, focal_plane_y, 0)))

ez_focal = np.array(ez_focal)
intensity = np.abs(ez_focal) ** 2
intensity_norm = intensity / np.max(intensity)

xs_sample = np.linspace(-lens_diameter / 2, lens_diameter / 2, nx_sample)

# 焦斑半宽 (FWHM)
half_max = 1.0 / 2
above_hm = intensity_norm >= half_max
fwhm_indices = np.where(above_hm)[0]
if len(fwhm_indices) >= 2:
    fwhm = xs_sample[fwhm_indices[-1]] - xs_sample[fwhm_indices[0]]
else:
    fwhm = np.nan

# 衍射极限: λ/(2·NA)
NA = lens_diameter / 2 / focal_length
diffraction_limit = wavelength / (2 * NA)

# 聚焦效率
# 积分焦平面能量
focal_power = np.sum(intensity) * (xs_sample[1] - xs_sample[0])

# --- 结果输出 ---
print(f"=== Metalens Performance ===")
print(f"波长:         {wavelength:.3f} μm")
print(f"NA:           {NA:.3f}")
print(f"FWHM:         {fwhm:.3f} μm")
print(f"衍射极限:     {diffraction_limit:.3f} μm")
print(f"Strehl 近似:  {diffraction_limit / fwhm:.2f}" if not np.isnan(fwhm) else "N/A")

# --- 可视化 ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# 结构图
eps = sim.get_array(component=mp.Dielectric)
axes[0].imshow(eps.transpose(), interpolation='spline36', cmap='binary',
               extent=[-cell_x/2, cell_x/2, -cell_y/2, cell_y/2])
axes[0].axhline(y=focal_plane_y, color='red', linestyle='--',
                label='Focal Plane')
axes[0].set_title('Metalens Structure')
axes[0].set_xlabel('x (μm)')
axes[0].set_ylabel('y (μm)')
axes[0].legend()

# 焦斑剖面
axes[1].plot(xs_sample, intensity_norm, linewidth=1.5)
axes[1].axhline(y=0.5, color='gray', linestyle=':', alpha=0.5,
                label='Half Maximum')
if not np.isnan(fwhm):
    axes[1].axvspan(-fwhm / 2, fwhm / 2, alpha=0.1, color='red',
                    label=f'FWHM = {fwhm:.3f} μm')
axes[1].set_xlabel('x (μm)')
axes[1].set_ylabel('Normalized Intensity')
axes[1].set_title('Focal Spot Profile')
axes[1].legend()
axes[1].grid(True, alpha=0.3)
axes[1].set_xlim(-lens_diameter / 2, lens_diameter / 2)

plt.tight_layout()
plt.savefig('metalens.png', dpi=150)
plt.show()
