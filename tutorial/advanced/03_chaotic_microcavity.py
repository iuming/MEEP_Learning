"""
Advanced 03: 混沌光学微腔 — WG 模式与射线动力学

仿真变形微盘腔中的光学回音壁模式 (WGM)。
通过改变微盘形状研究混沌射线动力学对模式的影响。

应用：低阈值激光器、高灵敏度传感器、定向发射
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# --- 参数 ---
resolution = 80
n_cavity  = 3.3    # 腔折射率 (如 GaAs)
n_clad    = 1.0    # 空气
wavelength = 1.0
freq       = 1.0 / wavelength
fwidth     = 0.4 * freq

# 微腔参数
radius_basic = 2.0  # 基圆半径
deformation  = 0.1  # 变形参数 (四极形变)
sharpness    = 3.0  # 边界锐度 (越大边界越清晰)

dpml = 1.0
pad  = 1.0
cell_size = 2 * (radius_basic + deformation + dpml + pad)
cell = mp.Vector3(cell_size, cell_size, 0)

# --- 变形微盘形状 ---
def deformed_disk_r(theta, R, eps):
    """
    四极形变: r(θ) = R * (1 + ε * cos(2θ))
    """
    return R * (1 + eps * np.cos(2 * theta))


# 使用 MaterialGrid 定义光滑形状
n_grid = 400
rho = np.ones((n_grid, n_grid)) * n_clad ** 2
xs = np.linspace(-cell_size / 2, cell_size / 2, n_grid)
ys = np.linspace(-cell_size / 2, cell_size / 2, n_grid)

X, Y = np.meshgrid(xs, ys, indexing='ij')
R_grid = np.sqrt(X ** 2 + Y ** 2)
Theta_grid = np.arctan2(Y, X)

# 变形边界
boundary = deformed_disk_r(Theta_grid, radius_basic, deformation)

# 软化边界 (sigmoid 过渡)
delta = 0.1  # 过渡宽度
# sigmoid 软化: 1/(1+exp(-(boundary - r)/delta))
mask = 1.0 / (1.0 + np.exp(-(boundary - R_grid) / delta))
mask_smoothed = mask

# 介电常数分布
epsilon_grid = n_clad ** 2 + (n_cavity ** 2 - n_clad ** 2) * mask_smoothed

# 由于 MaterialGrid 需要 MaterialGrid 对象，这里直接用几何体近似
# 使用极坐标定义轮廓
n_vertices = 180
theta_vals = np.linspace(0, 2 * np.pi, n_vertices, endpoint=False)
r_vals = deformed_disk_r(theta_vals, radius_basic, deformation)

vertices = []
for t, r in zip(theta_vals, r_vals):
    vertices.append(mp.Vector3(r * np.cos(t), r * np.sin(t), 0))

# 使用 prism 定义多边形腔
geometry = [
    mp.Prism(
        vertices,
        height=mp.inf,
        material=mp.Medium(epsilon=n_cavity ** 2),
    )
]

# --- 寻找 WGM ---
# 使用带限脉冲激发并分析谐振模式
sources = [mp.Source(
    mp.GaussianSource(frequency=freq, fwidth=fwidth),
    component=mp.Ez,
    center=mp.Vector3(radius_basic * 0.8, 0, 0),
    size=mp.Vector3(0, 0.5, 0),
)]

# 仿真
sim = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    sources=sources,
    boundary_layers=[mp.PML(dpml)],
    resolution=resolution,
    default_material=mp.Medium(epsilon=n_clad ** 2),
)

# Harminv 分析 - 在腔内放置探针
probe_pos = mp.Vector3(radius_basic * 0.3, 0, 0)
harm = mp.Harminv(mp.Ez, probe_pos, freq, fwidth)

sim.run(mp.after_sources(harm), until_after_sources=500)

# --- 分析模式 ---
modes = harm.modes

print("=== 变形微盘腔 WGM 分析 ===")
print(f"基圆半径: R = {radius_basic} μm")
print(f"形变参数: ε = {deformation}")
print(f"折射率: n = {n_cavity}")
print()
print(f"发现 {len(modes)} 个模式:")

wg_modes = []
for i, m in enumerate(modes):
    Q = m.Q
    f = m.freq
    wl = 1.0 / f
    amp = m.amp
    err = m.err

    # 判断是否为 WGM (高 Q)
    is_wgm = Q > 100

    print(f"  Mode {i + 1}: λ = {wl:.3f} μm, Q = {Q:.0f}, "
          f"amp = {abs(amp):.4f} {'★ WGM' if is_wgm else ''}")

    if is_wgm:
        wg_modes.append(m)

# 选择最高 Q 的模式可视化
if wg_modes:
    best_mode = max(wg_modes, key=lambda m: m.Q)

    # 重新运行在谐振频率处
    sim_res = mp.Simulation(
        cell_size=cell,
        geometry=geometry,
        sources=[mp.Source(
            mp.GaussianSource(frequency=best_mode.freq, fwidth=0.01),
            component=mp.Ez,
            center=mp.Vector3(radius_basic * 0.8, 0, 0),
        )],
        boundary_layers=[mp.PML(dpml)],
        resolution=resolution,
        default_material=mp.Medium(epsilon=n_clad ** 2),
    )

    sim_res.run(until=300)

    # 获取稳态场分布
    eps_data = sim_res.get_array(component=mp.Dielectric)
    ez_data  = sim_res.get_array(component=mp.Ez)

    # --- 可视化 ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 14))

    # 介质结构
    im0 = axes[0, 0].imshow(eps_data.transpose(),
                            interpolation='spline36', cmap='binary',
                            extent=[-cell_size/2, cell_size/2,
                                    -cell_size/2, cell_size/2],
                            origin='lower')
    # 绘制理想变形边界
    theta_circle = np.linspace(0, 2 * np.pi, 200)
    boundary_x = deformed_disk_r(theta_circle, radius_basic, deformation) * \
                 np.cos(theta_circle)
    boundary_y = deformed_disk_r(theta_circle, radius_basic, deformation) * \
                 np.sin(theta_circle)
    axes[0, 0].plot(boundary_x, boundary_y, 'r-', linewidth=1,
                    label=f'Deformed boundary (ε={deformation})')
    axes[0, 0].set_title('Deformed Microdisk Cavity')
    axes[0, 0].set_xlabel('x (μm)')
    axes[0, 0].set_ylabel('y (μm)')
    axes[0, 0].legend()

    # 场分布 (强度)
    im1 = axes[0, 1].imshow(abs(ez_data.transpose()) ** 0.3,
                            interpolation='spline36', cmap='hot',
                            extent=[-cell_size/2, cell_size/2,
                                    -cell_size/2, cell_size/2],
                            origin='lower')
    axes[0, 1].set_title(f'WGM: λ = {1.0 / best_mode.freq:.3f} μm, '
                         f'Q = {best_mode.Q:.0f}')
    axes[0, 1].set_xlabel('x (μm)')
    axes[0, 1].set_ylabel('y (μm)')

    # 模式质量
    axes[1, 0].bar(range(len(modes)), [m.Q for m in modes],
                   color=['#e74c3c' if m.Q > 100 else 'gray' for m in modes])
    axes[1, 0].set_xlabel('Mode Index')
    axes[1, 0].set_ylabel('Quality Factor Q')
    axes[1, 0].set_title('Mode Quality Factors')
    axes[1, 0].set_yscale('log')
    axes[1, 0].grid(True, alpha=0.3, axis='y')
    axes[1, 0].axhline(y=100, color='gray', linestyle=':', alpha=0.5,
                       label='WGM threshold (Q=100)')
    axes[1, 0].legend()

    # 模式波长分布
    wavelengths = [1.0 / m.freq for m in modes]
    qs = [m.Q for m in modes]
    axes[1, 1].scatter(wavelengths, qs, c=qs, cmap='plasma',
                       s=100, edgecolors='black', linewidth=0.5)
    axes[1, 1].set_xlabel('Wavelength (μm)')
    axes[1, 1].set_ylabel('Q Factor')
    axes[1, 1].set_title('Mode Spectrum')
    axes[1, 1].set_yscale('log')
    axes[1, 1].grid(True, alpha=0.3)

    # 标记自由光谱范围 (FSR)
    if len(wg_modes) >= 2:
        wg_wls = sorted([1.0 / m.freq for m in wg_modes])
        fsrs = np.diff(wg_wls)
        avg_fsr = np.mean(fsrs)
        print(f"\n平均 FSR: {avg_fsr:.4f} μm")
        print(f"近似 FSR = λ²/(2πR n_g) ≈ {wavelength ** 2 / (2 * np.pi * radius_basic * n_cavity):.4f} μm")

        for wl in wg_wls:
            axes[1, 1].axvline(x=wl, color='red', linestyle='--',
                              alpha=0.2, linewidth=0.5)

    plt.suptitle('Chaotic Microcavity: Deformed Whispering-Gallery Modes',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('chaotic_microcavity.png', dpi=150)
    plt.show()

    print(f"\n最佳模式: λ = {1.0 / best_mode.freq:.3f} μm, "
          f"Q = {best_mode.Q:.0f}")
    print(f"应用: 低阈值激光器、定向发射、混沌辅助光传输")
else:
    print("未检测到高 Q WGM 模式，尝试增大变形参数或扫描频率范围")

print("\n=== 混沌腔物理 ===")
print("四极形变打破旋转对称性 → 部分射线轨道变为混沌")
print("回音壁模式仍然存在但在边界上形成焦散")
print("混沌辅助隧穿可导致定向发射")
