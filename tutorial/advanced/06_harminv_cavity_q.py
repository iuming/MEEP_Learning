"""
Advanced 06: Harminv 提取谐振腔频率与 Q 值

Harminv 会对时间信号做谐波反演，估计共振频率、衰减率和品质因子 Q。
它非常适合微环、光子晶体缺陷腔、Fabry-Perot 腔等高 Q 结构的初步分析。

模型：2D 介质圆盘腔。用短脉冲激发后，在腔边缘记录 Ez 衰减信号，
再由 Harminv 自动寻找指定频带内的共振模式。
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# ------------------------- 参数 -------------------------
resolution = 30
sx = sy = 10
dpml = 1.0
radius = 1.2
index = 3.4
fcen = 0.32                  # 搜索中心频率
df = 0.18                    # 搜索带宽

cell = mp.Vector3(sx, sy, 0)
geometry = [mp.Cylinder(
    radius=radius,
    height=mp.inf,
    center=mp.Vector3(),
    material=mp.Medium(index=index),
)]

# 源放在圆盘边缘附近，更容易激发 whispering-gallery-like 模式。
source_point = mp.Vector3(radius * 0.8, 0)
sources = [mp.Source(
    src=mp.GaussianSource(frequency=fcen, fwidth=df),
    component=mp.Ez,
    center=source_point,
)]

sim = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    sources=sources,
    boundary_layers=[mp.PML(dpml)],
    resolution=resolution,
)

probe_point = mp.Vector3(radius * 0.9, 0)
harminv = mp.Harminv(mp.Ez, probe_point, fcen, df)
trace = []

def record_after_source(sim):
    """源关闭后记录自由衰减信号，用于画图检查是否真的在衰减。"""
    trace.append(sim.get_field_point(mp.Ez, probe_point))

# after_sources(harminv) 告诉 MEEP 在源结束后把时间信号交给 Harminv。
sim.run(
    mp.after_sources(harminv),
    mp.after_sources(mp.at_every(1.0, record_after_source)),
    until_after_sources=450,
)

print("\nHarminv 找到的模式：")
print("  freq        Q          |amp|       decay")
for m in harminv.modes:
    print(f"  {m.freq:9.6f}  {m.Q:9.1f}  {abs(m.amp):9.3e}  {m.decay:9.3e}")

# 选择振幅最大的模式作为主模式，方便新手快速读结果。
if harminv.modes:
    main_mode = max(harminv.modes, key=lambda m: abs(m.amp))
    print(f"\n主模式估计：频率={main_mode.freq:.6f}, Q={main_mode.Q:.1f}, 波长={1/main_mode.freq:.3f}")
else:
    print("未找到模式：可尝试增大 df、延长 until_after_sources 或改变源位置。")

# 画介电常数和最终场分布。
eps = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
ez = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Ez)

fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
axes[0].plot(np.real(trace))
axes[0].set_title("源关闭后的 Ez 衰减信号")
axes[0].set_xlabel("采样点")
axes[0].set_ylabel("Ez")
axes[0].grid(alpha=0.3)

im0 = axes[1].imshow(eps.T, origin="lower", cmap="gray", extent=[-sx/2, sx/2, -sy/2, sy/2])
axes[1].set_title("圆盘腔介电常数")
axes[1].set_xlabel("x")
axes[1].set_ylabel("y")
fig.colorbar(im0, ax=axes[1], fraction=0.046)

im1 = axes[2].imshow(np.real(ez).T, origin="lower", cmap="RdBu", extent=[-sx/2, sx/2, -sy/2, sy/2])
axes[2].set_title("最终 Ez 场")
axes[2].set_xlabel("x")
axes[2].set_ylabel("y")
fig.colorbar(im1, ax=axes[2], fraction=0.046)

plt.tight_layout()
plt.savefig("harminv_cavity_q.png", dpi=160)
print("图像已保存到 harminv_cavity_q.png")
plt.show()
