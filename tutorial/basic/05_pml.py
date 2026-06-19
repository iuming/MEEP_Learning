"""
Example 05: PML 吸收边界条件

演示 Perfectly Matched Layer (PML) 的吸收效果。
对比有无 PML 的反射，验证 PML 的反射抑制能力。
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

resolution = 20
freq = 1.0
fwidth = 0.5
cell_size = mp.Vector3(10, 0, 0)

sources = [mp.Source(
    src=mp.GaussianSource(frequency=freq, fwidth=fwidth),
    component=mp.Ez,
    center=mp.Vector3(-3, 0, 0)
)]

# --- 方案1: 有 PML ---
sim_pml = mp.Simulation(
    cell_size=cell_size,
    sources=sources,
    boundary_layers=[mp.PML(thickness=1.0)],
    resolution=resolution,
)

pml_data = []
sim_pml.run(
    mp.at_every(0.5, lambda s: pml_data.append(
        s.get_field_point(mp.Ez, mp.Vector3(3, 0, 0))
    )),
    until=15
)

# --- 方案2: 无 PML（金属边界） ---
sim_nopml = mp.Simulation(
    cell_size=cell_size,
    sources=sources,
    boundary_layers=[],
    resolution=resolution,
)

nopml_data = []
sim_nopml.run(
    mp.at_every(0.5, lambda s: nopml_data.append(
        s.get_field_point(mp.Ez, mp.Vector3(3, 0, 0))
    )),
    until=15
)

# --- 绘图 ---
fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

axes[0].plot(np.array(pml_data), linewidth=1)
axes[0].set_title('With PML — 脉冲到达后无反射')
axes[0].set_ylabel('Ez @ x=3')
axes[0].grid(True, alpha=0.3)

axes[1].plot(np.array(nopml_data), linewidth=1)
axes[1].set_title('Without PML (金属边界) — 明显反射回波')
axes[1].set_xlabel('Time Step (×0.5)')
axes[1].set_ylabel('Ez @ x=3')
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('pml_comparison.png', dpi=150)
plt.show()

print("PML 对比完成！")
print("■ 带 PML: 电磁波被完全吸收，无反射")
print("■ 无 PML: 金属边界产生强烈反射回波")
