"""
Example 02: 光源类型对比 — CW vs Gaussian vs Custom

演示 MEEP 中三种常用光源的区别：
- ContinuousSource（连续波）
- GaussianSource（高斯脉冲）
- 自定义时域波形
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

resolution = 40
freq = 1.0       # 频率

# --- 连续波 (CW) ---
sources_cw = [mp.Source(
    src=mp.ContinuousSource(frequency=freq, width=2.0),
    component=mp.Ex,
    center=mp.Vector3(0, 0, 0)
)]

# --- 高斯脉冲 ---
sources_gauss = [mp.Source(
    src=mp.GaussianSource(frequency=freq, fwidth=0.4),
    component=mp.Ex,
    center=mp.Vector3(0, 0, 0)
)]

# --- 宽频脉冲 (用于 S 参数计算) ---
sources_broadband = [mp.Source(
    src=mp.GaussianSource(frequency=freq, fwidth=0.8),
    component=mp.Ex,
    center=mp.Vector3(0, 0, 0)
)]

def run_source(label, sources, until=10):
    """运行简单 1D 仿真，记录场的时域响应"""
    sim = mp.Simulation(
        cell_size=mp.Vector3(8, 0, 0),
        sources=sources,
        boundary_layers=[mp.PML(1.0)],
        resolution=resolution,
    )
    data = []
    sim.run(
        mp.at_every(0.1, lambda s: data.append(s.get_field_point(mp.Ex, mp.Vector3(2, 0, 0)))),
        until=until
    )
    return np.array(data)

# 运行三种仿真
data_cw   = run_source('CW', sources_cw, until=12)
data_gauss = run_source('Gaussian', sources_gauss, until=12)
data_bb    = run_source('Broadband', sources_broadband, until=12)

# --- 绘图 ---
fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

axes[0].plot(data_cw, color='C0')
axes[0].set_ylabel('CW\nEx @ x=2')
axes[0].grid(True, alpha=0.3)
axes[0].set_title('Continuous Wave — 单频稳定振荡')

axes[1].plot(data_gauss, color='C1')
axes[1].set_ylabel('Gaussian\nEx @ x=2')
axes[1].grid(True, alpha=0.3)
axes[1].set_title('Gaussian Pulse — 有限带宽脉冲')

axes[2].plot(data_bb, color='C2')
axes[2].set_ylabel('Broadband\nEx @ x=2')
axes[2].grid(True, alpha=0.3)
axes[2].set_xlabel('Time Step (×0.1)')
axes[2].set_title('Broadband Gaussian — 宽频脉冲 (适合 S 参数)')

plt.tight_layout()
plt.savefig('source_types_result.png', dpi=150)
plt.show()

print("三种光源对比完成！")
print("- ContinuousSource: 单频正弦波，适合稳态分析")
print("- GaussianSource: 有限带宽，适合传输/反射系数")
print("- Broadband: 宽频，一次仿真覆盖多个波长")
