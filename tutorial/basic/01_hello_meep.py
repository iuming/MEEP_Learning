"""
Example 01: Hello MEEP — 平面波在真空中传播

最简单的 MEEP 仿真：一个高斯脉冲在真空中传播，
观察电场 Ey 随时间的变化。
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# --- 仿真参数 ---
resolution = 20      # 每微米的网格点数
cell_size  = mp.Vector3(10, 0, 0)  # 1D 仿真区域: 10 μm
freq       = 0.15    # 中心频率 (1/μm)
fwidth     = 0.1     # 频率宽度

# --- 光源 ---
# 高斯脉冲放置在区域左侧
sources = [mp.Source(
    src=mp.GaussianSource(frequency=freq, fwidth=fwidth),
    component=mp.Ey,
    center=mp.Vector3(-4, 0, 0)
)]

# --- PML 吸收边界 ---
pml_layers = [mp.PML(thickness=1.0)]

# --- 创建仿真对象 ---
sim = mp.Simulation(
    cell_size=cell_size,
    sources=sources,
    boundary_layers=pml_layers,
    resolution=resolution,
)

# --- 记录场数据 ---
# 在区域中心放置一个点监视器，记录每个时间步的电场 Ey
ey_data = []

def record_ey(sim):
    """每个时间步调用，记录 Ey 在 x=0 处的值"""
    ey_data.append(sim.get_field_point(mp.Ey, mp.Vector3(0, 0, 0)))

# 运行仿真，直到光源能量降至设定的阈值以下
sim.run(
    mp.at_every(1.0, record_ey),  # 每个单位时间记录一次
    until=200                     # 最大运行时间
)

print(f"仿真完成，共记录 {len(ey_data)} 个时间点")

# --- 绘制结果 ---
ey_array = np.array(ey_data)

plt.figure(figsize=(10, 4))
plt.plot(ey_array)
plt.xlabel('Time Step')
plt.ylabel('Electric Field Ey')
plt.title('Hello MEEP: Ey at x=0 vs. Time')
plt.grid(True, alpha=0.3)
plt.tight_layout()

# 保存图像
plt.savefig('hello_meep_result.png', dpi=150)
print("图像已保存到 hello_meep_result.png")

plt.show()
