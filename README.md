# MEEP Learning ⚡

从零开始学习 [MEEP](https://meep.readthedocs.io/) —— MIT 开发的时域有限差分（FDTD）电磁仿真工具。

## 什么是 MEEP？

MEEP（MIT Electromagnetic Equation Propagation）是一个免费开源的 FDTD 电磁仿真软件，支持：

- 1D / 2D / 3D 电磁场仿真
- 任意介电常数和磁导率分布
- 色散材料（Drude-Lorentz 模型）
- PML 吸收边界条件
- 频域和时域分析
- Python 和 Scheme 两种接口

## 目录结构

```
MEEP_Learning/
├── tutorial/              # 入门教程（按难度递进）
│   ├── basic/             # 基础篇：光源、结构、监视器
│   ├── intermediate/      # 进阶篇（待添加）
│   └── advanced/          # 高级篇（待添加）
├── examples/              # 独立仿真示例
├── projects/              # 综合项目练习
└── README.md
```

## 安装

```bash
# Conda（推荐）
conda create -n mp -c conda-forge pymeep pymeep-extras
conda activate mp

# 或 pip
pip install meep
```

## 学习路线

### 🟢 基础篇 `tutorial/basic/`

| 编号 | 文件 | 内容 |
|------|------|------|
| 01 | `hello_meep` | 第一个仿真：平面波透过真空 |
| 02 | `source_types` | 连续波、高斯脉冲、宽频源 |
| 03 | `dielectric_slab` | 介质平板传输与反射（含解析验证） |
| 04 | `bragg_grating` | 布拉格光栅与光子带隙 |
| 05 | `pml` | PML 吸收边界条件对比 |
| 06 | `band_structure` | 光子晶体能带结构计算 |
| 07 | `waveguide` | 介质波导模式分析与场可视化 |
| 08 | `ring_resonator` | 环形谐振腔传输谱 |
| 09 | `waveguide_bend` | 弯曲波导传播损耗 |
| 10 | `grating_coupler` | 亚波长光栅耦合器设计 |

### 🟡 进阶篇 `tutorial/intermediate/`

| 编号 | 文件 | 内容 |
|------|------|------|
| 01 | `mmi_coupler` | 1×2 多模干涉耦合器设计 |
| 02 | `photonic_crystal_wg` | 光子晶体波导与慢光效应 |
| 03 | `plasmonic_antenna` | 金纳米棒 LSPR 仿真（Drude 色散模型） |
| 04 | `directional_coupler` | 定向耦合器与模式拍长 |
| 05 | `nonlinear_kerr` | Kerr 非线性（SPM 光谱展宽） |

### 🔴 高级篇 `tutorial/advanced/`

> 即将添加：拓扑优化、逆设计、量子发射体

### 🛠 实用工具 `examples/`

| 文件 | 内容 |
|------|------|
| `visualization_tools.py` | 2D/矢量场/动画可视化函数库 |
| `s_parameter_scanner.py` | 自动化 S 参数扫描框架 |
| `near_to_far_field.py` | 天线远场辐射方向图计算 |

### 🚀 综合项目 `projects/`

| 项目 | 内容 |
|------|------|
| `wdm_demux.py` | 4 通道 DWDM 波分复用器（级联微环） |
| `metalens.py` | 介质超表面透镜（TiO2 纳米柱） |

## 快速开始

```python
import meep as mp

# 仿真区域
cell = mp.Vector3(16, 8, 0)

# 几何体
geometry = [mp.Block(
    mp.Vector3(1, 1e20, 1e20),
    center=mp.Vector3(-5, 0),
    material=mp.Medium(epsilon=12)
)]

# 光源
sources = [mp.Source(
    mp.GaussianSource(frequency=0.4, fwidth=0.1),
    component=mp.Ey,
    center=mp.Vector3(-7, 0)
)]

# PML
pml_layers = [mp.PML(1.0)]

# 运行
sim = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    sources=sources,
    boundary_layers=pml_layers,
    resolution=20
)

sim.run(until=200)
```

## 资源

- [MEEP 官方文档](https://meep.readthedocs.io/)
- [MEEP GitHub](https://github.com/NanoComp/meep)
- [FDTD 入门（英文）](https://en.wikipedia.org/wiki/Finite-difference_time-domain_method)
- [计算电磁学教程 (中文)](https://www.bilibili.com/video/BV1KJ411W7zp)

## License

MIT
