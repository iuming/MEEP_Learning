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
| 03 | `dielectric_slab` | 电介质平板传输与反射 |
| 04 | `bragg_grating` | 布拉格光栅与带隙 |
| 05 | `pml_absorption` | PML 吸收边界条件 |
| 06 | `band_structure` | 光子晶体能带结构计算 |
| 07 | `waveguide` | 介质波导模式分析 |
| 08 | `ring_resonator` | 环形谐振腔 |
| 09 | `field_monitors` | 场监视器与可视化 |
| 10 | `s_parameters` | S 参数计算 |

### 🟡 进阶篇 `tutorial/intermediate/`

> 即将添加：非线性材料、各向异性介质、片上器件、光栅耦合器

### 🔴 高级篇 `tutorial/advanced/`

> 即将添加：拓扑优化、逆设计、量子发射体、近场-远场变换

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
