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
├── tutorial/              # 23 个教程（按难度递进）
│   ├── basic/             # 11 基础篇：光源、波导、光栅、能带、对称性
│   ├── intermediate/      #  6 进阶篇：MMI、光子晶体、等离子体、非线性、模式分解
│   └── advanced/          #  6 高级篇：拓扑优化、腔QED、混沌、光镊、超连续谱、Harminv
├── examples/              #  5 实用工具与框架
├── docs/                  # 学习路线与概念说明
├── cheatsheets/           # API 速查表
├── projects/              #  2 综合项目
├── RUNNING.md              # 运行说明与环境检查
├── environment.yml         # Conda 环境
├── requirements.txt        # pip 依赖参考
└── README.md
```

## 安装

```bash
# Conda（推荐）
conda env create -f environment.yml
conda activate meep-learning

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
| 11 | `symmetry_and_monitors` | 对称性、点监视器与二维场切片 |

### 🟡 进阶篇 `tutorial/intermediate/`

| 编号 | 文件 | 内容 |
|------|------|------|
| 01 | `mmi_coupler` | 1×2 多模干涉耦合器设计 |
| 02 | `photonic_crystal_wg` | 光子晶体波导与慢光效应 |
| 03 | `plasmonic_antenna` | 金纳米棒 LSPR 仿真（Drude 色散模型） |
| 04 | `directional_coupler` | 定向耦合器与模式拍长 |
| 05 | `nonlinear_kerr` | Kerr 非线性（SPM 光谱展宽） |
| 06 | `eigenmode_decomposition` | 本征模源与波导端口模式分解 |

### 🔴 高级篇 `tutorial/advanced/`

| 编号 | 文件 | 内容 |
|------|------|------|
| 01 | `topology_optimization` | 伴随法拓扑优化 → 逆向设计 1×2 分束器 |
| 02 | `purcell_effect` | 腔 QED — 量子发射体 Purcell 因子与自发辐射增强 |
| 03 | `chaotic_microcavity` | 变形微盘腔 WG 模式与混沌射线动力学 |
| 04 | `optical_tweezers` | 光力效应 — Maxwell 应力张量光阱力仿真 |
| 05 | `supercontinuum` | PCF 超连续谱 — 非线性脉冲传播与光谱展宽 |
| 06 | `harminv_cavity_q` | Harminv 提取谐振腔频率与 Q 值 |

### 🛠 实用工具 `examples/`

| 文件 | 内容 |
|------|------|
| `visualization_tools.py` | 2D/矢量场/动画可视化函数库 |
| `s_parameter_scanner.py` | 自动化 S 参数扫描框架 |
| `near_to_far_field.py` | 天线远场辐射方向图计算 |
| `flux_normalization_reflection.py` | 空结构归一化与反射/透射谱模板 |
| `mode_decomposition_ports.py` | 二端口波导模式分解 S 参数模板 |

### 🚀 综合项目 `projects/`

| 项目 | 内容 |
|------|------|
| `wdm_demux.py` | 4 通道 DWDM 波分复用器（级联微环） |
| `metalens.py` | 介质超表面透镜（TiO2 纳米柱） |


## Star History

[![Star History Chart](https://api.star-history.com/svg?repos=iuming/MEEP_Learning&type=Date)](https://www.star-history.com/#iuming/MEEP_Learning&Date)

## 📊 内容总览

| 难度 | 数量 | 涵盖主题 |
|------|------|----------|
| 🟢 基础 | 11 | 光源、边界、波导、光栅、能带、谐振腔、监视器 |
| 🟡 进阶 | 6 | MMI、光子晶体、等离激元、定向耦合、Kerr 非线性、模式分解 |
| 🔴 高级 | 6 | 拓扑优化、腔 QED、混沌腔、光镊、超连续谱、Harminv |
| 🛠 工具 | 5 | 可视化库、S 参数扫描、近远场变换、通量归一化、端口模板 |
| 🚀 项目 | 2 | DWDM 波分复用、超表面透镜 |
| 📚 文档 | 2 | 学习地图、Python API 速查表 |
| **合计** | **30** | |

## 快速开始

```bash
# 克隆仓库
git clone git@github.com:iuming/MEEP_Learning.git
cd MEEP_Learning

# 运行第一个示例
python tutorial/basic/01_hello_meep.py

# 更多运行细节见 RUNNING.md
```

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

# 运行
sim = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    sources=sources,
    boundary_layers=[mp.PML(1.0)],
    resolution=20
)
sim.run(until=200)
```

## 🔗 学习路径建议

```
基础 01-05  →  掌握 MEEP 基本工作流
基础 06-11  →  片上集成光子学核心器件与监视器技巧
进阶 01-06  →  现代纳米光子学设计与端口模式分析
工具        →  建立自己的仿真工具箱
高级 01-06  →  前沿研究方向与方法
项目        →  独立完成完整器件设计
```

## 贡献

欢迎提交 PR！新增教程请参考已有文件的代码风格：
- 中文注释说明关键步骤
- 自包含、可直接运行
- 包含结果分析与物理讨论

## 资源

- [MEEP 官方文档](https://meep.readthedocs.io/)
- [MEEP GitHub](https://github.com/NanoComp/meep)
- [FDTD 入门](https://en.wikipedia.org/wiki/Finite-difference_time-domain_method)

## License

MIT
