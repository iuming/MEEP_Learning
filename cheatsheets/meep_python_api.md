# MEEP Python API 速查表

本速查表面向本仓库脚本，强调最常用对象和容易踩坑的参数。

## 1. 最小仿真骨架

```python
import meep as mp

cell = mp.Vector3(16, 8, 0)       # z=0 表示 2D
pml_layers = [mp.PML(1.0)]        # PML 厚度
geometry = []                     # 几何体列表
sources = [mp.Source(
    mp.GaussianSource(frequency=0.5, fwidth=0.2),
    component=mp.Ez,
    center=mp.Vector3(-5, 0),
)]

sim = mp.Simulation(
    cell_size=cell,
    boundary_layers=pml_layers,
    geometry=geometry,
    sources=sources,
    resolution=20,
)
sim.run(until=100)
```

## 2. 常用几何体

```python
# 长方体/矩形
mp.Block(size=mp.Vector3(4, 1, mp.inf), center=mp.Vector3(), material=mp.Medium(epsilon=12))

# 圆柱/圆盘；2D 中是圆
mp.Cylinder(radius=1.0, height=mp.inf, center=mp.Vector3(), material=mp.Medium(index=3.45))

# 球体；3D 才常用
mp.Sphere(radius=0.5, center=mp.Vector3(), material=mp.Medium(epsilon=2.25))

# 棱柱；vertices 使用 xy 平面坐标
mp.Prism(vertices=[mp.Vector3(-1,0), mp.Vector3(1,0), mp.Vector3(0,1)], height=mp.inf, material=mp.Medium(index=2))
```

提示：`mp.inf` 常用于 2D 几何的无限厚度；不要把结构直接贴到 PML 上。

## 3. 光源选择

| 光源 | 用途 | 关键参数 |
|---|---|---|
| `mp.ContinuousSource` | 单频稳态场图 | `frequency`, `width` |
| `mp.GaussianSource` | 宽频谱、透射/反射谱 | `frequency`, `fwidth` |
| `mp.EigenModeSource` | 注入波导本征模式 | `eig_band`, `eig_parity`, `size` |


示例：

```python
mp.EigenModeSource(
    src=mp.GaussianSource(frequency=0.645, fwidth=0.08),
    center=mp.Vector3(-6, 0),
    size=mp.Vector3(0, 3),          # 覆盖整个波导截面
    direction=mp.X,
    eig_band=1,
    eig_parity=mp.ODD_Z,            # 2D Ez 常用 ODD_Z；需按模型确认
)
```

## 4. 监视器与数据提取

### 点场记录

```python
values = []
def record(sim):
    values.append(sim.get_field_point(mp.Ez, mp.Vector3(0, 0)))
sim.run(mp.at_every(1.0, record), until=200)
```

### 通量监视器

```python
fr = mp.FluxRegion(center=mp.Vector3(5, 0), size=mp.Vector3(0, 4))
flux = sim.add_flux(fcen, df, nfreq, fr)
sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(5,0), 1e-8))
spectrum = mp.get_fluxes(flux)
freqs = mp.get_flux_freqs(flux)
```

### 空结构归一化

```python
# 第一次：空结构
incident_flux_data = sim.get_flux_data(refl_flux)
incident = mp.get_fluxes(tran_flux)

# 第二次：有结构，反射监视器要减去入射场
sim.load_minus_flux_data(refl_flux, incident_flux_data)
```

## 5. 对称性

```python
symmetries = [mp.Mirror(mp.Y, phase=+1)]  # 偶对称
sim = mp.Simulation(..., symmetries=symmetries)
```

- `phase=+1`：镜像后场分量同号。
- `phase=-1`：镜像后场分量反号。
- 对称性只在结构、源、边界、监视器都满足时可用；不确定时先不用。

## 6. Harminv 共振分析

```python
h = mp.Harminv(mp.Ez, mp.Vector3(0, 0), fcen, df)
sim.run(mp.after_sources(h), until_after_sources=300)
for mode in h.modes:
    print(mode.freq, mode.Q, abs(mode.amp))
```

适合高 Q 腔、环形谐振器、缺陷腔。仿真时间越长，频率和 Q 越可靠。

## 7. 常见报错定位

- `fields do not decay`：仿真时间不够、高 Q 共振、PML 反射或源仍在激发。
- 结果随 `resolution` 大幅变化：网格太粗或几何台阶误差明显。
- 反射率为负/大于 1：通量方向、归一化或监视器位置有问题。
- 模式源激发失败：`eig_parity`、截面尺寸或频率不支持目标模式。
