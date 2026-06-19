# MEEP 学习地图：从入门到高级

本文把仓库中的脚本串成一条可执行的学习路线，并补充每一阶段应掌握的 MEEP 功能。所有示例默认从仓库根目录运行。

## 0. 先理解 MEEP 的基本单位

MEEP 使用无量纲单位：通常把长度单位设为 `1 μm`，则频率单位为 `1/μm`，真空光速 `c=1`。

- 波长：`wavelength = 1 / frequency`
- 折射率：`n = sqrt(epsilon)`（非磁性材料）
- 时间：光在 `1 μm` 中传播所需时间为 1 个 MEEP 时间单位
- 分辨率：`resolution=20` 表示每个长度单位 20 个网格点

经验规则：每个介质中的最短波长至少 8–12 个网格点：

```text
resolution >= 10 * n_max / wavelength_min
```

## 1. 入门：会搭一个能跑的仿真

建议顺序：

1. `tutorial/basic/01_hello_meep.py`：仿真区域、光源、PML、时间推进
2. `tutorial/basic/02_source_types.py`：连续源与脉冲源的差异
3. `tutorial/basic/05_pml.py`：为什么边界会反射，PML 厚度如何选
4. `tutorial/basic/11_symmetry_and_monitors.py`：对称性、切片输出、点监视器

你应该能回答：

- `cell_size`、`resolution`、`boundary_layers` 分别控制什么？
- 为什么仿真区域不能无限大？PML 的作用是什么？
- 什么时候可以用 `mp.Mirror` 对称性加速？

## 2. 基础器件：从场图到透射率

建议顺序：

1. `tutorial/basic/03_dielectric_slab.py`：通量监视器、反射率/透射率
2. `examples/flux_normalization_reflection.py`：入射场归一化、减去背景场
3. `tutorial/basic/07_waveguide.py`：波导模式与场分布
4. `tutorial/basic/08_ring_resonator.py`：频谱、共振、品质因子直觉

关键习惯：

- 先跑“空结构”得到入射功率，再跑“有结构”得到散射功率。
- 监视器离结构太近会采到近场振荡；通常要留出一段传播距离。
- 对宽频结果做收敛检查：提高 `resolution`、增大 PML、延长 `until_after_sources`。

## 3. 进阶：模式、S 参数与远场

建议顺序：

1. `tutorial/intermediate/06_eigenmode_decomposition.py`：本征模源、模式分解、端口功率
2. `examples/mode_decomposition_ports.py`：二端口波导的模式 S 参数模板
3. `examples/s_parameter_scanner.py`：参数扫描框架
4. `examples/near_to_far_field.py`：近远场变换与辐射方向图

这一阶段要把“看场图”升级为“可复用测量流程”：

- 端口定义：位置、尺寸、频率范围、模式编号
- 归一化：入射模式功率为 1，输出为相对透射/反射
- 扫描：把单次仿真封装成函数或类，保存 CSV/图像

## 4. 高级：谐振、色散、非线性、优化

建议顺序：

1. `tutorial/advanced/06_harminv_cavity_q.py`：Harminv 提取共振频率和 Q 值
2. `tutorial/intermediate/03_plasmonic_antenna.py`：Drude/Lorentz 色散材料
3. `tutorial/intermediate/05_nonlinear_kerr.py`：Kerr 非线性
4. `tutorial/advanced/01_topology_optimization.py`：逆向设计思想

高级仿真的共同风险：

- 物理时间不够长会导致频谱/Q 值不准。
- 高 Q 腔需要更细网格、更小 Courant 因子或更严格收敛。
- 色散/非线性参数必须确认单位和稳定性；先用低分辨率做 sanity check。

## 5. 每个新模型的检查清单

1. **几何检查**：画出介电常数分布或场图，确认结构位置正确。
2. **边界检查**：PML 与结构/源/监视器保持足够距离。
3. **源检查**：中心频率、带宽、偏振分量符合目标模式。
4. **归一化检查**：有无空结构参考仿真？功率符号是否正确？
5. **收敛检查**：至少改变一次 `resolution` 或 PML 厚度比较结果。
6. **保存结果**：脚本输出图片/CSV，并在文件名里写清楚含义。
