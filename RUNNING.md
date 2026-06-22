# 运行说明

## 推荐环境

MEEP 的 Python 包在 conda-forge 维护得最好，推荐使用 `environment.yml`：

```bash
conda env create -f environment.yml
conda activate meep-learning
```

如果使用已有环境，也可参考：

```bash
conda install -c conda-forge pymeep pymeep-extras numpy scipy matplotlib h5py
```

部分进阶工具会额外用到 `pandas`（参数扫描表格）和 `autograd`（拓扑优化示例）。如果不是用本仓库的 `environment.yml` 创建环境，请一并安装：

```bash
conda install -c conda-forge pandas autograd
```

## 运行教程

所有脚本均为自包含示例，可以从仓库根目录运行：

```bash
python tutorial/basic/01_hello_meep.py
python examples/near_to_far_field.py
python projects/wdm_demux.py
```


## 新增学习材料

```bash
# 阅读学习路线和 API 速查表
less docs/learning_map.md
less cheatsheets/meep_python_api.md

# 运行新增示例
python tutorial/basic/11_symmetry_and_monitors.py
python tutorial/intermediate/06_eigenmode_decomposition.py
python tutorial/advanced/06_harminv_cavity_q.py
python examples/flux_normalization_reflection.py
python examples/mode_decomposition_ports.py
```

新增脚本会在当前目录保存对应 PNG 图像。部分模式分解/Harminv 示例计算量比基础 1D 示例更大，建议先保持默认分辨率，确认流程后再提高精度。

## 输出文件

仿真可能生成 HDF5、图片、CSV、视频等结果文件，这些已在 `.gitignore` 中排除，避免污染仓库。

## 快速检查

```bash
python -m py_compile $(find tutorial examples projects -name '*.py')
```

这个检查只验证 Python 语法，不会启动完整 FDTD 仿真；适合在提交前快速确认脚本没有语法错误。
