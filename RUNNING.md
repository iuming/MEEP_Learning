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

## 运行教程

所有脚本均为自包含示例，可以从仓库根目录运行：

```bash
python tutorial/basic/01_hello_meep.py
python examples/near_to_far_field.py
python projects/wdm_demux.py
```

## 输出文件

仿真可能生成 HDF5、图片、CSV、视频等结果文件，这些已在 `.gitignore` 中排除，避免污染仓库。

## 快速检查

```bash
python -m py_compile $(find tutorial examples projects -name '*.py')
```
