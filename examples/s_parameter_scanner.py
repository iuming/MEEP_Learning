"""
Examples: S 参数扫描工具

自动扫描几何参数并计算 S 参数。
用于器件优化和参数扫描。
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt
from itertools import product


class SParameterScanner:
    """
    S 参数扫描器
    支持多参数扫描：几何尺寸、材料、波长等
    """

    def __init__(self, param_ranges, resolution=30, dpml=1.0):
        """
        param_ranges: dict, e.g. {'width': [0.3, 0.5, 0.7], 'length': [5, 10, 15]}
        """
        self.param_ranges = param_ranges
        self.resolution = resolution
        self.dpml = dpml
        self.results = []

    def create_simulation(self, params):
        """
        需要子类重写此方法：根据参数返回 mp.Simulation 对象。
        示例实现见下方 run()。
        """
        raise NotImplementedError("请重写 create_simulation 方法")

    def add_monitors(self, sim):
        """
        需要子类重写：添加 S 参数监视器。
        返回 dict: {'S11': flux_obj, 'S21': flux_obj, ...}
        """
        raise NotImplementedError("请重写 add_monitors 方法")

    def run_scan(self, verbose=True):
        """
        执行完整参数扫描。
        """
        keys = list(self.param_ranges.keys())
        values_list = list(self.param_ranges.values())

        total = np.prod([len(v) for v in values_list])

        idx = 0
        for combo in product(*values_list):
            params = dict(zip(keys, combo))

            if verbose:
                print(f"[{idx + 1}/{total}] {params}")

            try:
                sim = self.create_simulation(params)
                monitors = self.add_monitors(sim)
                sim.run(until=300)

                # 提取 S 参数
                result = {'params': params}
                for name, flux in monitors.items():
                    result[name] = np.array(mp.get_fluxes(flux))[0]

                self.results.append(result)

            except Exception as e:
                print(f"  ⚠ 扫描失败: {e}")

            idx += 1

        if verbose:
            print(f"\n扫描完成，共 {len(self.results)} 个有效结果")

        return self.results

    def to_dataframe(self):
        """转换为 pandas DataFrame"""
        try:
            import pandas as pd
            records = []
            for r in self.results:
                record = {**r['params']}
                for k, v in r.items():
                    if k != 'params':
                        record[k] = v
                records.append(record)
            return pd.DataFrame(records)
        except ImportError:
            print("需要安装 pandas: pip install pandas")
            return self.results

    def plot_2d(self, x_param, y_param, z_key='S21',
                xlabel=None, ylabel=None, title=None, save_path=None):
        """
        绘制 2D 参数扫描热力图。
        """
        if xlabel is None:
            xlabel = x_param
        if ylabel is None:
            ylabel = y_param

        # 构建网格
        x_vals = self.param_ranges[x_param]
        y_vals = self.param_ranges[y_param]

        Z = np.zeros((len(y_vals), len(x_vals)))

        for r in self.results:
            xi = list(x_vals).index(r['params'][x_param])
            yi = list(y_vals).index(r['params'][y_param])
            Z[yi, xi] = r[z_key]

        plt.figure(figsize=(8, 6))
        plt.imshow(Z, aspect='auto', origin='lower', cmap='viridis',
                   extent=[x_vals[0], x_vals[-1], y_vals[0], y_vals[-1]])
        plt.colorbar(label=z_key)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if title:
            plt.title(title)
        else:
            plt.title(f'{z_key} vs {x_param} and {y_param}')

        if save_path:
            plt.savefig(save_path, dpi=150)
        plt.show()
        return Z


# ============================================================
# 演示：介质平板厚度扫描
# ============================================================
class DielectricSlabScanner(SParameterScanner):
    """
    介质平板反射率/透射率 vs 厚度
    """

    def __init__(self, thicknesses, n_slab=3.0, freq=0.5):
        super().__init__({'thickness': thicknesses})
        self.n_slab = n_slab
        self.freq = freq
        self.fwidth = 0.2 * freq

    def create_simulation(self, params):
        t = params['thickness']
        cell = mp.Vector3(2 * self.dpml + 2 + t + 2, 0, 0)

        geometry = [mp.Block(
            size=mp.Vector3(t, mp.inf, mp.inf),
            center=mp.Vector3(0, 0, 0),
            material=mp.Medium(epsilon=self.n_slab ** 2),
        )]

        sources = [mp.Source(
            mp.GaussianSource(frequency=self.freq, fwidth=self.fwidth),
            component=mp.Ey,
            center=mp.Vector3(-self.dpml - 0.5, 0, 0),
        )]

        return mp.Simulation(
            cell_size=cell,
            geometry=geometry,
            sources=sources,
            boundary_layers=[mp.PML(self.dpml)],
            resolution=self.resolution,
        )

    def add_monitors(self, sim):
        refl_fr = mp.FluxRegion(
            center=mp.Vector3(-self.dpml - 0.2, 0, 0),
            size=mp.Vector3(0, mp.inf, mp.inf),
        )
        tran_fr = mp.FluxRegion(
            center=mp.Vector3(self.dpml + 0.2, 0, 0),
            size=mp.Vector3(0, mp.inf, mp.inf),
        )
        return {
            'R': sim.add_flux(self.freq, self.fwidth, 1, refl_fr),
            'T': sim.add_flux(self.freq, self.fwidth, 1, tran_fr),
        }


if __name__ == '__main__':
    wavelengths = np.arange(0.2, 2.0, 0.1)
    scanner = DielectricSlabScanner(thicknesses=wavelengths)
    results = scanner.run_scan()

    thicknesses = [r['params']['thickness'] for r in results]
    T_values = [r['T'] for r in results]
    R_values = [r['R'] for r in results]

    plt.figure(figsize=(10, 5))
    plt.plot(thicknesses, T_values, 'o-', label='Transmission', markersize=4)
    plt.plot(thicknesses, R_values, 's-', label='Reflection', markersize=4)
    plt.xlabel('Slab Thickness (μm)')
    plt.ylabel('Normalized Power')
    plt.title('Fabry-Pérot Resonances')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig('sparameter_scan_demo.png', dpi=150)
    plt.show()

    print(f"扫描了 {len(results)} 个厚度值")
    print(f"最大透射: {max(T_values):.3f} @ thickness = {thicknesses[np.argmax(T_values)]:.1f} μm")
