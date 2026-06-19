"""
Example: 通量归一化与反射率计算模板

很多新手第一次算反射率时会忘记“空结构参考仿真”。正确流程通常是：
1. 空结构：记录入射到反射监视器的场数据，以及透射监视器的入射功率；
2. 有结构：加载负的入射场数据，只保留真正的反射散射场；
3. 反射率 R = -反射通量 / 入射通量，透射率 T = 透射通量 / 入射通量。

本例为 1D 介质平板，速度快，适合作为所有频域通量计算的模板。
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

resolution = 40
dpml = 1.0
slab_thickness = 1.0
slab_index = 3.0
fcen = 0.5
fwidth = 0.4
nfreq = 101
sx = 12

refl_pt = mp.Vector3(-2.5)
tran_pt = mp.Vector3(2.5)
src_pt = mp.Vector3(-4.0)


def make_sim(with_slab):
    """根据开关创建空结构或介质平板仿真。"""
    geometry = []
    if with_slab:
        geometry = [mp.Block(
            size=mp.Vector3(slab_thickness, mp.inf, mp.inf),
            center=mp.Vector3(),
            material=mp.Medium(index=slab_index),
        )]

    sources = [mp.Source(
        src=mp.GaussianSource(frequency=fcen, fwidth=fwidth),
        component=mp.Ey,
        center=src_pt,
    )]

    return mp.Simulation(
        cell_size=mp.Vector3(sx, 0, 0),
        geometry=geometry,
        sources=sources,
        boundary_layers=[mp.PML(dpml)],
        resolution=resolution,
    )

# ------------------------- 1. 空结构参考 -------------------------
sim_ref = make_sim(with_slab=False)
refl_ref = sim_ref.add_flux(fcen, fwidth, nfreq, mp.FluxRegion(center=refl_pt))
tran_ref = sim_ref.add_flux(fcen, fwidth, nfreq, mp.FluxRegion(center=tran_pt))
sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ey, tran_pt, 1e-9))

incident_flux = np.array(mp.get_fluxes(tran_ref))
incident_refl_data = sim_ref.get_flux_data(refl_ref)
freqs = np.array(mp.get_flux_freqs(tran_ref))

# ------------------------- 2. 有结构仿真 -------------------------
sim = make_sim(with_slab=True)
refl = sim.add_flux(fcen, fwidth, nfreq, mp.FluxRegion(center=refl_pt))
tran = sim.add_flux(fcen, fwidth, nfreq, mp.FluxRegion(center=tran_pt))

# 关键步骤：在反射监视器中减去空结构的入射场，只留下向左传播的散射场。
sim.load_minus_flux_data(refl, incident_refl_data)
sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ey, tran_pt, 1e-9))

reflected_flux = np.array(mp.get_fluxes(refl))
transmitted_flux = np.array(mp.get_fluxes(tran))

R = -reflected_flux / incident_flux
T = transmitted_flux / incident_flux
wavelengths = 1 / freqs

print(f"能量守恒检查：R+T 的平均值 = {np.mean(R + T):.4f}")

plt.figure(figsize=(7, 4))
plt.plot(wavelengths, R, label="反射率 R")
plt.plot(wavelengths, T, label="透射率 T")
plt.plot(wavelengths, R + T, "k--", label="R+T")
plt.gca().invert_xaxis()
plt.xlabel("波长 / μm（若长度单位为 μm）")
plt.ylabel("功率比例")
plt.title("介质平板的通量归一化反射/透射谱")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("flux_normalization_reflection.png", dpi=160)
print("图像已保存到 flux_normalization_reflection.png")
plt.show()
