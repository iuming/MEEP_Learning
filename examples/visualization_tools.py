"""
Examples: 二维场可视化工具

通用函数集合，用于可视化 MEEP 仿真结果：
- 场分布 2D 彩色图
- 矢量场图 (quiver)
- 动画生成 (时间演化)
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def plot_field_2d(sim, component=mp.Ez, cmap='RdBu',
                  title='Field Distribution', save_path=None):
    """
    绘制 2D 场分布的彩色图。
    """
    field = sim.get_array(component=component)
    size = sim.cell_size

    plt.figure(figsize=(8, 6))
    plt.imshow(field.transpose(), interpolation='spline36', cmap=cmap,
               extent=[-size.x/2, size.x/2, -size.y/2, size.y/2],
               origin='lower')
    plt.colorbar(label=f'|{component.component}|')
    plt.title(title)
    plt.xlabel('x (μm)')
    plt.ylabel('y (μm)')

    if save_path:
        plt.savefig(save_path, dpi=150)
    plt.show()


def plot_structure(sim, cmap='binary', title='Dielectric Structure',
                   save_path=None):
    """
    绘制介质结构图。
    """
    eps = sim.get_array(component=mp.Dielectric)
    size = sim.cell_size

    plt.figure(figsize=(8, 6))
    plt.imshow(eps.transpose(), interpolation='spline36', cmap=cmap,
               extent=[-size.x/2, size.x/2, -size.y/2, size.y/2],
               origin='lower')
    plt.title(title)
    plt.xlabel('x (μm)')
    plt.ylabel('y (μm)')

    if save_path:
        plt.savefig(save_path, dpi=150)
    plt.show()


def plot_vector_field(sim, scale=1.0, title='Vector Field',
                      save_path=None):
    """
    绘制 Ex, Ey 矢量场。
    """
    size = sim.cell_size
    shape = (int(size.x * sim.resolution), int(size.y * sim.resolution))

    ex = sim.get_array(component=mp.Ex)
    ey = sim.get_array(component=mp.Ey)

    x = np.linspace(-size.x/2, size.x/2, ex.shape[0])
    y = np.linspace(-size.y/2, size.y/2, ex.shape[1])
    X, Y = np.meshgrid(x, y, indexing='ij')

    magnitude = np.sqrt(ex**2 + ey**2)

    plt.figure(figsize=(8, 6))
    plt.streamplot(X, Y, ex.T, ey.T, color=magnitude.T,
                   cmap='plasma', density=2, linewidth=1)
    plt.colorbar(label='|E|')
    plt.title(title)
    plt.xlabel('x (μm)')
    plt.ylabel('y (μm)')

    if save_path:
        plt.savefig(save_path, dpi=150)
    plt.show()


def animate_field(sim, component=mp.Ez, nframes=100,
                  filename='field_animation.gif', dt=0.5):
    """
    生成场动画 GIF。
    注意：需要先运行仿真并在每个时间步保存场数据。
    """
    frames = []

    def record_frame(sim):
        frames.append(sim.get_array(component=component))

    # 运行并记录
    for i in range(nframes):
        sim.run(until=i * dt)
        record_frame(sim)

    size = sim.cell_size
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(frames[0].transpose(), cmap='RdBu',
                   extent=[-size.x/2, size.x/2, -size.y/2, size.y/2],
                   origin='lower', animated=True)
    plt.colorbar(im, ax=ax)
    ax.set_xlabel('x (μm)')
    ax.set_ylabel('y (μm)')

    def update(frame_idx):
        im.set_array(frames[frame_idx].transpose())
        ax.set_title(f'{component.component} @ t = {frame_idx * dt:.1f}')
        return [im]

    anim = FuncAnimation(fig, update, frames=len(frames),
                         interval=100, blit=True)
    anim.save(filename, writer='pillow', fps=10)
    plt.close()
    print(f"Animation saved to {filename}")

    return frames


# ============================================================
# 快速演示
# ============================================================
if __name__ == '__main__':
    # 创建一个简单仿真来演示可视化工具
    sim = mp.Simulation(
        cell_size=mp.Vector3(8, 8, 0),
        geometry=[mp.Cylinder(radius=2.0,
                               material=mp.Medium(epsilon=12))],
        sources=[mp.Source(
            mp.GaussianSource(frequency=0.3, fwidth=0.1),
            component=mp.Ez,
            center=mp.Vector3(-3, 0, 0),
        )],
        resolution=20,
        boundary_layers=[mp.PML(1.0)],
    )
    sim.run(until=100)

    plot_structure(sim, save_path='example_structure.png')
    plot_field_2d(sim, component=mp.Ez, save_path='example_field.png')
    print("可视化示例已保存。")
