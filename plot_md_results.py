#!/usr/bin/env python3
"""
绘制MD模拟结果图表
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# 设置绘图风格
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (15, 10)
plt.rcParams['font.size'] = 12

def read_xvg(filename):
    """读取GROMACS xvg文件"""
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    data.append([float(parts[0]), float(parts[1])])
                except ValueError:
                    continue
    return np.array(data)

# 创建输出目录
output_dir = Path("/home/tony/colabfold/results/figures")
output_dir.mkdir(exist_ok=True)

# 分析野生型和突变体
proteins = {
    'WT': '/home/tony/colabfold/results/trajectories/p53_DNA_binding_domain_WT',
    'R273H': '/home/tony/colabfold/results/trajectories/p53_DNA_binding_domain_R273H'
}

# 创建子图
fig, axes = plt.subplots(2, 2, figsize=(15, 12))

# 1. 温度
ax = axes[0, 0]
for name, path in proteins.items():
    energy_file = Path(path) / "energy.xvg"
    if energy_file.exists():
        data = read_xvg(energy_file)
        if len(data) > 0 and data.shape[1] >= 3:  # 时间, 温度
            ax.plot(data[:, 0], data[:, 2], label=f'{name}', linewidth=2)
ax.set_xlabel('Time (ps)')
ax.set_ylabel('Temperature (K)')
ax.set_title('Temperature During Production MD')
ax.legend()
ax.axhline(y=310, color='r', linestyle='--', label='Target (310 K)')

# 2. 势能
ax = axes[0, 1]
for name, path in proteins.items():
    energy_file = Path(path) / "energy.xvg"
    if energy_file.exists():
        data = read_xvg(energy_file)
        if len(data) > 0 and data.shape[1] >= 2:
            ax.plot(data[:, 0], data[:, 1], label=f'{name}', linewidth=2)
ax.set_xlabel('Time (ps)')
ax.set_ylabel('Potential Energy (kJ/mol)')
ax.set_title('Potential Energy During Production MD')
ax.legend()

# 3. RMSD (NPT阶段)
ax = axes[1, 0]
for name, path in proteins.items():
    rmsd_file = Path(path) / "rmsd_npt.xvg"
    if rmsd_file.exists():
        data = read_xvg(rmsd_file)
        if len(data) > 0:
            ax.plot(data[:, 0], data[:, 1], label=f'{name}', linewidth=2)
ax.set_xlabel('Time (ps)')
ax.set_ylabel('RMSD (nm)')
ax.set_title('RMSD During NPT Equilibration')
ax.legend()

# 4. 回转半径 (NPT阶段)
ax = axes[1, 1]
for name, path in proteins.items():
    rg_file = Path(path) / "rg_npt.xvg"
    if rg_file.exists():
        data = read_xvg(rg_file)
        if len(data) > 0:
            ax.plot(data[:, 0], data[:, 1], label=f'{name}', linewidth=2)
ax.set_xlabel('Time (ps)')
ax.set_ylabel('Radius of Gyration (nm)')
ax.set_title('Radius of Gyration During NPT Equilibration')
ax.legend()

plt.tight_layout()
plt.savefig(output_dir / 'md_analysis_comparison.png', dpi=300, bbox_inches='tight')
print(f"✅ 图表已保存: {output_dir / 'md_analysis_comparison.png'}")

# 生成第二个图：密度和压力
fig, axes = plt.subplots(1, 2, figsize=(15, 5))

# 密度
ax = axes[0]
for name, path in proteins.items():
    energy_file = Path(path) / "energy.xvg"
    if energy_file.exists():
        data = read_xvg(energy_file)
        if len(data) > 0 and data.shape[1] >= 6:
            ax.plot(data[:, 0], data[:, 5], label=f'{name}', linewidth=2)
ax.set_xlabel('Time (ps)')
ax.set_ylabel('Density (kg/m³)')
ax.set_title('Density During Production MD')
ax.legend()
ax.axhline(y=1000, color='r', linestyle='--', label='Water (1000 kg/m³)')

# 压力
ax = axes[1]
for name, path in proteins.items():
    energy_file = Path(path) / "energy.xvg"
    if energy_file.exists():
        data = read_xvg(energy_file)
        if len(data) > 0 and data.shape[1] >= 4:
            ax.plot(data[:, 0], data[:, 3], label=f'{name}', linewidth=2)
ax.set_xlabel('Time (ps)')
ax.set_ylabel('Pressure (bar)')
ax.set_title('Pressure During Production MD')
ax.legend()
ax.axhline(y=1.0, color='r', linestyle='--', label='Target (1 bar)')

plt.tight_layout()
plt.savefig(output_dir / 'md_density_pressure.png', dpi=300, bbox_inches='tight')
print(f"✅ 图表已保存: {output_dir / 'md_density_pressure.png'}")

print("\n✅ 所有图表生成完成！")
