#!/usr/bin/env python3
"""
分子动力学模拟模块 - 使用GROMACS进行MD模拟
作者: iFlow CLI
目标: p53 DNA结合域野生型与R273H突变体的MD模拟
"""

import os
import sys
import argparse
import subprocess
import yaml
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.PDB import PDBParser
import MDAnalysis as mda
from MDAnalysis.analysis import rms

class MDSimulator:
    def __init__(self, config_path="config/protein_config.yaml"):
        """初始化MD模拟器"""
        self.config = self._load_config(config_path)
        self.setup_directories()
        self.check_gromacs()
        
    def _load_config(self, config_path):
        """加载配置文件"""
        with open(config_path, 'r', encoding='utf-8') as f:
            return yaml.safe_load(f)
    
    def setup_directories(self):
        """创建输出目录"""
        base_dir = Path(self.config['output']['base_dir'])
        self.trajectory_dir = base_dir / self.config['output']['trajectory_dir']
        self.analysis_dir = base_dir / self.config['output']['analysis_dir']
        
        self.trajectory_dir.mkdir(parents=True, exist_ok=True)
        self.analysis_dir.mkdir(parents=True, exist_ok=True)
    
    def check_gromacs(self):
        """检查GROMACS是否可用"""
        try:
            result = subprocess.run(['gmx', '-version'],
                                  capture_output=True, text=True, check=True)
            print("GROMACS可用")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError, PermissionError):
            print("警告: GROMACS未安装或无法访问")
            print("如需完整MD模拟功能，请安装GROMACS: conda install -c bioconda gromacs")
            return False
    
    def run_gmx_command(self, cmd, input_data=None, cwd=None):
        """执行GROMACS命令"""
        try:
            process = subprocess.Popen(cmd, stdin=subprocess.PIPE, 
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                     text=True, cwd=cwd)
            
            stdout, stderr = process.communicate(input=input_data)
            
            if process.returncode != 0:
                print(f"GROMACS命令失败: {' '.join(cmd)}")
                print(f"错误信息: {stderr}")
                return False
            
            return True
        except Exception as e:
            print(f"执行GROMACS命令时出错: {e}")
            return False
    
    def prepare_topology(self, pdb_file, output_dir):
        """准备拓扑结构"""
        print(f"准备拓扑结构: {pdb_file}")

        work_dir = Path(output_dir).resolve()
        work_dir.mkdir(exist_ok=True)

        # 转换为绝对路径
        pdb_file_abs = Path(pdb_file).resolve()

        # 1. pdb2gmx - 生成拓扑文件
        cmd = ['gmx', 'pdb2gmx',
               '-f', str(pdb_file_abs),
               '-o', str(work_dir / 'processed.gro'),
               '-p', str(work_dir / 'topol.top'),
               '-ff', self.config['md_simulation']['force_field'],
               '-water', self.config['md_simulation']['water_model']]

        # 自动选择力场参数
        input_data = "1\n"  # 选择第一个力场选项

        if not self.run_gmx_command(cmd, input_data=input_data):
            return None

        return work_dir / 'processed.gro'
    
    def define_box_and_solvate(self, gro_file, work_dir):
        """定义盒子并溶剂化"""
        print("定义模拟盒子和溶剂化...")

        # 转换为绝对路径
        gro_file_abs = Path(gro_file).resolve()

        # 1. 定义盒子
        cmd = ['gmx', 'editconf',
               '-f', str(gro_file_abs),
               '-o', str(work_dir / 'boxed.gro'),
               '-c', '-d', str(self.config['md_simulation']['box_distance']),
               '-bt', self.config['md_simulation']['box_type']]

        if not self.run_gmx_command(cmd):
            return None

        # 2. 溶剂化
        cmd = ['gmx', 'solvate',
               '-cp', str(work_dir / 'boxed.gro'),
               '-cs', 'spc216.gro',
               '-o', str(work_dir / 'solvated.gro'),
               '-p', str(work_dir / 'topol.top')]

        if not self.run_gmx_command(cmd):
            return None

        return work_dir / 'solvated.gro'
    
    def add_ions(self, gro_file, work_dir):
        """添加离子平衡电荷"""
        print("添加离子...")

        work_dir = Path(work_dir).resolve()
        gro_file_abs = Path(gro_file).resolve()

        # 1. 生成.tpr文件用于添加离子
        cmd = ['gmx', 'grompp',
               '-f', str(work_dir / 'ions.mdp'),
               '-c', str(gro_file_abs),
               '-p', str(work_dir / 'topol.top'),
               '-o', str(work_dir / 'ions.tpr'),
               '-maxwarn', '1']  # 允许净电荷警告

        if not self.run_gmx_command(cmd):
            return None

        # 2. 添加离子
        cmd = ['gmx', 'genion',
               '-s', str(work_dir / 'ions.tpr'),
               '-o', str(work_dir / 'ionized.gro'),
               '-p', str(work_dir / 'topol.top'),
               '-pname', 'NA',
               '-nname', 'CL',
               '-neutral', '-conc', str(self.config['md_simulation']['ion_concentration'])]

        # 选择溶剂组
        input_data = "13\n"  # 通常SOL组是13

        if not self.run_gmx_command(cmd, input_data=input_data):
            return None

        return work_dir / 'ionized.gro'
    
    def create_mdp_files(self, work_dir):
        """创建MD参数文件"""
        print("创建MD参数文件...")
        
        mdp_configs = {
            'ions.mdp': self._get_ions_mdp(),
            'em.mdp': self._get_em_mdp(),
            'nvt.mdp': self._get_nvt_mdp(),
            'npt.mdp': self._get_npt_mdp(),
            'md.mdp': self._get_md_mdp()
        }
        
        for filename, content in mdp_configs.items():
            mdp_file = work_dir / filename
            with open(mdp_file, 'w') as f:
                f.write(content)
    
    def _get_ions_mdp(self):
        """离子化参数文件"""
        return """
; ions.mdp - used to generate a topology with a single ion
integrator  = steep        ; algorithm (steep = steepest descent minimization)
emtol       = 1000.0      ; stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01        ; initial step size
nsteps      = 50          ; maximum number of (minimization) steps

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1         ; frequency to update the neighbor list and long-range forces
cutoff-scheme   = Verlet    ; Buffered neighbor list
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = PME       ; Treatment of long-range electrostatic interactions
rcoulomb        = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw            = 1.0       ; short-range van der Waals cutoff (in nm)
pbc             = xyz       ; Periodic Boundary Conditions
"""
    
    def _get_em_mdp(self):
        """能量最小化参数文件"""
        return f"""
; em.mdp - Energy minimization
integrator  = steep        ; algorithm
emtol       = {self.config['md_simulation']['em_max_force']}     ; stop minimization when max force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; initial step size
nsteps      = {self.config['md_simulation']['em_nsteps']}        ; maximum steps

; Neighbor searching
nstlist         = 1         ; frequency to update neighbor list
cutoff-scheme   = Verlet
ns_type         = grid
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
"""
    
    def _get_nvt_mdp(self):
        """NVT平衡参数文件"""
        return f"""
; nvt.mdp - NVT equilibration
title                   = NVT Equilibration

integrator              = md
dt                      = {self.config['md_simulation']['nvt_time_step']}
nsteps                  = {int(self.config['md_simulation']['nvt_time'] / self.config['md_simulation']['nvt_time_step'])}

nstxout                 = 1000
nstvout                 = 1000
nstenergy               = 1000
nstlog                  = 1000

continuation            = no
constraint_algorithm    = Lincs
constraints             = h-bonds

cutoff-scheme           = Verlet
nstlist                 = 10
rcoulomb                = 1.0
rvdw                    = 1.0
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16

tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1     0.1
ref_t                   = {self.config['md_simulation']['nvt_temperature']} {self.config['md_simulation']['nvt_temperature']}

pcoupl                  = no
pbc                     = xyz
"""
    
    def _get_npt_mdp(self):
        """NPT平衡参数文件"""
        return f"""
; npt.mdp - NPT equilibration
title                   = NPT Equilibration

integrator              = md
dt                      = {self.config['md_simulation']['npt_time_step']}
nsteps                  = {int(self.config['md_simulation']['npt_time'] / self.config['md_simulation']['npt_time_step'])}

nstxout                 = 1000
nstvout                 = 1000
nstenergy               = 1000
nstlog                  = 1000

continuation            = yes
constraint_algorithm    = Lincs
constraints             = h-bonds

cutoff-scheme           = Verlet
nstlist                 = 10
rcoulomb                = 1.0
rvdw                    = 1.0
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16

tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1     0.1
ref_t                   = {self.config['md_simulation']['npt_temperature']} {self.config['md_simulation']['npt_temperature']}

pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = {self.config['md_simulation']['npt_pressure']}
compressibility         = 4.5e-5

pbc                     = xyz
"""

    def _get_md_mdp(self):
        """生产MD参数文件"""
        return f"""
; md.mdp - Production MD
title                   = Production MD

integrator              = md
dt                      = {self.config['md_simulation']['production_time_step']}
nsteps                  = {int(self.config['md_simulation']['production_time'] * 1000 / self.config['md_simulation']['production_time_step'])}

nstxout                 = 0
nstvout                 = 0
nstenergy               = {int(self.config['md_simulation']['save_interval'] * 1000 / self.config['md_simulation']['production_time_step'])}
nstlog                  = {int(self.config['md_simulation']['save_interval'] * 1000 / self.config['md_simulation']['production_time_step'])}

continuation            = yes
constraint_algorithm    = Lincs
constraints             = h-bonds

cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 10
rcoulomb                = 1.0
rvdw                    = 1.0
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16

tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1     0.1
ref_t                   = {self.config['md_simulation']['production_temperature']} {self.config['md_simulation']['production_temperature']}

pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = {self.config['md_simulation']['production_pressure']}
compressibility         = 4.5e-5

pbc                     = xyz
"""
    
    def energy_minimization(self, work_dir):
        """能量最小化"""
        print("执行能量最小化...")
        
        work_dir = Path(work_dir).resolve()
        
        # 1. 生成.tpr文件
        cmd = ['gmx', 'grompp',
               '-f', str(work_dir / 'em.mdp'),
               '-c', str(work_dir / 'ionized.gro'),
               '-p', str(work_dir / 'topol.top'),
               '-o', str(work_dir / 'em.tpr')]
        
        if not self.run_gmx_command(cmd):
            return False
        
        # 2. 运行能量最小化（使用CUDA GPU加速）
        cmd = ['gmx', 'mdrun', '-nb', 'gpu', '-ntmpi', '1', '-ntomp', '8', '-v', '-deffnm', 'em']

        if not self.run_gmx_command(cmd, cwd=str(work_dir)):
            return False
            
        return True
    
    def nvt_equilibration(self, work_dir):
        """NVT平衡"""
        print("执行NVT平衡...")

        work_dir = Path(work_dir).resolve()

        # 1. 生成.tpr文件
        cmd = ['gmx', 'grompp',
               '-f', str(work_dir / 'nvt.mdp'),
               '-c', str(work_dir / 'em.gro'),
               '-r', str(work_dir / 'em.gro'),
               '-p', str(work_dir / 'topol.top'),
               '-o', str(work_dir / 'nvt.tpr')]

        if not self.run_gmx_command(cmd):
            return False

        # 2. 运行NVT模拟（使用CUDA GPU加速）
        cmd = ['gmx', 'mdrun', '-nb', 'gpu', '-ntmpi', '1', '-ntomp', '8', '-v', '-deffnm', 'nvt']

        if not self.run_gmx_command(cmd, cwd=str(work_dir)):
            return False

        return True
    
    def npt_equilibration(self, work_dir):
        """NPT平衡"""
        print("执行NPT平衡...")

        work_dir = Path(work_dir).resolve()

        # 1. 生成.tpr文件
        cmd = ['gmx', 'grompp',
               '-f', str(work_dir / 'npt.mdp'),
               '-c', str(work_dir / 'nvt.gro'),
               '-r', str(work_dir / 'nvt.gro'),
               '-t', str(work_dir / 'nvt.cpt'),
               '-p', str(work_dir / 'topol.top'),
               '-o', str(work_dir / 'npt.tpr')]

        if not self.run_gmx_command(cmd):
            return False

        # 2. 运行NPT模拟（使用CUDA GPU加速）
        cmd = ['gmx', 'mdrun', '-nb', 'gpu', '-ntmpi', '1', '-ntomp', '8', '-v', '-deffnm', 'npt']

        if not self.run_gmx_command(cmd, cwd=str(work_dir)):
            return False

        return True

    def production_md(self, work_dir):
        """生产MD模拟"""
        print("执行生产MD模拟...")

        work_dir = Path(work_dir).resolve()

        # 1. 生成.tpr文件
        cmd = ['gmx', 'grompp',
               '-f', str(work_dir / 'md.mdp'),
               '-c', str(work_dir / 'npt.gro'),
               '-t', str(work_dir / 'npt.cpt'),
               '-p', str(work_dir / 'topol.top'),
               '-o', str(work_dir / 'md.tpr')]

        if not self.run_gmx_command(cmd):
            return False

        # 2. 运行生产MD模拟（使用CUDA GPU加速）
        cmd = ['gmx', 'mdrun', '-nb', 'gpu', '-ntmpi', '1', '-ntomp', '8', '-v', '-deffnm', 'md']

        if not self.run_gmx_command(cmd, cwd=str(work_dir)):
            return False

        return True
    
    def run_complete_md(self, pdb_file, protein_name):
        """运行完整的MD模拟流程"""
        print(f"开始 {protein_name} 的MD模拟...")
        
        # 创建工作目录
        work_dir = self.trajectory_dir / protein_name
        work_dir.mkdir(exist_ok=True)
        
        try:
            # 1. 准备拓扑
            processed_gro = self.prepare_topology(pdb_file, work_dir)
            if not processed_gro:
                return False
            
            # 2. 创建参数文件
            self.create_mdp_files(work_dir)
            
            # 3. 定义盒子和溶剂化
            solvated_gro = self.define_box_and_solvate(processed_gro, work_dir)
            if not solvated_gro:
                return False
            
            # 4. 添加离子
            ionized_gro = self.add_ions(solvated_gro, work_dir)
            if not ionized_gro:
                return False
            
            # 5. 能量最小化
            if not self.energy_minimization(work_dir):
                return False
            
            # 6. NVT平衡
            if not self.nvt_equilibration(work_dir):
                return False
            
            # 7. NPT平衡
            if not self.npt_equilibration(work_dir):
                return False
            
            # 8. 生产MD模拟
            if not self.production_md(work_dir):
                return False
            
            print(f"{protein_name} MD模拟完成!")
            return True
            
        except Exception as e:
            print(f"{protein_name} MD模拟失败: {e}")
            return False

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description="分子动力学模拟")
    parser.add_argument("--config", default="config/protein_config.yaml", help="配置文件路径")
    parser.add_argument("--pdb", help="PDB文件路径")
    parser.add_argument("--name", help="蛋白名称")
    
    args = parser.parse_args()
    
    # 初始化MD模拟器
    simulator = MDSimulator(args.config)
    
    if not simulator.check_gromacs():
        sys.exit(1)
    
    if args.pdb and args.name:
        # 运行单个蛋白的MD模拟
        success = simulator.run_complete_md(args.pdb, args.name)
        if success:
            print("MD模拟完成!")
        else:
            print("MD模拟失败!")
            sys.exit(1)
    else:
        print("请提供PDB文件和蛋白名称")
        sys.exit(1)

if __name__ == "__main__":
    main()
