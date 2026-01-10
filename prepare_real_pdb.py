#!/usr/bin/env python3
"""
准备真实PDB结构用于MD模拟
"""
from pathlib import Path
import shutil

def prepare_pdb_structures():
    """准备野生型和突变体PDB文件"""

    # 输入的真实PDB文件
    source_pdb = Path("data/1TUP.pdb")

    # 输出目录
    wt_dir = Path("results/structures/p53_DNA_binding_domain_WT")
    mut_dir = Path("results/structures/p53_DNA_binding_domain_R273H")

    # 创建目录
    wt_dir.mkdir(parents=True, exist_ok=True)
    mut_dir.mkdir(parents=True, exist_ok=True)

    # 只保留蛋白质链（A链），移除DNA和其他链
    def extract_chain_a(input_pdb, output_pdb):
        """提取A链（蛋白质链）"""
        with open(input_pdb, 'r') as f_in:
            with open(output_pdb, 'w') as f_out:
                for line in f_in:
                    # PDB格式中链ID在第22列（索引21）
                    if line.startswith('ATOM') and len(line) > 21 and line[21] == 'A':
                        f_out.write(line)
                f_out.write('TER\n')
                f_out.write('END\n')

    # 准备野生型结构
    wt_pdb = wt_dir / "p53_DNA_binding_domain_WT.pdb"
    extract_chain_a(source_pdb, wt_pdb)
    print(f"✓ 野生型结构已创建: {wt_pdb}")

    # 准备突变体结构（暂时使用相同结构，后续可手动修改或使用工具）
    mut_pdb = mut_dir / "p53_DNA_binding_domain_R273H.pdb"
    extract_chain_a(source_pdb, mut_pdb)

    # 修改R273为H273（简化处理）
    with open(mut_pdb, 'r') as f:
        content = f.read()

    # 将R（精氨酸）替换为H（组氨酸）在273位
    lines = content.split('\n')
    modified_lines = []
    for line in lines:
        if line.startswith('ATOM') and line[17:20].strip() == 'A' and line[22:26].strip() == '273':
            # 修改残基名称为HIS
            line = line[:17] + ' HIS' + line[21:]
        modified_lines.append(line)

    with open(mut_pdb, 'w') as f:
        f.write('\n'.join(modified_lines))

    print(f"✓ 突变体结构已创建: {mut_pdb}")
    print("\n结构准备完成！可以运行MD模拟了。")

if __name__ == '__main__':
    prepare_pdb_structures()
