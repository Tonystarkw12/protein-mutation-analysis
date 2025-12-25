#!/bin/bash

# 蛋白质突变影响分析项目环境配置脚本
# 作者: iFlow CLI
# 用途: 一键配置项目运行环境

set -e  # 遇到错误立即退出

echo "=========================================="
echo "蛋白质突变影响分析项目环境配置"
echo "=========================================="

# 检查conda是否可用
if ! command -v conda &> /dev/null; then
    echo "错误: conda未安装或未激活"
    echo "请先安装Miniconda或Anaconda"
    exit 1
fi

echo "✓ Conda检查通过"

# 检查Python版本
python_version=$(python3 --version 2>&1 | awk '{print $2}' | cut -d. -f1,2)
required_version="3.8"

if [ "$(printf '%s\n' "$required_version" "$python_version" | sort -V | head -n1)" = "$required_version" ]; then
    echo "✓ Python版本检查通过: $python_version"
else
    echo "错误: Python版本过低，需要3.8+，当前版本: $python_version"
    exit 1
fi

# 检查GPU可用性
if command -v nvidia-smi &> /dev/null; then
    echo "✓ GPU可用:"
    nvidia-smi --query-gpu=name,memory.total --format=csv,noheader,nounits | head -1
else
    echo "⚠ 警告: 未检测到GPU，将使用CPU运行（速度较慢）"
fi

# 创建conda环境
echo "创建conda环境..."
if conda env list | grep -q "protein_mutation_analysis"; then
    echo "环境已存在，跳过创建"
else
    conda env create -f environment.yml
    echo "✓ Conda环境创建完成"
fi

# 激活环境
echo "激活conda环境..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate protein_mutation_analysis

# 安装Python依赖
echo "安装Python依赖..."
pip install --upgrade pip
pip install -r requirements.txt

# 检查GROMACS安装
echo "检查GROMACS..."
if command -v gmx &> /dev/null; then
    echo "✓ GROMACS已安装:"
    gmx -version | head -1
else
    echo "⚠ GROMACS未安装，尝试安装..."
    conda install -c bioconda gromacs -y
    
    if command -v gmx &> /dev/null; then
        echo "✓ GROMACS安装成功"
    else
        echo "❌ GROMACS安装失败，请手动安装"
        echo "运行: conda install -c bioconda gromacs"
    fi
fi

# 测试Python导入
echo "测试Python依赖导入..."
python3 -c "
import sys
packages = ['yaml', 'Bio', 'pandas', 'numpy', 'matplotlib', 'seaborn']
failed = []
for pkg in packages:
    try:
        __import__(pkg)
        print(f'✓ {pkg}')
    except ImportError:
        failed.append(pkg)
        print(f'❌ {pkg}')

if failed:
    print(f'\\n导入失败的包: {failed}')
    print('请运行: pip install -r requirements.txt')
    sys.exit(1)
else:
    print('\\n✓ 所有Python依赖检查通过')
"

# 创建结果目录
echo "创建结果目录..."
mkdir -p results/{structures,trajectories,analysis,figures}

# 设置权限
chmod +x main.py
chmod +x src/*.py

echo "=========================================="
echo "环境配置完成!"
echo "=========================================="
echo ""
echo "下一步操作:"
echo "1. 激活环境: conda activate protein_mutation_analysis"
echo "2. 运行分析: python main.py"
echo "3. 查看结果: cat results/final_report.md"
echo ""
echo "如需帮助，请查看 README.md"