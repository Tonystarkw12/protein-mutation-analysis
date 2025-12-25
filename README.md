# 蛋白质突变影响分析项目

## 项目概述

本项目实现了一个完整的蛋白质突变影响计算验证流程，以p53 DNA结合域R273H突变为例，展示了从结构预测到分子动力学模拟再到结构分析的全过程。

**目标**: 量化单点突变对蛋白质结构稳定性与功能的影响  
**时间**: 1-2天完成  
**适配**: 蛋白质结构预测、分子动力学模拟、蛋白质设计

## 项目特点

- ✅ **完整流程**: 序列→结构预测→MD模拟→结构分析→结果报告
- ✅ **课题组适配**: 紧扣蛋白质结构预测、分子动力学核心技术
- ✅ **可复现性**: 标准化配置文件，一键运行
- ✅ **科学价值**: 量化突变影响，为实验提供计算指导

## 技术栈

- **结构预测**: ColabFold (AlphaFold2)
- **分子动力学**: GROMACS 2024
- **结构分析**: MDAnalysis, PyMOL, py3Dmol
- **可视化**: matplotlib, seaborn, plotly
- **编程语言**: Python 3.11+

## 环境要求

### 硬件要求
- GPU: 8GB+ 显存 (推荐RTX 4060及以上)
- 内存: 16GB+ RAM
- 存储: 50GB+ 可用空间

### 软件依赖
```bash
# 创建conda环境
conda env create -f environment.yml

# 激活环境
conda activate protein_mutation_analysis

# 安装Python依赖
pip install -r requirements.txt
```

## 快速开始

### 1. 环境配置
```bash
# 克隆项目
git clone <repository_url>
cd colabfold

# 创建环境
conda env create -f environment.yml
conda activate protein_mutation_analysis

# 安装依赖
pip install -r requirements.txt
```

### 2. 运行完整分析
```bash
# 运行完整流程（推荐）
python main.py

# 分步运行
python main.py --step 1  # 结构预测
python main.py --step 2  # MD模拟
python main.py --step 3  # 结构分析
python main.py --step 4  # 生成报告
```

### 3. 查看结果
```bash
# 查看最终报告
cat results/final_report.md

# 查看可视化结果
ls results/figures/
```

## 项目结构

```
colabfold/
├── config/
│   └── protein_config.yaml      # 项目配置文件
├── data/
│   └── p53_sequences.fasta      # 蛋白质序列文件
├── src/
│   ├── structure_prediction.py  # 结构预测模块
│   ├── md_simulation.py         # MD模拟模块
│   └── structure_analysis.py    # 结构分析模块
├── results/                     # 结果输出目录
│   ├── structures/              # 预测结构
│   ├── trajectories/            # MD轨迹
│   ├── analysis/               # 分析数据
│   └── figures/                # 图表可视化
├── main.py                      # 主执行脚本
├── requirements.txt             # Python依赖
├── environment.yml              # Conda环境配置
└── README.md                    # 项目说明
```

## 核心模块说明

### 1. 结构预测模块 (`structure_prediction.py`)
- 使用ColabFold进行AlphaFold2预测
- 支持野生型和突变体并行预测
- 自动质量评估和可视化

### 2. MD模拟模块 (`md_simulation.py`)
- 完整GROMACS工作流程
- 包含能量最小化、NVT/NPT平衡、生产MD
- 自动参数配置和错误处理

### 3. 结构分析模块 (`structure_analysis.py`)
- RMSD、回转半径、氢键、SASA分析
- 交互式可视化和统计报告
- 野生型与突变体对比分析

## 配置说明

主要配置在 `config/protein_config.yaml`:

```yaml
protein:
  name: "p53_DNA_binding_domain"
  pdb_id: "1TUP"
  
mutation:
  position: 273
  wild_type: "R"
  mutant: "H"

md_simulation:
  production_time: 10  # ns
  temperature: 310      # K
```

## 预期结果

### 1. 结构预测
- 野生型和突变体3D结构
- 置信度评分 (pLDDT)
- 结构叠加比较

### 2. MD模拟
- 10ns轨迹文件
- 能量收敛曲线
- 温度/压力平衡验证

### 3. 结构分析
- RMSD时间序列
- 回转半径分布
- 氢键网络变化
- 溶剂可及表面积

### 4. 科学结论
- 突变对结构稳定性的定量影响
- 构象变化机制分析
- 功能影响预测

## 时间预估

| 步骤 | 预计时间 | 说明 |
|------|----------|------|
| 环境配置 | 2-3小时 | 依赖安装和配置 |
| 结构预测 | 2-3小时 | GPU加速 |
| MD模拟 | 6-8小时 | 10ns模拟 |
| 结构分析 | 1-2小时 | 数据处理和可视化 |
| 报告生成 | 0.5小时 | 整理结果 |

**总计**: 1-2天

## 课题组适配要点

### 1. 工具替换
- 可使用课题组GDFold2替代ColabFold
- 可集成课题组专用分析流程
- 支持课题组标准力场参数

### 2. 深度拓展
- 增加更多突变体对比
- 使用不同力场交叉验证
- 集成自由能计算方法

### 3. 标准化
- 遵循课题组数据格式标准
- 支持课题组计算集群
- 集成现有分析流程

## 故障排除

### 常见问题

1. **GROMACS安装失败**
   ```bash
   conda install -c bioconda gromacs
   ```

2. **GPU内存不足**
   - 减少batch_size
   - 缩短MD模拟时间
   - 使用CPU模式

3. **依赖包冲突**
   ```bash
   pip install --upgrade pip
   pip install -r requirements.txt --force-reinstall
   ```

## 贡献指南

欢迎提交Issue和Pull Request来改进项目。

## 许可证

本项目采用MIT许可证。

## 联系方式

如有问题，请通过以下方式联系：
- 项目Issues
- 邮件联系 zhou-zh23@mails.tsinghua.edu.cn
