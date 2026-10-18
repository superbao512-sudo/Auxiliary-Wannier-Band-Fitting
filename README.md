# Auxiliary Wannier Band Fitting 1.0

辅助 Wannier 能带拟合程序用于从 VASP、QE 或 EPW 的能带、PDOS 和 Wannier 拟合能带文件中自动提取信息，给出 Wannier90 初始投影、窗口参数和拟合质量判断。程序只读取已有计算输出，不调用 VASP、QE、EPW 或 Wannier90。

## 1. 运行环境

推荐在 Linux 终端中运行：

```bash
python AWBF.py
```

Python 依赖：

- `numpy`
- `matplotlib`

程序路径全部使用 `pathlib.Path` 管理，默认工作目录是程序所在目录下的 `user_work/`。

## 2. 工作目录结构

第一次使用时运行功能 5 初始化；也可以直接使用命令行初始化：

```bash
python AWBF.py --task init
```

初始化后会生成：

```text
user_work/
├── input.dat
├── inputs-file/
│   ├── VASP/
│   ├── QE/
│   └── EPW/
└── outputs-file/
    ├── VASP/
    ├── QE/
    └── EPW/
```

说明：

- `user_work/input.dat` 是用户主配置文件。
- `user_work/inputs-file/<backend>/` 放输入文件。
- `user_work/outputs-file/<backend>/` 保存输出结果。
- 功能 5 如果发现 `input.dat` 已存在，会保留原文件，只补齐目录。
- `user_work/input.dat` 的默认内容来自程序内置模板；在开发目录中，该模板与项目根目录下的参考 `input.dat` 保持一致。

## 3. 启动界面

运行 `python AWBF.py` 后先选择语言：

```text
1) 中文
2) English
```

随后进入功能菜单：

```text
1) 费米面轨道贡献
2) 拟合能带建议
3) 拟合质量分析
4) 拟合扫描/拟合检查
5) 仅初始化目录和 input.dat
0) 退出
```

除功能 5 和退出外，程序会要求选择后端：

```text
1) VASP
2) QE
3) EPW
```

## 4. 输入文件

### VASP

常用文件：

- `BAND.dat`：DFT 能带。
- `OUTCAR`：读取费米能级。
- `PDOS_*.dat`：功能一使用的 PDOS 文件。
- `*_band.dat`：Wannier 拟合能带，功能三使用。
- `KLABELS`：可选，高对称点标签。

### QE

常用文件：

- `band.dat.gnu`：DFT 能带。
- `scf.out`：读取费米能级。
- `*.pdos.pdos_atm*`：功能一使用的 PDOS 文件。
- `*_band.dat`：Wannier 拟合能带，功能三使用。
- `band.out`、`bands.out` 或 `fit.in`：可选，高对称点信息。

### EPW

常用文件：

- `band.dat.gnu`：DFT 能带。如果 EPW 目录缺失，程序可回退读取 `inputs-file/QE/` 中的 QE 能带。
- `scf.out`：读取费米能级。如果 EPW 目录缺失，程序可回退读取 `inputs-file/QE/`。
- `*_band.dat`：EPW/Wannier 拟合能带。

## 5. 功能说明

### 功能 1：费米面轨道贡献

读取结构、费米能级和 PDOS 文件，在 `fermi_interval` 指定的窗口内积分轨道贡献。

主要输出：

- `pdos_integrated_detail.csv`
- `pdos_integrated_summary.csv`
- `projection_suggestion.md`
- `dft_pdos_plot.png/.pdf`
- 轨道贡献柱状图

基本原理：

1. 读取费米能级 `E_F`。
2. 将 PDOS 能量统一到绝对能量标尺。
3. 在 `[E_F + low, E_F + high]` 范围内积分。
4. 按元素和轨道壳层累计贡献。
5. 根据贡献阈值生成 Wannier90 `begin projections ... end projections` 初始建议。

### 功能 2：拟合能带建议

读取 DFT 能带，统计每条能带的能量范围和是否穿过费米能级，然后推荐：

- `fit_start` / `fit_end`
- `num_wann`
- `num_bands`
- `dis_win_min` / `dis_win_max`
- `dis_froz_min` / `dis_froz_max`

主要输出：

- `dft_band.csv`
- `band_suggestion.md`
- `dft_band_plot.png/.pdf`
- `dft_band_with_window.png/.pdf`

功能二现在可以独立运行，不再强制自动运行功能一。

运行逻辑：

- 如果 `outputs-file/<backend>/projection_suggestion.md` 已存在，功能二自动读取功能一给出的 `num_wann` 和 projections。
- 如果没有功能一输出，则读取 `input.dat` 功能二区域中的 `num_band`。
- 若 `num_band = auto` 且没有功能一输出，功能二会提示用户先运行功能一或手动填写 `num_band`。

`num_band` 示例：

```ini
num_band = 8
```

这里的 `num_band` 用作功能二手动参考的目标能带数量。程序会围绕费米能级或空带间隙选择对应数量的能带，并据此给出窗口建议。

### 功能 3：拟合质量分析

读取 DFT 能带和 Wannier 拟合能带，比较：

```text
error = E_wannier_absolute - E_dft_absolute
```

主要输出：

- `fit_quality_report.md`
- `dft_vs_wannier_band.png/.pdf`
- `dft_vs_wannier_error.png/.pdf`

报告包含：

- 实际读取的 DFT 文件和 Wannier 文件。
- 费米能级及其来源。
- 横坐标对齐因子。
- Overall RMSE 和最大误差。
- Near-Fermi RMSE 和最大误差。
- 费米面交叉数量是否一致。
- 每条匹配能带的误差。
- 最终判断：`good`、`warning` 或 `bad`。

`metric_focus` 的含义：

- `metric_focus = fermi`：最终质量主要依据费米附近误差和费米交叉一致性；全能量范围误差作为提示。
- `metric_focus = full`：最终质量同时严格检查全窗口误差。

### 功能 4：拟合扫描/拟合检查

当前为预留功能入口，尚未完整实现。

## 6. 能量零点规则

程序所有能带图和拟合质量分析都采用未扣除费米能的绝对能量标尺。

规则：

- QE DFT 能带 `band.dat.gnu`：直接使用原始能量。
- QE/EPW/Wannier `*_band.dat`：直接使用原始能量。
- VASP `BAND.dat`：后处理文件通常已扣除费米能，因此程序读取后执行：

```text
E_vasp_absolute = E_BAND.dat + E_F
```

- VASP Wannier `*_band.dat`：直接使用原始能量。
- 图中如果显示费米能级，画的是真实 `E_F` 位置，不是 `E = 0`。

## 7. 关键配置参数

### 全局

```ini
backend = VASP
input_root = ./inputs-file
output_root = ./outputs-file
task = all
```

这些路径相对于 `user_work/`。

### 功能一

```ini
fermi_interval = -2 2
projection_threshold_1 = 0.85
projection_threshold_2 = 0.95
orbital_anisotropy_ratio = 3.0
```

### 功能二

```ini
num_band = auto
min_gap_threshold = 0.05
extra_bands_for_diswin = 8
above_fermi_ratio = 0.30
force_rule = auto
```

`num_band` 仅在功能一输出不存在时作为手动目标能带数量。

### 功能三

```ini
fit_fermi_energy = auto
metric_focus = fermi
near_fermi_window = 1
max_allowed_rmse = 0.05
max_allowed_near_fermi_error = 0.03
```

`fit_fermi_energy = auto` 时，功能三优先读取功能一/二输出中的费米能级；若不存在，再从 `OUTCAR` 或 `scf.out` 回退读取。

## 8. 推荐工作流程

### 初次使用

```bash
python AWBF.py
```

选择：

```text
5) 仅初始化目录和 input.dat
```

然后把计算文件放入对应目录。

### 完整推荐流程

1. 运行功能 1，得到轨道贡献和 projections 建议。
2. 运行功能 2，读取功能一输出并生成窗口建议。
3. 根据建议设置 Wannier90/EPW。
4. 运行 Wannier90/EPW 得到 `*_band.dat`。
5. 运行功能 3 检查拟合质量。

### 只做能带窗口建议

如果不想运行功能一：

1. 在 `input.dat` 中设置：

```ini
num_band = 8
```

2. 直接运行功能 2。

程序会用该数值作为目标能带数量生成窗口建议。

## 9. 输出判断解释

RMSE：

```text
sqrt(mean((E_wannier - E_dft)^2))
```

最大误差：

```text
max(abs(E_wannier - E_dft))
```

Near-Fermi 误差只在 `E_F +/- near_fermi_window` 内统计。Overall 误差在当前 metric 窗口内统计。

当 `metric_focus = fermi` 时，即使 Overall RMSE 较大，只要费米附近误差满足阈值且费米交叉一致，最终仍可判断为 `good`，但报告会保留全局误差提示。

## 10. 常见问题

### 未检测到 user_work

先运行功能 5 初始化。

### 功能二提示 num_band=auto

说明没有找到功能一输出。解决办法：

- 先运行功能一；或
- 在 `input.dat` 中手动设置 `num_band`。

### 功能三缺少 Wannier 能带

请将 `*_band.dat` 放入对应的 `inputs-file/<backend>/` 目录。

### EPW 自动选错 Wannier 文件

如果目录中存在多个 `*_band.dat`，程序优先选择最近修改的文件。建议只保留当前要分析的拟合能带文件，或手动指定路径。

## 11. 打包内容

发布压缩包只包含：

- `AWBF.py`
- `src/`
- `example/`
- `README.md`

`user_work/` 是用户运行后生成的工作目录，不放入发布压缩包。
