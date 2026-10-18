# Band Fitting Suggestion

## 1. 读取情况

- 后端：VASP
- 输入目录：C:\Users\Admin\Documents\wannier自动拟合\Auxiliary Wannier Band Fitting 1.0\example\inputs-file\VASP
- DFT band 文件：C:\Users\Admin\Documents\wannier自动拟合\Auxiliary Wannier Band Fitting 1.0\example\inputs-file\VASP\BAND.dat
- scf.out 或 OUTCAR：C:\Users\Admin\Documents\wannier自动拟合\Auxiliary Wannier Band Fitting 1.0\example\inputs-file\VASP\OUTCAR
- 能量单位：eV
- Fermi energy：5.9512 eV
- 能带数量：16
- k 点数量：180
- Spin mode：auto
- 元素和原子数：{'Mg': 1, 'B': 3}

## 2. Wannier90 初始建议

### 主要建议

- projections：B:p, B:s, Mg:s
- num_wann：13
- num_bands：16
- fit bands：1 - 13
- dis_win_min：-8.2
- dis_win_max：23.5
- dis_froz_min：-4.2
- dis_froz_max：15.4
- 使用规则：规则二

Reason:
- 围绕费米能级选择目标能带并保留一定比例的费米能级以上能带。
- 窗口值为未扣除 Fermi energy 的绝对能量。

### 建议二

- projections：B:p, B:s, Mg:s, Mg:p
- num_wann：16
- num_bands：16
- fit bands：1 - 16
- dis_win_min：-8.2
- dis_win_max：23.5
- dis_froz_min：-4.2
- dis_froz_max：19.6
- 使用规则：规则二

Reason:
- 围绕费米能级选择目标能带并保留一定比例的费米能级以上能带。
- 窗口值为未扣除 Fermi energy 的绝对能量。

### 建议三

- projections：B:pz, B:s, Mg:s, Mg:p
- num_wann：10
- num_bands：16
- fit bands：1 - 10
- dis_win_min：-8.2
- dis_win_max：23.5
- dis_froz_min：-4.2
- dis_froz_max：12.1
- 使用规则：规则二

Reason:
- 围绕费米能级选择目标能带并保留一定比例的费米能级以上能带。
- 窗口值为未扣除 Fermi energy 的绝对能量。

## 3. Warnings

- 未找到合适空带区域，已使用规则二。

## 4. Notes

- 以上窗口值默认是未扣除 Fermi energy 的绝对能量。
- 如果 Wannier90 流程使用相对费米能级的能量窗口，需要手动减去 E_F。
- 本程序只给出初始建议，最终仍需结合 Wannier 拟合能带图人工检查。
