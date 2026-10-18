# Projection Suggestion

## 1. 读取情况

- Backend: VASP
- Input directory: C:\Users\Admin\Documents\wannier自动拟合\Auxiliary Wannier Band Fitting 1.0\example\inputs-file\VASP
- Structure file: C:\Users\Admin\Documents\wannier自动拟合\Auxiliary Wannier Band Fitting 1.0\example\inputs-file\VASP\POSCAR
- PDOS files: 2 files
- Fermi energy: 5.9512 eV
- Fermi interval: [3.951200, 7.951200] eV
- Spin mode: soc_noncollinear
- SOC: True
- Formula: MgB3
- Elements and atom numbers: {'Mg': 1, 'B': 3}

## 2. PDOS contribution near Fermi level

| Rank | Element | Orbital | Weight | Percent |
|---:|---|---|---:|---:|
| 1 | B | p | 1.47619 | 66.34% |
| 2 | B | s | 0.385025 | 17.30% |
| 3 | Mg | s | 0.163627 | 7.35% |
| 4 | Mg | p | 0.102697 | 4.62% |
| 5 | Mg | d | 0.0811944 | 3.65% |
| 6 | Mg | dx2 | 0.0164992 | 0.74% |
| 7 | B | d | 0 | 0.00% |
| 8 | B | dx2 | 0 | 0.00% |

## 3. Recommended projections

### Recommended projections 1: Priority recommendation

Fermi_energy = 5.9512 eV

num_wann = 13

begin projections
B:p
B:s
Mg:s
end projections

Reason:
- 该方案使用费米窗口内最主要的元素-轨道壳层。
- 当前累计贡献约为 90.99%。

### Recommended projections 2: Expanded recommendation

Fermi_energy = 5.9512 eV

num_wann = 16

begin projections
B:p
B:s
Mg:s
Mg:p
end projections

Reason:
- 该方案包含更多杂化轨道，对复杂能带更稳健。
- 当前累计贡献约为 95.61%。

### Recommended projections 3: Orbital-resolved recommendation

Fermi_energy = 5.9512 eV

num_wann = 10

begin projections
B:pz
B:s
Mg:s
Mg:p
end projections

Reason:
- B:pz 明显强于同壳层其他分轨道。
- 当前累计贡献约为 95.61%。

## 4. Notes

- The fermi_interval is defined relative to the Fermi energy.
- All plotted energies are absolute energies and are not shifted by Fermi energy.
- The recommended num_wann is counted from orbital degeneracy.
- For SOC or noncollinear calculations, please manually check whether spinor-related settings require further adjustment.
- 注意：检测到 SOC 或 noncollinear 计算。当前 num_wann 按轨道数统计，请自行检查是否需要考虑 spinor 自由度。
