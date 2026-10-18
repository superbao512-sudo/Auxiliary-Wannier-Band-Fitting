#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""通用配置、路径、命令行和流程调度工具。

本程序只读取已有 VASP/QE/EPW 输出文件，给出 Wannier90 初始 projections、
num_wann、dis_win 和 dis_froz 建议；不会调用任何第一性原理或 Wannier 程序。
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

try:
    import numpy as np
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:  # pragma: no cover - 顶层给出清楚错误，避免 traceback 淹没原因
    print("Error: 需要可用的 numpy 和 matplotlib。当前 Python 导入失败：%s" % exc, file=sys.stderr)
    print("建议使用 Codex bundled Python 或重新安装科学计算库。", file=sys.stderr)
    raise SystemExit(1)


from .default_inputs import DEFAULT_INPUTS_DAT
from . import workspace

APP_ROOT = Path(__file__).resolve().parents[1]
ORBITAL_DEGENERACY = {
    "s": 1,
    "p": 3,
    "d": 5,
    "f": 7,
    "sp3": 4,
    "px": 1,
    "py": 1,
    "pz": 1,
    "dxy": 1,
    "dyz": 1,
    "dz2": 1,
    "dxz": 1,
    "dx2-y2": 1,
}
FLOAT_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
SUB_ORBITAL_GROUP = {
    "s": "s",
    "px": "p",
    "py": "p",
    "pz": "p",
    "p": "p",
    "dxy": "d",
    "dyz": "d",
    "dz2": "d",
    "dxz": "d",
    "dx2-y2": "d",
    "d": "d",
    "f": "f",
}
DEFAULT_INPUT = DEFAULT_INPUTS_DAT


@dataclass
class Context:
    backend: str
    task: str
    cfg: Dict[str, str]
    input_dir: Path
    output_dir: Path
    warnings: List[str]


@dataclass
class StructureInfo:
    structure_file: Optional[Path]
    elements: Dict[str, int]
    formula: str


@dataclass
class PdosResult:
    fermi_energy: float
    e_min: float
    e_max: float
    structure: StructureInfo
    spin_mode: str
    soc: bool
    detail_rows: List[Dict[str, object]]
    summary_rows: List[Dict[str, object]]
    energy_curves: Dict[str, Tuple[np.ndarray, np.ndarray]]
    recommendations: List[Dict[str, object]]


@dataclass
class BandData:
    path: Path
    k_distances: List[float]
    bands: List[List[float]]
    high_symmetry: List[Tuple[float, str]]


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def fmt_fermi(value: float) -> str:
    return str(value)


def trunc1(value: float) -> str:
    scaled = math.trunc(float(value) * 10.0) / 10.0
    if scaled == 0:
        scaled = 0.0
    return f"{scaled:.1f}"


def cfg_float(cfg: Dict[str, str], key: str, default: float) -> float:
    value = cfg.get(key, "auto").strip()
    if value.lower() == "auto":
        return default
    try:
        return float(value)
    except ValueError:
        return default


def cfg_bool(cfg: Dict[str, str], key: str, default: bool) -> bool:
    value = cfg.get(key, "").strip().lower()
    if not value:
        return default
    return value in {"1", "yes", "true", "on", "y"}


def strip_comment(line: str) -> str:
    return line.split("#", 1)[0].strip()


def parse_input(path: Path) -> Dict[str, str]:
    if not path.exists():
        raise FileNotFoundError(f"缺少配置文件：{path}")
    cfg: Dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = strip_comment(raw)
        if not line or "=" not in line:
            continue
        key, value = line.split("=", 1)
        cfg[key.strip().lower()] = value.strip()
    return cfg


def resolve_path(text: str, base: Path) -> Path:
    path = Path(text).expanduser()
    return path if path.is_absolute() else (base / path)


def backend_input_dir(cfg: Dict[str, str], backend: str) -> Path:
    root = resolve_path(cfg.get("input_root", "./inputs-file"), workspace.USER_WORK)
    primary = root / backend
    return primary


def backend_output_dir(cfg: Dict[str, str], backend: str) -> Path:
    root = resolve_path(cfg.get("output_root", "./outputs-file"), workspace.USER_WORK)
    out = root / backend
    out.mkdir(parents=True, exist_ok=True)
    return out


def make_context(args: argparse.Namespace) -> Context:
    workspace.require_workspace()
    cfg = parse_input(workspace.config_path())
    backend = (args.backend or cfg.get("backend", "VASP")).upper()
    task = (args.task or cfg.get("task", "all")).lower()
    if backend not in {"VASP", "QE", "EPW"}:
        raise ValueError("backend 只能是 VASP、QE 或 EPW。")
    if task not in {"projection", "band", "all"}:
        raise ValueError("task 只能是 projection、band 或 all。")
    return Context(
        backend=backend,
        task=task,
        cfg=cfg,
        input_dir=backend_input_dir(cfg, backend),
        output_dir=backend_output_dir(cfg, backend),
        warnings=[],
    )


def find_first(directory: Path, names: Sequence[str]) -> Optional[Path]:
    for name in names:
        path = directory / name
        if path.exists():
            return path
    return None


def parse_vasp_structure(directory: Path) -> StructureInfo:
    path = find_first(directory, ["POSCAR", "CONTCAR"])
    if not path:
        return StructureInfo(None, {}, "unknown")
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    elements: Dict[str, int] = {}
    if len(lines) >= 7:
        names = lines[5].split()
        counts = lines[6].split()
        if all(re.match(r"^[A-Za-z][A-Za-z0-9_+-]*$", item) for item in names):
            for name, count in zip(names, counts):
                if count.isdigit():
                    elements[name] = elements.get(name, 0) + int(count)
    return StructureInfo(path, elements, formula_from_elements(elements))


def parse_qe_structure(directory: Path) -> StructureInfo:
    path = find_first(directory, ["scf.in", "relax.in", "vc-relax.in"])
    if not path:
        return StructureInfo(None, {}, "unknown")
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    elements: Dict[str, int] = {}
    in_positions = False
    for raw in lines:
        line = strip_comment(raw)
        if not line:
            if in_positions:
                break
            continue
        upper = line.upper()
        if upper.startswith("ATOMIC_POSITIONS"):
            in_positions = True
            continue
        if in_positions:
            if re.match(r"^[A-Z_]+", upper) and not re.match(r"^[A-Z][a-z]?\s+", line):
                break
            symbol = line.split()[0]
            if re.match(r"^[A-Z][a-z]?$", symbol):
                elements[symbol] = elements.get(symbol, 0) + 1
    return StructureInfo(path, elements, formula_from_elements(elements))


def formula_from_elements(elements: Dict[str, int]) -> str:
    if not elements:
        return "unknown"
    chunks = []
    for element, count in elements.items():
        chunks.append(element if count == 1 else f"{element}{count}")
    return "".join(chunks)


def read_fermi(ctx: Context) -> float:
    manual = ctx.cfg.get("fermi_energy", "auto").strip()
    if manual.lower() != "auto":
        try:
            return float(manual.split()[0].replace(",", ""))
        except ValueError:
            ctx.warnings.append("fermi_energy 手动值无法解析，改为从输出文件读取。")
    if ctx.backend == "VASP":
        outcar = find_first(ctx.input_dir, ["OUTCAR"])
        if not outcar:
            raise FileNotFoundError("缺少 OUTCAR，无法读取 VASP Fermi energy。")
        value = None
        for line in outcar.read_text(encoding="utf-8", errors="ignore").splitlines():
            match = re.search(r"E-fermi\s*:\s*(%s)" % FLOAT_RE, line, re.I)
            if match:
                value = float(match.group(1).split()[0])
            match = re.search(r"\bEFERMI\s*=\s*(%s)" % FLOAT_RE, line, re.I)
            if match and value is None:
                value = float(match.group(1).split()[0])
        if value is None:
            raise ValueError("OUTCAR 中未找到 E-fermi。")
        return value
    scf_out = find_first(ctx.input_dir, ["scf.out", "pw.out"])
    if not scf_out and ctx.backend == "EPW":
        scf_out = find_first(workspace.backend_input_dir("QE"), ["scf.out", "pw.out"])
    if not scf_out:
        raise FileNotFoundError("缺少 scf.out，无法读取 QE/EPW Fermi energy。")
    for line in scf_out.read_text(encoding="utf-8", errors="ignore").splitlines():
        match = re.search(r"Fermi energy is\s+(%s)" % FLOAT_RE, line, re.I)
        if match:
            return float(match.group(1))
    raise ValueError("scf.out 中未找到 Fermi energy。")


def detect_spin(ctx: Context) -> Tuple[str, bool]:
    spin_mode = ctx.cfg.get("spin_mode", "auto").lower()
    if spin_mode != "auto":
        return spin_mode, "soc" in spin_mode or "noncol" in spin_mode
    text = ""
    for name in ["OUTCAR", "scf.in", "scf.out"]:
        path = ctx.input_dir / name
        if path.exists():
            text += "\n" + path.read_text(encoding="utf-8", errors="ignore").lower()
    soc = any(key in text for key in ["lsorbit", "lspinorb", "noncolin", "non-collinear"])
    if soc:
        return "soc_noncollinear", True
    if any(key in text for key in ["ispin  =      2", "ispin = 2", "nspin=2", "nspin = 2"]):
        return "spin_polarized", False
    return "nonspin", False


def is_relative_energy(energies: np.ndarray, fermi: float) -> bool:
    if len(energies) == 0:
        return False
    lo, hi = float(np.nanmin(energies)), float(np.nanmax(energies))
    return not (lo <= fermi <= hi) and (lo <= 0.0 <= hi)


def convert_to_absolute_energy(
    code_type: str,
    file_role: str,
    energy: np.ndarray,
    efermi: Optional[float] = None,
) -> np.ndarray:
    """Convert band energies to the absolute, non-EF-subtracted scale."""
    if code_type.lower() == "vasp" and file_role.lower() == "dft":
        if efermi is None:
            raise ValueError("VASP BAND.dat requires efermi to add back the Fermi energy.")
        return energy + efermi
    return energy


def read_numeric_table(path: Path) -> Tuple[List[str], np.ndarray]:
    header: List[str] = []
    rows: List[List[float]] = []
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            tokens = line.lstrip("#").split()
            if tokens and any(token.lower().startswith(("energy", "e(")) for token in tokens):
                header = tokens
            continue
        parts = line.split()
        try:
            rows.append([float(part.replace("D", "E")) for part in parts])
        except ValueError:
            continue
    if not rows:
        return header, np.empty((0, 0))
    width = max(len(row) for row in rows)
    clean = [row for row in rows if len(row) == width]
    return header, np.array(clean, dtype=float)


def integrate_vasp_pdos(ctx: Context, fermi: float) -> Tuple[List[Dict[str, object]], Dict[str, Tuple[np.ndarray, np.ndarray]]]:
    rows: List[Dict[str, object]] = []
    curves: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    for path in sorted(ctx.input_dir.glob(ctx.cfg.get("vasp_pdos", "PDOS_*.dat"))):
        header, data = read_numeric_table(path)
        if data.size == 0 or data.shape[1] < 2:
            ctx.warnings.append(f"PDOS 文件为空或格式无法识别：{path.name}")
            continue
        element = parse_element_from_name(path.name) or "X"
        energies = data[:, 0]
        if is_relative_energy(energies, fermi):
            energies = energies + fermi
        column_names = normalize_vasp_columns(header, data.shape[1])
        total_curve = np.zeros_like(energies)
        for index, orbital in enumerate(column_names[1:], start=1):
            if orbital.lower() in {"tot", "total"}:
                continue
            values = np.maximum(data[:, index], 0.0)
            total_curve += values
            rows.append(make_pdos_row(ctx.backend, element, "", orbital, "none", energies, values, fermi, ctx))
        curves[element] = (energies, total_curve)
    return rows, curves


def normalize_vasp_columns(header: Sequence[str], width: int) -> List[str]:
    if len(header) >= width:
        names = list(header[:width])
        names[0] = "Energy"
        return names
    default = ["Energy", "s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "dx2-y2", "tot"]
    if width <= len(default):
        return default[:width]
    return ["Energy"] + [f"orb{i}" for i in range(1, width)]


def parse_element_from_name(name: str) -> Optional[str]:
    match = re.search(r"PDOS[_-]?([A-Z][a-z]?)", name)
    if match:
        return match.group(1)
    match = re.search(r"\(([A-Z][a-z]?)\)", name)
    if match:
        return match.group(1)
    match = re.match(r"([A-Z][a-z]?)\.", name)
    return match.group(1) if match else None


def integrate_qe_pdos(ctx: Context, fermi: float) -> Tuple[List[Dict[str, object]], Dict[str, Tuple[np.ndarray, np.ndarray]]]:
    directory = ctx.input_dir
    if ctx.backend == "EPW" and not list(directory.glob("*.pdos.pdos_atm*")):
        fallback = workspace.backend_input_dir("QE")
        if fallback.exists():
            directory = fallback
            ctx.warnings.append("EPW 目录未找到 PDOS，已回退读取 QE 输入目录。")
    rows: List[Dict[str, object]] = []
    curve_acc: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    for path in sorted(directory.glob("*.pdos.pdos_atm*")):
        parsed = re.search(r"atm#(\d+)\(([^)]+)\).*?\(([spdf])\)", path.name)
        if not parsed:
            ctx.warnings.append(f"QE PDOS 文件名无法识别：{path.name}")
            continue
        atom_index, element, shell = parsed.group(1), parsed.group(2), parsed.group(3)
        _header, data = read_numeric_table(path)
        if data.size == 0 or data.shape[1] < 2:
            ctx.warnings.append(f"PDOS 文件为空或格式无法识别：{path.name}")
            continue
        energies = data[:, 0]
        if is_relative_energy(energies, fermi):
            energies = energies + fermi
        orbital_columns = qe_orbital_columns(shell, data.shape[1] - 1)
        total_curve = np.zeros_like(energies)
        for orbital, column in orbital_columns:
            values = np.maximum(data[:, column], 0.0)
            total_curve += values
            rows.append(make_pdos_row(ctx.backend, element, atom_index, orbital, "none", energies, values, fermi, ctx))
        if element in curve_acc:
            old_e, old_v = curve_acc[element]
            if len(old_v) == len(total_curve):
                curve_acc[element] = (old_e, old_v + total_curve)
        else:
            curve_acc[element] = (energies, total_curve)
    return rows, curve_acc


def qe_orbital_columns(shell: str, ncols_after_energy: int) -> List[Tuple[str, int]]:
    # QE projwfc 输出常见为 E, ldos, 分轨道...；这里跳过 ldos，按文档规定 pz,px,py / dz2,dxz,dyz,dx2-y2,dxy。
    if shell == "s":
        return [("s", ncols_after_energy)] if ncols_after_energy >= 1 else []
    if shell == "p":
        names = ["pz", "px", "py"]
        start = max(1, ncols_after_energy - 2)
        return [(name, start + i) for i, name in enumerate(names) if start + i <= ncols_after_energy]
    if shell == "d":
        names = ["dz2", "dxz", "dyz", "dx2-y2", "dxy"]
        start = max(1, ncols_after_energy - 4)
        return [(name, start + i) for i, name in enumerate(names) if start + i <= ncols_after_energy]
    return [(shell, ncols_after_energy)]


def make_pdos_row(
    backend: str,
    element: str,
    atom_index: object,
    orbital: str,
    spin_channel: str,
    energies: np.ndarray,
    values: np.ndarray,
    fermi: float,
    ctx: Context,
) -> Dict[str, object]:
    low, high = fermi_interval(ctx, fermi)
    mask = (energies >= low) & (energies <= high)
    if np.count_nonzero(mask) >= 2:
        weight = float(np.trapezoid(values[mask], energies[mask]))
    else:
        weight = 0.0
    return {
        "backend": backend,
        "element": element,
        "atom_index": atom_index,
        "orbital": orbital,
        "orbital_group": SUB_ORBITAL_GROUP.get(orbital, orbital),
        "spin_channel": spin_channel,
        "integral_weight": max(weight, 0.0),
        "percent": 0.0,
    }


def fermi_interval(ctx: Context, fermi: float) -> Tuple[float, float]:
    parts = ctx.cfg.get("fermi_interval", "-2 2").replace(",", " ").split()
    lo, hi = (float(parts[0]), float(parts[1])) if len(parts) >= 2 else (-2.0, 2.0)
    return fermi + lo, fermi + hi


def aggregate_pdos(rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    grouped: Dict[Tuple[str, str], float] = {}
    total = sum(float(row["integral_weight"]) for row in rows)
    for row in rows:
        key = (str(row["element"]), str(row["orbital_group"]))
        grouped[key] = grouped.get(key, 0.0) + float(row["integral_weight"])
    summary = [
        {
            "backend": rows[0]["backend"] if rows else "",
            "element": element,
            "orbital_group": orbital,
            "integral_weight": weight,
            "percent": 100.0 * weight / total if total > 0 else 0.0,
        }
        for (element, orbital), weight in grouped.items()
    ]
    summary.sort(key=lambda item: float(item["integral_weight"]), reverse=True)
    for row in rows:
        row["percent"] = 100.0 * float(row["integral_weight"]) / total if total > 0 else 0.0
    return summary


def count_num_wann(projections: Sequence[str], elements: Dict[str, int]) -> int:
    total = 0
    for projection in projections:
        if ":" not in projection:
            continue
        element, orbital_text = projection.split(":", 1)
        atom_count = elements.get(element.strip(), 1)
        orbitals = [item.strip() for item in orbital_text.split(",") if item.strip()]
        total += atom_count * sum(ORBITAL_DEGENERACY.get(orb, 1) for orb in orbitals)
    return total


def choose_until(summary: Sequence[Dict[str, object]], threshold: float) -> Tuple[List[str], float]:
    projections: List[str] = []
    cumulative = 0.0
    for row in summary:
        projections.append(f"{row['element']}:{row['orbital_group']}")
        cumulative += float(row["percent"]) / 100.0
        if cumulative >= threshold:
            break
    if not projections and summary:
        row = summary[0]
        projections.append(f"{row['element']}:{row['orbital_group']}")
        cumulative = float(row["percent"]) / 100.0
    return projections, cumulative


def anisotropy_projection(
    detail_rows: Sequence[Dict[str, object]],
    base: Sequence[str],
    ratio: float,
) -> Tuple[Optional[List[str]], List[str]]:
    by_element_group: Dict[Tuple[str, str], List[Dict[str, object]]] = {}
    for row in detail_rows:
        group = str(row["orbital_group"])
        orbital = str(row["orbital"])
        if group in {"p", "d"} and orbital != group:
            by_element_group.setdefault((str(row["element"]), group), []).append(row)
    replacements: Dict[str, List[str]] = {}
    reasons: List[str] = []
    base_set = set(base)
    for (element, group), rows in by_element_group.items():
        projection_key = f"{element}:{group}"
        if projection_key not in base_set:
            continue
        by_orbital: Dict[str, float] = {}
        for row in rows:
            by_orbital[str(row["orbital"])] = by_orbital.get(str(row["orbital"]), 0.0) + float(row["integral_weight"])
        nonzero = [value for value in by_orbital.values() if value > 0]
        if len(nonzero) >= 2 and max(nonzero) / min(nonzero) >= ratio:
            dominant = max(by_orbital, key=by_orbital.get)
            replacements[projection_key] = [f"{element}:{dominant}"]
            reasons.append(f"{element}:{dominant} 明显强于同壳层其他分轨道。")
    if not replacements:
        return None, []
    resolved: List[str] = []
    for projection in base:
        resolved.extend(replacements.get(projection, [projection]))
    return resolved, reasons


def build_recommendations(ctx: Context, result: PdosResult) -> List[Dict[str, object]]:
    threshold1 = float(ctx.cfg.get("projection_threshold_1", "0.85"))
    threshold2 = float(ctx.cfg.get("projection_threshold_2", "0.95"))
    projections1, cumulative1 = choose_until(result.summary_rows, threshold1)
    projections2, cumulative2 = choose_until(result.summary_rows, threshold2)
    recs = [
        {
            "index": 1,
            "title": "Priority recommendation",
            "projections": projections1,
            "num_wann": count_num_wann(projections1, result.structure.elements),
            "cumulative": cumulative1,
            "reason": "该方案使用费米窗口内最主要的元素-轨道壳层。",
        },
        {
            "index": 2,
            "title": "Expanded recommendation",
            "projections": projections2,
            "num_wann": count_num_wann(projections2, result.structure.elements),
            "cumulative": cumulative2,
            "reason": "该方案包含更多杂化轨道，对复杂能带更稳健。",
        },
    ]
    resolved, reasons = anisotropy_projection(
        result.detail_rows,
        projections2,
        float(ctx.cfg.get("orbital_anisotropy_ratio", "3.0")),
    )
    if resolved:
        recs.append(
            {
                "index": 3,
                "title": "Orbital-resolved recommendation",
                "projections": resolved,
                "num_wann": count_num_wann(resolved, result.structure.elements),
                "cumulative": cumulative2,
                "reason": "；".join(reasons),
            }
        )
    return recs


def run_projection(ctx: Context) -> PdosResult:
    fermi = read_fermi(ctx)
    structure = parse_vasp_structure(ctx.input_dir) if ctx.backend == "VASP" else parse_qe_structure(ctx.input_dir)
    spin_mode, soc = detect_spin(ctx)
    if ctx.backend == "VASP":
        detail_rows, curves = integrate_vasp_pdos(ctx, fermi)
    else:
        detail_rows, curves = integrate_qe_pdos(ctx, fermi)
    summary_rows = aggregate_pdos(detail_rows)
    e_min, e_max = fermi_interval(ctx, fermi)
    result = PdosResult(
        fermi_energy=fermi,
        e_min=e_min,
        e_max=e_max,
        structure=structure,
        spin_mode=spin_mode,
        soc=soc,
        detail_rows=detail_rows,
        summary_rows=summary_rows,
        energy_curves=curves,
        recommendations=[],
    )
    result.recommendations = build_recommendations(ctx, result)
    write_projection_outputs(ctx, result)
    return result


def write_csv(path: Path, rows: Sequence[Dict[str, object]], columns: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns))
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in columns})


def write_projection_outputs(ctx: Context, result: PdosResult) -> None:
    write_csv(
        ctx.output_dir / "pdos_integrated_detail.csv",
        result.detail_rows,
        ["backend", "element", "atom_index", "orbital", "spin_channel", "integral_weight", "percent"],
    )
    write_csv(
        ctx.output_dir / "pdos_integrated_summary.csv",
        result.summary_rows,
        ["backend", "element", "orbital_group", "integral_weight", "percent"],
    )
    write_text(ctx.output_dir / "projection_suggestion.md", projection_markdown(ctx, result))
    plot_pdos(ctx, result)
    plot_bar(
        ctx.output_dir / "pdos_orbital_bar1",
        [
            (f"{row['element']}_{row['orbital']}", float(row["percent"]))
            for row in sorted(result.detail_rows, key=lambda r: float(r["percent"]), reverse=True)
        ],
        "Detailed orbital contribution near Fermi level",
    )
    plot_bar(
        ctx.output_dir / "pdos_orbital_bar2",
        [
            (f"{row['element']}_{row['orbital_group']}", float(row["percent"]))
            for row in result.summary_rows
        ],
        "Shell contribution near Fermi level",
    )


def projection_markdown(ctx: Context, result: PdosResult) -> str:
    pdos_files = sorted(path.name for path in ctx.input_dir.glob("PDOS_*.dat"))
    if ctx.backend != "VASP":
        pdos_files = sorted(path.name for path in ctx.input_dir.glob("*.pdos.pdos_atm*"))
    if not pdos_files and ctx.backend == "EPW":
        pdos_files = sorted(path.name for path in workspace.backend_input_dir("QE").glob("*.pdos.pdos_atm*"))
    lines = [
        "# Projection Suggestion",
        "",
        "## 1. 读取情况",
        "",
        f"- Backend: {ctx.backend}",
        f"- Input directory: {ctx.input_dir}",
        f"- Structure file: {result.structure.structure_file or 'not found'}",
        f"- PDOS files: {len(pdos_files)} files",
        f"- Fermi energy: {fmt_fermi(result.fermi_energy)} eV",
        f"- Fermi interval: [{result.e_min:.6f}, {result.e_max:.6f}] eV",
        f"- Spin mode: {result.spin_mode}",
        f"- SOC: {result.soc}",
        f"- Formula: {result.structure.formula}",
        f"- Elements and atom numbers: {result.structure.elements or 'unknown'}",
        "",
        "## 2. PDOS contribution near Fermi level",
        "",
        "| Rank | Element | Orbital | Weight | Percent |",
        "|---:|---|---|---:|---:|",
    ]
    for idx, row in enumerate(result.summary_rows, start=1):
        lines.append(
            f"| {idx} | {row['element']} | {row['orbital_group']} | "
            f"{float(row['integral_weight']):.6g} | {float(row['percent']):.2f}% |"
        )
    if not result.summary_rows:
        lines.append("| 1 | - | - | 0 | 0.00% |")
        lines += [
            "",
            "Warning: The total PDOS weight in the selected Fermi interval is zero.",
            "Please check the PDOS files or fermi_interval setting.",
        ]
    lines += ["", "## 3. Recommended projections", ""]
    for rec in result.recommendations:
        lines += [
            f"### Recommended projections {rec['index']}: {rec['title']}",
            "",
            f"Fermi_energy = {fmt_fermi(result.fermi_energy)} eV",
            "",
            f"num_wann = {rec['num_wann']}",
            "",
            "begin projections",
            *[str(item) for item in rec["projections"]],
            "end projections",
            "",
            "Reason:",
            f"- {rec['reason']}",
            f"- 当前累计贡献约为 {100.0 * float(rec['cumulative']):.2f}%。",
            "",
        ]
    if len(result.recommendations) < 3:
        lines += [
            "### Recommended projections 3: Orbital-resolved recommendation",
            "",
            "No strong orbital anisotropy was detected under the current threshold.",
            "",
        ]
    lines += [
        "## 4. Notes",
        "",
        "- The fermi_interval is defined relative to the Fermi energy.",
        "- All plotted energies are absolute energies and are not shifted by Fermi energy.",
        "- The recommended num_wann is counted from orbital degeneracy.",
        "- For SOC or noncollinear calculations, please manually check whether spinor-related settings require further adjustment.",
    ]
    if result.soc:
        lines.append("- 注意：检测到 SOC 或 noncollinear 计算。当前 num_wann 按轨道数统计，请自行检查是否需要考虑 spinor 自由度。")
    if ctx.warnings:
        lines += ["", "## Warnings", ""]
        lines += [f"- {warning}" for warning in ctx.warnings]
    return "\n".join(lines) + "\n"


def setup_plot() -> None:
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["axes.unicode_minus"] = False


def save_fig(path_base: Path, dpi: int = 300) -> None:
    path_base.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path_base.with_suffix(".png"), dpi=dpi, facecolor="white")
    plt.savefig(path_base.with_suffix(".pdf"), dpi=dpi, facecolor="white")
    plt.close()


def plot_pdos(ctx: Context, result: PdosResult) -> None:
    setup_plot()
    plt.figure(figsize=(7.2, 4.6))
    for label, (energies, values) in result.energy_curves.items():
        plt.plot(energies, values, label=label, lw=1.2)
    plt.axvline(result.fermi_energy, color="k", ls="--", lw=1.0, label="Fermi energy")
    plt.xlabel("Energy (eV)")
    plt.ylabel("PDOS")
    plt.title("DFT PDOS")
    if result.energy_curves:
        plt.legend(fontsize=8)
    save_fig(ctx.output_dir / "dft_pdos_plot")


def plot_bar(path_base: Path, data: Sequence[Tuple[str, float]], title: str) -> None:
    setup_plot()
    data = [(label, value) for label, value in data if value > 0]
    if not data:
        data = [("none", 0.0)]
    labels, values = zip(*data[:40])
    width = max(7.2, min(16.0, 0.35 * len(labels) + 5.0))
    plt.figure(figsize=(width, 4.8))
    bars = plt.bar(range(len(labels)), values, color="#4C78A8")
    plt.ylabel("Percent (%)")
    plt.title(title)
    plt.xticks(range(len(labels)), labels, rotation=55, ha="right", fontsize=8)
    for bar, value in zip(bars, values):
        plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{value:.1f}%", ha="center", va="bottom", fontsize=7)
    save_fig(path_base)


def parse_projection_file(path: Path) -> Tuple[float, List[Dict[str, object]]]:
    if not path.exists():
        raise FileNotFoundError("projection_suggestion.md not found and function 1 failed.")
    text = path.read_text(encoding="utf-8")
    fermi_match = re.search(r"Fermi_energy\s*=\s*([-+0-9. Ee]+)", text)
    if not fermi_match:
        fermi_match = re.search(r"Fermi energy:\s*([-+0-9. Ee]+)", text)
    fermi = float(fermi_match.group(1).split()[0]) if fermi_match else math.nan
    recs: List[Dict[str, object]] = []
    for block in re.split(r"### Recommended projections ", text)[1:]:
        head = block.splitlines()[0].strip()
        index = int(re.match(r"(\d+)", head).group(1)) if re.match(r"(\d+)", head) else len(recs) + 1
        num_match = re.search(r"num_wann\s*=\s*(\d+)", block)
        proj_match = re.search(r"begin projections\s*(.*?)\s*end projections", block, re.S)
        if not num_match or not proj_match:
            continue
        projections = [line.strip() for line in proj_match.group(1).splitlines() if line.strip()]
        recs.append({"index": index, "num_wann": int(num_match.group(1)), "projections": projections})
    return fermi, recs


def manual_num_band_recommendations(ctx: Context) -> List[Dict[str, object]]:
    value = ctx.cfg.get("num_band", "auto").strip()
    if value.lower() == "auto":
        raise ValueError(
            "功能二未找到 projection_suggestion.md，且 input.dat 中 num_band=auto。"
            "请先运行功能一，或在功能二区域手动设置 num_band 为整数。"
        )
    try:
        num_band = int(float(value))
    except ValueError as exc:
        raise ValueError(f"input.dat 中 num_band 无法解析为整数：{value}") from exc
    if num_band <= 0:
        raise ValueError("input.dat 中 num_band 必须为正整数。")
    ctx.warnings.append("未读取到功能一 projection_suggestion.md，功能二使用 input.dat 中的 num_band 手动生成能带窗口建议。")
    return [
        {
            "index": 1,
            "num_wann": num_band,
            "projections": [f"manual_num_band = {num_band}"],
        }
    ]


def read_band(ctx: Context, fermi: float) -> BandData:
    if ctx.backend == "VASP":
        path = resolve_auto_file(ctx, ctx.cfg.get("vasp_band_file", "auto"), ["BAND.dat", "Band.dat", "band.dat"])
        bands, kdist = parse_band_file(path)
        labels = parse_klabels(resolve_auto_file(ctx, ctx.cfg.get("vasp_klabels_file", "auto"), ["KLABELS"], required=False))
    else:
        directory = ctx.input_dir
        path = resolve_auto_file(ctx, ctx.cfg.get("qe_band_gnu", "auto"), ["band.dat.gnu"], required=False)
        if path is None and ctx.backend == "EPW":
            fallback = workspace.backend_input_dir("QE")
            path = find_first(fallback, ["band.dat.gnu"])
            directory = fallback
            ctx.warnings.append("EPW 目录未找到 band.dat.gnu，已回退读取 QE 输入目录。")
        if path is None:
            raise FileNotFoundError("缺少 band.dat.gnu。")
        bands, kdist = parse_band_file(path)
        labels = parse_klabels(find_first(directory, ["band.out", "bands.out", "fit.in"]))
    band_array = np.array(bands, dtype=float)
    band_array = convert_to_absolute_energy(ctx.backend, "dft", band_array, fermi)
    bands = band_array.tolist()
    return BandData(path=path, k_distances=kdist, bands=bands, high_symmetry=labels)


def read_eigenval_energy_range(path: Optional[Path]) -> Optional[Tuple[float, float]]:
    if not path or not path.exists():
        return None
    values: List[float] = []
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        parts = raw.split()
        if len(parts) < 3:
            continue
        try:
            band_index = int(parts[0])
            energy = float(parts[1].replace("D", "E"))
        except ValueError:
            continue
        if band_index >= 1:
            values.append(energy)
    if not values:
        return None
    return min(values), max(values)


def detect_vasp_relative_band(ctx: Context, flat: np.ndarray, fermi: float, fallback: bool) -> bool:
    eigenval = resolve_auto_file(ctx, ctx.cfg.get("vasp_eigenval_file", "auto"), ["EIGENVAL"], required=False)
    eigen_range = read_eigenval_energy_range(eigenval)
    if not eigen_range:
        return fallback
    raw_range = (float(np.nanmin(flat)), float(np.nanmax(flat)))
    shifted_range = (raw_range[0] + fermi, raw_range[1] + fermi)
    raw_error = abs(raw_range[0] - eigen_range[0]) + abs(raw_range[1] - eigen_range[1])
    shifted_error = abs(shifted_range[0] - eigen_range[0]) + abs(shifted_range[1] - eigen_range[1])
    if shifted_error + 1e-4 < raw_error:
        return True
    if raw_error + 1e-4 < shifted_error:
        return False
    return fallback


def resolve_auto_file(ctx: Context, value: str, names: Sequence[str], required: bool = True) -> Optional[Path]:
    if value and value.lower() != "auto":
        path = resolve_path(value, ctx.input_dir)
        if path.exists():
            return path
    path = find_first(ctx.input_dir, names)
    if required and path is None:
        raise FileNotFoundError("缺少文件：%s" % " / ".join(names))
    return path


def parse_band_file(path: Path) -> Tuple[List[List[float]], List[float]]:
    blocks: List[List[Tuple[float, float]]] = []
    current: List[Tuple[float, float]] = []
    multicol_rows: List[List[float]] = []
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.strip()
        if not line:
            if current:
                blocks.append(current)
                current = []
            continue
        if line.startswith("#"):
            if "Band-Index" in line and current:
                blocks.append(current)
                current = []
            continue
        parts = line.split()
        try:
            nums = [float(part.replace("D", "E")) for part in parts]
        except ValueError:
            continue
        if len(nums) == 2:
            current.append((nums[0], nums[1]))
        elif len(nums) > 2:
            multicol_rows.append(nums)
    if current:
        blocks.append(current)
    if multicol_rows:
        arr = np.array(multicol_rows, dtype=float)
        kdist = arr[:, 0].tolist()
        return [arr[:, i].tolist() for i in range(1, arr.shape[1])], kdist
    if not blocks:
        raise ValueError(f"能带文件无法解析：{path}")
    kdist = [point[0] for point in blocks[0]]
    return [[point[1] for point in block] for block in blocks], kdist


def parse_klabels(path: Optional[Path]) -> List[Tuple[float, str]]:
    if not path or not path.exists():
        return []
    labels: List[Tuple[float, str]] = []
    highsym_count = 0
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "high-symmetry point" in line and "x coordinate" in line:
            match = re.search(r"x coordinate\s+(%s)" % FLOAT_RE, line)
            if match:
                highsym_count += 1
                labels.append((float(match.group(1)), f"K{highsym_count}"))
            continue
        parts = line.split()
        if len(parts) != 2:
            continue
        if parts[0].lower() in {"k-label", "label", "coordinate"}:
            continue
        candidates = [(parts[0], parts[-1]), (parts[-1], parts[0])]
        for label, coord in candidates:
            try:
                x = float(coord)
            except ValueError:
                continue
            if not re.match(r"^[A-Za-zΓ|_][A-Za-z0-9Γ|_.-]*$", label):
                continue
            label = label.replace("GAMMA", "Γ").replace("Gamma", "Γ")
            labels.append((x, label))
            break
    unique: List[Tuple[float, str]] = []
    for x, label in labels:
        if not any(abs(x - old_x) < 1e-5 and label == old_label for old_x, old_label in unique):
            unique.append((x, label))
    return unique


def band_stats(ctx: Context, data: BandData, fermi: float) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for idx, band in enumerate(data.bands, start=1):
        emin, emax = min(band), max(band)
        rows.append(
            {
                "backend": ctx.backend,
                "spin_channel": "none",
                "band_index": idx,
                "energy_min_eV": emin,
                "energy_max_eV": emax,
                "band_width_eV": emax - emin,
                "cross_fermi": str(emin <= fermi <= emax).lower(),
            }
        )
    return rows


def find_gap_below(stats: Sequence[Dict[str, object]], fermi: float, threshold: float) -> Optional[Tuple[int, float]]:
    best: Optional[Tuple[int, float, float]] = None
    ordered = sorted(stats, key=lambda row: int(row["band_index"]))
    for left, right in zip(ordered[:-1], ordered[1:]):
        left_max = float(left["energy_max_eV"])
        right_min = float(right["energy_min_eV"])
        if left_max < right_min < fermi and right_min - left_max >= threshold:
            gap = 0.5 * (left_max + right_min)
            distance = fermi - gap
            candidate = (int(right["band_index"]), gap, distance)
            if best is None or distance < best[2]:
                best = candidate
    return (best[0], best[1]) if best else None


def window_for_recommendation(
    ctx: Context,
    stats: Sequence[Dict[str, object]],
    rec: Dict[str, object],
    fermi: float,
) -> Dict[str, object]:
    nbands = len(stats)
    m = max(1, min(int(rec["num_wann"]), nbands))
    extra = max(0, int(float(ctx.cfg.get("extra_bands_for_diswin", "8"))))
    force_rule = ctx.cfg.get("force_rule", "auto").lower()
    threshold = float(ctx.cfg.get("min_gap_threshold", "0.05"))
    cross = [int(row["band_index"]) for row in stats if str(row["cross_fermi"]) == "true"]
    gap = find_gap_below(stats, fermi, threshold)
    use_rule1 = gap is not None and force_rule != "rule2"
    if force_rule == "rule1" and gap is None:
        ctx.warnings.append("force_rule=rule1，但未找到满足阈值的费米能级下方空带区域，已回退规则二。")
        use_rule1 = False
    if use_rule1 and gap:
        start, egap = gap
        end = min(nbands, start + m - 1)
        dis_froz_max_idx = min(nbands, start + m)
        dis_win_max_idx = min(nbands, start + m + extra)
        covered = set(range(start, end + 1))
        warnings = []
        if any(index not in covered for index in cross):
            warnings.append("The current num_wann does not include all bands crossing the Fermi level.")
        return {
            "projection_index": rec["index"],
            "projections": rec["projections"],
            "num_wann": rec["num_wann"],
            "num_bands": nbands,
            "fit_start": start,
            "fit_end": end,
            "dis_win_min": egap,
            "dis_win_max": float(stats[dis_win_max_idx - 1]["energy_max_eV"]),
            "dis_froz_min": egap,
            "dis_froz_max": float(stats[dis_froz_max_idx - 1]["energy_min_eV"]),
            "rule": "规则一",
            "warnings": warnings,
            "reason": "费米能级下方存在满足阈值的空带区域。",
        }
    above_count = max(1, int(round(m * float(ctx.cfg.get("above_fermi_ratio", "0.30")))))
    first_above = next((int(row["band_index"]) for row in stats if float(row["energy_min_eV"]) > fermi), nbands)
    start = max(1, min(nbands - m + 1, first_above - (m - above_count)))
    if cross:
        start = min(start, min(cross))
        start = max(1, min(start, max(cross) - m + 1))
    end = min(nbands, start + m - 1)
    lo_idx = max(1, start - extra)
    hi_idx = min(nbands, end + extra)
    froz_hi_idx = min(nbands, end + 1)
    warnings = []
    if len(cross) > int(rec["num_wann"]):
        warnings.append("Critical Warning: The number of bands crossing the Fermi level is larger than num_wann.")
    if gap is None:
        warnings.append("未找到合适空带区域，已使用规则二。")
    return {
        "projection_index": rec["index"],
        "projections": rec["projections"],
        "num_wann": rec["num_wann"],
        "num_bands": nbands,
        "fit_start": start,
        "fit_end": end,
        "dis_win_min": float(stats[lo_idx - 1]["energy_min_eV"]),
        "dis_win_max": float(stats[hi_idx - 1]["energy_max_eV"]),
        "dis_froz_min": float(stats[start - 1]["energy_max_eV"]),
        "dis_froz_max": float(stats[froz_hi_idx - 1]["energy_min_eV"]),
        "rule": "规则二",
        "warnings": warnings,
        "reason": "围绕费米能级选择目标能带并保留一定比例的费米能级以上能带。",
    }


def run_band(ctx: Context, projection_result: Optional[PdosResult] = None) -> None:
    projection_file = ctx.output_dir / "projection_suggestion.md"
    if projection_result is not None:
        fermi = projection_result.fermi_energy
        recs = [
            {
                "index": rec["index"],
                "num_wann": rec["num_wann"],
                "projections": rec["projections"],
            }
            for rec in projection_result.recommendations
        ]
    elif projection_file.exists():
        fermi, recs = parse_projection_file(projection_file)
    else:
        fermi = read_fermi(ctx)
        recs = manual_num_band_recommendations(ctx)
    data = read_band(ctx, fermi)
    stats = band_stats(ctx, data, fermi)
    write_csv(
        ctx.output_dir / "dft_band.csv",
        stats,
        ["backend", "spin_channel", "band_index", "energy_min_eV", "energy_max_eV", "band_width_eV", "cross_fermi"],
    )
    windows = [window_for_recommendation(ctx, stats, rec, fermi) for rec in recs]
    write_text(ctx.output_dir / "band_suggestion.md", band_markdown(ctx, data, stats, windows, fermi))
    plot_band(ctx.output_dir / "dft_band_plot", data, fermi, None, ctx.cfg)
    plot_band_windows(ctx.output_dir / "dft_band_with_window", data, fermi, windows, ctx.cfg)


def band_markdown(ctx: Context, data: BandData, stats: Sequence[Dict[str, object]], windows: Sequence[Dict[str, object]], fermi: float) -> str:
    all_warnings = list(ctx.warnings)
    for window in windows:
        all_warnings.extend(str(item) for item in window.get("warnings", []))
    structure = parse_vasp_structure(ctx.input_dir) if ctx.backend == "VASP" else parse_qe_structure(ctx.input_dir)
    lines = [
        "# Band Fitting Suggestion",
        "",
        "## 1. 读取情况",
        "",
        f"- 后端：{ctx.backend}",
        f"- 输入目录：{ctx.input_dir}",
        f"- DFT band 文件：{data.path}",
        f"- scf.out 或 OUTCAR：{find_first(ctx.input_dir, ['scf.out', 'OUTCAR']) or 'not found'}",
        "- 能量单位：eV",
        f"- Fermi energy：{fmt_fermi(fermi)} eV",
        f"- 能带数量：{len(data.bands)}",
        f"- k 点数量：{len(data.k_distances)}",
        f"- Spin mode：{ctx.cfg.get('spin_mode', 'auto')}",
        f"- 元素和原子数：{structure.elements or 'unknown'}",
        "",
        "## 2. Wannier90 初始建议",
        "",
    ]
    titles = ["主要建议", "建议二", "建议三"]
    for i, window in enumerate(windows):
        title = titles[i] if i < len(titles) else f"建议 {i + 1}"
        lines += [
            f"### {title}",
            "",
            f"- projections：{', '.join(str(x) for x in window['projections'])}",
            f"- num_wann：{window['num_wann']}",
            f"- num_bands：{window['num_bands']}",
            f"- fit bands：{window['fit_start']} - {window['fit_end']}",
            f"- dis_win_min：{trunc1(float(window['dis_win_min']))}",
            f"- dis_win_max：{trunc1(float(window['dis_win_max']))}",
            f"- dis_froz_min：{trunc1(float(window['dis_froz_min']))}",
            f"- dis_froz_max：{trunc1(float(window['dis_froz_max']))}",
            f"- 使用规则：{window['rule']}",
            "",
            "Reason:",
            f"- {window['reason']}",
            "- 窗口值为未扣除 Fermi energy 的绝对能量。",
            "",
        ]
    lines += ["## 3. Warnings", ""]
    if all_warnings:
        lines += [f"- {warning}" for warning in dict.fromkeys(all_warnings)]
    else:
        lines.append("- 未发现需要特别提示的警告。")
    lines += [
        "",
        "## 4. Notes",
        "",
        "- 以上窗口值默认是未扣除 Fermi energy 的绝对能量。",
        "- 如果 Wannier90 流程使用相对费米能级的能量窗口，需要手动减去 E_F。",
        "- 本程序只给出初始建议，最终仍需结合 Wannier 拟合能带图人工检查。",
    ]
    return "\n".join(lines) + "\n"


def band_plot_options(cfg: Dict[str, str], bands: Sequence[Sequence[float]]) -> Dict[str, object]:
    flat = [float(value) for band in bands for value in band]
    auto_ymin = min(flat) if flat else -1.0
    auto_ymax = max(flat) if flat else 1.0
    return {
        "ymin": cfg_float(cfg, "band_plot_ymin", auto_ymin),
        "ymax": cfg_float(cfg, "band_plot_ymax", auto_ymax),
        "width": cfg_float(cfg, "band_plot_width", 7.2),
        "height": cfg_float(cfg, "band_plot_height", 4.8),
        "linewidth": cfg_float(cfg, "band_plot_linewidth", 0.65),
        "dpi": int(cfg_float(cfg, "band_plot_dpi", 300)),
        "show_klabels": cfg_bool(cfg, "band_plot_show_klabels", True),
        "show_grid": cfg_bool(cfg, "band_plot_show_grid", True),
    }


def plot_band(path_base: Path, data: BandData, fermi: float, window: Optional[Dict[str, object]], cfg: Optional[Dict[str, str]] = None) -> None:
    setup_plot()
    opts = band_plot_options(cfg or {}, data.bands)
    plt.figure(figsize=(float(opts["width"]), float(opts["height"])))
    for band in data.bands:
        plt.plot(data.k_distances[: len(band)], band, color="#222222", lw=float(opts["linewidth"]))
    plt.axhline(fermi, color="red", ls="--", lw=1.0, label="Fermi energy")
    if window:
        for key, color in [("dis_win_min", "#4C78A8"), ("dis_win_max", "#4C78A8"), ("dis_froz_min", "#F58518"), ("dis_froz_max", "#F58518")]:
            plt.axhline(float(window[key]), color=color, ls=":", lw=1.0)
            plt.text(data.k_distances[-1], float(window[key]), f" {key}={trunc1(float(window[key]))}", fontsize=7, va="center")
    if data.high_symmetry and bool(opts["show_klabels"]):
        xs, labels = zip(*data.high_symmetry)
        plt.xticks(xs, labels)
        if bool(opts["show_grid"]):
            for x in xs:
                plt.axvline(x, color="#BBBBBB", lw=0.5)
    plt.ylim(float(opts["ymin"]), float(opts["ymax"]))
    plt.xlabel("k-distance")
    plt.ylabel("Energy (eV)")
    plt.title("DFT band")
    save_fig(path_base, int(opts["dpi"]))


def plot_band_windows(path_base: Path, data: BandData, fermi: float, windows: Sequence[Dict[str, object]], cfg: Optional[Dict[str, str]] = None) -> None:
    if not windows:
        plot_band(path_base, data, fermi, None, cfg)
        return
    setup_plot()
    opts = band_plot_options(cfg or {}, data.bands)
    n = len(windows)
    cols = 2 if n > 1 else 1
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(float(opts["width"]) * cols, float(opts["height"]) * rows), squeeze=False)
    for ax, window in zip(axes.ravel(), windows):
        for band in data.bands:
            ax.plot(data.k_distances[: len(band)], band, color="#222222", lw=float(opts["linewidth"]))
        ax.axhline(fermi, color="red", ls="--", lw=0.9)
        ax.axhspan(float(window["dis_win_min"]), float(window["dis_win_max"]), color="#4C78A8", alpha=0.12)
        ax.axhspan(float(window["dis_froz_min"]), float(window["dis_froz_max"]), color="#F58518", alpha=0.18)
        ax.set_title(f"{', '.join(str(x) for x in window['projections'])}; num_wann={window['num_wann']}", fontsize=10)
        ax.set_xlabel("k-distance")
        ax.set_ylabel("Energy (eV)")
        ax.set_ylim(float(opts["ymin"]), float(opts["ymax"]))
        ax.text(0.01, 0.98, f"win [{trunc1(float(window['dis_win_min']))}, {trunc1(float(window['dis_win_max']))}]\nfroz [{trunc1(float(window['dis_froz_min']))}, {trunc1(float(window['dis_froz_max']))}]", transform=ax.transAxes, va="top", fontsize=8)
        if data.high_symmetry and bool(opts["show_klabels"]):
            xs, labels = zip(*data.high_symmetry)
            ax.set_xticks(xs, labels)
            if bool(opts["show_grid"]):
                for x in xs:
                    ax.axvline(x, color="#BBBBBB", lw=0.5)
    for ax in axes.ravel()[len(windows) :]:
        ax.axis("off")
    plt.tight_layout()
    plt.savefig(path_base.with_suffix(".png"), dpi=int(opts["dpi"]), facecolor="white")
    plt.savefig(path_base.with_suffix(".pdf"), dpi=int(opts["dpi"]), facecolor="white")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Auxiliary Wannier Band Fitting")
    parser.add_argument("--task", choices=["projection", "band", "all"], default=None)
    parser.add_argument("--backend", choices=["VASP", "QE", "EPW", "vasp", "qe", "epw"], default=None)
    return parser.parse_args()


def main() -> int:
    try:
        ctx = make_context(parse_args())
        projection_result = None
        if ctx.task in {"projection", "all"}:
            projection_result = run_projection(ctx)
        if ctx.task in {"band", "all"}:
            run_band(ctx, projection_result)
        if ctx.warnings:
            write_text(ctx.output_dir / "warnings.log", "\n".join(dict.fromkeys(ctx.warnings)) + "\n")
        print(f"Backend: {ctx.backend}")
        print(f"Task: {ctx.task}")
        print(f"Input directory: {ctx.input_dir}")
        print(f"Output directory: {ctx.output_dir}")
        print("Done.")
        return 0
    except Exception as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
