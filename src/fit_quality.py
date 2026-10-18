#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""功能三：DFT 与 Wannier 拟合能带质量分析。"""

from __future__ import annotations

import argparse
import configparser
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:
    print(f"Error: matplotlib 无法导入，不能生成图像：{exc}", file=sys.stderr)
    raise SystemExit(1)


APP_ROOT = Path(__file__).resolve().parents[1]
from .default_inputs import DEFAULT_INPUTS_DAT
from . import workspace
FLOAT_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
DEFAULT_INPUTS = DEFAULT_INPUTS_DAT


@dataclass
class BandSet:
    path: Path
    x: np.ndarray
    bands: np.ndarray


@dataclass
class FitResult:
    backend: str
    dft: BandSet
    wannier: BandSet
    align_factor: float
    align_method: str
    y_shift: float
    fermi_energy: float
    fermi_source: str
    energy_note: str
    metric_energy_min: float
    metric_energy_max: float
    near_fermi_window: float
    quality_basis: str
    matched_dft_indices: List[int]
    matched_dft: np.ndarray
    matched_wannier: np.ndarray
    errors: np.ndarray
    overall_rmse: float
    overall_max: float
    near_rmse: float
    near_max: float
    dft_cross_count: int
    wannier_cross_count: int
    crossing_consistent: bool
    per_band: List[Dict[str, object]]
    quality: str
    suggestions: List[str]
    warnings: List[str]


def bool_cfg(config: configparser.ConfigParser, section: str, key: str, default: bool) -> bool:
    if not config.has_option(section, key):
        return default
    return config.get(section, key).split("#", 1)[0].strip().lower() in {"1", "true", "yes", "on", "y"}


def float_cfg(config: configparser.ConfigParser, section: str, key: str, default: float) -> float:
    if not config.has_option(section, key):
        return default
    value = config.get(section, key).split("#", 1)[0].strip()
    if value.lower() == "auto":
        return default
    try:
        return float(value)
    except ValueError:
        return default


def optional_float_cfg(config: configparser.ConfigParser, section: str, key: str) -> Optional[float]:
    if not config.has_option(section, key):
        return None
    value = config.get(section, key).split("#", 1)[0].strip()
    if value.lower() == "auto" or not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def text_cfg(config: configparser.ConfigParser, section: str, key: str, default: str) -> str:
    if not config.has_option(section, key):
        return default
    return config.get(section, key).split("#", 1)[0].strip()


def ensure_input_dat(path: Path) -> None:
    if not path.exists():
        path.write_text(DEFAULT_INPUTS, encoding="utf-8", newline="\n")
        return
    text = path.read_text(encoding="utf-8", errors="ignore")
    additions: List[str] = []
    if "[read dat]" not in text:
        additions.append(DEFAULT_INPUTS.split("[fit_check]", 1)[0].rstrip())
    if "[fit_check]" not in text:
        additions.append("[fit_check]" + DEFAULT_INPUTS.split("[fit_check]", 1)[1].rstrip())
    if additions:
        text = text.rstrip() + "\n\n" + "\n\n".join(additions) + "\n"
    if additions:
        path.write_text(text, encoding="utf-8")


def load_config() -> configparser.ConfigParser:
    workspace.require_workspace()
    path = workspace.config_path()
    parser = configparser.ConfigParser(inline_comment_prefixes=("#",))
    parser.optionxform = str
    text = path.read_text(encoding="utf-8", errors="ignore")
    section_start = text.find("[read dat]")
    if section_start >= 0:
        parser.read_string(text[section_start:])
    else:
        parser.read_string(DEFAULT_INPUTS)
    return parser


def candidate_dirs(backend: str) -> List[Path]:
    backend = backend.upper()
    dirs = [workspace.backend_input_dir(backend)]
    if backend == "EPW":
        dirs.append(workspace.backend_input_dir("QE"))
    return [path for path in dirs if path.exists()]


def resolve_project_file(path_text: str) -> Path:
    path = Path(path_text).expanduser()
    if not path.is_absolute():
        path = workspace.USER_WORK / path
    resolved = path.resolve()
    try:
        resolved.relative_to(workspace.USER_WORK.resolve())
    except ValueError as exc:
        raise ValueError(f"不允许读取 user_work 目录以外的文件：{resolved}") from exc
    return resolved


def find_dft_band(backend: str, user_path: Optional[str]) -> Path:
    if user_path and user_path.lower() != "auto":
        path = resolve_project_file(user_path)
        if path.exists():
            return path
        raise FileNotFoundError(f"DFT 能带文件不存在：{path}")
    names = ["BAND.dat"] if backend.upper() == "VASP" else ["band.dat.gnu"]
    for directory in candidate_dirs(backend):
        for name in names:
            path = directory / name
            if path.exists():
                return path.resolve()
    raise FileNotFoundError(f"未找到 {backend} DFT 本征能带文件：{', '.join(names)}")


def find_wannier_band(backend: str, user_path: str) -> Path:
    if user_path and user_path.lower() != "auto":
        path = resolve_project_file(user_path)
        if path.exists():
            return path
        raise FileNotFoundError(f"Wannier 拟合能带文件不存在：{path}")
    matches: List[Path] = []
    for directory in candidate_dirs(backend):
        matches.extend(path for path in directory.glob("*_band.dat") if path.is_file())
    if not matches:
        raise FileNotFoundError("未找到 Wannier 拟合能带文件：*_band.dat")
    matches = sorted(set(matches), key=lambda path: path.stat().st_mtime, reverse=True)
    if len(matches) > 1:
        print("发现多个 Wannier 拟合能带候选文件：")
        for path in matches:
            print(f"  - {path}")
        print(f"当前自动选择最近修改文件：{matches[0]}")
    return matches[0].resolve()


def read_fermi_energy(backend: str, dft_path: Path) -> float:
    search_dirs = [dft_path.parent, *candidate_dirs(backend)]
    if backend.upper() == "VASP":
        for directory in search_dirs:
            path = directory / "OUTCAR"
            if not path.exists():
                continue
            value = None
            for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
                if "E-fermi" in line:
                    try:
                        value = float(line.split(":", 1)[1].split()[0])
                    except (IndexError, ValueError):
                        pass
            if value is not None:
                return value
    else:
        for directory in search_dirs:
            path = directory / "scf.out"
            if not path.exists():
                continue
            for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
                if "Fermi energy is" in line:
                    try:
                        return float(line.split("Fermi energy is", 1)[1].split()[0])
                    except (IndexError, ValueError):
                        pass
    return 0.0


def parse_fermi_from_markdown(path: Path) -> Optional[float]:
    if not path.exists():
        return None
    text = path.read_text(encoding="utf-8", errors="ignore")
    patterns = [
        r"Fermi_energy\s*=\s*(%s)" % FLOAT_RE,
        r"Fermi\s+energy\s*[:：]\s*(%s)" % FLOAT_RE,
    ]
    for pattern in patterns:
        match = re.search(pattern, text, re.I)
        if match:
            try:
                return float(match.group(1).replace("D", "E"))
            except ValueError:
                continue
    return None


def resolve_fit_fermi_energy(config: configparser.ConfigParser, backend: str, dft_path: Path) -> Tuple[float, str]:
    value_text = text_cfg(config, "fit_check", "fit_fermi_energy", "auto")
    if value_text.lower() != "auto":
        try:
            return float(value_text.replace("D", "E")), "manual fit_check.fit_fermi_energy"
        except ValueError as exc:
            raise ValueError(f"fit_fermi_energy 手动值无法解析：{value_text}") from exc

    output_dir = workspace.backend_output_dir(backend.upper())
    candidates = [
        output_dir / "projection_suggestion.md",
        output_dir / "band_suggestion.md",
    ]
    existing = sorted((path for path in candidates if path.exists()), key=lambda item: item.stat().st_mtime, reverse=True)
    for path in existing:
        fermi = parse_fermi_from_markdown(path)
        if fermi is not None:
            return fermi, f"auto from {path.name}"

    fallback = read_fermi_energy(backend, dft_path)
    return fallback, "auto fallback from OUTCAR/scf.out because function 1/2 output was not found"


def read_vasp_formula(directory: Path) -> Tuple[str, List[str]]:
    path = directory / "POSCAR"
    if not path.exists():
        path = directory / "CONTCAR"
    if not path.exists():
        return "unknown", []
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    if len(lines) < 7:
        return "unknown", []
    elements = lines[5].split()
    counts = lines[6].split()
    formula = []
    clean_elements: List[str] = []
    for element, count_text in zip(elements, counts):
        try:
            count = int(float(count_text))
        except ValueError:
            continue
        clean_elements.append(element)
        formula.append(element if count == 1 else f"{element}{count}")
    return "".join(formula) if formula else "unknown", clean_elements


def build_input_warnings(backend: str, dft_path: Path, wannier_path: Path) -> List[str]:
    warnings: List[str] = []
    if backend.upper() == "VASP":
        formula, elements = read_vasp_formula(dft_path.parent)
        name = wannier_path.stem.lower()
        name_tokens = {token for token in re.split(r"[^a-z0-9]+", name) if token}
        element_match = any(element.lower() in name_tokens for element in elements)
        if formula != "unknown" and formula.lower() not in name and not element_match:
            warnings.append(
                f"VASP 输入结构为 {formula}，但当前 Wannier 文件名为 {wannier_path.name}，可能不是同一体系。"
            )
    return warnings


def split_two_column_blocks(path: Path) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    blocks: List[List[Tuple[float, float]]] = []
    current: List[Tuple[float, float]] = []
    last_x: Optional[float] = None
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    has_band_index = any("Band-Index" in line for line in lines)
    for raw in lines:
        line = raw.strip()
        if not line:
            if current and not has_band_index:
                blocks.append(current)
                current = []
                last_x = None
            continue
        if line.startswith("#"):
            if "Band-Index" in line and current:
                blocks.append(current)
                current = []
                last_x = None
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            x = float(parts[0].replace("D", "E"))
            y = float(parts[1].replace("D", "E"))
        except ValueError:
            continue
        if not has_band_index and last_x is not None and x < last_x - 1e-10 and current:
            blocks.append(current)
            current = []
        current.append((x, y))
        last_x = x
    if current:
        blocks.append(current)
    xs = [np.array([point[0] for point in block], dtype=float) for block in blocks if block]
    ys = [np.array([point[1] for point in block], dtype=float) for block in blocks if block]
    if xs:
        base_x = xs[0]
        aligned_xs: List[np.ndarray] = []
        aligned_ys: List[np.ndarray] = []
        for x_values, y_values in zip(xs, ys):
            if len(x_values) == len(base_x) and np.allclose(x_values[::-1], base_x, rtol=1e-6, atol=1e-8):
                aligned_xs.append(x_values[::-1])
                aligned_ys.append(y_values[::-1])
            else:
                aligned_xs.append(x_values)
                aligned_ys.append(y_values)
        xs, ys = aligned_xs, aligned_ys
    return xs, ys


def read_band_file(path: Path) -> BandSet:
    text = path.read_text(encoding="utf-8", errors="ignore")
    if "Band-Index" in text:
        xs, ys = split_two_column_blocks(path)
        if not ys:
            raise ValueError(f"无法解析二维能带数据：{path}")
        min_points = min(len(y) for y in ys)
        return BandSet(path=path, x=xs[0][:min_points], bands=np.array([y[:min_points] for y in ys], dtype=float))
    rows: List[List[float]] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        try:
            rows.append([float(part.replace("D", "E")) for part in line.split()])
        except ValueError:
            continue
    if not rows:
        raise ValueError(f"能带文件没有可读取数值：{path}")
    widths = {len(row) for row in rows}
    if max(widths) > 2 and len(widths) == 1:
        arr = np.array(rows, dtype=float)
        return BandSet(path=path, x=arr[:, 0], bands=arr[:, 1:].T)
    xs, ys = split_two_column_blocks(path)
    if not ys:
        raise ValueError(f"无法解析二维能带数据：{path}")
    min_points = min(len(y) for y in ys)
    return BandSet(path=path, x=xs[0][:min_points], bands=np.array([y[:min_points] for y in ys], dtype=float))


def convert_to_absolute_energy(
    code_type: str,
    file_role: str,
    energy: np.ndarray,
    efermi: Optional[float] = None,
) -> np.ndarray:
    """
    Convert all band energies to the absolute, non-EF-subtracted scale.

    Rule:
        QE DFT: already absolute, do nothing.
        QE Wannier: already absolute, do nothing.
        EPW DFT: already absolute, do nothing.
        EPW Wannier: already absolute, do nothing.
        VASP DFT BAND.dat: EF was subtracted by post-processing, so add EF back.
        VASP Wannier: already absolute, do nothing.
    """
    if code_type.lower() == "vasp" and file_role.lower() == "dft":
        if efermi is None:
            raise ValueError("VASP BAND.dat requires efermi to add back the Fermi energy.")
        return energy + efermi
    return energy


def normalize_fit_energies(
    backend: str,
    dft_bands: np.ndarray,
    wannier_bands: np.ndarray,
    fermi: float,
) -> Tuple[np.ndarray, np.ndarray, str]:
    """Convert DFT and Wannier bands to absolute, non-EF-subtracted energies."""
    backend_upper = backend.upper()
    dft_absolute = convert_to_absolute_energy(backend_upper, "dft", dft_bands, fermi)
    wannier_absolute = convert_to_absolute_energy(backend_upper, "wannier", wannier_bands, fermi)
    if backend_upper == "VASP":
        note = "VASP DFT BAND.dat is shifted by +E_F; VASP Wannier bands are used directly."
    else:
        note = f"{backend_upper} DFT and Wannier bands are used directly without any E_F subtraction."
    return dft_absolute, wannier_absolute, note


def align_x(dft: BandSet, wannier: BandSet, config: configparser.ConfigParser) -> Tuple[np.ndarray, float, str]:
    manual = text_cfg(config, "fit_check", "align_x", "auto")
    method = text_cfg(config, "fit_check", "align_x_method", "endpoint").lower()
    if manual.lower() != "auto":
        factor = float(manual)
        return wannier.x * factor, factor, "manual"
    dft_span = float(dft.x[-1] - dft.x[0]) if len(dft.x) > 1 else 1.0
    wannier_span = float(wannier.x[-1] - wannier.x[0]) if len(wannier.x) > 1 else 1.0
    factor = dft_span / wannier_span if abs(wannier_span) > 1e-12 else 1.0
    if method == "least_squares" and len(dft.x) == len(wannier.x):
        denom = float(np.dot(wannier.x, wannier.x))
        if abs(denom) > 1e-12:
            factor = float(np.dot(wannier.x, dft.x) / denom)
    return wannier.x * factor, factor, method


def interpolate_bands(source_x: np.ndarray, source_bands: np.ndarray, target_x: np.ndarray) -> np.ndarray:
    return np.array([np.interp(target_x, source_x, band) for band in source_bands], dtype=float)


def rmse(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    return float(np.sqrt(np.mean(values * values))) if values.size else math.nan


def max_abs(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    return float(np.max(np.abs(values))) if values.size else math.nan


def match_dft_bands(dft_bands: np.ndarray, wannier_bands: np.ndarray, energy_min: float, energy_max: float) -> Tuple[np.ndarray, np.ndarray, List[int]]:
    selected: List[np.ndarray] = []
    selected_wannier: List[np.ndarray] = []
    selected_indices: List[int] = []
    used: set[int] = set()
    for wannier_band in wannier_bands:
        best_idx: Optional[int] = None
        best_score = math.inf
        for idx, dft_band in enumerate(dft_bands):
            if idx in used:
                continue
            mask = ((dft_band >= energy_min) & (dft_band <= energy_max)) | ((wannier_band >= energy_min) & (wannier_band <= energy_max))
            if not np.any(mask):
                mask = np.ones_like(dft_band, dtype=bool)
            score = rmse(wannier_band[mask] - dft_band[mask])
            if score < best_score:
                best_score = score
                best_idx = idx
        if best_idx is not None:
            used.add(best_idx)
            selected.append(dft_bands[best_idx])
            selected_wannier.append(wannier_band)
            selected_indices.append(best_idx + 1)
    if not selected:
        raise ValueError("无法在 DFT 能带中匹配 Wannier 目标能带。")
    return np.array(selected), np.array(selected_wannier), selected_indices


def crossing_count(bands: np.ndarray, fermi: float, tolerance: float) -> int:
    return sum(1 for band in bands if np.nanmin(band) <= fermi + tolerance and np.nanmax(band) >= fermi - tolerance)


def classify_quality(
    overall_rmse: float,
    overall_max: float,
    near_rmse: float,
    near_max: float,
    crossing_consistent: bool,
    config: configparser.ConfigParser,
) -> Tuple[str, str, List[str]]:
    max_rmse = float_cfg(config, "fit_check", "max_allowed_rmse", 0.05)
    max_near = float_cfg(config, "fit_check", "max_allowed_near_fermi_error", 0.03)
    focus = text_cfg(config, "fit_check", "metric_focus", "fermi").lower()
    suggestions: List[str] = []

    if focus == "fermi":
        basis = "metric_focus=fermi: final quality is based on near-Fermi RMSE/max error and Fermi crossing consistency."
        near_ok = near_rmse < max_near and near_max < max_near
        if crossing_consistent and near_ok:
            if overall_rmse >= max_rmse:
                suggestions.append("全能量范围 RMSE 超过阈值，但 metric_focus=fermi 时不作为最终 bad 判据；如需全窗口严格检查，请设置 metric_focus=full。")
            return "good", basis, suggestions
        if not crossing_consistent:
            suggestions.append("费米面交叉数量不一致：请增加/调整 projection，或检查 Wannier 能带文件是否对应同一路径。")
        if not near_ok:
            suggestions.append("费米附近误差偏大：优先检查 dis_froz 窗口和费米面附近关键轨道。")
        if not crossing_consistent or near_rmse >= 2 * max_near or near_max >= 2 * max_near:
            return "bad", basis, suggestions
        return "warning", basis, suggestions

    basis = "metric_focus=full/selected: final quality is based on full metric-window RMSE plus near-Fermi checks."
    if crossing_consistent and overall_rmse < max_rmse and near_rmse < max_near and near_max < max_near:
        return "good", basis, suggestions
    if not crossing_consistent:
        suggestions.append("费米面交叉数量不一致：请增加/调整 projection，或检查 Wannier 能带文件是否对应同一路径。")
    if near_rmse >= max_near or near_max >= max_near:
        suggestions.append("费米附近误差偏大：优先检查 dis_froz 窗口和费米面附近关键轨道。")
    if overall_rmse >= max_rmse:
        suggestions.append("整体 RMSE 超过阈值：建议增大 dis_win 或加入更多杂化轨道。")
    if not crossing_consistent or overall_rmse >= 2 * max_rmse or near_rmse >= 2 * max_near or near_max >= 2 * max_near:
        return "bad", basis, suggestions
    return "warning", basis, suggestions


def run_fit_quality(backend: str, config: configparser.ConfigParser, dft_arg: Optional[str], wannier_arg: Optional[str]) -> FitResult:
    backend = backend.upper()
    read_key = {"VASP": "vasp_fit_band", "QE": "qe_fit_band", "EPW": "epw_fit_band"}[backend]
    dft_path = find_dft_band(backend, dft_arg)
    wannier_path = find_wannier_band(backend, wannier_arg or text_cfg(config, "read dat", read_key, "auto"))
    input_warnings = build_input_warnings(backend, dft_path, wannier_path)
    dft = read_band_file(dft_path)
    wannier = read_band_file(wannier_path)
    fermi, fermi_source = resolve_fit_fermi_energy(config, backend, dft_path)
    dft_bands, wannier_bands, energy_note = normalize_fit_energies(backend, dft.bands, wannier.bands, fermi)
    dft = BandSet(dft.path, dft.x, dft_bands)
    wannier = BandSet(wannier.path, wannier.x, wannier_bands)
    aligned_x, factor, method = align_x(dft, wannier, config)
    nband = min(dft.bands.shape[0], wannier.bands.shape[0])
    wannier_interp_all = interpolate_bands(aligned_x, wannier.bands[:nband], dft.x)
    finite_for_metric = np.concatenate([dft.bands.ravel(), wannier_interp_all.ravel()])
    finite_for_metric = finite_for_metric[np.isfinite(finite_for_metric)]
    metric_min = optional_float_cfg(config, "fit_check", "metric_energy_min")
    metric_max = optional_float_cfg(config, "fit_check", "metric_energy_max")
    if finite_for_metric.size:
        if metric_min is None:
            metric_min = float(np.nanmin(finite_for_metric))
        if metric_max is None:
            metric_max = float(np.nanmax(finite_for_metric))
    else:
        metric_min = -math.inf if metric_min is None else metric_min
        metric_max = math.inf if metric_max is None else metric_max
    dft_match, wannier_interp, matched_indices = match_dft_bands(dft.bands, wannier_interp_all, metric_min, metric_max)
    y_shift = 0.0
    if bool_cfg(config, "fit_check", "align_y", False):
        y_shift = float(np.nanmean(dft_match - wannier_interp))
        wannier_interp = wannier_interp + y_shift
    errors = wannier_interp - dft_match
    metric_mask = (dft_match >= metric_min) & (dft_match <= metric_max)
    near_window = max(float_cfg(config, "fit_check", "near_fermi_window", 2.0), 1e-12)
    near_mask = np.abs(dft_match - fermi) <= near_window
    overall_values = errors[metric_mask] if np.any(metric_mask) else errors
    near_values = errors[near_mask] if np.any(near_mask) else overall_values
    tolerance = float_cfg(config, "fit_check", "fermi_tolerance", 0.02)
    dft_cross = crossing_count(dft_match, fermi, tolerance)
    wannier_cross = crossing_count(wannier_interp, fermi, tolerance)
    crossing_consistent = dft_cross == wannier_cross if bool_cfg(config, "fit_check", "check_band_crossing_count", True) else True
    per_band = []
    for idx in range(errors.shape[0]):
        band_mask = metric_mask[idx] if idx < metric_mask.shape[0] else np.ones_like(errors[idx], dtype=bool)
        band_errors = errors[idx][band_mask] if np.any(band_mask) else errors[idx]
        per_band.append(
            {
                "band_index": matched_indices[idx] if idx < len(matched_indices) else idx + 1,
                "rmse": rmse(band_errors),
                "max_error": max_abs(band_errors),
                "dft_cross_fermi": bool(np.nanmin(dft_match[idx]) <= fermi + tolerance and np.nanmax(dft_match[idx]) >= fermi - tolerance),
                "wannier_cross_fermi": bool(np.nanmin(wannier_interp[idx]) <= fermi + tolerance and np.nanmax(wannier_interp[idx]) >= fermi - tolerance),
            }
        )
    overall = rmse(overall_values)
    overall_abs_max = max_abs(overall_values)
    near = rmse(near_values)
    near_abs_max = max_abs(near_values)
    quality, quality_basis, suggestions = classify_quality(overall, overall_abs_max, near, near_abs_max, crossing_consistent, config)
    suggestions.extend(input_warnings)
    return FitResult(
        backend=backend,
        dft=dft,
        wannier=wannier,
        align_factor=factor,
        align_method=method,
        y_shift=y_shift,
        fermi_energy=fermi,
        fermi_source=fermi_source,
        energy_note=energy_note,
        metric_energy_min=metric_min,
        metric_energy_max=metric_max,
        near_fermi_window=near_window,
        quality_basis=quality_basis,
        matched_dft_indices=matched_indices,
        matched_dft=dft_match,
        matched_wannier=wannier_interp,
        errors=errors,
        overall_rmse=overall,
        overall_max=overall_abs_max,
        near_rmse=near,
        near_max=near_abs_max,
        dft_cross_count=dft_cross,
        wannier_cross_count=wannier_cross,
        crossing_consistent=crossing_consistent,
        per_band=per_band,
        quality=quality,
        suggestions=suggestions,
        warnings=input_warnings,
    )


def save_current_figure(base: Path, config: configparser.ConfigParser) -> List[Path]:
    dpi = int(float_cfg(config, "fit_check", "plot_dpi", 300))
    outputs: List[Path] = []
    if bool_cfg(config, "fit_check", "plot_format_png", True):
        path = base.with_suffix(".png")
        plt.savefig(path, dpi=dpi, facecolor="white", bbox_inches="tight")
        outputs.append(path.resolve())
    if bool_cfg(config, "fit_check", "plot_format_pdf", True):
        path = base.with_suffix(".pdf")
        plt.savefig(path, dpi=dpi, facecolor="white", bbox_inches="tight")
        outputs.append(path.resolve())
    plt.close()
    return outputs


def plot_fit_bands(result: FitResult, config: configparser.ConfigParser, output_dir: Path) -> List[Path]:
    plt.figure(figsize=(7.5, 5.0))
    for band in result.matched_dft:
        plt.plot(result.dft.x, band, color="#222222", lw=0.8, alpha=0.85)
    for band in result.matched_wannier:
        plt.plot(result.dft.x, band, color="#D55E00", lw=0.7, ls="--", alpha=0.85)
    plt.plot([], [], color="#222222", label="DFT")
    plt.plot([], [], color="#D55E00", ls="--", label="Wannier")
    if bool_cfg(config, "fit_check", "plot_show_fermi", True):
        plt.axhline(result.fermi_energy, color="red", lw=0.9, ls=":", label="Fermi level")
    ymin = optional_float_cfg(config, "fit_check", "plot_energy_min")
    ymax = optional_float_cfg(config, "fit_check", "plot_energy_max")
    if ymin is None or ymax is None:
        values = np.concatenate([result.matched_dft.ravel(), result.matched_wannier.ravel()])
        finite = values[np.isfinite(values)]
        if finite.size:
            span = max(float(np.nanmax(finite) - np.nanmin(finite)), 1.0)
            if ymin is None:
                ymin = float(np.nanmin(finite) - 0.05 * span)
            if ymax is None:
                ymax = float(np.nanmax(finite) + 0.05 * span)
    if ymin is not None and ymax is not None:
        plt.ylim(ymin, ymax)
    plt.xlabel("Aligned k-distance")
    plt.ylabel("Energy / eV")
    plt.title(f"DFT vs Wannier band ({result.quality})")
    plt.legend()
    return save_current_figure(output_dir / "dft_vs_wannier_band", config)


def plot_fit_errors(result: FitResult, config: configparser.ConfigParser, output_dir: Path) -> List[Path]:
    plt.figure(figsize=(7.5, 4.6))
    for idx, error in enumerate(result.errors, start=1):
        plt.plot(result.dft.x, error, lw=0.7, alpha=0.8, label=f"band {idx}" if idx <= 8 else None)
    if bool_cfg(config, "fit_check", "plot_show_error_zero", True):
        plt.axhline(0.0, color="black", lw=0.9, ls=":")
    plt.xlabel("Aligned k-distance")
    plt.ylabel("E_wannier - E_dft (eV)")
    plt.title("DFT vs Wannier fitting error")
    if result.errors.shape[0] <= 8:
        plt.legend(fontsize=8)
    return save_current_figure(output_dir / "dft_vs_wannier_error", config)


def write_fit_report(result: FitResult, config: configparser.ConfigParser, output_dir: Path) -> Path:
    lines = [
        "# Fit Quality Report",
        "",
        "## 1. 读取文件",
        "",
        f"- Backend: {result.backend}",
        f"- DFT file: {result.dft.path}",
        f"- Wannier file: {result.wannier.path}",
        f"- Fermi energy: {result.fermi_energy:.10g} eV",
        f"- Fermi energy source: {result.fermi_source}",
        f"- Energy handling: {result.energy_note}",
        "",
        "## Energy reference convention",
        "",
        "All energies in this analysis are converted to the absolute, non-EF-subtracted scale.",
        "",
        "- QE DFT band energies are used directly.",
        "- QE/EPW/Wannier band energies are used directly.",
        "- VASP BAND.dat energies are generated by post-processing and have already been shifted by -EF, so EF is added back before comparison.",
        "- No global EF subtraction is applied in plotting or fitting-quality analysis.",
        "",
        "## 2. 横坐标与能量对齐",
        "",
        f"- align_x_method: {result.align_method}",
        f"- align_x_factor: {result.align_factor:.10g}",
        f"- align_y_shift: {result.y_shift:.10g} eV",
        f"- Matched DFT band indices: {', '.join(str(item) for item in result.matched_dft_indices)}",
        "",
        "## 3. 误差统计",
        "",
        f"- metric_energy_min: {result.metric_energy_min:.6g} eV",
        f"- metric_energy_max: {result.metric_energy_max:.6g} eV",
        f"- near_fermi_window: E_F +/- {result.near_fermi_window:.6g} eV",
        f"- Overall RMSE: {result.overall_rmse:.6g} eV",
        f"- Overall max error: {result.overall_max:.6g} eV",
        f"- Near-Fermi RMSE: {result.near_rmse:.6g} eV",
        f"- Near-Fermi max error: {result.near_max:.6g} eV",
        f"- DFT Fermi crossing bands: {result.dft_cross_count}",
        f"- Wannier Fermi crossing bands: {result.wannier_cross_count}",
        f"- Fermi crossing consistent: {result.crossing_consistent}",
        f"- Quality basis: {result.quality_basis}",
        "",
        "## 4. 每条关键能带误差",
        "",
        "| Band | RMSE (eV) | Max error (eV) | DFT cross E_F | Wannier cross E_F |",
        "|---:|---:|---:|---|---|",
    ]
    for row in sorted(result.per_band, key=lambda item: float(item["rmse"]), reverse=True):
        lines.append(
            f"| {row['band_index']} | {float(row['rmse']):.6g} | {float(row['max_error']):.6g} | "
            f"{row['dft_cross_fermi']} | {row['wannier_cross_fermi']} |"
        )
    lines += ["", "## 5. 最终判断", "", f"**{result.quality}**", ""]
    if result.warnings:
        lines += ["## 6. Warnings", ""]
        lines.extend(f"- {item}" for item in result.warnings)
        lines.append("")
    if result.suggestions:
        lines += ["## 7. 修改建议", ""]
        lines.extend(f"- {item}" for item in result.suggestions)
    else:
        lines.append("当前阈值下拟合质量较好。")
    path = output_dir / "fit_quality_report.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path.resolve()


def print_summary(result: FitResult, outputs: Sequence[Path]) -> None:
    print(f"读取的 DFT 文件: {result.dft.path}")
    print(f"读取的 Wannier 文件: {result.wannier.path}")
    print(f"费米能级: {result.fermi_energy:.10g} eV")
    print(f"费米能级来源: {result.fermi_source}")
    print(f"横坐标对齐因子: {result.align_factor:.10g}")
    print(f"RMSE: {result.overall_rmse:.6g} eV")
    print(f"最大误差: {result.overall_max:.6g} eV")
    print(f"费米面交叉是否一致: {result.crossing_consistent}")
    print(f"最终判断: {result.quality}")
    print("输出文件路径:")
    for path in outputs:
        print(f"  {path}")


def run_feature_three(args: argparse.Namespace) -> int:
    config = load_config()
    if not bool_cfg(config, "fit_check", "enable", True):
        print("fit_check.enable = false，已跳过拟合质量分析。")
        return 0
    output_dir = workspace.backend_output_dir(args.backend.upper())
    output_dir.mkdir(parents=True, exist_ok=True)
    try:
        result = run_fit_quality(args.backend, config, args.dft_band, args.wannier_band)
    except Exception as exc:
        print(f"错误：{exc}")
        return 1
    outputs: List[Path] = [write_fit_report(result, config, output_dir)]
    outputs.extend(plot_fit_bands(result, config, output_dir))
    outputs.extend(plot_fit_errors(result, config, output_dir))
    print_summary(result, outputs)
    return 0


def menu() -> int:
    text = """功能选择：

    1) 费米面轨道贡献
    2) 拟合能带建议
    3) 拟合质量分析
    4) 拟合扫描/拟合检查
    5) 仅初始化目录和 input.dat
    6) 退出
       请输入功能编号："""
    try:
        choice = input(text).strip()
    except EOFError:
        print("\n未收到输入，已退出。")
        return 0
    if choice == "3":
        backend = input("请输入后端 VASP/QE/EPW（默认 QE）：").strip() or "QE"
        args = argparse.Namespace(backend=backend.upper(), dft_band=None, wannier_band=None)
        return run_feature_three(args)
    if choice == "5":
        workspace.init_workspace()
        print(f"已初始化 user_work 和 input.dat: {workspace.USER_WORK}")
        return 0
    if choice == "6":
        print("已退出。")
        return 0
    print("当前阶段仅完整实现功能三；请选择 3 进行拟合质量分析。")
    return 0

