#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Command-line entry point for Auxiliary Wannier Band Fitting."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence

from src import fit_quality
from src import utils as core
from src import workspace


def parse_entry_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Auxiliary Wannier Band Fitting")
    parser.add_argument(
        "--task",
        choices=["projection", "band", "all", "fit_quality", "fit_check", "init", "menu"],
        default="menu",
    )
    parser.add_argument("--backend", choices=["VASP", "QE", "EPW", "vasp", "qe", "epw"], default="QE")
    parser.add_argument("--dft-band", default=None, help="Feature 3: DFT band file; default is auto.")
    parser.add_argument("--wannier-band", default=None, help="Feature 3: Wannier *_band.dat file; default is auto.")
    return parser.parse_args(argv)


def choose_language() -> str:
    text = """请选择中文还是英文模式：
  1) 中文
  2) English
请输入编号："""
    try:
        choice = input(text).strip()
    except EOFError:
        return "zh"
    return "en" if choice == "2" else "zh"


def main_banner(lang: str) -> str:
    if lang == "en":
        return """================================================================================
                    Auxiliary Wannier Band Fitting
--------------------------------------------------------------------------------
Features:
1. Read DFT/Wannier band files.
2. Analyze band min/max values and Fermi crossings.
3. Suggest Wannier90 projections from PDOS/projections.
4. Suggest num_wann, num_bands, dis_win and dis_froz.
5. Check Wannier fitting quality and output good/warning/bad.
6. Optionally generate fitting scan scripts.
================================================================================
Function menu:
  1) Fermi-level orbital contribution
  2) Band fitting suggestion
  3) Fitting quality analysis
  4) Fitting scan/check
  5) Initialize user_work and input.dat only
  0) Exit
Enter function number:"""
    return """================================================================================
                    辅助 Wannier 能带拟合
             Auxiliary Wannier Band Fitting
--------------------------------------------------------------------------------
功能：
1. 读取 DFT/Wannier 能带文件。
2. 统计每条能带 min/max 和费米面交叉情况。
3. 根据 PDOS/投影给出 Wannier90 projections 建议。
4. 自动推荐 num_wann、num_bands、dis_win、dis_froz。
5. 检查 Wannier 拟合质量并输出 good/warning/bad。
6. 可选生成窗口扫描脚本。
================================================================================
功能选择：
  1) 费米面轨道贡献
  2) 拟合能带建议
  3) 拟合质量分析
  4) 拟合扫描/拟合检查
  5) 仅初始化目录和 input.dat
  0) 退出
请输入功能编号："""


def select_backend(lang: str) -> Optional[str]:
    text = (
        "请选择计算软件/数据来源：\n  1) VASP\n  2) QE\n  3) EPW\n请输入编号："
        if lang == "zh"
        else "Select code/data source:\n  1) VASP\n  2) QE\n  3) EPW\nEnter number:"
    )
    try:
        choice = input(text).strip()
    except EOFError:
        return None
    return {"1": "VASP", "2": "QE", "3": "EPW"}.get(choice)


def print_uninitialized(lang: str) -> None:
    if lang == "en":
        print(
            "Incomplete user_work workspace was not detected.\n"
            "Please choose feature 5 first: initialize user_work and input.dat.\n"
            "Then place calculation files under inputs-file/VASP, inputs-file/QE or inputs-file/EPW."
        )
        return
    print(
        "未检测到完整 user_work 工作目录。\n"
        "请先选择功能 5：仅初始化目录和 input.dat。\n"
        "初始化后，将计算文件放入 inputs-file/VASP、inputs-file/QE 或 inputs-file/EPW。"
    )


def initialize_workspace(lang: str) -> None:
    created_inputs = workspace.init_workspace()
    if lang == "en":
        print("Created/checked the following workspace:")
        print()
        print("user_work/")
        print("|-- input.dat%s" % ("" if created_inputs else "  (exists, kept)"))
        print("|-- inputs-file/")
        print("|   |-- VASP/")
        print("|   |-- QE/")
        print("|   `-- EPW/")
        print("`-- outputs-file/")
        print("    |-- VASP/")
        print("    |-- QE/")
        print("    `-- EPW/")
    else:
        print("已构建以下目录结构：")
        print()
        print("user_work/")
        print("├── input.dat%s" % ("" if created_inputs else "  (已存在，已保留)"))
        print("├── inputs-file/")
        print("│   ├── VASP/")
        print("│   ├── QE/")
        print("│   └── EPW/")
        print("└── outputs-file/")
        print("    ├── VASP/")
        print("    ├── QE/")
        print("    └── EPW/")
    print()
    if lang == "en":
        print("Place calculation files in:")
    else:
        print("请将对应计算软件的能带、费米能、投影或 PDOS 文件放入：")
    for backend in workspace.BACKENDS:
        print(f"  {workspace.backend_input_dir(backend).relative_to(workspace.USER_WORK).as_posix()}/")


def show_run_log(backend: str, lang: str) -> None:
    config = workspace.config_path()
    input_dir = workspace.backend_input_dir(backend)
    output_dir = workspace.backend_output_dir(backend)
    if lang == "en":
        print(f"Reading {config.relative_to(workspace.APP_ROOT).as_posix()} ...")
        print(f"Scanning {input_dir.relative_to(workspace.USER_WORK).as_posix()} ...")
    else:
        print(f"正在读取 {config.relative_to(workspace.APP_ROOT).as_posix()} ...")
        print(f"正在扫描 {input_dir.relative_to(workspace.USER_WORK).as_posix()} ...")
    files = sorted(path.name for path in input_dir.iterdir() if path.is_file()) if input_dir.exists() else []
    expected = {
        "VASP": ["BAND.dat", "OUTCAR", "PDOS_*.dat", "*_band.dat"],
        "QE": ["band.dat.gnu", "scf.out", "*.pdos.pdos_atm*", "*_band.dat"],
        "EPW": ["band.dat.gnu", "scf.out", "*_band.dat"],
    }[backend]
    print(("已读取：" if lang == "zh" else "Found:") + (" " + ", ".join(files) if files else " none"))
    missing = []
    for item in expected:
        if "*" in item:
            if not list(input_dir.glob(item)):
                missing.append(item)
        elif item not in files:
            missing.append(item)
    print(("缺失：" if lang == "zh" else "Missing:") + (" " + ", ".join(missing) if missing else " none"))
    output_dir.mkdir(parents=True, exist_ok=True)


def run_core_task(task: str, backend: str) -> int:
    sys.argv = [sys.argv[0], "--task", task, "--backend", backend.upper()]
    return core.main()


def run_selected(choice: str, backend: str) -> int:
    if choice == "1":
        return run_core_task("projection", backend)
    if choice == "2":
        return run_core_task("band", backend)
    if choice == "3":
        args = argparse.Namespace(backend=backend, dft_band=None, wannier_band=None)
        return fit_quality.run_feature_three(args)
    print("功能 4 暂未实现。")
    return 0


def menu() -> int:
    lang = choose_language()
    while True:
        try:
            choice = input(main_banner(lang)).strip()
        except EOFError:
            print("\n未收到输入，已退出。" if lang == "zh" else "\nNo input received. Exit.")
            return 0
        if choice == "0":
            print("已退出。" if lang == "zh" else "Exit.")
            return 0
        if choice == "5":
            initialize_workspace(lang)
            continue
        if choice not in {"1", "2", "3", "4"}:
            print("请输入有效功能编号。" if lang == "zh" else "Please enter a valid function number.")
            continue
        if not workspace.is_initialized():
            print_uninitialized(lang)
            continue
        backend = select_backend(lang)
        if backend is None:
            print("后端选择无效。" if lang == "zh" else "Invalid backend selection.")
            continue
        show_run_log(backend, lang)
        try:
            run_selected(choice, backend)
        except Exception as exc:
            print(f"错误：{exc}" if lang == "zh" else f"Error: {exc}")
        print()


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_entry_args(argv)
    args.backend = args.backend.upper()
    if args.task == "menu":
        return menu()
    if args.task == "init":
        initialize_workspace("zh")
        return 0
    if not workspace.is_initialized():
        print_uninitialized("zh")
        return 1
    try:
        if args.task in {"fit_quality", "fit_check"}:
            return fit_quality.run_feature_three(args)
        return run_core_task(args.task, args.backend)
    except Exception as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
