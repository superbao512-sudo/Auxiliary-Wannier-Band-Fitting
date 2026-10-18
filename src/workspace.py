#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Workspace paths and initialization helpers."""

from __future__ import annotations

from pathlib import Path
from typing import List

from .default_inputs import DEFAULT_INPUTS_DAT

APP_ROOT = Path(__file__).resolve().parents[1]
USER_WORK = APP_ROOT / "user_work"
INPUT_DAT = USER_WORK / "input.dat"
REFERENCE_INPUT_DAT = APP_ROOT / "input.dat"
INPUTS_FILE = USER_WORK / "inputs-file"
OUTPUTS_FILE = USER_WORK / "outputs-file"
BACKENDS = ("VASP", "QE", "EPW")


def backend_input_dir(backend: str) -> Path:
    return INPUTS_FILE / backend.upper()


def backend_output_dir(backend: str) -> Path:
    return OUTPUTS_FILE / backend.upper()


def expected_dirs() -> List[Path]:
    dirs = [USER_WORK, INPUTS_FILE, OUTPUTS_FILE]
    for backend in BACKENDS:
        dirs.append(backend_input_dir(backend))
        dirs.append(backend_output_dir(backend))
    return dirs


def is_initialized() -> bool:
    return INPUT_DAT.exists() and all(path.exists() for path in expected_dirs())


def init_workspace(overwrite_inputs: bool = False) -> bool:
    for path in expected_dirs():
        path.mkdir(parents=True, exist_ok=True)
    created = False
    if overwrite_inputs or not INPUT_DAT.exists():
        INPUT_DAT.write_text(DEFAULT_INPUTS_DAT, encoding="utf-8", newline="\n")
        created = True
    return created


def config_path() -> Path:
    if INPUT_DAT.exists():
        return INPUT_DAT
    if REFERENCE_INPUT_DAT.exists():
        return REFERENCE_INPUT_DAT
    return INPUT_DAT


def require_workspace() -> None:
    if not is_initialized():
        raise FileNotFoundError(
            "未检测到完整 user_work 工作目录。请先选择功能 5：仅初始化目录和 input.dat。"
        )
