#!/usr/bin/env python3
"""Static-only AST inventory for the System A Stage A audit.

Audited Python files are parsed as text.  They are never imported or executed.
All generated output stays below the dedicated audit directory.
"""

from __future__ import annotations

import ast
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


TARGET_ROOT_REL = Path("次元の生成構造/第9論文_フェルミオンの生成構造")
ENV_REL = TARGET_ROOT_REL / "対照実験_波束収縮_実行環境_v1"
AUDIT_ROOT_REL = TARGET_ROOT_REL / "codex_systemA_scattering_audit_v1"

FILES = (
    (
        "system_a_primary",
        "inside_target",
        ENV_REL / "20260715" / "run_system_A_localization_exchange_R_sweep_preliminary_v1.py",
    ),
    (
        "system_a_source",
        "inside_target",
        ENV_REL / "20260713" / "run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py",
    ),
    (
        "system_a_instrumented",
        "inside_target",
        ENV_REL / "20260715" / "run_system_A_localization_exchange_R_sweep_instrumented_v1.py",
    ),
    (
        "system_b_local_compare",
        "inside_target",
        ENV_REL / "20260715" / "run_system_B_gray_cat_metastable_R_sweep_preliminary_v1.py",
    ),
    (
        "parity_suite_wrapper",
        "inside_target",
        ENV_REL / "parity_suite_v1" / "run_parity_suite_v1.py",
    ),
    (
        "production_dump_wrapper",
        "inside_target",
        ENV_REL / "production_dump_v1" / "run_production_dump_v1.py",
    ),
    (
        "fullkernel_wrapper",
        "inside_target",
        ENV_REL / "additional_fullkernel_wavelength_v1" / "run_fullkernel_wavelength_experiment_v1.py",
    ),
    (
        "external_system_a_copy",
        "outside_readonly",
        Path("時間軸Q軸とフェルミオンの生成構造")
        / "検証_対照実験"
        / "N3_有限位数共鳴_対照"
        / "finite_order_resonance_v1"
        / "src"
        / "run_system_A_localization_exchange_R_sweep_preliminary_v1.py",
    ),
    (
        "external_system_b_candidate",
        "outside_readonly",
        Path("時間軸Q軸とフェルミオンの生成構造")
        / "検証_対照実験"
        / "N3_有限位数共鳴_対照"
        / "finite_order_resonance_v1"
        / "src"
        / "run_minimal_system_B_gray_direct_check_v5.py",
    ),
    (
        "external_phase5_candidate",
        "outside_readonly",
        Path("時間軸Q軸とフェルミオンの生成構造")
        / "検証_対照実験"
        / "N3_有限位数共鳴_対照"
        / "finite_order_resonance_v1"
        / "src"
        / "phase5_eigenphase_resonance_v2.py",
    ),
    (
        "external_roots_candidate",
        "outside_readonly",
        Path("時間軸Q軸とフェルミオンの生成構造")
        / "検証_対照実験"
        / "N3_有限位数共鳴_対照"
        / "finite_order_resonance_v1"
        / "src"
        / "run_two_physical_roots_multiprecision_v1.py",
    ),
)


def dotted_name(node: ast.AST) -> str | None:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        prefix = dotted_name(node.value)
        return f"{prefix}.{node.attr}" if prefix else node.attr
    return None


def literal_summary(node: ast.AST) -> Any:
    try:
        value = ast.literal_eval(node)
    except (ValueError, TypeError):
        return None
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if isinstance(value, tuple) and len(value) <= 12:
        return list(value)
    return None


def function_record(node: ast.FunctionDef | ast.AsyncFunctionDef) -> dict[str, Any]:
    calls: list[dict[str, Any]] = []
    names_read: set[str] = set()
    names_written: set[str] = set()
    for child in ast.walk(node):
        if isinstance(child, ast.Call):
            name = dotted_name(child.func)
            calls.append({"name": name or "<dynamic>", "line": child.lineno})
        elif isinstance(child, ast.Name):
            if isinstance(child.ctx, ast.Load):
                names_read.add(child.id)
            elif isinstance(child.ctx, (ast.Store, ast.Del)):
                names_written.add(child.id)
    return {
        "name": node.name,
        "line_start": node.lineno,
        "line_end": node.end_lineno,
        "calls": sorted(calls, key=lambda item: (item["line"], item["name"])),
        "names_read": sorted(names_read),
        "names_written": sorted(names_written),
    }


def inspect_file(repo_root: Path, role: str, scope: str, relative_path: Path) -> dict[str, Any]:
    path = repo_root / relative_path
    text = path.read_text(encoding="utf-8")
    tree = ast.parse(text, filename=relative_path.as_posix())
    imports: list[dict[str, Any]] = []
    assignments: list[dict[str, Any]] = []
    functions: list[dict[str, Any]] = []
    classes: list[dict[str, Any]] = []
    top_level_calls: list[dict[str, Any]] = []

    for node in tree.body:
        if isinstance(node, ast.Import):
            imports.append(
                {
                    "line": node.lineno,
                    "kind": "import",
                    "names": [alias.name for alias in node.names],
                }
            )
        elif isinstance(node, ast.ImportFrom):
            imports.append(
                {
                    "line": node.lineno,
                    "kind": "from",
                    "module": node.module,
                    "names": [alias.name for alias in node.names],
                }
            )
        elif isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            value = node.value
            assignments.append(
                {
                    "line": node.lineno,
                    "targets": [dotted_name(target) or ast.dump(target) for target in targets],
                    "literal": literal_summary(value) if value is not None else None,
                    "expression": ast.unparse(value) if value is not None else None,
                }
            )
            if isinstance(value, ast.Call):
                top_level_calls.append(
                    {
                        "line": node.lineno,
                        "name": dotted_name(value.func) or "<dynamic>",
                        "expression": ast.unparse(value),
                    }
                )
        elif isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
            top_level_calls.append(
                {
                    "line": node.lineno,
                    "name": dotted_name(node.value.func) or "<dynamic>",
                    "expression": ast.unparse(node.value),
                }
            )
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            functions.append(function_record(node))
        elif isinstance(node, ast.ClassDef):
            classes.append(
                {
                    "name": node.name,
                    "line_start": node.lineno,
                    "line_end": node.end_lineno,
                    "methods": [
                        function_record(child)
                        for child in node.body
                        if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef))
                    ],
                }
            )

    return {
        "role": role,
        "scope": scope,
        "path": relative_path.as_posix(),
        "line_count": len(text.splitlines()),
        "imports": imports,
        "top_level_assignments": assignments,
        "top_level_calls": top_level_calls,
        "functions": functions,
        "classes": classes,
    }


def main() -> int:
    repo_root = Path(__file__).resolve().parents[4]
    if not (repo_root / ".git").exists():
        raise RuntimeError(f"repository root not found: {repo_root}")
    output_path = repo_root / AUDIT_ROOT_REL / "logs" / "static_code_inventory.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": "codex_systemA_static_code_inventory_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "execution_policy": "AST parse only; audited modules were not imported or executed",
        "files": [
            inspect_file(repo_root, role, scope, relative_path)
            for role, scope, relative_path in FILES
        ],
    }
    output_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
