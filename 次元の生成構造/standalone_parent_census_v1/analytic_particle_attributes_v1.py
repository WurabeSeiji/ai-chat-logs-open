#!/usr/bin/env python3
"""分解能 N だけから、定常倍音族の粒子属性候補を解析的に数える。

物理入力は分解能 N だけである。H、M、seed、DFT、数値波形は使わない。

規則
----
* lambda0 = 2*pi/N。
* 基底波長は q*lambda0 (q | N)。
* 基底 q*lambda0 に載る倍音波長は d*lambda0 (d | q, d < q)。
  倍音倍率は r=q/d。
* 各選択成分の位相は 0 度または 180 度。
* 現行の万能相互作用演算器と同じく A（基底）と B（倍音束）を個別に
  単位規格化する。等振幅 B 束では R=P_f/(P_f+P_b)=O/[2(O+E)]。
* theta=asin(sqrt(R)) とし、有限位数になる場合だけ閉鎖住所 m/n を読む。
  正準読出しで可能なのは Niven 点 R=0,1/4,1/2（R<=1/2）である。
  R=0 は不変な自由基底で住所なし、R=1/4 は (m,n)=(1,3)、
  R=1/2 は (m,n)=(1,4)。モデル電荷量は sin^2(pi*m/n)。
* 粒子・反粒子は有限住所 +m/n と -m/n の共役方向対として数える。
* 内部倍音が純奇数なら 2:1 被覆（半整数スピン型）、純偶数なら
  1:1 被覆（整数スピン型）、混合なら単一被覆型に分類しない。

反射率は等振幅での解析値である。振幅が与えられた場合は、個数比ではなく
奇数・偶数成分のパワー比から再計算する。
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class HarmonicMode:
    wavelength_multiple: int
    harmonic_order: int
    parity: str


@dataclass(frozen=True)
class AttributeFamily:
    family_id: str
    resolution: int
    base_wavelength_multiple: int
    base_wavelength: str
    available_odd_harmonic_orders: tuple[int, ...]
    available_even_harmonic_orders: tuple[int, ...]
    selected_odd_count: int
    selected_even_count: int
    harmonic_support_count: int
    phase_count_per_support: int
    waveform_count: int
    reflection_fraction: str
    reflection_value: float
    theta_over_pi: str
    finite_order_lock: bool
    real_rotation_period: int | None
    closure_address: str
    conjugate_address: str
    charge_magnitude: float | None
    particle_antiparticle: str
    conjugate_direction_count: int
    winding_cover_spin: str
    full_attribute_state_count: int


def divisors(n: int) -> list[int]:
    if n < 1:
        raise ValueError("N must be a positive integer")
    low: list[int] = []
    high: list[int] = []
    root = math.isqrt(n)
    for d in range(1, root + 1):
        if n % d == 0:
            low.append(d)
            if d * d != n:
                high.append(n // d)
    return low + high[::-1]


def harmonic_modes(base_q: int) -> list[HarmonicMode]:
    modes: list[HarmonicMode] = []
    for d in divisors(base_q):
        if d == base_q:
            continue
        order = base_q // d
        modes.append(
            HarmonicMode(
                wavelength_multiple=d,
                harmonic_order=order,
                parity="奇数" if order % 2 else "偶数",
            )
        )
    return modes


def reflection(odd_count: int, even_count: int) -> tuple[str, float]:
    # A（基底）と B（倍音束）は個別に単位規格化される。B 内の F 側比率
    # O/(O+E) は AB 全体では半分になる。
    total_harmonics = odd_count + even_count
    if total_harmonics == 0:
        return "0", 0.0
    value = Fraction(odd_count, 2 * total_harmonics)
    return str(value), float(value)


def address_and_charge(
    reflection_fraction: str,
) -> tuple[str, bool, int | None, str, str, float | None, int, str]:
    """正準読出しの R から有限位数住所を解析する。"""

    r = Fraction(reflection_fraction)
    if r == 0:
        return "0", False, None, "なし", "なし", 0.0, 1, "閉鎖共役対なし（自由基底）"
    niven_theta: dict[Fraction, Fraction] = {
        Fraction(1, 4): Fraction(1, 6),
        Fraction(1, 2): Fraction(1, 4),
        Fraction(3, 4): Fraction(1, 3),
        Fraction(1, 1): Fraction(1, 2),
    }
    if r not in niven_theta:
        theta = math.asin(math.sqrt(float(r))) / math.pi
        return (
            f"{theta:.12g}（無理数）",
            False,
            None,
            "なし（有限位数非成立）",
            "なし",
            None,
            1,
            "閉鎖共役対なし",
        )
    theta_fraction = niven_theta[r]
    address_fraction = Fraction(1, 2) - theta_fraction
    m, order = address_fraction.numerator, address_fraction.denominator
    period = (2 * theta_fraction.denominator) // math.gcd(
        theta_fraction.numerator, 2 * theta_fraction.denominator
    )
    theta_text = f"{theta_fraction.numerator}/{theta_fraction.denominator}"
    if m == 0:
        return theta_text, True, period, "0/1", "0/1", 0.0, 1, "自己共役"
    charge = math.sin(math.pi * m / order) ** 2
    positive = f"+{m}/{order}"
    negative_numerator = order - m
    negative = f"-{m}/{order} (={negative_numerator}/{order})"
    return (
        theta_text,
        True,
        period,
        positive,
        negative,
        charge,
        2,
        "粒子・反粒子の共役方向対（+住所 / -住所）",
    )


def spin_readout(odd_count: int, even_count: int) -> str:
    if odd_count == 0 and even_count == 0:
        return "内部倍音なし（スピン型未分類）"
    if odd_count > 0 and even_count == 0:
        return "奇数巻数のみ: 2:1被覆（半整数スピン型）"
    if even_count > 0 and odd_count == 0:
        return "偶数巻数のみ: 1:1被覆（整数スピン型）"
    return "奇偶巻数混合: 混合被覆（単一スピン型に未分類）"


def build_families(n: int) -> tuple[list[AttributeFamily], list[dict[str, object]]]:
    families: list[AttributeFamily] = []
    mode_rows: list[dict[str, object]] = []
    serial = 0
    for base_q in divisors(n):
        modes = harmonic_modes(base_q)
        odd_orders = tuple(m.harmonic_order for m in modes if m.parity == "奇数")
        even_orders = tuple(m.harmonic_order for m in modes if m.parity == "偶数")
        for mode in modes:
            mode_rows.append(
                {
                    "resolution": n,
                    "base_wavelength_multiple": base_q,
                    "base_wavelength": f"{base_q}lambda0",
                    "harmonic_wavelength_multiple": mode.wavelength_multiple,
                    "harmonic_wavelength": f"{mode.wavelength_multiple}lambda0",
                    "harmonic_order": mode.harmonic_order,
                    "parity": mode.parity,
                }
            )

        for odd_count in range(len(odd_orders) + 1):
            for even_count in range(len(even_orders) + 1):
                serial += 1
                selected_count = odd_count + even_count
                support_count = math.comb(len(odd_orders), odd_count) * math.comb(
                    len(even_orders), even_count
                )
                phase_count = 2 ** (selected_count + 1)  # 基底も 0/180 の二択
                waveform_count = support_count * phase_count
                r_fraction, r_value = reflection(odd_count, even_count)
                (
                    theta_over_pi,
                    finite_order,
                    rotation_period,
                    address,
                    conjugate,
                    charge,
                    direction_count,
                    pair_label,
                ) = address_and_charge(r_fraction)
                families.append(
                    AttributeFamily(
                        family_id=f"N{n}-F{serial:03d}",
                        resolution=n,
                        base_wavelength_multiple=base_q,
                        base_wavelength=f"{base_q}lambda0",
                        available_odd_harmonic_orders=odd_orders,
                        available_even_harmonic_orders=even_orders,
                        selected_odd_count=odd_count,
                        selected_even_count=even_count,
                        harmonic_support_count=support_count,
                        phase_count_per_support=phase_count,
                        waveform_count=waveform_count,
                        reflection_fraction=r_fraction,
                        reflection_value=r_value,
                        theta_over_pi=theta_over_pi,
                        finite_order_lock=finite_order,
                        real_rotation_period=rotation_period,
                        closure_address=address,
                        conjugate_address=conjugate,
                        charge_magnitude=charge,
                        particle_antiparticle=pair_label,
                        conjugate_direction_count=direction_count,
                        winding_cover_spin=spin_readout(odd_count, even_count),
                        full_attribute_state_count=waveform_count * direction_count,
                    )
                )
    return families, mode_rows


def expected_waveform_count(n: int) -> int:
    return sum(2 * (3 ** (len(divisors(q)) - 1)) for q in divisors(n))


def exact_waveforms(n: int) -> list[dict[str, object]]:
    """小さい N の検算用。位相まで展開した個別波形を返す。"""

    from itertools import product

    rows: list[dict[str, object]] = []
    serial = 0
    for base_q in divisors(n):
        modes = harmonic_modes(base_q)
        # 各倍音は 0=不在、1=0度、2=180度。基底は必ず存在する。
        for base_phase in (0, 180):
            for choices in product((0, 1, 2), repeat=len(modes)):
                selected: list[tuple[HarmonicMode, int]] = []
                for mode, choice in zip(modes, choices):
                    if choice:
                        selected.append((mode, 0 if choice == 1 else 180))
                odd_count = sum(mode.harmonic_order % 2 for mode, _ in selected)
                even_count = len(selected) - odd_count
                r_fraction, r_value = reflection(odd_count, even_count)
                (
                    theta_over_pi,
                    finite_order,
                    rotation_period,
                    address,
                    conjugate,
                    charge,
                    direction_count,
                    pair_label,
                ) = address_and_charge(r_fraction)
                serial += 1
                rows.append(
                    {
                        "wave_id": f"N{n}-W{serial:04d}",
                        "resolution": n,
                        "base_wave": f"{base_q}lambda0({base_phase}deg)",
                        "harmonics": "+".join(
                            f"{mode.wavelength_multiple}lambda0({phase}deg);order={mode.harmonic_order}"
                            for mode, phase in selected
                        )
                        or "なし",
                        "phase_composition": ", ".join(
                            [f"base:{base_phase}deg"]
                            + [
                                f"order{mode.harmonic_order}:{phase}deg"
                                for mode, phase in selected
                            ]
                        ),
                        "reflection_fraction": r_fraction,
                        "reflection_value": r_value,
                        "theta_over_pi": theta_over_pi,
                        "finite_order_lock": finite_order,
                        "real_rotation_period": rotation_period,
                        "closure_address": address,
                        "conjugate_address": conjugate,
                        "charge_magnitude": charge,
                        "particle_antiparticle": pair_label,
                        "conjugate_direction_count": direction_count,
                        "winding_cover_spin": spin_readout(odd_count, even_count),
                    }
                )
    return rows


def write_csv(path: Path, rows: Iterable[dict[str, object]]) -> None:
    materialized = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not materialized:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(materialized[0]))
        writer.writeheader()
        writer.writerows(materialized)


def family_csv_row(row: AttributeFamily) -> dict[str, object]:
    return {
        "族ID": row.family_id,
        "分解能N": row.resolution,
        "基底波長": row.base_wavelength,
        "利用可能な奇数倍音次数": ",".join(map(str, row.available_odd_harmonic_orders)) or "なし",
        "利用可能な偶数倍音次数": ",".join(map(str, row.available_even_harmonic_orders)) or "なし",
        "選択する奇数倍音数": row.selected_odd_count,
        "選択する偶数倍音数": row.selected_even_count,
        "波長・倍音構成数": row.harmonic_support_count,
        "各構成の位相数": row.phase_count_per_support,
        "位相込み波形数": row.waveform_count,
        "反射率R（等振幅）": row.reflection_fraction,
        "theta/pi": row.theta_over_pi,
        "有限位数閉鎖": "成立" if row.finite_order_lock else "不成立",
        "実回転周期": row.real_rotation_period if row.real_rotation_period else "—",
        "閉鎖住所": row.closure_address,
        "共役住所": row.conjugate_address,
        "住所電荷量sin2(pi*m/n)": (
            f"{row.charge_magnitude:.12g}" if row.charge_magnitude is not None else "—"
        ),
        "粒子・反粒子": row.particle_antiparticle,
        "巻数・被覆から読むスピン": row.winding_cover_spin,
        "共役方向を区別した候補状態数": row.full_attribute_state_count,
    }


def markdown_table_a(rows: list[AttributeFamily]) -> list[str]:
    out = [
        "| 族ID | 基底波長 | 選択する奇数倍音数 | 選択する偶数倍音数 | 反射率 R | 波長構成数 | 各構成の位相数 | 位相込み波形数 |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        out.append(
            f"| {r.family_id} | {r.base_wavelength} | {r.selected_odd_count} | "
            f"{r.selected_even_count} | {r.reflection_fraction} | "
            f"{r.harmonic_support_count} | {r.phase_count_per_support} | {r.waveform_count} |"
        )
    return out


def markdown_mode_catalog(n: int) -> list[str]:
    out = [
        "| 基底波長 | 許容する倍音（波長; 倍率; 偶奇） |",
        "|---:|---|",
    ]
    for base_q in divisors(n):
        modes = harmonic_modes(base_q)
        description = ", ".join(
            f"{m.wavelength_multiple}lambda0; r={m.harmonic_order}; {m.parity}"
            for m in modes
        ) or "なし"
        out.append(f"| {base_q}lambda0 | {description} |")
    return out


def markdown_exact_waveforms(rows: list[dict[str, object]]) -> list[str]:
    out = [
        "| 波形ID | 基底波 | 倍音 | 反射率 R | 閉鎖住所と電荷量 | 共役方向 | スピン型 |",
        "|---|---|---|---:|---|---|---|",
    ]
    for row in rows:
        charge = (
            f"{float(row['charge_magnitude']):.9f}"
            if row["charge_magnitude"] is not None
            else "—"
        )
        out.append(
            f"| {row['wave_id']} | {row['base_wave']} | {row['harmonics']} | "
            f"{row['reflection_fraction']} | {row['closure_address']}; {charge} | "
            f"{row['particle_antiparticle']} | "
            f"{row['winding_cover_spin']} |"
        )
    return out


def markdown_table_b(rows: list[AttributeFamily]) -> list[str]:
    out = [
        "| 族ID | theta/pi | 閉鎖住所 | 住所から読む電荷量 | 粒子・反粒子 | 巻数・被覆から読むスピン | 方向込み候補数 |",
        "|---|---:|---|---:|---|---|---:|",
    ]
    for r in rows:
        charge = f"{r.charge_magnitude:.9f}" if r.charge_magnitude is not None else "—"
        out.append(
            f"| {r.family_id} | {r.theta_over_pi} | "
            f"{r.closure_address} ↔ {r.conjugate_address} | {charge} | "
            f"{r.particle_antiparticle} | "
            f"{r.winding_cover_spin} | {r.full_attribute_state_count} |"
        )
    return out


def write_resolution_report(output_dir: Path, n: int) -> dict[str, object]:
    families, mode_rows = build_families(n)
    wave_count = sum(row.waveform_count for row in families)
    formula_count = expected_waveform_count(n)
    if wave_count != formula_count:
        raise AssertionError((n, wave_count, formula_count))
    attribute_count = sum(row.full_attribute_state_count for row in families)
    locked_waveform_count = sum(
        row.waveform_count for row in families if row.finite_order_lock and row.reflection_value > 0
    )
    locked_particle_state_count = sum(
        row.full_attribute_state_count
        for row in families
        if row.finite_order_lock and row.reflection_value > 0
    )
    n_dir = output_dir / f"N{n}"
    n_dir.mkdir(parents=True, exist_ok=True)
    write_csv(n_dir / "粒子属性族一覧.csv", (family_csv_row(row) for row in families))
    write_csv(n_dir / "許容倍音一覧.csv", mode_rows)

    # N=5 は8波形なので、位相を含む全波形も検算表として出す。
    exact_rows: list[dict[str, object]] = []
    if formula_count <= 100:
        exact_rows = exact_waveforms(n)
        if len(exact_rows) != formula_count:
            raise AssertionError((n, len(exact_rows), formula_count))
        write_csv(n_dir / "全波形一覧.csv", exact_rows)

    lines = [
        f"# N={n} 解析的粒子属性表",
        "",
        f"- $\\lambda_0=2\\pi/{n}$",
        f"- 許容基底波長 $q\\lambda_0$: {', '.join(str(q) for q in divisors(n))} ($q\\mid N$)",
        f"- 解析族数: **{len(families)}**",
        f"- 位相を区別した波形数: **{wave_count}**",
        f"- 共役方向まで区別した候補状態数: **{attribute_count}**",
        f"- 非自明な有限位数閉鎖を持つ波形数: **{locked_waveform_count}**",
        f"- その粒子・反粒子方向まで区別した状態数: **{locked_particle_state_count}**",
        "",
        "反射率 $R$ は、現行の万能相互作用演算器と同じく A（基底）と B（倍音束）を"
        "個別に単位規格化した等振幅値 $R=O/[2(O+E)]$ である。"
        "ここで $O$ は奇数倍音数、$E$ は偶数倍音数である。"
        "振幅を導入した場合は個数比ではなくパワー比で再計算する。",
        "閉鎖住所は静的な波長から仮定せず、$R=\\sin^2\\theta$ が厳密有限位数を作る"
        "場合だけ $m/n=1/2-\\theta/\\pi$ から表示する。正準読出しの非自明な"
        " Niven 点は $R=1/4,1/2$ である。",
        "",
        "## 許容する波長と倍音",
        "",
        *markdown_mode_catalog(n),
        "",
        "## A. 波長・倍音・位相・反射率",
        "",
        *markdown_table_a(families),
        "",
        "## B. 閉鎖住所・電荷・共役方向・スピン",
        "",
        *markdown_table_b(families),
        "",
        *(
            ["## C. 位相まで展開した全波形", "", *markdown_exact_waveforms(exact_rows), ""]
            if exact_rows
            else []
        ),
        "## 数え上げ式",
        "",
        "基底 $q\\lambda_0$ ごとに、基底位相が2通り、各内部倍音が"
        "「不在・0度・180度」の3通りなので、",
        "",
        "$$",
        f"W({n})=\\sum_{{q\\mid {n}}}2\\,3^{{\\tau(q)-1}}={wave_count}.",
        "$$",
        "",
    ]
    (n_dir / "粒子属性表.md").write_text("\n".join(lines), encoding="utf-8")
    return {
        "resolution": n,
        "lambda0": f"2*pi/{n}",
        "allowed_base_wavelength_multiples": divisors(n),
        "analytic_family_count": len(families),
        "phase_distinguished_waveform_count": wave_count,
        "full_attribute_state_count": attribute_count,
        "nontrivial_finite_order_waveform_count": locked_waveform_count,
        "nontrivial_finite_order_particle_antiparticle_state_count": locked_particle_state_count,
        "families": [asdict(row) for row in families],
        "exact_waveforms_written": bool(exact_rows),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("N", nargs="+", type=int, help="分解能 N（複数指定可）")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("analytic_particle_tables_v1"),
        help="出力先（物理パラメータではない）",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if any(n < 1 for n in args.N):
        raise SystemExit("N must be a positive integer")
    args.output.mkdir(parents=True, exist_ok=True)
    reports = [write_resolution_report(args.output, n) for n in args.N]
    (args.output / "summary.json").write_text(
        json.dumps(reports, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    summary_lines = [
        "# 解析的粒子属性表 v1",
        "",
        "| N | 約数（許容基底波長数） | 解析族数 | 位相込み波形数 | 非自明有限閉鎖波形数 | 粒子・反粒子方向込み閉鎖状態数 |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for report in reports:
        summary_lines.append(
            f"| {report['resolution']} | {len(report['allowed_base_wavelength_multiples'])} | "
            f"{report['analytic_family_count']} | {report['phase_distinguished_waveform_count']} | "
            f"{report['nontrivial_finite_order_waveform_count']} | "
            f"{report['nontrivial_finite_order_particle_antiparticle_state_count']} |"
        )
    summary_lines.extend(
        [
            "",
            "物理入力は分解能 $N$ だけである。$H$、$M=N(N-1)/2$、seed、DFTは使わない。",
            "",
            "各 N の詳細は `N*/粒子属性表.md`、表計算用データは `N*/粒子属性族一覧.csv` にある。",
            "",
        ]
    )
    (args.output / "summary.md").write_text("\n".join(summary_lines), encoding="utf-8")
    for report in reports:
        print(
            f"N={report['resolution']}: families={report['analytic_family_count']}, "
            f"waveforms={report['phase_distinguished_waveform_count']}, "
            f"full_states={report['full_attribute_state_count']}"
        )


if __name__ == "__main__":
    main()
