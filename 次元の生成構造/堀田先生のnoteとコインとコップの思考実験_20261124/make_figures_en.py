"""論文で使う 17 枚の図（figures/ の SVG）の中の日本語を英語に置き換えて、figures_en/ に書き出す。

- 図形と数値には触れない。置き換えるのは文字列だけである。
- 置き換えた後に日本語が一文字でも残っていれば、エラーで止まる。
- 使い方：python3 make_figures_en.py
"""
import pathlib
import re

HERE = pathlib.Path(__file__).parent
FIGS = [
    "fig_cup_coin_schematic_covered.svg", "fig_cup_coin_schematic_beside.svg",
    "fig_cup_coin_s1_phase90.svg", "fig_cup_coin_s1_phase60.svg",
    "cross_term/fig_cup_coin_s1_cross_phase090.svg",
    "parity/fig_cup_s1_parity_even.svg", "parity/fig_cup_s1_parity_odd.svg", "parity/fig_cup_s1_parity_both.svg",
    "open_system/fig_cup_s1_all_harmonics.svg", "open_system/fig_cup_cut_open_270.svg",
    "open_system/fig_line_fixed_ends_equal16.svg", "open_system/fig_line_fixed_ends_fejer16.svg",
    "well_cross/fig_cup_well_cross_phase090.svg", "well_cross/fig_cup_well_cross_phase075.svg",
    "well_cross/fig_cup_well_cross_phase060.svg", "well_cross/fig_cup_well_cross_summary.svg",
    "well_cross_derived_R/fig_cup_well_cross_derivedR_phase060.svg",
]

# 上から順に当てる。長い言い回しを先に、短い語を後に置く。
RULES = [
    ("R² の次元なので、R で一度割って長さにして描く:", "It has the dimension R², so it is divided by R once and drawn as a length:"),
    ("重みを掛けると波打ちは消えるが、正の側に片寄った長い裾 (黒の点線) は残る。中央での値",
     "With the weights the ripples vanish, but a long tail biased to the positive side (black dotted line) remains. Value at the middle"),
    ("山の先は、正の側に片寄ったまま波打って小さくなり、反対の端で 0 になる。最小は",
     "Beyond the peak the waveform ripples, staying biased to the positive side, decays, and becomes 0 at the opposite end. Minimum"),
    (re.compile(r"合成波形は切り口で (-[\d.]+ R) と小さいが、それは奇数と偶数が打ち消し合った結果であって、節ではない"),
     r"The sum waveform is small at the cut (\1), but that is the result of odd and even harmonics cancelling; it is not a node"),
    ("°)　最小", "°), minimum"), ("°　円の半径 = 波から計算した R", "°; radius of the circle = R computed from the waves"),
    ("切り口 (270°) では、どの倍音も 0 にならない。値は奇数倍音が −a R w、偶数倍音が +a R w で、振れが最大の位置にある",
     "At the cut (270°) no harmonic is 0. The values are −a R w for odd harmonics and +a R w for even harmonics; all are at maximum swing"),
    ("270° では、どの倍音も振れが最大 (奇数は内向き、偶数は外向き)。節になる倍音は一つもない",
     "At 270° every harmonic is at maximum swing (odd inward, even outward). No harmonic has a node there"),
    ("この波形は、偶数倍音だけの波形と奇数倍音だけの波形の和に等しい。270° では両者が打ち消し合う",
     "This waveform equals the sum of the even-only and the odd-only waveforms. At 270° the two cancel"),
    ("270° にも、90° と同じ向きで同じ高さの山が立つ", "At 270° too, a peak stands with the same direction and height as at 90°"),
    ("270° には、90° と逆向きで同じ大きさの凹みが立つ", "At 270°, a dip stands with the opposite direction and the same size as at 90°"),
    ("270° で切り開いて直線にした図", "Cut open at 270° and laid out as a line"),
    ("コインに DC 成分が無いと、コインを一周で足した値が 0 なので、覆われていないコインは C に出ない (波形には出る)。",
     "Without a constant term the coin sums to 0 over one turn, so an uncovered coin does not appear in C (it appears in the waveform)."),
    ("横はコップの中心位相。コインは 90° に固定 (赤の縦の破線)。井戸の幅", "Horizontal axis: center phase of the cup. The coin is fixed at 90° (red vertical dashed line). Well width"),
    ("描画: 合成波形は r = R + f(ψ)。倍音の点線は小さいので、半径方向に 5 倍して r = R + 5 h(ψ) で描いた",
     "Drawn as: sum waveform r = R + f(ψ). The dotted harmonics are small, so they are magnified 5 times radially: r = R + 5 h(ψ)"),
    ("描画: r = R(Δ) + 値。波形・交差項の数と縮尺は前の図と同じ。変えたのは円の半径だけ。",
     "Drawn as: r = R(Δ) + value. Numbers and scale of the waveforms and the cross term are as in the previous figure. Only the radius of the circle is changed."),
    ("描画: 中心からの距離", "Drawn as: distance from the center"),
    ("山は直線の端に寄る: 最大", "The peak sits near the end of the line: maximum"),
    (re.compile(r" が x/R = ([\d.]+) \(端から全長の ([\d.]+%)\)"), r" at x/R = \1 (\2 of the full length from the end)"),
    ("両端を固定した長さ L = 2πR の直線に、定常波 1〜16 番を足した波形", "Sum of standing waves 1–16 on a line of length L = 2πR with both ends fixed"),
    ("両端では全ての定常波が 0", "All standing waves are 0 at both ends"),
    ("参考: コインに m = 0 の項 (a R / 2) を足した場合の C", "Reference: C when the m = 0 term (a R / 2) is added to the coin"),
    ("(コインはきのうのまま。DC 成分なし)", "(coin as on the previous day; no constant term)"),
    ("半径 R(Δ) の円。", "Circle of radius R(Δ). "),
    ("(倍音のノルム。直流は含めない)", "(norm of the harmonics; the constant term is not included)"),
    ("和の R² =", "R² of the sum ="),
    (re.compile(r"。Δ = (-?\d+°) で "), r". At Δ = \1, "),
    ("コインの R₁", "R₁ of the coin"), ("コップの R₂", "R₂ of the cup"), ("比 R₁/R₂", "ratio R₁/R₂"),
    ("2. 交差項を円周全体で足した値 C (単位 R²)", "2. Cross term summed over the whole circle, C (in units of R²)"),
    ("円周全体で足した値", "Value summed over the whole circle"),
    ("1. 交差項の波形 (その場所での値。単位 R²)", "1. Cross-term waveform (local value; in units of R²)"),
    ("円周の上の位相 ψ (切り口から一周)", "Phase ψ on the circle (one turn from the cut)"),
    ("直線の上の位置", "Position on the line"),
    ("コップをコインの隣に置いた状態", "The cup placed beside the coin"),
    ("コップをコインの上に被せた状態", "The cup placed over the coin"),
    ("コップの倍音を一本ずつ描いた図", "The harmonics of the cup drawn one by one"),
    ("コップ (井戸型) の位置と交差項の関係", "Position of the cup (well-shaped) and the cross term"),
    ("コップ (井戸型) の中心位相", "Center phase of the cup (well-shaped)"),
    ("円の半径 = 波から計算した R", "radius of the circle = R computed from the waves"),
    ("コップの中心位相", "Center phase of the cup"),
    ("コップの合成波形", "Sum waveform of the cup"),
    ("コインの中心でのコップの値", "Value of the cup at the coin center"),
    ("コインの中心 (90°) での交差項", "Cross term at the coin center (90°)"),
    ("コインの中心 (90°) での値", "Value at the coin center (90°)"),
    ("そこでのコップの値", "value of the cup there"),
    ("コインの山", "coin peak"),
    ("交差項の波形の最大", "Maximum of the cross-term waveform"),
    ("波形の最大", "waveform maximum"), ("最小", "minimum"),
    ("合成波形 f(ψ) は m = 1…23 の和", "The sum waveform f(ψ) is the sum over m = 1…23"),
    ("合成波形 (R 単位)", "Sum waveform (in units of R)"),
    ("合成波形", "Sum waveform"),
    ("定常波 一本ずつ (R 単位)", "Each standing wave (in units of R)"),
    ("倍音 一本ずつ (R 単位)", "Each harmonic (in units of R)"),
    ("奇数番の定常波 (一本ずつ)", "Odd-numbered standing waves (each)"),
    ("偶数番の定常波 (一本ずつ)", "Even-numbered standing waves (each)"),
    ("偶数倍音と奇数倍音の両方", "Both even and odd harmonics"),
    ("奇数倍音だけ", "Odd harmonics only"), ("偶数倍音だけ", "Even harmonics only"),
    ("奇数倍音", "Odd harmonics"), ("偶数倍音", "Even harmonics"),
    ("(一本ずつ)", "(each)"),
    (re.compile(r"、Σ は (m = [^。]+?) について。"), r", Σ over \1. "),
    ("等振幅で足した波形 (比較)", "Sum with equal amplitudes (for comparison)"),
    ("裾の曲線", "Tail curve"),
    ("等振幅", "equal amplitudes"),
    ("テーブル (円周、半径 R)", "Table (circle, radius R)"),
    ("テーブル", "Table"),
    ("交差項 2 × コイン × コップ", "Cross term 2 × coin × cup"),
    ("f(コイン)", "f(coin)"), ("g(コップ)", "g(cup)"), ("f(コップ)", "f(cup)"),
    ("交差項:", "Cross term:"),
    ("固定端", "fixed end"), ("切り口", "cut"),
    ("井戸の幅 W、深さ h", "well width W, depth h"),
    ("井戸型。幅", "well-shaped; width"), ("井戸型", "well-shaped"),
    ("深さ", "depth"),
    ("フェイエール重み", "Fejér weights"),
    ("一本あたりの振幅", "amplitude per harmonic"),
    ("山の高さ", "peak height"),
    ("中心", "center"), ("凸。", "convex; "), ("凸", "convex"),
    (re.compile(r"倍音 (\d+)〜(\d+)"), r"harmonics \1–\2"),
    ("倍音:", "Harmonics:"),
    ("波形:", "Waveform:"),
    ("コップ:", "Cup:"), ("コイン:", "Coin:"),
    ("コップ", "Cup"), ("コイン", "Coin"),
    ("〜", "–"), ("、", ", "), ("。", ". "), ("　　", ";  "), ("　", "  "),
]

JP = re.compile(r"[　-ヿ一-鿿＀-￯]")


def translate(s):
    for a, b in RULES:
        s = a.sub(b, s) if hasattr(a, "sub") else s.replace(a, b)
    return s


TEXT = re.compile(r"<text ([^>]*)>(.*?)</text>", re.S)


def fit(svg):
    """英語にして長くなった行が右端からはみ出すときは、その行の文字だけを小さくする（図形と数値は変えない）。"""
    m = re.search(r'viewBox="[\d.\-]+ [\d.\-]+ ([\d.]+) [\d.]+"', svg) or re.search(r'<svg[^>]*\swidth="([\d.]+)', svg)
    width = float(m.group(1))

    def one(t):
        attrs, inner = t.group(1), t.group(2)
        if 'text-anchor="start"' not in attrs or "transform" in attrs:
            return t.group(0)
        x = float(re.search(r'\bx="([\d.\-]+)"', attrs).group(1))
        size = float(re.search(r'font-size="([\d.]+)"', attrs).group(1))
        plain = re.sub(r"<[^>]+>", "", inner)
        est = len(plain) * size * 0.54
        avail = width - x - 8
        if est <= avail:
            return t.group(0)
        new = max(8.0, round(size * avail / est, 1))
        return "<text %s>%s</text>" % (attrs.replace('font-size="%s"' % re.search(r'font-size="([\d.]+)"', attrs).group(1), 'font-size="%s"' % new), inner)

    return TEXT.sub(one, svg)


def main():
    for rel in FIGS:
        src = HERE / "figures" / rel
        text = fit(translate(src.read_text(encoding="utf-8")))
        left = sorted(set(m.group(0) for m in JP.finditer(text)))
        if left:
            lines = [l for l in text.splitlines() if JP.search(l)]
            raise SystemExit("日本語が残っている: %s %s\n%s" % (rel, "".join(left), "\n".join(lines[:5])))
        dst = HERE / "figures_en" / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(text, encoding="utf-8")
        print("ok", rel)


if __name__ == "__main__":
    main()
