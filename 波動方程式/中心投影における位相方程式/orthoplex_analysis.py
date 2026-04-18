"""
5次元正軸体（cross-polytope）の数値解析プログラム
===================================================
W6論文「5次元正軸体の配向構造から導かれるスピンの幾何学的分類」の
数値的検証と、折り紙帰着に向けた探索を行う。
"""

import sys
import numpy as np
from itertools import combinations
from collections import Counter, defaultdict
import time


def log(msg=""):
    print(msg, flush=True)


# ============================================================
# 1. 基本構造
# ============================================================

def build_orthoplex_5d():
    vertices = []
    for i in range(5):
        v_pos = np.zeros(5); v_pos[i] = 1.0
        v_neg = np.zeros(5); v_neg[i] = -1.0
        vertices.append(v_pos)
        vertices.append(v_neg)
    vertices = np.array(vertices)

    edge_set = set()
    edges = []
    for i in range(10):
        for j in range(i+1, 10):
            if i // 2 != j // 2:
                edges.append((i, j))
                edge_set.add((i, j))

    faces = []
    for i in range(10):
        for j in range(i+1, 10):
            for k in range(j+1, 10):
                if len({i//2, j//2, k//2}) == 3:
                    if (i,j) in edge_set and (i,k) in edge_set and (j,k) in edge_set:
                        faces.append((i, j, k))

    cells = []
    for combo in combinations(range(10), 4):
        if len({v//2 for v in combo}) == 4:
            if all((a,b) in edge_set or (b,a) in edge_set for a,b in combinations(combo, 2)):
                cells.append(combo)

    facets_4d = []
    for combo in combinations(range(10), 5):
        if len({v//2 for v in combo}) == 5:
            if all((a,b) in edge_set or (b,a) in edge_set for a,b in combinations(combo, 2)):
                facets_4d.append(combo)

    # 辺→インデックスの逆引き
    edge_index = {}
    for idx, (i, j) in enumerate(edges):
        edge_index[(i,j)] = idx
        edge_index[(j,i)] = idx

    return vertices, edges, faces, cells, facets_4d, edge_index


def print_basic_structure(vertices, edges, faces, cells, facets_4d):
    from math import comb
    log("=" * 60)
    log("1. 5次元正軸体の基本構造")
    log("=" * 60)

    theoretical = [
        ('頂点 (0-胞)', 2**1 * comb(5,1), len(vertices)),
        ('辺 (1-胞)',   2**2 * comb(5,2), len(edges)),
        ('面 (2-胞)',   2**3 * comb(5,3), len(faces)),
        ('胞 (3-胞)',   2**4 * comb(5,4), len(cells)),
        ('超胞 (4-胞)', 2**5 * comb(5,5), len(facets_4d)),
    ]
    for name, t, a in theoretical:
        log(f"  {name}: {a}  (理論値: {t})  {'✓' if t == a else '✗'}")

    degree = Counter()
    for i, j in edges:
        degree[i] += 1; degree[j] += 1
    log(f"  各頂点の次数: {set(degree.values()).pop()}  (理論値: 8)")
    log()


# ============================================================
# 2. 辺の型分類（S₃対称性）
# ============================================================

def classify_edges(vertices, edges):
    log("=" * 60)
    log("2. 辺の型分類（S₃対称性）")
    log("=" * 60)

    edge_types = defaultdict(list)
    for idx, (i, j) in enumerate(edges):
        ax_i, ax_j = i // 2, j // 2
        if ax_i < 3 and ax_j < 3:
            etype = "空間-空間"
        elif (ax_i < 3 and ax_j == 3) or (ax_i == 3 and ax_j < 3):
            etype = "空間-第4軸"
        elif (ax_i < 3 and ax_j == 4) or (ax_i == 4 and ax_j < 3):
            etype = "空間-第5軸"
        elif {ax_i, ax_j} == {3, 4}:
            etype = "第4軸-第5軸"
        else:
            etype = "other"
        edge_types[etype].append(idx)

    for etype in ["空間-空間", "空間-第4軸", "空間-第5軸", "第4軸-第5軸"]:
        log(f"  {etype}: {len(edge_types[etype])}本")

    counts = [len(edge_types[t]) for t in ["空間-空間", "空間-第4軸", "空間-第5軸", "第4軸-第5軸"]]
    log(f"  合計: {sum(counts)} = {' + '.join(map(str, counts))}")
    log()
    return edge_types


# ============================================================
# 3. 中心投影
# ============================================================

def analyze_projection(vertices, edges, edge_types):
    log("=" * 60)
    log("3. 中心投影の効果")
    log("=" * 60)

    proj = vertices[:, :4]

    log(f"  +x5 {vertices[8]} → {proj[8]}")
    log(f"  -x5 {vertices[9]} → {proj[9]}")
    log(f"  *** +x5 と -x5 は射影後に原点に潰れる ***")
    log()

    log("  射影後の辺の長さ:")
    for etype in ["空間-空間", "空間-第4軸", "空間-第5軸", "第4軸-第5軸"]:
        lengths = set()
        for idx in edge_types[etype]:
            i, j = edges[idx]
            d = np.linalg.norm(proj[i] - proj[j])
            lengths.add(round(d, 6))
        log(f"    {etype}: {lengths}")
    log()

    parallel_groups = defaultdict(list)
    for idx, (i, j) in enumerate(edges):
        direction = proj[j] - proj[i]
        norm = np.linalg.norm(direction)
        if norm < 1e-10:
            parallel_groups["退化"].append(idx)
            continue
        direction = direction / norm
        for d in range(4):
            if abs(direction[d]) > 1e-10:
                if direction[d] < 0:
                    direction = -direction
                break
        key = tuple(np.round(direction, 6))
        parallel_groups[key].append(idx)

    n_deg = len(parallel_groups.get("退化", []))
    dirs = {k: v for k, v in parallel_groups.items() if k != "退化"}
    log(f"  射影後の方向グループ数: {len(dirs)}")
    log(f"  退化辺（長さ0）: {n_deg}本")
    for key, group in sorted(dirs.items(), key=lambda x: -len(x[1])):
        log(f"    方向 {key}: {len(group)}本")
    log()
    return proj


# ============================================================
# 4. 面の配向と射影後の位置
# ============================================================

def analyze_face_orientations(vertices, faces, proj):
    log("=" * 60)
    log("4. 面の配向と射影後の位置")
    log("=" * 60)

    inside = outside = degenerate = 0
    for face in faces:
        i, j, k = face
        v1 = proj[j] - proj[i]
        v2 = proj[k] - proj[i]
        area = np.sqrt(max(0, np.linalg.norm(v1)**2 * np.linalg.norm(v2)**2 - np.dot(v1, v2)**2))
        if area < 1e-10:
            degenerate += 1
        else:
            dist = np.linalg.norm((proj[i] + proj[j] + proj[k]) / 3)
            if dist > 0.3:
                outside += 1
            else:
                inside += 1

    log(f"  射影後の面の分類:")
    log(f"    退化面（面積≈0）: {degenerate}")
    log(f"    外側の面: {outside}")
    log(f"    内側の面: {inside}")

    face_axis_types = Counter()
    for face in faces:
        axes = tuple(sorted({v//2 for v in face}))
        n_sp = sum(1 for a in axes if a < 3)
        key = f"空間{n_sp}_軸4{'有' if 3 in axes else '無'}_軸5{'有' if 4 in axes else '無'}"
        face_axis_types[key] += 1

    log(f"\n  面の軸構成別分類:")
    for key, count in sorted(face_axis_types.items()):
        log(f"    {key}: {count}面")
    log()


# ============================================================
# 5. 結合行列
# ============================================================

def analyze_coupling_matrix():
    log("=" * 60)
    log("5. 結合行列 C = diag(s₁,...,s₅) による変換")
    log("=" * 60)

    log("  非空間軸（x4,x5）の符号反転パターン:")
    log(f"  {'パターン':>12} {'非空間反転数':>8} {'det':>5}")
    log("  " + "-" * 35)

    for s_space in [(1,1,1), (-1,-1,-1)]:
        for s4 in [1, -1]:
            for s5 in [1, -1]:
                signs = s_space + (s4, s5)
                sign_str = ''.join('+' if s > 0 else '-' for s in signs)
                n_flip = sum(1 for s in (s4, s5) if s < 0)
                det = 1
                for s in signs: det *= s
                log(f"  ({sign_str}): {n_flip:>8} {'+'if det>0 else '-':>5}")
    log()


# ============================================================
# 6. スピン状態の数値検証
# ============================================================

def verify_spin_states():
    from math import comb
    log("=" * 60)
    log("6. スピン状態数の検証（W6 命題3.1, 4.1）")
    log("=" * 60)

    log("  k-ベクトルの独立配向数 C(5,k):")
    total = 0
    for k in range(6):
        c = comb(5, k); total += c
        log(f"    k={k}: C(5,{k}) = {c}")
    log(f"    合計: {total} = 2^5 = {2**5}")

    log("\n  各スピンの状態数:")
    total_states = 0
    for s in [0, 0.5, 1, 1.5, 2, 2.5]:
        n = 2 * int(s) + 1 if s == int(s) else 2 * int(s + 0.5)
        total_states += n
        log(f"    スピン {s}: {n} 状態")
    log(f"    合計: {total_states}")

    log(f"\n  数論的対応:")
    log(f"    21 = T(6) = C(7,2) = {comb(7,2)} ✓")
    log(f"    21 = dim SO(7) = {7*6//2} ✓")
    log(f"    32 = 2^5 = dim Cl(5) ✓")
    log(f"    10 = dim SO(5) = {5*4//2} ✓")
    log()


# ============================================================
# 7. 前川型制約の完全列挙（z3）
# ============================================================

def build_face_edge_map(edges, faces, edge_index):
    """各面を構成する3辺のインデックスを事前計算"""
    face_edges = []
    for face in faces:
        i, j, k = face
        fe = [edge_index[(i,j)], edge_index[(i,k)], edge_index[(j,k)]]
        face_edges.append(fe)
    return face_edges


def maekawa_z3(vertices, edges, faces, edge_index, timeout_sec=300):
    log("=" * 60)
    log("7. 前川型制約の完全列挙（z3 SAT solver）")
    log("=" * 60)

    try:
        import z3
    except ImportError:
        log("  z3未インストール。スキップ。")
        log()
        return None

    start = time.time()

    n_edges = len(edges)
    edge_vars = [z3.Bool(f'e{i}') for i in range(n_edges)]

    solver = z3.Solver()
    solver.set("timeout", timeout_sec * 1000)

    # 前川制約: 各頂点で接続辺の符号和 = ±2
    # Bool: True=+1(山), False=-1(谷)
    # 8辺中 山5谷3 → 和=+2、山3谷5 → 和=-2
    vertex_edges = defaultdict(list)
    for idx, (i, j) in enumerate(edges):
        vertex_edges[i].append(idx)
        vertex_edges[j].append(idx)

    for v_idx in range(len(vertices)):
        connected = vertex_edges[v_idx]
        # 山の数 = 5 or 3
        mountain_count = z3.Sum([z3.If(edge_vars[e], 1, 0) for e in connected])
        solver.add(z3.Or(mountain_count == 5, mountain_count == 3))

    log(f"  変数数: {n_edges} (Bool)")
    log(f"  頂点制約数: {len(vertices)}")
    log(f"  各頂点の次数: 8 → 山5谷3 or 山3谷5")
    log(f"  タイムアウト: {timeout_sec}秒")
    log()

    solutions = []
    max_solutions = 100000
    check_interval = 100

    while len(solutions) < max_solutions:
        result = solver.check()
        if result != z3.sat:
            if result == z3.unknown:
                log(f"  タイムアウト（{len(solutions)}解まで発見）")
            break
        model = solver.model()
        sol = tuple(1 if z3.is_true(model[v]) else -1 for v in edge_vars)
        solutions.append(sol)
        solver.add(z3.Or([v if not z3.is_true(model[v]) else z3.Not(v) for v in edge_vars]))

        if len(solutions) % check_interval == 0:
            elapsed = time.time() - start
            log(f"  ... {len(solutions)}解発見 ({elapsed:.1f}秒)")

    elapsed = time.time() - start
    log(f"\n  前川制約の解の総数: {len(solutions)}")
    if len(solutions) >= max_solutions:
        log(f"  （上限 {max_solutions} に達したため打ち切り）")
    log(f"  計算時間: {elapsed:.2f}秒")

    if solutions:
        plus_counts = [sum(1 for s in sol if s == 1) for sol in solutions]
        log(f"\n  山(+1)の数の分布:")
        for count, freq in sorted(Counter(plus_counts).items()):
            log(f"    山{count}谷{40-count}: {freq}解")

        face_edges = build_face_edge_map(edges, faces, edge_index)
        log(f"\n  面の符号積の分析（最初の解）:")
        sol = solutions[0]
        pos = neg = 0
        for fe in face_edges:
            p = sol[fe[0]] * sol[fe[1]] * sol[fe[2]]
            if p > 0: pos += 1
            else: neg += 1
        log(f"    正の面: {pos}, 負の面: {neg}")

        log(f"\n  超八面体群による同値類の推定:")
        log(f"    群の位数: 2^5 × 5! = 3840")
        log(f"    解の数 / 群の位数 ≈ {len(solutions) / 3840:.1f} 同値類")

    log()
    return solutions


# ============================================================
# 8. 面整合条件（z3）
# ============================================================

def face_consistency_z3(edges, faces, edge_index, timeout_sec=120):
    log("=" * 60)
    log("8. 面整合条件の探索（z3）")
    log("=" * 60)

    try:
        import z3
    except ImportError:
        log("  z3未インストール。スキップ。")
        log()
        return

    face_edges = build_face_edge_map(edges, faces, edge_index)
    n_edges = len(edges)
    start = time.time()

    for target_label, target_val in [("全面符号積=+1", True), ("全面符号積=-1", False)]:
        edge_vars = [z3.Bool(f'e{i}') for i in range(n_edges)]
        solver = z3.Solver()
        solver.set("timeout", timeout_sec * 1000)

        # 3辺の積 = +1 ⟺ 偶数個がFalse
        # 3辺の積 = -1 ⟺ 奇数個がFalse
        for fe in face_edges:
            a, b, c = edge_vars[fe[0]], edge_vars[fe[1]], edge_vars[fe[2]]
            if target_val:  # 積=+1: 0個または2個がFalse
                solver.add(z3.Or(
                    z3.And(a, b, c),
                    z3.And(z3.Not(a), z3.Not(b), c),
                    z3.And(z3.Not(a), b, z3.Not(c)),
                    z3.And(a, z3.Not(b), z3.Not(c)),
                ))
            else:  # 積=-1: 1個または3個がFalse
                solver.add(z3.Or(
                    z3.And(z3.Not(a), z3.Not(b), z3.Not(c)),
                    z3.And(a, b, z3.Not(c)),
                    z3.And(a, z3.Not(b), c),
                    z3.And(z3.Not(a), b, c),
                ))

        result = solver.check()
        if result == z3.sat:
            model = solver.model()
            plus = sum(1 for v in edge_vars if z3.is_true(model[v]))
            log(f"  {target_label}: 解が存在 (山{plus}谷{n_edges-plus})")

            # 解の数を数える
            count = 0
            while solver.check() == z3.sat and count < 10000:
                model = solver.model()
                solver.add(z3.Or([v if not z3.is_true(model[v]) else z3.Not(v) for v in edge_vars]))
                count += 1
                if count % 100 == 0:
                    log(f"    ... {count}解")
            log(f"    解の総数: {count}{'以上' if count >= 10000 else ''}")
        elif result == z3.unknown:
            log(f"  {target_label}: タイムアウト")
        else:
            log(f"  {target_label}: 解なし")

    elapsed = time.time() - start
    log(f"  計算時間: {elapsed:.2f}秒")
    log()


# ============================================================
# 9. 複合制約（前川 + 面整合）
# ============================================================

def combined_z3(vertices, edges, faces, edge_index, timeout_sec=120):
    log("=" * 60)
    log("9. 複合制約（前川 + 面整合）")
    log("=" * 60)

    try:
        import z3
    except ImportError:
        log("  z3未インストール。")
        log()
        return

    n_edges = len(edges)
    edge_vars = [z3.Bool(f'e{i}') for i in range(n_edges)]
    solver = z3.Solver()
    solver.set("timeout", timeout_sec * 1000)

    start = time.time()

    # 前川制約
    vertex_edges = defaultdict(list)
    for idx, (i, j) in enumerate(edges):
        vertex_edges[i].append(idx)
        vertex_edges[j].append(idx)

    for v_idx in range(len(vertices)):
        connected = vertex_edges[v_idx]
        mountain_count = z3.Sum([z3.If(edge_vars[e], 1, 0) for e in connected])
        solver.add(z3.Or(mountain_count == 5, mountain_count == 3))

    # 面整合: 全面の符号積 = +1
    face_edges = build_face_edge_map(edges, faces, edge_index)
    for fe in face_edges:
        a, b, c = edge_vars[fe[0]], edge_vars[fe[1]], edge_vars[fe[2]]
        solver.add(z3.Or(
            z3.And(a, b, c),
            z3.And(z3.Not(a), z3.Not(b), c),
            z3.And(z3.Not(a), b, z3.Not(c)),
            z3.And(a, z3.Not(b), z3.Not(c)),
        ))

    result = solver.check()
    elapsed = time.time() - start

    if result == z3.sat:
        model = solver.model()
        plus = sum(1 for v in edge_vars if z3.is_true(model[v]))
        log(f"  両方満たす解が存在: 山{plus}谷{n_edges-plus}")
        count = 0
        while solver.check() == z3.sat and count < 10000:
            model = solver.model()
            solver.add(z3.Or([v if not z3.is_true(model[v]) else z3.Not(v) for v in edge_vars]))
            count += 1
            if count % 100 == 0:
                log(f"    ... {count}解")
        log(f"  解の総数: {count}{'以上' if count >= 10000 else ''}")
    elif result == z3.unknown:
        log(f"  タイムアウト")
    else:
        log(f"  前川 + 面整合を同時に満たす解は存在しない")

    log(f"  計算時間: {elapsed:.2f}秒")
    log()


# ============================================================
# メイン
# ============================================================

def main():
    log("╔" + "═" * 58 + "╗")
    log("║  5次元正軸体（cross-polytope）数値解析           ║")
    log("║  W6論文の検証と折り紙帰着の探索                  ║")
    log("╚" + "═" * 58 + "╝")
    log()

    vertices, edges, faces, cells, facets_4d, edge_index = build_orthoplex_5d()

    print_basic_structure(vertices, edges, faces, cells, facets_4d)
    edge_types = classify_edges(vertices, edges)
    proj = analyze_projection(vertices, edges, edge_types)
    analyze_face_orientations(vertices, faces, proj)
    analyze_coupling_matrix()
    verify_spin_states()

    log(">>> z3 SAT solver による制約解析を開始します <<<")
    log()
    maekawa_z3(vertices, edges, faces, edge_index, timeout_sec=300)
    face_consistency_z3(edges, faces, edge_index, timeout_sec=120)
    combined_z3(vertices, edges, faces, edge_index, timeout_sec=120)

    log("=" * 60)
    log("解析完了")
    log("=" * 60)


if __name__ == "__main__":
    main()
