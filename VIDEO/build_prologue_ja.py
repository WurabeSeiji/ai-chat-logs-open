#!/usr/bin/env python3
"""
「創世記・第2章 〜カードと雛壇の2500年〜」プロローグ（日本語版）動画生成スクリプト
- サブピクセル補間（Pillow box crop）による振動ゼロのスムーズなズーム＆パン
- 曲全体の尺（約178秒）に合わせてテロップとシーン切り替えを再配分
"""

import glob
import os
import subprocess
import unicodedata
from PIL import Image, ImageDraw, ImageFont

VIDEO_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_MP4 = os.path.join(VIDEO_DIR, "01_プロローグ_日本語.mp4")

FONT_MINCHO = "/System/Library/Fonts/ヒラギノ明朝 ProN.ttc"
FONT_GOTHIC = "/System/Library/Fonts/ヒラギノ角ゴシック W6.ttc"

W, H = 1920, 1080
FPS = 30
TOTAL_SEC = 178.42
TOTAL_FRAMES = int(round(TOTAL_SEC * FPS))


def find_asset(prefix: str, ext: str) -> str:
    for p in glob.glob(os.path.join(VIDEO_DIR, f"*{ext}")):
        name = unicodedata.normalize("NFC", os.path.basename(p))
        if name.startswith(prefix):
            return p
    raise FileNotFoundError(f"Asset not found: {prefix}*{ext}")


def load_font(path: str, size: int) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(path, size)


def ease_in_out(t: float) -> float:
    t = max(0.0, min(1.0, t))
    return t * t * (3.0 - 2.0 * t)


def draw_text_with_shadow(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    font: ImageFont.FreeTypeFont,
    fill=(255, 255, 255, 255),
    shadow_color=(0, 0, 0, 210),
    stroke_width: int = 3,
    stroke_fill=(10, 12, 22, 230),
    anchor: str = "mm",
):
    x, y = xy
    draw.text((x + 3, y + 4), text, font=font, fill=shadow_color, anchor=anchor)
    draw.text(
        (x, y),
        text,
        font=font,
        fill=fill,
        stroke_width=stroke_width,
        stroke_fill=stroke_fill,
        anchor=anchor,
    )


def create_overlay_image(
    lines: list[str],
    section_label: str = "",
    badge: dict | None = None,
    title_card: tuple[str, str] | None = None,
) -> Image.Image:
    img = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)

    if lines:
        bar_top = 845
        bar_bottom = 1045
        overlay_bar = Image.new("RGBA", (W, H), (0, 0, 0, 0))
        bd = ImageDraw.Draw(overlay_bar)
        bd.rounded_rectangle(
            [(140, bar_top), (W - 140, bar_bottom)],
            radius=22,
            fill=(8, 12, 26, 175),
            outline=(235, 200, 120, 120),
            width=2,
        )
        img = Image.alpha_composite(img, overlay_bar)
        draw = ImageDraw.Draw(img)

        if section_label:
            f_sec = load_font(FONT_GOTHIC, 22)
            draw_text_with_shadow(
                draw,
                (175, bar_top + 26),
                section_label,
                f_sec,
                fill=(245, 215, 140, 240),
                stroke_width=2,
                anchor="lm",
            )

        f_main = load_font(FONT_MINCHO, 44 if len(lines) <= 2 else 38)
        total_lines = len(lines)
        line_h = 58 if total_lines <= 2 else 50
        center_y = (bar_top + bar_bottom) // 2 + (10 if section_label else 0)
        start_y = center_y - (total_lines - 1) * line_h // 2

        for idx, line in enumerate(lines):
            y = start_y + idx * line_h
            draw_text_with_shadow(
                draw,
                (W // 2, y),
                line,
                f_main,
                fill=(255, 255, 255, 255),
                stroke_width=4,
                anchor="mm",
            )

    if title_card:
        main_title, sub_title = title_card
        box = Image.new("RGBA", (W, H), (0, 0, 0, 0))
        bd = ImageDraw.Draw(box)
        bd.rounded_rectangle(
            [(260, 220), (W - 260, 470)],
            radius=28,
            fill=(10, 14, 30, 185),
            outline=(245, 210, 130, 190),
            width=3,
        )
        img = Image.alpha_composite(img, box)
        draw = ImageDraw.Draw(img)
        f_t1 = load_font(FONT_MINCHO, 54)
        f_t2 = load_font(FONT_GOTHIC, 34)
        draw_text_with_shadow(
            draw,
            (W // 2, 305),
            main_title,
            f_t1,
            fill=(255, 240, 195, 255),
            stroke_width=4,
            anchor="mm",
        )
        draw_text_with_shadow(
            draw,
            (W // 2, 395),
            sub_title,
            f_t2,
            fill=(220, 235, 255, 255),
            stroke_width=3,
            anchor="mm",
        )

    if badge:
        kind = badge.get("kind")
        if kind == "symbols":
            items = badge["items"]
            card_w = 240
            gap = 36
            total_w = len(items) * card_w + (len(items) - 1) * gap
            start_x = (W - total_w) // 2
            y_top = 90
            layer = Image.new("RGBA", (W, H), (0, 0, 0, 0))
            ld = ImageDraw.Draw(layer)
            for idx, (_, _, col) in enumerate(items):
                x0 = start_x + idx * (card_w + gap)
                ld.rounded_rectangle(
                    [(x0, y_top), (x0 + card_w, y_top + 175)],
                    radius=22,
                    fill=(12, 16, 36, 195),
                    outline=col,
                    width=3,
                )
            img = Image.alpha_composite(img, layer)
            draw = ImageDraw.Draw(img)
            f_sym = load_font(FONT_MINCHO, 82)
            f_lbl = load_font(FONT_GOTHIC, 24)
            for idx, (sym, lbl, col) in enumerate(items):
                cx = start_x + idx * (card_w + gap) + card_w // 2
                draw_text_with_shadow(
                    draw, (cx, y_top + 68), sym, f_sym, fill=col, stroke_width=4
                )
                draw_text_with_shadow(
                    draw,
                    (cx, y_top + 138),
                    lbl,
                    f_lbl,
                    fill=(240, 245, 255, 245),
                    stroke_width=2,
                )

        elif kind == "banner":
            big = badge["big"]
            small = badge.get("small", "")
            col = badge.get("color", (255, 225, 130, 255))
            layer = Image.new("RGBA", (W, H), (0, 0, 0, 0))
            ld = ImageDraw.Draw(layer)
            ld.rounded_rectangle(
                [(440, 85), (W - 440, 255)],
                radius=26,
                fill=(10, 14, 32, 195),
                outline=col,
                width=3,
            )
            img = Image.alpha_composite(img, layer)
            draw = ImageDraw.Draw(img)
            f_big = load_font(FONT_MINCHO, 64)
            f_small = load_font(FONT_GOTHIC, 28)
            draw_text_with_shadow(
                draw, (W // 2, 150 if small else 170), big, f_big, fill=col, stroke_width=4
            )
            if small:
                draw_text_with_shadow(
                    draw,
                    (W // 2, 215),
                    small,
                    f_small,
                    fill=(235, 242, 255, 250),
                    stroke_width=2,
                )

    return img


def build_cues() -> list[dict]:
    # 元の構成案 (0:00〜2:15 = 135s) を曲尺 (~168s + エンドカード) に合わせてスケーリング
    scale = 1.24
    raw = [
        (0.0, 5.0, "Intro", ["初めに、名前のないカードがあった。"], None,
         ("創世記・第2章 〜カードと雛壇の2500年〜", "プロローグ：神様は足し算ができない")),
        (5.0, 10.0, "Intro", ["カードには数字も模様も書かれておらず、", "裏も表もなかった。"], None, None),
        (10.0, 15.0, "Intro", ["1枚だけでは比べる相手がいないので、", "神様は「2」から始めることにした。"], None, None),
        (15.0, 20.0, "Verse 1", ["実は、神様は計算が大の苦手だった。"], None, None),
        (20.0, 25.0, "Verse 1", ["カードを数えることも、", "面倒な足し算もできない。"], None, None),
        (25.0, 32.0, "Verse 1", ["「ゼロ」や「負の数」「実数」などという", "ややこしいものは、まだ発明すらしていなかった。"], None, None),
        (32.0, 37.0, "Pre-Chorus", ["その代わり、神様には無限の根気があり…"], None, None),
        (37.0, 40.0, "Pre-Chorus", ["「かたまりを作って丸ごとコピペする」"], None, None),
        (40.0, 45.0, "Pre-Chorus", ["この単純な遊び（指数）だけは、", "天才的に得意だった。"], None, None),
        (45.0, 50.0, "Chorus 1", ["神様は2枚のカードをコピペし、", "そのかたまりをまた束ねてコピペした。"],
         {"kind": "banner", "big": "2 → 4 → 8 → 16 ...", "small": "かたまりを丸ごとコピペ（指数）"}, None),
        (50.0, 55.0, "Chorus 1", ["1回コピペするごとに山は倍々に膨れ上がり、", "大小さまざまな山ができた。"],
         {"kind": "banner", "big": "× 2  × 2  × 2 ...", "small": "倍々に膨れ上がるカードの山"}, None),
        (55.0, 60.0, "Chorus 1", ["全体の枚数は数えられなかったが、", "隣り合う山を比べることだけはできた。"], None, None),
        (60.0, 65.0, "Chorus 1", ["そこに初めて、", "宇宙最初の言葉である「有理比」が生まれた。"],
         {"kind": "banner", "big": "「有理比」の誕生", "small": "隣り合う山と山を比べる"}, None),
        (65.0, 72.0, "Verse 2", ["やがて、真っ直ぐコピペを続けるだけでは", "散らかってしまうことに気づいた神様は…"], None, None),
        (72.0, 80.0, "Verse 2", ["コピペの向きを少しひねって、", "ぐるっと輪を描くように並べる遊びを思いついた。"], None, None),
        (80.0, 85.0, "Chorus 2", ["すると不思議なことに、何も計算していないのに", "かたまりの増殖の勢いから e が立ち上がり…"],
         {"kind": "symbols", "items": [("e", "増殖の自然定規", (255, 215, 110, 255))]}, None),
        (85.0, 90.0, "Chorus 2", ["ひねりの直角方向から i が顔を出し…"],
         {"kind": "symbols", "items": [
             ("e", "増殖の自然定規", (255, 215, 110, 255)),
             ("i", "ひねりの直角方向", (120, 225, 255, 255)),
         ]}, None),
        (90.0, 95.0, "Chorus 2", ["一周して最初の山とぴったり噛み合う周期から", "π が勝手に浮かび上がった。"],
         {"kind": "symbols", "items": [
             ("e", "増殖の自然定規", (255, 215, 110, 255)),
             ("i", "ひねりの直角方向", (120, 225, 255, 255)),
             ("π", "一周の周期", (255, 165, 200, 255)),
         ]}, None),
        (95.0, 102.0, "Bridge", ["そして、半周ずれたかたまり同士が正面からぶつかり、", "きれいさっぱり打ち消し合った瞬間…"],
         {"kind": "banner", "big": "e^(iπ) + 1 = 0", "small": "半周の衝突と打ち消し合い", "color": (160, 235, 255, 255)}, None),
        (102.0, 105.0, "Bridge", ["神様は手を叩いて喜んだ。"],
         {"kind": "banner", "big": "0  ＆  負の数", "small": "計算なしで生まれた奇跡", "color": (255, 225, 130, 255)}, None),
        (105.0, 110.0, "Bridge", ["そこで初めて", "「ゼロ」と「負の数」が発明されたのである。"],
         {"kind": "banner", "big": "0  ＆  負の数", "small": "ここに誕生", "color": (255, 225, 130, 255)}, None),
        (110.0, 120.0, "Outro", ["気まぐれに200回ほどコピペを繰り返したころ、", "カードの山は 10⁶⁰ 枚ほどの規模に膨れ上がっていた。"],
         {"kind": "banner", "big": "2²⁰⁰ ≈ 10⁶⁰ 枚", "small": "途方もない宇宙のスケール"}, None),
        (120.0, 133.0, "Outro", ["神様は自分の作った山の大きさに満足し、", "そのまま遊び続けていた。"], None, None),
    ]

    cues = []
    for st, en, sec, lines, badge, title in raw:
        cues.append({
            "start": st * scale,
            "end": en * scale,
            "section": sec,
            "lines": lines,
            "badge": badge,
            "title_card": title,
        })

    # 133 * 1.24 = 164.92s → 残り165〜174.5sでエンドカード
    cues.append({
        "start": 165.0,
        "end": 175.0,
        "section": "",
        "lines": [],
        "badge": None,
        "title_card": ("プロローグ・了", "つづく ─ 第1幕：いきなり核心を突いた男・ピタゴラス"),
    })
    return cues


def load_video_frames(mp4_path: str) -> tuple[list[Image.Image], int, int]:
    """ffmpeg を使って動画の全フレームを RGB Image のリストとして読み込む。"""
    probe = subprocess.check_output(
        [
            "ffprobe", "-v", "error",
            "-select_streams", "v:0",
            "-show_entries", "stream=width,height",
            "-of", "csv=p=0",
            mp4_path,
        ],
        text=True,
    ).strip()
    vw, vh = map(int, probe.split(","))
    raw = subprocess.check_output(
        [
            "ffmpeg", "-v", "error",
            "-i", mp4_path,
            "-f", "rawvideo",
            "-pix_fmt", "rgb24",
            "-",
        ]
    )
    frame_bytes = vw * vh * 3
    n_frames = len(raw) // frame_bytes
    frames = [
        Image.frombytes("RGB", (vw, vh), raw[i * frame_bytes : (i + 1) * frame_bytes])
        for i in range(n_frames)
    ]
    return frames, vw, vh


def render_slow_video_frame(
    frames: list[Image.Image],
    progress: float,
    z0: float = 1.00,
    z1: float = 1.08,
) -> Image.Image:
    """
    動画フレームをスロー再生（前後フレーム補間あり）しつつ 1920x1080 に合わせてズーム描画する。
    """
    n = len(frames)
    p = max(0.0, min(1.0, progress))
    pos = p * (n - 1)
    idx0 = int(pos)
    idx1 = min(n - 1, idx0 + 1)
    frac = pos - idx0

    if frac > 1e-4 and idx0 != idx1:
        frame = Image.blend(frames[idx0], frames[idx1], frac)
    else:
        frame = frames[idx0]

    iw, ih = frame.size
    # 16:9 に合わせてクロップ基準サイズを計算
    target_ar = W / H
    if iw / ih > target_ar:
        base_h = float(ih)
        base_w = ih * target_ar
    else:
        base_w = float(iw)
        base_h = iw / target_ar

    z = z0 + (z1 - z0) * p
    crop_w = base_w / z
    crop_h = base_h / z
    x0 = (iw - crop_w) / 2.0
    y0 = (ih - crop_h) / 2.0

    return frame.resize(
        (W, H),
        resample=Image.Resampling.BILINEAR,
        box=(x0, y0, x0 + crop_w, y0 + crop_h),
    )


def render_pan_zoom(
    img: Image.Image,
    progress: float,
    z0: float,
    z1: float,
    cx0: float,
    cy0: float,
    cx1: float,
    cy1: float,
) -> Image.Image:
    """
    サブピクセル精度の box 指定で、振動ゼロのスムーズなズーム＆パンを描画する。
    cx, cy は正規化中心座標 (0.5, 0.5 が画像中央)
    """
    iw, ih = img.size
    p = max(0.0, min(1.0, progress))
    z = z0 + (z1 - z0) * p
    cx = cx0 + (cx1 - cx0) * p
    cy = cy0 + (cy1 - cy0) * p

    crop_w = iw / z
    crop_h = ih / z

    x0 = cx * iw - crop_w / 2.0
    y0 = cy * ih - crop_h / 2.0
    x0 = max(0.0, min(iw - crop_w, x0))
    y0 = max(0.0, min(ih - crop_h, y0))

    return img.resize(
        (W, H),
        resample=Image.Resampling.BILINEAR,
        box=(x0, y0, x0 + crop_w, y0 + crop_h),
    )


def load_scene_source(prefix: str):
    """同名プレフィックスの .mp4 があればフレームリストを、なければ .jpeg を読み込む。"""
    try:
        mp4_path = find_asset(prefix, ".mp4")
        print(f"Loading video (slow-motion fit): {os.path.basename(mp4_path)}")
        frames, _, _ = load_video_frames(mp4_path)
        print(f"  -> Loaded {len(frames)} frames.")
        return frames
    except FileNotFoundError:
        jpg_path = find_asset(prefix, ".jpeg")
        print(f"Loading image: {os.path.basename(jpg_path)}")
        return Image.open(jpg_path).convert("RGB")


def main():
    src1 = load_scene_source("01) 天国")
    src2 = load_scene_source("02) 神様")
    src3 = load_scene_source("03) 神様")
    src4 = load_scene_source("04) 神様")
    audio = find_asset("01) プロローグ（日本語)", ".mp3")

    # シーン定義（曲尺 178.42s に合わせた割り当て）
    # (start, end, source, z0, z1, cx0, cy0, cx1, cy1)
    xfade_dur = 1.6
    scenes = [
        # Intro (0〜12s, +1.6s xfade)
        (0.0,   12.0 + xfade_dur, src1, 1.00, 1.08, 0.50, 0.50, 0.50, 0.48),
        # Verse 1 前半 (12〜30s)
        (12.0,  30.0 + xfade_dur, src2, 1.02, 1.15, 0.46, 0.50, 0.54, 0.48),
        # Verse 1 後半〜Pre-Chorus〜Chorus 1 (30〜82s)
        (30.0,  82.0 + xfade_dur, src3, 1.00, 1.18, 0.48, 0.52, 0.52, 0.48),
        # Verse 2〜Chorus 2〜Bridge (82〜138s)
        (82.0,  138.0 + xfade_dur, src4, 1.04, 1.20, 0.44, 0.55, 0.52, 0.50),
        # Outro (138〜178.42s) 引き構図へ
        (138.0, TOTAL_SEC + 1.0,  src4, 1.20, 1.00, 0.52, 0.50, 0.50, 0.50),
    ]

    cues = build_cues()
    for cue in cues:
        cue["img"] = create_overlay_image(
            lines=cue.get("lines", []),
            section_label=cue.get("section", ""),
            badge=cue.get("badge"),
            title_card=cue.get("title_card"),
        )

    cmd = [
        "ffmpeg", "-y",
        "-f", "rawvideo",
        "-vcodec", "rawvideo",
        "-pix_fmt", "rgb24",
        "-s", f"{W}x{H}",
        "-r", str(FPS),
        "-i", "-",
        "-i", audio,
        "-c:v", "libx264",
        "-preset", "medium",
        "-crf", "18",
        "-pix_fmt", "yuv420p",
        "-c:a", "aac",
        "-b:a", "192k",
        "-shortest",
        OUT_MP4,
    ]

    print(f"Rendering {TOTAL_FRAMES} frames ({TOTAL_SEC:.2f}s @ {FPS}fps)...")
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE)
    assert proc.stdin is not None

    black = Image.new("RGB", (W, H), (0, 0, 0))
    cue_idx = 0

    for frame_idx in range(TOTAL_FRAMES):
        t = frame_idx / FPS

        # アクティブなシーン（クロスフェード考慮で最大2つ）
        active = []
        for s_idx, (st, en, src, z0, z1, cx0, cy0, cx1, cy1) in enumerate(scenes):
            if st <= t < en:
                prog = (t - st) / max(0.001, en - st)
                if isinstance(src, list):
                    frame_img = render_slow_video_frame(src, prog, z0, z1)
                else:
                    frame_img = render_pan_zoom(src, prog, z0, z1, cx0, cy0, cx1, cy1)
                active.append((s_idx, st, en, frame_img))

        if len(active) == 1:
            canvas = active[0][3]
        elif len(active) >= 2:
            # 後ろのシーンへクロスフェード
            _, _, _, img_prev = active[0]
            _, st_next, _, img_next = active[1]
            alpha = ease_in_out((t - st_next) / xfade_dur)
            canvas = Image.blend(img_prev, img_next, alpha)
        else:
            canvas = black.copy()

        # テロップオーバーレイ
        while cue_idx < len(cues) and cues[cue_idx]["end"] <= t:
            cue_idx += 1
        if cue_idx < len(cues):
            c = cues[cue_idx]
            if c["start"] <= t < c["end"]:
                ov = c["img"]
                canvas.paste(ov, (0, 0), ov)

        # フェードイン・フェードアウト
        if t < 1.2:
            canvas = Image.blend(black, canvas, ease_in_out(t / 1.2))
        elif t > TOTAL_SEC - 3.5:
            rem = max(0.0, TOTAL_SEC - t)
            canvas = Image.blend(black, canvas, ease_in_out(rem / 3.5))

        proc.stdin.write(canvas.tobytes())

    proc.stdin.close()
    ret = proc.wait()
    if ret != 0:
        raise RuntimeError(f"ffmpeg failed with exit code {ret}")
    print("Done! Output:", OUT_MP4)


if __name__ == "__main__":
    main()
