#!/usr/bin/env python3
"""
「創世記・第2章 〜カードと雛壇の2500年〜」第3幕〜エピローグ（日本語版）動画生成スクリプト
- キーフレーム静止画 (30〜39) にサブピクセル補間ズーム＆パン＋クロスフェードを適用
- 楽曲 BGM (30) 第3幕：神様のズルい方程式.mp3, 107.96s) に同期
- テロップ帯＋上部要約プレート＋OP/EDタイトルカードを重畳
"""

import glob
import os
import subprocess
import unicodedata
from PIL import Image, ImageDraw, ImageFont

VIDEO_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_MP4 = os.path.join(VIDEO_DIR, "30_第3幕_日本語.mp4")

FONT_MINCHO = "/System/Library/Fonts/ヒラギノ明朝 ProN.ttc"
FONT_GOTHIC = "/System/Library/Fonts/ヒラギノ角ゴシック W6.ttc"

W, H = 1920, 1080
FPS = 30
TOTAL_SEC = 107.96
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
    shadow_color=(0, 0, 0, 220),
    stroke_width: int = 3,
    stroke_fill=(10, 12, 22, 235),
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


def fit_font_size(draw: ImageDraw.ImageDraw, lines: list[str], font_path: str, max_w: int, base_size: int) -> ImageFont.FreeTypeFont:
    size = base_size
    while size > 24:
        f = load_font(font_path, size)
        if all(draw.textlength(line, font=f) <= max_w for line in lines):
            return f
        size -= 2
    return load_font(font_path, 24)


def create_overlay_image(
    lines: list[str],
    section_label: str = "",
    plate: tuple[str, str] | None = None,
    title_card: tuple[str, str] | None = None,
) -> Image.Image:
    img = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)

    # 下部テロップ帯
    if lines:
        bar_top = 868
        bar_bottom = 1054
        overlay_bar = Image.new("RGBA", (W, H), (0, 0, 0, 0))
        bd = ImageDraw.Draw(overlay_bar)
        bd.rounded_rectangle(
            [(100, bar_top), (W - 100, bar_bottom)],
            radius=22,
            fill=(8, 12, 24, 232),
            outline=(235, 200, 120, 155),
            width=2,
        )
        img = Image.alpha_composite(img, overlay_bar)
        draw = ImageDraw.Draw(img)

        if section_label:
            f_sec = load_font(FONT_GOTHIC, 22)
            draw_text_with_shadow(
                draw,
                (136, bar_top + 24),
                section_label,
                f_sec,
                fill=(245, 215, 140, 245),
                stroke_width=2,
                anchor="lm",
            )

        base_sz = 40 if len(lines) <= 2 else 34
        f_main = fit_font_size(draw, lines, FONT_MINCHO, W - 280, base_sz)
        total_lines = len(lines)
        line_h = 54 if total_lines <= 2 else 46
        center_y = (bar_top + bar_bottom) // 2 + (8 if section_label else 0)
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

    # 上部要約プレート
    if plate:
        p_title, p_sub = plate
        box = Image.new("RGBA", (W, H), (0, 0, 0, 0))
        bd = ImageDraw.Draw(box)
        bd.rounded_rectangle(
            [(460, 34), (W - 460, 162)],
            radius=20,
            fill=(10, 14, 30, 205),
            outline=(240, 205, 125, 185),
            width=2,
        )
        img = Image.alpha_composite(img, box)
        draw = ImageDraw.Draw(img)
        f_p1 = load_font(FONT_MINCHO, 36)
        f_p2 = load_font(FONT_GOTHIC, 24)
        draw_text_with_shadow(
            draw,
            (W // 2, 76),
            p_title,
            f_p1,
            fill=(255, 232, 155, 255),
            stroke_width=3,
            anchor="mm",
        )
        draw_text_with_shadow(
            draw,
            (W // 2, 128),
            p_sub,
            f_p2,
            fill=(230, 240, 255, 250),
            stroke_width=2,
            anchor="mm",
        )

    # オープニング／エンディング タイトルカード
    if title_card:
        main_title, sub_title = title_card
        box = Image.new("RGBA", (W, H), (0, 0, 0, 0))
        bd = ImageDraw.Draw(box)
        bd.rounded_rectangle(
            [(210, 210), (W - 210, 460)],
            radius=28,
            fill=(10, 14, 30, 205),
            outline=(245, 210, 130, 205),
            width=3,
        )
        img = Image.alpha_composite(img, box)
        draw = ImageDraw.Draw(img)
        f_t1 = fit_font_size(draw, [main_title], FONT_MINCHO, W - 480, 50)
        f_t2 = fit_font_size(draw, [sub_title], FONT_GOTHIC, W - 480, 32)
        draw_text_with_shadow(
            draw,
            (W // 2, 295),
            main_title,
            f_t1,
            fill=(255, 240, 195, 255),
            stroke_width=4,
            anchor="mm",
        )
        draw_text_with_shadow(
            draw,
            (W // 2, 385),
            sub_title,
            f_t2,
            fill=(220, 235, 255, 255),
            stroke_width=3,
            anchor="mm",
        )

    return img


def build_cues() -> list[dict]:
    gacha_plate = ("加速器ガチャと標準模型", "U(1)×SU(2)×SU(3) と 19パラメータ")
    god_plate = ("たった1枚のカード遊び", "コピペとひねりの断面")
    renorm_plate = ("連続体の罠と「究極のズル」", "∞ − ∞ ＝ finite（繰り込み）")
    born_plate = ("有限のカード枚数なら発散しない", "決定論的カード vs Born則の cos² 波")
    epilogue_plate = ("エピローグ：原点回帰", "八元数・複素三角形・ピタゴラスの比")

    cues = [
        # 01: Intro (0.0 - 6.0s) / Title Card (0.0 - 4.0s)
        (0.0, 4.0, "Intro",
         ["時代は一気に現代へ。", "舞台は巨大な地下トンネルの加速器と、"],
         None,
         ("創世記・第2章 〜カードと雛壇の2500年〜", "第3幕：雛壇のガチャ芸人と「究極のズル」〜エピローグ")),
        (4.0, 6.0, "Intro",
         ["時代は一気に現代へ。", "舞台は巨大な地下トンネルの加速器と、"],
         None, None),

        # 02: Intro (6.0 - 11.0s)
        (6.0, 11.0, "Intro",
         ["黒板がずらりと並ぶスタジオの「雛壇」へと変わる。", "そこに座るのは、理論物理学者という名の雛壇芸人たち。"],
         None, None),

        # 03-04: Intro Break & Verse 1 (11.0 - 21.0s)
        (11.0, 16.0, "Intro Break",
         ["さあ始まるよ、運命のガチャタイム！", "神様が気まぐれに、カードの山をちょっとひねった。"],
         None, None),
        (16.0, 21.0, "Verse 1",
         ["その瞬間、巨大加速器から眩い粒子が噴き出し、", "天使も白衣の研究者も、雛壇の学者たちも大騒ぎ！"],
         gacha_plate, None),

        # 05-09: Verse 1 / Pre-Chorus 1 / Chorus 1 (21.0 - 51.0s)
        (21.0, 29.0, "Verse 1",
         ["彼らは「コピペの重なり」の断面をガチャにかけ、", "出てくる破片を見ては黒板の前で毎回大騒ぎしている。"],
         gacha_plate, None),
        (29.0, 36.0, "Pre-Chorus 1",
         ["「うおおおっ！ 対称性の破れから", "U(1)×SU(2)×SU(3) が出たぞ！ 標準理論コンプリートだ！」"],
         gacha_plate, None),
        (36.0, 41.0, "Chorus 1",
         ["「いや待て、まだ重力のSSRが排出されてない！」"],
         gacha_plate, None),
        (41.0, 46.0, "Chorus 1",
         ["「確率が渋すぎる！ 次元を11個に増やして", "超ひも理論の10連ガチャを回せば出るはずだ！」"],
         gacha_plate, None),
        (46.0, 51.0, "Chorus 1",
         ["「なんで標準模型のパラメータが19個もあるんだ、", "神様の微調整が鬼畜すぎる！」"],
         gacha_plate, None),

        # 10-12: Post-Chorus & Verse 2 (51.0 - 63.0s)
        (51.0, 55.0, "Post-Chorus",
         ["汗だくで黒板に複雑な数式を書き殴る姿を、", "天界の神様と天使たちは笑い転げて見ていた。"],
         god_plate, None),
        (55.0, 59.0, "Post-Chorus",
         ["「ただコピペしてひねっただけの模様なのに、", "なんで19個に分解して勝手に苦しんでんの！？」"],
         god_plate, None),
        (59.0, 63.0, "Verse 2",
         ["さらに「宇宙は連続で滑らかだ」と思い込んだため、", "計算のあちこちで存在しないはずの「無限」が発散した。"],
         god_plate, None),

        # 13-15: Pre-Chorus 2 / 繰り込みとノーベル賞 (63.0 - 81.0s)
        (63.0, 69.0, "Pre-Chorus 2",
         ["そこで雛壇芸人たちは、禁断の技に手を染める。", "無限大から無限大を引く「繰り込み」という究極のズル！"],
         renorm_plate, None),
        (69.0, 75.0, "Pre-Chorus 2",
         ["しかもあまつさえ、そのズルを発明した者たちは", "誇らしげに金メダルを掲げて大喜びしているのだ。"],
         renorm_plate, None),
        (75.0, 81.0, "Pre-Chorus 2",
         ["「なんでそれで祝ってるの？」と首をかしげる白衣の天使。", "両脇の天使も、不思議そうな顔でお祝いのラッパを吹く。"],
         renorm_plate, None),

        # 16-17: Chorus 2 / ポップコーン爆笑 (81.0 - 91.0s)
        (81.0, 86.0, "Chorus 2",
         ["これにはたまらず、神様も天使たちも涙を流して笑い転げた。", "「いや、離散的なカードの枚数で止めとけば発散しないのに！」"],
         None, None),
        (86.0, 91.0, "Chorus 2",
         ["「自分たちで勝手に無限を作っといて、", "自分で消しゴムで消してドヤ顔するなよ！！」"],
         None, None),

        # 18-19: Bridge / 決定論 vs Born則 (91.0 - 98.0s)
        (91.0, 94.5, "Bridge",
         ["神様の手元では、一糸乱れぬ決定論的な手順通りに", "カードが整然と引かれているだけだった。"],
         born_plate, None),
        (94.5, 98.0, "Bridge",
         ["なのに人間たちはそれをランダムにサンプリングし、", "Born則の cos² 波という「確率」として読んでいた。"],
         born_plate, None),

        # 20: Epilogue (98.0 - 103.0s)
        (98.0, 103.0, "Epilogue",
         ["喧騒の講義室で、彼女が複雑な数式をすべて消し去ると、", "八元数と複素頂点の三角形、そして「3・4・5」が浮かび上がった！"],
         epilogue_plate, None),

        # 21-22: Outro (103.0 - 107.96s)
        (103.0, 105.6, "Outro",
         ["「……これ、無名のカードをコピペしただけの『比』だ！」", "天界の神様も天使たちも、身を乗り出して大喜びの拍手！"],
         None, None),
        (105.6, TOTAL_SEC, "Outro",
         ["神様は優しく微笑み、こう呟いた。", "──「……ようやく気づいたか」"],
         None,
         ("創世記・第2章 〜カードと雛壇の2500年〜　完", "── 本編「背景時空を仮定しない数論・代数からの時空創発」へ")),
    ]
    return [
        {"start": st, "end": en, "section": sec, "lines": lines, "plate": pl, "title_card": tc}
        for st, en, sec, lines, pl, tc in cues
    ]


def load_video_frames(mp4_path: str) -> list[Image.Image]:
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
    return [
        Image.frombytes("RGB", (vw, vh), raw[i * frame_bytes : (i + 1) * frame_bytes])
        for i in range(n_frames)
    ]


def load_scene_source(prefix: str):
    try:
        mp4_path = find_asset(prefix, ".mp4")
        print(f"Loading video (audio stripped): {os.path.basename(mp4_path)}")
        frames = load_video_frames(mp4_path)
        print(f"  -> Loaded {len(frames)} frames.")
        return frames
    except FileNotFoundError:
        jpg_path = find_asset(prefix, ".jpeg")
        print(f"Loading image: {os.path.basename(jpg_path)}")
        return Image.open(jpg_path).convert("RGB")


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
    iw, ih = img.size
    target_ar = W / H
    if iw / ih > target_ar:
        base_h = float(ih)
        base_w = ih * target_ar
    else:
        base_w = float(iw)
        base_h = iw / target_ar

    p = ease_in_out(max(0.0, min(1.0, progress)))
    z = z0 + (z1 - z0) * p
    cx = cx0 + (cx1 - cx0) * p
    cy = cy0 + (cy1 - cy0) * p

    crop_w = base_w / z
    crop_h = base_h / z
    x0 = max(0.0, min(iw - crop_w, cx * iw - crop_w / 2.0))
    y0 = max(0.0, min(ih - crop_h, cy * ih - crop_h / 2.0))

    return img.resize(
        (W, H),
        resample=Image.Resampling.BILINEAR,
        box=(x0, y0, x0 + crop_w, y0 + crop_h),
    )


def render_slow_video_frame(
    frames: list[Image.Image],
    progress: float,
    z0: float,
    z1: float,
    cx0: float,
    cy0: float,
    cx1: float,
    cy1: float,
) -> Image.Image:
    n = len(frames)
    p = max(0.0, min(1.0, progress))
    pos = p * (n - 1)
    idx0 = int(pos)
    idx1 = min(n - 1, idx0 + 1)
    frac = pos - idx0

    frame = Image.blend(frames[idx0], frames[idx1], frac) if (frac > 1e-4 and idx0 != idx1) else frames[idx0]
    return render_pan_zoom(frame, progress, z0, z1, cx0, cy0, cx1, cy1)


def main():
    src30 = load_scene_source("30)")
    src31 = load_scene_source("31)")
    src32 = load_scene_source("32)")
    src33 = load_scene_source("33)")
    src34 = load_scene_source("34)")
    src35 = load_scene_source("35)")
    src36 = load_scene_source("36)")
    src37 = load_scene_source("37)")
    src38 = load_scene_source("38)")
    src39 = load_scene_source("39)")
    bgm = find_asset("30)", ".mp3")

    xfade_dur = 0.85
    scenes = [
        (0.0,   6.0  + xfade_dur, src30, 1.00, 1.07, 0.50, 0.50, 0.56, 0.48),
        (6.0,   11.0 + xfade_dur, src31, 1.00, 1.07, 0.50, 0.50, 0.52, 0.48),
        (11.0,  21.0 + xfade_dur, src32, 1.00, 1.08, 0.45, 0.50, 0.56, 0.50),
        (21.0,  51.0 + xfade_dur, src33, 1.00, 1.09, 0.48, 0.50, 0.55, 0.47),
        (51.0,  63.0 + xfade_dur, src34, 1.00, 1.07, 0.46, 0.50, 0.53, 0.48),
        (63.0,  81.0 + xfade_dur, src35, 1.00, 1.07, 0.50, 0.50, 0.48, 0.52),
        (81.0,  91.0 + xfade_dur, src36, 1.00, 1.07, 0.48, 0.52, 0.45, 0.48),
        (91.0,  98.0 + xfade_dur, src37, 1.00, 1.08, 0.45, 0.55, 0.56, 0.46),
        (98.0,  103.0 + xfade_dur, src38, 1.01, 1.07, 0.50, 0.50, 0.52, 0.48),
        (103.0, TOTAL_SEC + 1.0,  src39, 1.00, 1.05, 0.48, 0.50, 0.50, 0.50),
    ]

    cues = build_cues()
    for cue in cues:
        cue["img"] = create_overlay_image(
            lines=cue.get("lines", []),
            section_label=cue.get("section", ""),
            plate=cue.get("plate"),
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
        "-i", bgm,
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

        active = []
        for s_idx, (st, en, src, z0, z1, cx0, cy0, cx1, cy1) in enumerate(scenes):
            if st <= t < en:
                prog = (t - st) / max(0.001, en - st)
                if isinstance(src, list):
                    frame_img = render_slow_video_frame(src, prog, z0, z1, cx0, cy0, cx1, cy1)
                else:
                    frame_img = render_pan_zoom(src, prog, z0, z1, cx0, cy0, cx1, cy1)
                active.append((s_idx, st, en, frame_img))

        if len(active) == 1:
            canvas = active[0][3]
        elif len(active) >= 2:
            _, _, _, img_prev = active[0]
            _, st_next, _, img_next = active[1]
            alpha = ease_in_out((t - st_next) / xfade_dur)
            canvas = Image.blend(img_prev, img_next, alpha)
        else:
            canvas = black.copy()

        while cue_idx < len(cues) and cues[cue_idx]["end"] <= t:
            cue_idx += 1
        if cue_idx < len(cues):
            c = cues[cue_idx]
            if c["start"] <= t < c["end"]:
                ov = c["img"]
                canvas.paste(ov, (0, 0), ov)

        if t < 0.8:
            canvas = Image.blend(black, canvas, ease_in_out(t / 0.8))
        elif t > TOTAL_SEC - 2.0:
            rem = max(0.0, TOTAL_SEC - t)
            canvas = Image.blend(black, canvas, ease_in_out(rem / 2.0))

        proc.stdin.write(canvas.tobytes())

    proc.stdin.close()
    ret = proc.wait()
    if ret != 0:
        raise RuntimeError(f"ffmpeg failed with exit code {ret}")
    print("Done! Output:", OUT_MP4)


if __name__ == "__main__":
    main()
