#!/usr/bin/env python3
"""
「創世記・第2章 〜カードと雛壇の2500年〜」第2幕（日本語版）動画生成スクリプト
- 動画クリップ (21, 23, 25) の音声はカットし、BGM (20) 第2幕：優秀すぎる直系の弟子たち.mp3) に差し替え
- サブピクセル補間による振動ゼロのズーム＆パン＋スロー再生フレーム補間
- 元画像の焼き込みテロップを覆う高不透明度テロップ帯＋上部要約プレートを配置
"""

import glob
import os
import subprocess
import unicodedata
from PIL import Image, ImageDraw, ImageFont

VIDEO_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_MP4 = os.path.join(VIDEO_DIR, "20_第2幕_日本語.mp4")

FONT_MINCHO = "/System/Library/Fonts/ヒラギノ明朝 ProN.ttc"
FONT_GOTHIC = "/System/Library/Fonts/ヒラギノ角ゴシック W6.ttc"

W, H = 1920, 1080
FPS = 30
TOTAL_SEC = 109.84
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


def create_overlay_image(
    lines: list[str],
    section_label: str = "",
    plate: tuple[str, str] | None = None,
    title_card: tuple[str, str] | None = None,
) -> Image.Image:
    img = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)

    # 下部テロップ帯（元画像の焼き込みテロップを覆える高さ＆濃さ）
    if lines:
        bar_top = 868
        bar_bottom = 1054
        overlay_bar = Image.new("RGBA", (W, H), (0, 0, 0, 0))
        bd = ImageDraw.Draw(overlay_bar)
        bd.rounded_rectangle(
            [(120, bar_top), (W - 120, bar_bottom)],
            radius=22,
            fill=(8, 12, 24, 255),
            outline=(235, 200, 120, 145),
            width=2,
        )
        img = Image.alpha_composite(img, overlay_bar)
        draw = ImageDraw.Draw(img)

        if section_label:
            f_sec = load_font(FONT_GOTHIC, 22)
            draw_text_with_shadow(
                draw,
                (155, bar_top + 24),
                section_label,
                f_sec,
                fill=(245, 215, 140, 245),
                stroke_width=2,
                anchor="lm",
            )

        f_main = load_font(FONT_MINCHO, 42 if len(lines) <= 2 else 36)
        total_lines = len(lines)
        line_h = 56 if total_lines <= 2 else 48
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
            [(480, 34), (W - 480, 162)],
            radius=20,
            fill=(10, 14, 30, 195),
            outline=(240, 205, 125, 180),
            width=2,
        )
        img = Image.alpha_composite(img, box)
        draw = ImageDraw.Draw(img)
        f_p1 = load_font(FONT_MINCHO, 38)
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
            [(230, 210), (W - 230, 460)],
            radius=28,
            fill=(10, 14, 30, 195),
            outline=(245, 210, 130, 200),
            width=3,
        )
        img = Image.alpha_composite(img, box)
        draw = ImageDraw.Draw(img)
        f_t1 = load_font(FONT_MINCHO, 52)
        f_t2 = load_font(FONT_GOTHIC, 33)
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
    kepler_plate = ("天体の幾何学（e, i, π）", "ひねりの輪っかから比率を発見")
    newton_plate = ("微分積分という近似ツール", "コピペの段差を滑らかに翻訳")
    einstein_plate = ("光速度不変 ＆ 時空の曲がり", "投影軸のロックと密度の偏り")

    cues = [
        # 01: Intro (0.0 - 7.5s)
        (0.0, 7.5, "Intro",
         ["その後、ピタゴラスの落とした鍵を拾うべく、", "優秀な弟子たちが次々と登場する。"],
         None,
         ("創世記・第2章 〜カードと雛壇の2500年〜", "第2幕：優秀すぎる直系の弟子たち")),

        # 02-04: Verse 1 / ケプラー (7.5 - 30.0s)
        (7.5, 13.5, "Verse 1",
         ["ケプラーは、神様がカードをコピペするときに"],
         kepler_plate, None),
        (13.5, 21.5, "Verse 1",
         ["「ぐるっとひねった輪っか（e, i, π）」の", "幾何学的な美しさに気づき、"],
         kepler_plate, None),
        (21.5, 30.0, "Verse 1",
         ["天体の動きからその比率を掘り起こした。"],
         kepler_plate, None),

        # 05-06: Pre-Chorus / ニュートン (30.0 - 44.5s)
        (30.0, 36.5, "Pre-Chorus",
         ["次に現れたニュートンは、", "神様のコピペの積み重なりを"],
         newton_plate, None),
        (36.5, 44.5, "Pre-Chorus",
         ["「無限に細かく足し合わせる（微分積分）」という", "超便利な近似ツールに翻訳してしまった。"],
         newton_plate, None),

        # 07-08: Chorus 1 / 神様のツッコミ (44.5 - 60.0s)
        (44.5, 52.5, "Chorus 1",
         ["神様は「俺は足し算すら苦手なのに、", "勝手に滑らかな計算機にするなよ！」と"],
         None, None),
        (52.5, 60.0, "Chorus 1",
         ["苦笑いしつつも、その鮮やかな手並みに感心した。"],
         None, None),

        # 09-11: Verse 2 / アインシュタイン (60.0 - 82.0s)
        (60.0, 66.0, "Verse 2",
         ["そして真打ち、", "アインシュタインの登場である。"],
         None, None),
        (66.0, 74.0, "Verse 2",
         ["彼は、神様がひねりを固定した", "「投影軸のロック（光速度不変）」と、"],
         einstein_plate, None),
        (74.0, 82.0, "Verse 2",
         ["かたまりの密度の偏りである", "「時空の曲がり」を見事に数式化した。"],
         einstein_plate, None),

        # 12-13: Bridge (82.0 - 95.5s)
        (82.0, 89.5, "Bridge",
         ["アインシュタインが", "「神はサイコロを振らない」と名言を吐いた瞬間、"],
         None, None),
        (89.5, 95.5, "Bridge",
         ["天界の神様は腹を抱えて笑った。"],
         None, None),

        # 14-16: Outro (95.5 - 106.5s)
        (95.5, 99.5, "Outro",
         ["「いや、当たってるけど惜しい！」"],
         None, None),
        (99.5, 103.0, "Outro",
         ["「サイコロじゃなくて、", "ただ決定論的にカード配って"],
         None, None),
        (103.0, 106.5, "Outro",
         ["コピペしてるだけなんだけどな！」"],
         None, None),

        # End Card (106.5 - 109.5s)
        (106.5, 109.5, "", [], None,
         ("第2幕・了", "つづく ── 第3幕：雛壇のガチャ芸人と「究極のズル」")),
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


def render_slow_video_frame(
    frames: list[Image.Image],
    progress: float,
    z0: float = 1.00,
    z1: float = 1.06,
) -> Image.Image:
    n = len(frames)
    p = max(0.0, min(1.0, progress))
    pos = p * (n - 1)
    idx0 = int(pos)
    idx1 = min(n - 1, idx0 + 1)
    frac = pos - idx0

    frame = Image.blend(frames[idx0], frames[idx1], frac) if (frac > 1e-4 and idx0 != idx1) else frames[idx0]
    iw, ih = frame.size
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
    iw, ih = img.size
    p = max(0.0, min(1.0, progress))
    z = z0 + (z1 - z0) * p
    cx = cx0 + (cx1 - cx0) * p
    cy = cy0 + (cy1 - cy0) * p

    crop_w = iw / z
    crop_h = ih / z
    x0 = max(0.0, min(iw - crop_w, cx * iw - crop_w / 2.0))
    y0 = max(0.0, min(ih - crop_h, cy * ih - crop_h / 2.0))

    return img.resize(
        (W, H),
        resample=Image.Resampling.BILINEAR,
        box=(x0, y0, x0 + crop_w, y0 + crop_h),
    )


def main():
    src20 = load_scene_source("20) 第2幕：優秀すぎる直系の弟子たち")
    src21 = load_scene_source("21)")
    src22 = load_scene_source("22)")
    src23 = load_scene_source("23)")
    src24 = load_scene_source("24)")
    src25 = load_scene_source("25)")
    src26 = load_scene_source("26)")
    bgm = find_asset("20) 第2幕：優秀すぎる直系の弟子たち", ".mp3")

    xfade_dur = 1.1
    scenes = [
        # 20) キーイメージ (0.0 - 7.5s)
        (0.0,   7.5 + xfade_dur,  src20, 1.00, 1.08, 0.50, 0.50, 0.50, 0.48),
        # 21) ケプラー登場.mp4 (7.5 - 18.0s)
        (7.5,   18.0 + xfade_dur, src21, 1.00, 1.05, 0.50, 0.50, 0.50, 0.50),
        # 22) ケプラーの法則発見.jpeg (18.0 - 30.0s)
        (18.0,  30.0 + xfade_dur, src22, 1.00, 1.04, 0.50, 0.50, 0.50, 0.50),
        # 23) ニュートン登場.mp4 (30.0 - 44.5s)
        (30.0,  44.5 + xfade_dur, src23, 1.00, 1.04, 0.50, 0.50, 0.50, 0.50),
        # 24) 勝手に滑らかな計算機にするなよ.jpeg (44.5 - 60.0s)
        (44.5,  60.0 + xfade_dur, src24, 1.00, 1.04, 0.50, 0.50, 0.50, 0.50),
        # 25) アインシュタインとボルン、ボーア.mp4 (60.0 - 82.0s)
        (60.0,  82.0 + xfade_dur, src25, 1.00, 1.06, 0.50, 0.50, 0.50, 0.50),
        # 26) アインシュタインとボルン、ボーアを笑う神様.jpeg (82.0 - 109.84s)
        (82.0,  TOTAL_SEC + 1.0,  src26, 1.00, 1.04, 0.50, 0.50, 0.50, 0.50),
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
                    frame_img = render_slow_video_frame(src, prog, z0, z1)
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

        if t < 1.0:
            canvas = Image.blend(black, canvas, ease_in_out(t / 1.0))
        elif t > TOTAL_SEC - 2.5:
            rem = max(0.0, TOTAL_SEC - t)
            canvas = Image.blend(black, canvas, ease_in_out(rem / 2.5))

        proc.stdin.write(canvas.tobytes())

    proc.stdin.close()
    ret = proc.wait()
    if ret != 0:
        raise RuntimeError(f"ffmpeg failed with exit code {ret}")
    print("Done! Output:", OUT_MP4)


if __name__ == "__main__":
    main()
