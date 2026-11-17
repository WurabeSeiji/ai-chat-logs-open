#!/usr/bin/env python3
"""
「創世記・第2章 〜カードと雛壇の2500年〜」第1幕（日本語版）動画生成スクリプト
- 動画クリップの音声はカットし、BGM (10) 第1幕：万物は数.mp3) に差し替え
- サブピクセル補間による振動ゼロのズーム＆パン＋スロー再生フレーム補間
- BGM 尺（約69.25秒）に合わせたテロップ＆シーン構成
"""

import glob
import os
import subprocess
import unicodedata
from PIL import Image, ImageDraw, ImageFont

VIDEO_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_MP4 = os.path.join(VIDEO_DIR, "10_第1幕_日本語.mp4")

FONT_MINCHO = "/System/Library/Fonts/ヒラギノ明朝 ProN.ttc"
FONT_GOTHIC = "/System/Library/Fonts/ヒラギノ角ゴシック W6.ttc"

W, H = 1920, 1080
FPS = 30
TOTAL_SEC = 69.25
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
            fill=(8, 12, 26, 180),
            outline=(235, 200, 120, 125),
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
            [(240, 210), (W - 240, 460)],
            radius=28,
            fill=(10, 14, 30, 185),
            outline=(245, 210, 130, 190),
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
    cues = [
        # Intro (0.0 - 9.5s)
        (0.0, 3.5, "Intro", ["それからずっと後のこと。"],
         ("創世記・第2章 〜カードと雛壇の2500年〜", "第1幕：いきなり核心を突いた男・ピタゴラス")),
        (3.5, 6.5, "Intro", ["巨大なカードの山の内側に閉じ込められた", "小さな住人たちは、"], None),
        (6.5, 9.5, "Intro", ["自分たちの世界がどうやってできているのか", "悩み始めた。"], None),

        # Verse 1 (9.5 - 18.0s)
        (9.5, 13.8, "Verse 1", ["天界の特等席でポップコーンを片手に", "モニターを覗き込む神様を"], None),
        (13.8, 18.0, "Verse 1", ["一番ドキッとさせたのは、", "古代ギリシャのピタゴラスだった。"], None),

        # Pre-Chorus (18.0 - 25.5s)
        (18.0, 20.5, "Pre-Chorus", ["彼は世界をじっと見つめ、"], None),
        (20.5, 22.8, "Pre-Chorus", ["ある日突然叫んだ。"], None),
        (22.8, 25.5, "Pre-Chorus", ["「万物は数（整数の比）なり！」"], None),

        # Chorus 1 (25.5 - 40.0s)
        (25.5, 29.5, "Chorus 1", ["神様が足し算すらできずに", "ただ並べて作った「有理比」の正体に、"], None),
        (29.5, 33.5, "Chorus 1", ["人類で一番最初に、", "しかもほぼ一発で辿り着いたのだ。"], None),
        (33.5, 36.5, "Chorus 1", ["天界の神様は思わず身を乗り出した。"], None),
        (36.5, 40.0, "Chorus 1", ["「おいおい、いきなりネタバレかよ、", "こいつが一番偉いな！」"], None),

        # Verse 2 / Rap (40.0 - 49.5s)
        (40.0, 44.0, "Verse 2", ["ところが、神様が200回もコピペを重ねて", "山の粒が細かくなりすぎていたせいで……"], None),
        (44.0, 46.8, "Verse 2", ["ピタゴラスの弟子たちがパニックを起こしてしまう。"], None),
        (46.8, 49.5, "Verse 2", ["「対角線の比（√2）が整数比で割り切れない！」"], None),

        # Bridge (49.5 - 57.5s)
        (49.5, 53.0, "Bridge", ["せっかく正解の「有理比」を", "握りしめていたのに、"], None),
        (53.0, 57.5, "Bridge", ["彼らは「滑らかな連続体（無理数・実数）」という", "錯覚の方へ逃げ出してしまった。"], None),

        # Outro (57.5 - 66.0s)
        (57.5, 60.0, "Outro", ["神様はモニターの前でズッコケた。"], None),
        (60.0, 62.2, "Outro", ["「あーあ、惜しい！」"], None),
        (62.2, 66.0, "Outro", ["「ただ粒が細かいだけなのに、", "勝手に『連続体』なんて発明しちゃったよ！」"], None),

        # End Card (66.0 - 69.0s)
        (66.0, 69.0, "", [],
         ("第1幕・了", "つづく ─ 第2幕：優秀すぎる直系の弟子たち")),
    ]
    return [
        {"start": st, "end": en, "section": sec, "lines": lines, "title_card": tc}
        for st, en, sec, lines, tc in cues
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
    src10 = load_scene_source("10) 第1幕：ピタゴラスと天使たち")
    src11 = load_scene_source("11)")
    src12 = load_scene_source("12)")
    src13 = load_scene_source("13)")
    src14 = load_scene_source("14)")
    src15 = load_scene_source("15)")
    bgm = find_asset("10) 第1幕：万物は数", ".mp3")

    xfade_dur = 1.0
    scenes = [
        # 10) Intro (0.0 - 9.5s)
        (0.0,   9.5 + xfade_dur,  src10, 1.00, 1.06, 0.50, 0.50, 0.50, 0.48),
        # 11) 悩み考える古代の人々 (9.5 - 17.5s)
        (9.5,   17.5 + xfade_dur, src11, 1.00, 1.12, 0.50, 0.52, 0.50, 0.48),
        # 12) 万物は数なり (17.5 - 26.0s)
        (17.5,  26.0 + xfade_dur, src12, 1.02, 1.12, 0.48, 0.50, 0.54, 0.48),
        # 13) 身を乗り出して拍手する神様 (26.0 - 40.0s)
        (26.0,  40.0 + xfade_dur, src13, 1.00, 1.12, 0.50, 0.50, 0.48, 0.48),
        # 14) 弟子たちが頭を抱えるパニック.mp4 (40.0 - 57.5s)
        (40.0,  57.5 + xfade_dur, src14, 1.00, 1.06, 0.50, 0.50, 0.50, 0.50),
        # 15) モニターの前でズッコケる神様 (57.5 - 69.25s)
        (57.5,  TOTAL_SEC + 1.0,  src15, 1.02, 1.10, 0.50, 0.50, 0.48, 0.48),
    ]

    cues = build_cues()
    for cue in cues:
        cue["img"] = create_overlay_image(
            lines=cue.get("lines", []),
            section_label=cue.get("section", ""),
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
