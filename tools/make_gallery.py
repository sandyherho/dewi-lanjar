"""Build the documentation gallery: run every case and write light GIFs to docs/img.

Usage:
    python tools/make_gallery.py                      # run all cases, then shrink
    python tools/make_gallery.py --from-dir outputs   # only shrink existing runs

The full-resolution outputs stay in gallery_runs/ (or --from-dir); docs/img
receives copies resized to --width pixels with every --step-th frame.
"""

import argparse
from pathlib import Path

from PIL import Image, ImageSequence


def shrink(src: Path, dst: Path, width: int, max_frames: int) -> None:
    im = Image.open(src)
    frames = [f.convert("RGB") for f in ImageSequence.Iterator(im)]
    duration = im.info.get("duration", 55)
    step = max(1, -(-len(frames) // max_frames))
    frames = frames[::step]
    h = round(frames[0].height * width / frames[0].width)
    frames = [f.resize((width, h), Image.LANCZOS) for f in frames]
    sample = Image.new("RGB", (width, h * 3))
    for k, i in enumerate((0, len(frames) // 2, len(frames) - 1)):
        sample.paste(frames[i], (0, k * h))
    pal = sample.quantize(colors=255, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE)
    out = [f.quantize(palette=pal, dither=Image.Dither.NONE) for f in frames]
    out[0].save(dst, save_all=True, append_images=out[1:], duration=duration * step, loop=0,
                optimize=True)
    print(f"{dst.name}: {len(out)} frames, {dst.stat().st_size / 1e6:.1f} MB")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from-dir", default=None)
    ap.add_argument("--width", type=int, default=900)
    ap.add_argument("--max-frames", type=int, default=75)
    args = ap.parse_args()
    root = Path(__file__).resolve().parents[1]
    src = Path(args.from_dir) if args.from_dir else root / "gallery_runs"
    if not args.from_dir:
        from dewi_lanjar.cli import main as cli
        cli(["run-all", "--out", str(src), "--quiet"])
    img = root / "docs" / "img"
    img.mkdir(parents=True, exist_ok=True)
    for gif in sorted(src.glob("*/*.gif")):
        shrink(gif, img / gif.name, args.width, args.max_frames)


if __name__ == "__main__":
    main()
