"""Create a production-sized compressed TIFF from Figure4_production.tiff.

This helper performs no scientific analysis. It only uniformly downsamples the
already-rendered Figure 4 raster, preserves aspect ratio, writes 600-dpi
metadata, and uses TIFF LZW compression. The default 6000-pixel width remains
well above journal resolution requirements at typical two-column print widths.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from PIL import Image


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        default=Path("production/figure4/Figure4_production.tiff"),
    )
    parser.add_argument(
        "output",
        nargs="?",
        type=Path,
        default=Path("production/figure4/Figure4_production_600dpi_lzw.tif"),
    )
    parser.add_argument("--max-width", type=int, default=6000)
    parser.add_argument("--dpi", type=int, default=600)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    src = args.input.resolve()
    dst = args.output.resolve()

    if not src.exists():
        raise FileNotFoundError(src)

    with Image.open(src) as img:
        img.load()
        original_size = img.size

        if img.width > args.max_width:
            scale = args.max_width / img.width
            new_size = (
                args.max_width,
                max(1, round(img.height * scale)),
            )
            img = img.resize(new_size, Image.Resampling.LANCZOS)
        else:
            new_size = img.size

        if img.mode not in ("RGB", "L"):
            img = img.convert("RGB")

        dst.parent.mkdir(parents=True, exist_ok=True)
        img.save(
            dst,
            format="TIFF",
            compression="tiff_lzw",
            dpi=(args.dpi, args.dpi),
        )

    print(f"Input:  {src}")
    print(f"Output: {dst}")
    print(f"Pixels: {original_size[0]}x{original_size[1]} -> {new_size[0]}x{new_size[1]}")
    print(f"DPI metadata: {args.dpi} x {args.dpi}")
    print(f"Output size: {dst.stat().st_size / (1024 * 1024):.1f} MiB")
    print("Only uniform raster downsampling and TIFF LZW compression were applied.")


if __name__ == "__main__":
    main()
