"""Generate an explicitly simulated gel; no experimental image is distributed."""
import sys
from pathlib import Path
import numpy as np
from PIL import Image


def main(outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    img = np.zeros((600, 800), dtype=float)
    bands = [[(200, .60)], [(210, .45)], [(120, .7), (200, .7), (300, .7), (420, .7)],
             [(120, .6)], [(200, .6), (300, .5)], [], [(200, 1.5)], [(200, .06)]]
    for i, values in enumerate(bands):
        c = 88 * (i + 1)
        img[20:34, c - 14:c + 14] = .3
        for row, strength in values:
            for dr in range(-10, 11):
                img[row + dr, c - 16:c + 16] += strength * np.exp(-(dr / 4) ** 2)
    img += np.random.default_rng(4).normal(.01, .001, img.shape)
    Image.fromarray((np.clip(img, 0, 1) * 255).astype('uint8')).save(outdir / 'synthetic_gel.png')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).parent)
