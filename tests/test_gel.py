import numpy as np
import pytest
from PIL import Image
from apcs import gel


def synthetic_gel(n_lanes=8, h=600, w=800):
    img = np.zeros((h, w))
    pitch = w // (n_lanes + 1)
    for i in range(n_lanes):
        c = pitch * (i + 1)
        img[20:34, c - 14:c + 14] += 0.5                      # loading well
        rows = [120, 200, 300, 420] if i == 0 else [190 + 6 * i]
        for r in rows:
            for dr in range(-7, 8):
                img[r + dr, c - 16:c + 16] += 0.9 * np.exp(-(dr / 4.0) ** 2)
    img += np.random.default_rng(0).normal(0, 0.004, img.shape)
    return np.clip(img, 0, 1)


def test_lane_grid_is_regular():
    img = synthetic_gel()
    lanes = gel.find_lanes(img)
    centres = [(a + b) // 2 for a, b in lanes]
    assert len(lanes) >= 8
    d = np.diff(centres)
    assert d.std() < 2


def test_bands_found_and_well_excluded():
    img = synthetic_gel()
    lanes = gel.find_lanes(img)
    bands = gel.detect_bands(img, lanes[3])
    assert bands and all(b.row > 60 for b in bands)


def test_sizing_against_a_ladder():
    img = synthetic_gel()
    lanes = gel.find_lanes(img)
    lb = gel.detect_bands(img, lanes[0], max_bands=4, prom_frac=0.04)
    rows = sorted(b.row for b in lb)
    sizer = gel.calibrate(rows, [10000, 5000, 2000, 500])
    assert sizer(rows[0]) > sizer(rows[-1])
    assert 400 < sizer(rows[-1]) < 650


def test_uncalibrated_gel_cannot_approve_pooling(tmp_path):
    path = tmp_path / 'gel.png'
    Image.fromarray((synthetic_gel() * 240).astype('uint8')).save(path)
    result = gel.analyse(path, expected_bp=5000, crop=False)
    assert not result['calibrated']
    assert all(lane['verdict'] != 'OK' and lane['pool_uL'] == 0 for lane in result['lanes'])


def test_ladder_in_middle_keeps_sample_labels_in_order(tmp_path):
    img = synthetic_gel()
    # Use the same image but treat the first lane as a sample. The middle
    # ladder assignment is tested with mocked band detection/calibration.
    path = tmp_path / 'gel.png'
    Image.fromarray((img * 240).astype('uint8')).save(path)
    lanes = gel.find_lanes(img)
    labels = [f's{i}' for i in range(len(lanes) - 1)]
    from unittest.mock import patch
    with patch.object(gel, 'detect_bands', return_value=[gel.Band(120, 1, 1),
                                                       gel.Band(200, .8, .8),
                                                       gel.Band(300, .5, .5)]):
        result = gel.analyse(path, ladder_lane=2, ladder_sizes=[10000,5000,2000],
                             sample_labels=labels, crop=False)
    assert [x['label'] for x in result['lanes'] if x['role'] == 'sample'] == labels


def test_missing_ladder_band_and_bad_labels_are_rejected(tmp_path):
    path = tmp_path / 'gel.png'
    Image.fromarray((synthetic_gel() * 240).astype('uint8')).save(path)
    with pytest.raises(ValueError, match='ladder bands'):
        gel.analyse(path, ladder_lane=0, ladder_sizes=[10000,5000,2000,500,100], crop=False)
    with pytest.raises(ValueError, match='unique label'):
        gel.analyse(path, sample_labels=['a1'], crop=False)


def test_non_clipped_image_is_not_rescaled_into_a_clipped_one(tmp_path):
    path = tmp_path / 'dim.png'
    Image.fromarray((synthetic_gel() * 150).astype('uint8')).save(path)
    img = gel.load_image(path)
    assert img.max() < .7
    assert not any(b.saturated for lane in gel.find_lanes(img) for b in gel.detect_bands(img, lane))
