# Simulated gel example

Run through `demo/run_demo.py`. The image is generated with seed 4; the ladder is lane 2, with four bands assigned 10000, 5000, 2000 and 500 bp. Expected sample size is 5000 bp. The seven sample lanes demonstrate two passes, wrong size, multiple bands, an empty lane, image clipping and a volume-cap review. The geometry, intensities and band sizes are artificial.

The lane detector and calibration must be visually checked on real images. The demo does not establish sizing accuracy, concentration accuracy or sensitivity to weak secondary bands.
