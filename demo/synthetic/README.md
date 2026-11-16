# Synthetic sequencing example

Run `python demo/run_demo.py --out demo-output` from the repository root after installation. This generates random circular references with compatible ColE1 primer sites, omits the designed gap in the amplicon, adds indexed ends, and simulates substitutions, insertions, deletions, end loss and handle erosion. The seeded toy error model is not an empirical nanopore model.

`truth.csv` records the input source and molecule count for each loaded well. The runner checks all 16 outcomes. The report is an illustration, not experimental validation or a clinical/production release criterion. Real sequencing data and internal reference maps are not included.
