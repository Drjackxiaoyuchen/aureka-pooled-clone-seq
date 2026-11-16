# Contributing

Use Python 3.10 or newer. Install with `python -m pip install '.[dev]' build`, then run `python -m pytest -q`, `python demo/run_demo.py --out fresh-demo-output` and `python -m build`.

For a calling change, add a regression based on an independent expected biological outcome or a known ambiguous case. Do not adjust simulation labels just to make a changed caller pass. Preserve the difference between assigned, mapped and informative molecules, and between reference screening and sequence validation. Include independent molecule counts, ambiguity and reference-panel limitations in reports.

Submit small pull requests with the problem, resulting behaviour and validation. Do not commit patient data, private constructs, raw experimental reads, credentials or lab-identifying metadata without explicit publication rights. Public demos should use generated inputs and clearly label them as simulated.

The project uses Aureka purple #4200AC, mint #6FC0A9 and ink #33266C. Keep status labels readable without colour and preserve printable output. Aureka marks identify this project and are not covered by a grant of trademark rights.
