# Primer order data

CSV sheets are generated from the canonical package configuration. All sequences are 5′→3′. Regenerate CSV and XLSX with `apcs order --set cole1 -o orders` or `--set ampr`. Check template compatibility and supplier requirements before ordering. AmpR is a proposed design with unvalidated wet-lab performance. Purification and scale are suggestions, not an established service agreement.

For the preferred A1–A4, B1–B4, … H1–H4 layout, use [the 12-primer order list](ColE1_ORDER_32well.txt) or the [order workbook and well map](Aureka_ColE1_32well_Order.xlsx). The workbook includes 12 required oligos and two explicitly optional reverse primers. The default CLI export contains the 12 required primers; `--columns 6` includes all 14. Existing `order_cole1_primers.csv` is the complete 14-oligo catalog.
