"""Oligo order sheets for a primer set."""
from __future__ import annotations
import csv, os

NOTE_PAGE = ("Check exact template compatibility and assay performance before ordering. "
             "Purification and synthesis scale are suggestions to confirm with your supplier. "
             "The AmpR design has not been wet-lab validated.")


def rows_for(pset, cfg, scale="10 nmol", purification="PAGE"):
    out = []
    for name, seq in cfg.fwd_index.items():
        out.append({"name": f"{pset.name}-F-{name}", "direction": "forward", "index": name,
                    "index_seq": seq, "sequence": pset.primer_f(seq),
                    "length": len(pset.primer_f(seq)), "scale": scale,
                    "purification": purification, "note": ""})
    for name, seq in cfg.rev_index.items():
        out.append({"name": f"{pset.name}-R-{name}", "direction": "reverse", "index": name,
                    "index_seq": seq, "sequence": pset.primer_r(seq),
                    "length": len(pset.primer_r(seq)), "scale": scale,
                    "purification": purification, "note": ""})
    if out:
        out[0]["note"] = NOTE_PAGE
    return out


def write_csv(path, rows):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def write_vendor_xlsx(path, rows, title="Primer synthesis order"):
    """A sheet with the column headers Chinese oligo vendors use, plus an English sheet."""
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Font, Alignment
    except ImportError:
        raise RuntimeError("openpyxl is required: pip install openpyxl")
    wb = Workbook()
    ws = wb.active
    ws.title = "订购表"
    headers = ["序号", "引物名称", "序列 (5'-3')", "碱基数", "合成规模", "纯化方式", "备注"]
    ws.append(headers)
    for c in ws[1]:
        c.font = Font(bold=True)
    for i, r in enumerate(rows, 1):
        ws.append([i, r["name"], r["sequence"], r["length"], r["scale"], r["purification"], r["note"]])
    widths = [6, 20, 56, 8, 12, 12, 70]
    for i, wd in enumerate(widths, 1):
        ws.column_dimensions[ws.cell(row=1, column=i).column_letter].width = wd
    ws.column_dimensions["C"].width = 56
    for row in ws.iter_rows(min_row=2):
        row[2].font = Font(name="Menlo", size=10)
        row[6].alignment = Alignment(wrap_text=True, vertical="top")

    ws2 = wb.create_sheet("design")
    ws2.append(["field", "value"])
    ws2.append(["title", title])
    ws2.append(["structure", "handle + 6-nt index + arm"])
    ws2.append(["note", NOTE_PAGE])
    for c in ws2["A"]:
        c.font = Font(bold=True)
    ws2.column_dimensions["B"].width = 100
    for r in ws2.iter_rows(min_row=2):
        r[1].alignment = Alignment(wrap_text=True, vertical="top")
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    wb.save(path)
