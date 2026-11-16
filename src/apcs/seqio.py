"""Minimal sequence IO. No Biopython dependency for the hot paths."""
from __future__ import annotations
import gzip, io, re, os

COMP = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhdbn")

def rc(s: str) -> str:
    return s.translate(COMP)[::-1]

def _open(path):
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"))
    return open(path, "r")

def read_fastq(path):
    """Yield (name, seq, qual). Tolerates blank trailing lines."""
    with _open(path) as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            h = h.rstrip("\r\n")
            if not h:
                continue
            if not h.startswith("@"):
                raise ValueError(f"{path}: expected '@' at record start, got {h[:40]!r}")
            seq = fh.readline().rstrip("\r\n")
            plus = fh.readline().rstrip("\r\n")
            qual = fh.readline().rstrip("\r\n")
            if not h[1:].strip() or not seq or not plus.startswith("+") or len(seq) != len(qual):
                raise ValueError(f"{path}: malformed or truncated FASTQ record {h[:60]!r}")
            if any(ord(c) < 33 or ord(c) > 126 for c in qual):
                raise ValueError(f"{path}: invalid FASTQ quality in {h[:60]!r}")
            if re.search(r"[^ACGTRYMKBDHVNacgtrymkbdhvn]", seq):
                raise ValueError(f"{path}: invalid DNA sequence in {h[:60]!r}")
            yield h[1:].split()[0], seq.upper(), qual

def write_fastq(path, records):
    with open(path, "w") as fh:
        for name, seq, qual in records:
            fh.write(f"@{name}\n{seq}\n+\n{qual}\n")

def read_fasta(path):
    name, buf = None, []
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf).upper()

_ORIGIN = re.compile(r"\nORIGIN(.*)", re.S)
_LOCUS = re.compile(r"^LOCUS\s+(\S+)", re.M)

def read_genbank(path):
    """Return (name, sequence). Features are not parsed - we only need the sequence."""
    with _open(path) as fh:
        text = fh.read()
    if text.count("\nORIGIN") != 1:
        raise ValueError(f"{path}: expected exactly one GenBank record")
    m = _ORIGIN.search(text)
    if not m:
        raise ValueError(f"{path}: no ORIGIN block")
    origin = m.group(1).split("//", 1)[0]
    seq = re.sub(r"[\s\d]", "", origin).upper()
    n = _LOCUS.search(text)
    name = n.group(1) if n else os.path.splitext(os.path.basename(path))[0]
    return name, seq

def load_reference(path):
    """GenBank or FASTA -> (name, sequence)."""
    suffix = str(path).lower().removesuffix(".gz")
    if suffix.endswith((".gb", ".gbk", ".genbank", ".ape")):
        name, seq = read_genbank(path)
    else:
        recs = list(read_fasta(path))
        if len(recs) != 1:
            raise ValueError(f"{path}: expected exactly one FASTA record, got {len(recs)}")
        name, seq = recs[0]
    if not seq or re.search(r"[^ACGT]", seq):
        raise ValueError(f"{path}: reference must contain nonempty, unambiguous A/C/G/T DNA")
    return name, seq

def mismatches(a: str, b: str, limit: int) -> int:
    """Hamming distance, early-exit at limit+1."""
    if len(a) != len(b):
        return limit + 1
    n = 0
    for x, y in zip(a, b):
        if x != y:
            n += 1
            if n > limit:
                return n
    return n
