"""Conservative discrimination of indexed circular-plasmid reads.

Minimap2 locates molecules and aligns circular references. Edlib explicitly
scores every diagnostic-window allele. Each molecule contributes at most one
genotype vote. This reference screen does not call consensus or novel variants.
"""
from __future__ import annotations
import os
from collections import Counter, defaultdict
import edlib
import mappy
from .seqio import load_reference, rc

SIMILAR_FRAC = 0.80


def _label(path):
    return os.path.splitext(os.path.basename(str(path).removesuffix(".gz")))[0]


def _breadth(intervals, n):
    """Circular union; each interval has a start and an unwrapped end."""
    flat = []
    for start, end in intervals:
        span = end - start
        if span >= n:
            return n
        if span <= 0:
            continue
        start %= n
        end = start + span
        flat.append((start, min(end, n)))
        if end > n:
            flat.append((0, end - n))
    total, stop = 0, 0
    for start, end in sorted(flat):
        total += max(0, end - max(start, stop))
        stop = max(stop, end)
    return total


class RefPanel:
    def __init__(self, paths, preset="map-ont"):
        self.paths = list(paths)
        if not self.paths:
            raise ValueError("at least one reference is required")
        self.names, self.seqs, self.lengths, self.aligners = [], {}, {}, {}
        for path in self.paths:
            _, seq = load_reference(path)
            name = _label(path)
            if name in self.seqs:
                raise ValueError(f"duplicate reference filename label: {name}")
            if len(seq) < 200:
                raise ValueError("plasmid references must be at least 200 bp")
            self.names.append(name)
            self.seqs[name], self.lengths[name] = seq, len(seq)
            self.aligners[name] = mappy.Aligner(seq=seq + seq, preset=preset, best_n=20)
        self.primary = self.names[0]
        self.posmap, self.oriented = {}, {}
        self.similar, self.sites = self._diagnostics()
        other = [name for name in self.names if name not in self.similar]
        for i, a in enumerate(other):
            for b in other[i + 1:]:
                if any(h.mlen >= SIMILAR_FRAC * min(self.lengths[a], self.lengths[b])
                       for h in self.aligners[a].map(self.seqs[b])):
                    raise ValueError("multiple related backbone groups: analyse panels separately")
        self.clusters = self._windows()

    def _diagnostics(self):
        primary = self.seqs[self.primary]
        n = len(primary)
        similar, sites = {self.primary}, set()
        self.posmap[self.primary] = {p: p for p in range(n)}
        self.oriented[self.primary] = primary
        for name in self.names[1:]:
            original = self.seqs[name]
            candidates = []
            # Move an origin-spanning variant away from the query ends so
            # local alignment cannot silently clip the diagnostic sequence.
            for offset in (0, len(original) // 3, 2 * len(original) // 3):
                query = original[offset:] + original[:offset]
                candidates.extend((h, query) for h in self.aligners[self.primary].map(query))
            if candidates:
                hit, seq = max(candidates, key=lambda item: (item[0].q_en - item[0].q_st, item[0].mlen))
            else:
                hit, seq = None, original
            if hit is None or hit.mlen < SIMILAR_FRAC * min(n, len(seq)):
                continue
            if hit.q_en - hit.q_st < 0.90 * len(seq) or hit.r_en - hit.r_st < 0.90 * n:
                raise ValueError(f"{name}: partial backbone homology; use a separate panel")
            similar.add(name)
            oriented = seq if hit.strand > 0 else rc(seq)
            self.oriented[name] = oriented
            qp = hit.q_st if hit.strand > 0 else len(seq) - hit.q_en
            rp, mapping = hit.r_st, {}
            for length, op in hit.cigar:
                if op == 0:
                    for k in range(length):
                        p = (rp + k) % n
                        mapping[p] = qp + k
                        if primary[p] != oriented[qp + k]:
                            sites.add(p)
                    rp += length
                    qp += length
                elif op == 1:
                    sites.update({(rp - 1) % n, rp % n})
                    qp += length
                elif op == 2:
                    sites.update((rp + k) % n for k in range(length))
                    rp += length
            self.posmap[name] = mapping
        return similar, sorted(sites)

    def _windows(self, gap=200, pad=60):
        if not self.sites:
            return []
        groups = [[self.sites[0]]]
        for pos in self.sites[1:]:
            if pos - groups[-1][-1] <= gap:
                groups[-1].append(pos)
            else:
                groups.append([pos])
        n = self.lengths[self.primary]
        if len(groups) > 1 and groups[0][0] + n - groups[-1][-1] <= gap:
            groups[0] = groups[-1] + [p + n for p in groups[0]]
            groups.pop()
        clusters = []
        for group in groups:
            lo, hi = group[0] - pad, group[-1] + pad + 1
            if hi - lo > n:
                lo, hi = 0, n
            alleles = {}
            for name in sorted(self.similar):
                mapping = self.posmap[name]
                qs = [mapping[p % n] for p in range(lo, hi) if p % n in mapping]
                if len(qs) < 0.60 * (hi - lo):
                    continue
                seq = self.oriented[name]
                length = (qs[-1] - qs[0]) % len(seq) + 1
                start = qs[0] % len(seq)
                allele = (seq + seq)[start:start + length]
                if len(allele) >= 40:
                    alleles[name] = allele
            if len(alleles) == len(self.similar):
                clusters.append({"id": len(clusters), "primary_span": (lo, hi), "seqs": alleles})
        return clusters

    def hits(self, seq):
        results = {}
        for name, aligner in self.aligners.items():
            h = max(aligner.map(seq), key=lambda h: h.mlen, default=None)
            if h is None:
                continue
            n = self.lengths[name]
            results[name] = {"ref": name, "matched": h.mlen,
                             "identity": h.mlen / max(1, h.blen),
                             "r_st": h.r_st % n,
                             "r_en": h.r_st % n + min(n, h.r_en - h.r_st)}
        return results

    def best_hit(self, seq):
        return max(self.hits(seq).values(), key=lambda h: h["matched"], default=None)

    def diagnostic_vote(self, seq, min_margin=2, max_error=0.15):
        """Require a unique winner across near-complete covered windows."""
        totals = Counter({name: 0 for name in self.similar})
        covered = 0
        reverse = rc(seq)
        for cluster in self.clusters:
            scores = {}
            for name, allele in cluster["seqs"].items():
                aligns = [edlib.align(allele, s, mode="HW", task="locations",
                                     k=max(1, int(len(allele) * max_error)))
                          for s in (seq, reverse)]
                valid = [a for a in aligns if a["editDistance"] >= 0
                         and any(end - start + 1 >= 0.90 * len(allele)
                                 for start, end in a["locations"])]
                if valid:
                    scores[name] = min(a["editDistance"] for a in valid)
            if not scores:
                continue
            penalty = max(len(s) for s in cluster["seqs"].values())
            for name in totals:
                totals[name] += scores.get(name, penalty)
            covered += 1
        order = sorted(totals.items(), key=lambda item: (item[1], item[0]))
        if covered and len(order) > 1 and order[1][1] - order[0][1] >= min_margin:
            return order[0][0], covered
        return None, covered


def genotype_well(reads, panel: RefPanel, min_read_len=1000):
    votes = Counter()
    mapped = informative = ambiguous = unmapped = 0
    per_ref_hits = defaultdict(list)
    for _, seq, _ in reads:
        if len(seq) < min_read_len:
            continue
        hits = panel.hits(seq)
        if not hits:
            unmapped += 1
            continue
        mapped += 1
        for name, h in hits.items():
            per_ref_hits[name].append(h)
        best = max(hits.values(), key=lambda h: h["matched"])
        outsiders = [h for name, h in hits.items() if name not in panel.similar]
        if len(panel.names) == 1:
            vote = panel.primary
        elif best["ref"] in panel.similar and len(panel.similar) > 1:
            if outsiders and max(h["matched"] for h in outsiders) >= best["matched"] - 2:
                vote = None
            else:
                vote, _ = panel.diagnostic_vote(seq)
        else:
            ordered = sorted(hits.values(), key=lambda h: -h["matched"])
            vote = (best["ref"] if len(ordered) == 1
                    or ordered[0]["matched"] - ordered[1]["matched"] >= 2 else None)
        if vote is None:
            ambiguous += 1
        else:
            votes[vote] += 1
            informative += 1
    order = sorted(votes.items(), key=lambda item: (-item[1], item[0]))
    best_ref, best_n = order[0] if order else (None, 0)
    second_ref, second_n = order[1] if len(order) > 1 else (None, 0)
    winner_hits = per_ref_hits.get(best_ref, [])
    n = panel.lengths.get(best_ref, panel.lengths[panel.primary])
    intervals = [(h["r_st"], h["r_en"]) for h in winner_hits]
    return {
        "n_reads": mapped, "unmapped": unmapped, "n_ambiguous_reads": ambiguous,
        "mean_identity": round(sum(h["identity"] for h in winner_hits) / max(1, len(winner_hits)), 4),
        "breadth_pct": round(100 * _breadth(intervals, n) / n, 1),
        "mean_depth": round(sum(e - s for s, e in intervals) / n, 1),
        "n_informative_reads": informative,
        "n_site_calls": informative,  # compatibility: now counts molecules
        "read_votes": dict(votes), "ref_counts": dict(votes),
        "best_ref": best_ref, "best_n": best_n,
        "second_ref": second_ref, "second_n": second_n,
        "purity": round(best_n / informative, 4) if informative else 0.0,
        "has_alternatives": len(panel.names) > 1,
    }
