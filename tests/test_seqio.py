import gzip
import pytest
from apcs.seqio import read_fastq, load_reference, rc


@pytest.mark.parametrize("body", ["@x\nACGT\n+\nIII\n", "@x\nACGT\n-\nIIII\n",
                                  "@x\nACGT\n+\n", "@\nA\n+\nI\n", "@x\nAX\n+\nII\n"])
def test_corrupt_fastq_rejected(tmp_path, body):
    p = tmp_path / "reads.fastq"
    p.write_text(body)
    with pytest.raises(ValueError):
        list(read_fastq(p))


def test_crlf_and_gzip_fastq(tmp_path):
    p = tmp_path / "reads.fastq.gz"
    with gzip.open(p, "wt", newline="") as fh:
        fh.write("@x desc\r\nacgt\r\n+\r\nIIII\r\n\r\n")
    assert list(read_fastq(p)) == [("x", "ACGT", "IIII")]


@pytest.mark.parametrize("suffix,body", [("fa", ">x\nACNGT\n"),
                                          ("fa", ">x\nACGT\n>y\nTGCA\n"),
                                          ("gb", "LOCUS x\nORIGIN\n 1 acngt\n//\n")])
def test_uncertain_reference_cannot_silently_change_coordinates(tmp_path, suffix, body):
    p = tmp_path / f"ref.{suffix}"
    p.write_text(body)
    with pytest.raises(ValueError):
        load_reference(p)


def test_genbank_stops_at_record_terminator(tmp_path):
    p = tmp_path / "ref.gb"
    p.write_text("LOCUS x\nORIGIN\n 1 acgt\n//\nCOMMENT acgt\n")
    assert load_reference(p) == ("x", "ACGT")
    assert rc("ARYN") == "NRYT"
