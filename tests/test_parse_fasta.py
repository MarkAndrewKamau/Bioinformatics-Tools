import csv
import hashlib
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import parse_fasta


def read_tsv_rows(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def md5_of(sequence: str) -> str:
    return hashlib.md5(sequence.encode("utf-8")).hexdigest()


def test_fasta_multi_record(tmp_path):
    fasta_path = tmp_path / "multi.fa"
    fasta_path.write_text(
        ">seq1 first record\n"
        "ACGTACGT\n"
        ">seq2 second record\n"
        "GGCCAA\n",
        encoding="utf-8",
    )

    out_path = tmp_path / "summary.tsv"
    exit_code = parse_fasta.run([str(fasta_path), "-o", str(out_path)])

    assert exit_code == 0
    assert out_path.exists()

    rows = read_tsv_rows(out_path)
    assert len(rows) == 2

    by_id = {row["id"]: row for row in rows}
    assert by_id["seq1"]["length"] == "8"
    assert by_id["seq1"]["md5"] == md5_of("ACGTACGT")
    assert by_id["seq2"]["length"] == "6"
    assert by_id["seq2"]["md5"] == md5_of("GGCCAA")


def test_min_length_filter(tmp_path):
    fasta_path = tmp_path / "filter.fa"
    fasta_path.write_text(
        ">short\n"
        "ACG\n"
        ">long\n"
        "ACGTACGTAC\n",
        encoding="utf-8",
    )

    out_path = tmp_path / "summary.tsv"
    exit_code = parse_fasta.run(
        [str(fasta_path), "-o", str(out_path), "--min-length", "5"]
    )

    assert exit_code == 0
    rows = read_tsv_rows(out_path)
    assert len(rows) == 1
    assert rows[0]["id"] == "long"
    assert rows[0]["length"] == "10"


def test_genbank_parsing(tmp_path):
    gb_path = tmp_path / "sample.gb"
    record = SeqRecord(Seq("ATGCGCAT"), id="gb1", description="Example GenBank record")
    record.annotations["molecule_type"] = "DNA"
    SeqIO.write([record], gb_path, "genbank")

    out_path = tmp_path / "summary.tsv"
    exit_code = parse_fasta.run([str(gb_path), "-o", str(out_path)])

    assert exit_code == 0
    rows = read_tsv_rows(out_path)
    assert len(rows) == 1

    row = rows[0]
    assert row["id"] == "gb1"
    assert row["format"] == "genbank"
    assert row["mol_type"] == "DNA"
    assert row["gc_percent"] == "50.00"
    assert row["md5"] == md5_of("ATGCGCAT")


def test_bad_file_continues(tmp_path, capsys):
    bad_path = tmp_path / "broken.txt"
    bad_path.write_bytes(b"\x80\x81\x82")

    good_path = tmp_path / "valid.fa"
    good_path.write_text(
        ">ok\n"
        "ACGT\n",
        encoding="utf-8",
    )

    out_path = tmp_path / "summary.tsv"
    exit_code = parse_fasta.run([str(bad_path), str(good_path), "-o", str(out_path)])
    captured = capsys.readouterr()

    assert exit_code == 0
    assert "warning" in captured.err.lower()

    rows = read_tsv_rows(out_path)
    assert len(rows) == 1
    assert rows[0]["id"] == "ok"
