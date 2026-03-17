#!/usr/bin/env python3
"""Parse FASTA/GenBank files and write a summary TSV."""

import argparse
import csv
import hashlib
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

__version__ = "0.1.0"

FASTA_EXTENSIONS = {".fa", ".fasta", ".fas", ".fna", ".faa"}
GENBANK_EXTENSIONS = {".gb", ".gbk", ".genbank"}
SCAN_EXTENSIONS = FASTA_EXTENSIONS | GENBANK_EXTENSIONS

NUCLEOTIDE_STRICT = {"A", "C", "G", "T", "U"}
NUCLEOTIDE_IUPAC = {"A", "C", "G", "T", "U", "N", "R", "Y", "K", "M", "S", "W", "B", "D", "H", "V"}
PROTEIN_MARKERS = {"E", "F", "I", "L", "P", "Q", "Z", "J", "O", "X"}

TSV_HEADER = [
    "id",
    "source_file",
    "format",
    "description",
    "length",
    "mol_type",
    "gc_percent",
    "md5",
]


class ParseFileError(Exception):
    """Raised when a sequence file cannot be parsed as FASTA or GenBank."""


def warn(message: str) -> None:
    """Write a warning to stderr."""
    print(f"Warning: {message}", file=sys.stderr)


def detect_format(path: Path) -> Optional[str]:
    """Detect FASTA/GenBank format from file extension."""
    suffix = path.suffix.lower()
    if suffix in FASTA_EXTENSIONS:
        return "fasta"
    if suffix in GENBANK_EXTENSIONS:
        return "genbank"
    return None


def list_input_files(paths: Sequence[str]) -> List[Path]:
    """Return sorted unique files from input paths (files and directories)."""
    files = set()

    for raw_path in paths:
        path = Path(raw_path)
        if path.is_file():
            files.add(path)
            continue

        if path.is_dir():
            found_any = False
            for child in path.rglob("*"):
                if child.is_file() and child.suffix.lower() in SCAN_EXTENSIONS:
                    files.add(child)
                    found_any = True
            if not found_any:
                warn(f"No FASTA/GenBank files found in directory: {path}")
            continue

        warn(f"Input path not found: {path}")

    return sorted(files, key=lambda p: str(p))


def _parse_records_with_format(path: Path, file_format: str) -> Tuple[List[SeqRecord], Optional[Exception]]:
    """Attempt to parse all records for one file-format pair."""
    records: List[SeqRecord] = []
    try:
        with path.open("r", encoding="utf-8") as handle:
            for record in SeqIO.parse(handle, file_format):
                records.append(record)
    except Exception as exc:  # pylint: disable=broad-except
        return [], exc

    return records, None


def _candidate_formats(path: Path, format_hint: Optional[str]) -> List[str]:
    """Build parser attempts in priority order."""
    if format_hint:
        return [format_hint]

    guessed = detect_format(path)
    if guessed == "fasta":
        return ["fasta", "genbank"]
    if guessed == "genbank":
        return ["genbank", "fasta"]
    return ["fasta", "genbank"]


def parse_file_records(path: Path, format_hint: Optional[str] = None) -> Tuple[str, List[SeqRecord]]:
    """
    Parse one file and return (format_used, records).

    Raises:
        ParseFileError: if no records can be parsed.
    """
    parse_errors: List[Tuple[str, Exception]] = []
    no_record_formats: List[str] = []

    for candidate_format in _candidate_formats(path, format_hint):
        records, error = _parse_records_with_format(path, candidate_format)
        if error is not None:
            parse_errors.append((candidate_format, error))
            continue
        if records:
            return candidate_format, records
        no_record_formats.append(candidate_format)

    if parse_errors:
        fmt, error = parse_errors[-1]
        raise ParseFileError(f"{path}: parse failed as {fmt} ({error})")

    tried = ", ".join(no_record_formats) if no_record_formats else "fasta, genbank"
    raise ParseFileError(f"{path}: no parseable records found (tried {tried})")


def infer_mol_type(record: SeqRecord) -> str:
    """Infer molecule type: DNA, RNA, protein, or unknown."""
    molecule_type = record.annotations.get("molecule_type")
    if molecule_type:
        mol_text = str(molecule_type).strip().lower()
        if "dna" in mol_text:
            return "DNA"
        if "rna" in mol_text:
            return "RNA"
        if "protein" in mol_text or "peptide" in mol_text:
            return "protein"

    letters = [base.upper() for base in str(record.seq) if base.isalpha()]
    if not letters:
        return "unknown"

    unique_letters = set(letters)

    if unique_letters.issubset(NUCLEOTIDE_STRICT):
        has_t = "T" in unique_letters
        has_u = "U" in unique_letters
        if has_t and has_u:
            return "unknown"
        if has_u:
            return "RNA"
        return "DNA"

    if unique_letters.issubset(NUCLEOTIDE_IUPAC):
        return "unknown"

    if unique_letters & PROTEIN_MARKERS:
        return "protein"

    non_iupac = [base for base in letters if base not in NUCLEOTIDE_IUPAC]
    if non_iupac:
        return "protein"

    return "unknown"


def calculate_gc_percent(sequence: str, mol_type: str) -> str:
    """Calculate GC percentage for DNA/RNA, else return empty string."""
    if mol_type not in {"DNA", "RNA"}:
        return ""

    normalized = sequence.upper().replace("U", "T")
    a_count = normalized.count("A")
    c_count = normalized.count("C")
    g_count = normalized.count("G")
    t_count = normalized.count("T")
    denominator = a_count + c_count + g_count + t_count

    if denominator == 0:
        return ""

    gc_percent = 100.0 * (g_count + c_count) / denominator
    return f"{gc_percent:.2f}"


def sequence_md5(sequence: str) -> str:
    """Return lower-case md5 hex digest for the raw sequence."""
    return hashlib.md5(sequence.encode("utf-8")).hexdigest()


def summarize_record(record: SeqRecord, source_file: Path, file_format: str) -> Dict[str, object]:
    """Convert one SeqRecord into a TSV-ready summary dict."""
    sequence = str(record.seq)
    mol_type = infer_mol_type(record)
    return {
        "id": record.id,
        "source_file": str(source_file),
        "format": file_format,
        "description": record.description.strip(),
        "length": len(sequence),
        "mol_type": mol_type,
        "gc_percent": calculate_gc_percent(sequence, mol_type),
        "md5": sequence_md5(sequence),
    }


def write_summary_tsv(rows: List[Dict[str, object]], output_path: Path) -> None:
    """Write summary rows to TSV with fixed header."""
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TSV_HEADER, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Read FASTA/GenBank files and write a summary TSV."
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="One or more input files or directories.",
    )
    parser.add_argument(
        "-f",
        "--format",
        choices=["fasta", "genbank"],
        help="Input format. If omitted, auto-detect per file and fall back to both parsers.",
    )
    parser.add_argument(
        "-o",
        "--out",
        default="summary.tsv",
        help="Output TSV path (default: ./summary.tsv).",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=0,
        help="Filter out records shorter than N bases/amino acids (default: 0).",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    args = parser.parse_args(argv)
    if args.min_length < 0:
        parser.error("--min-length must be >= 0")
    return args


def run(argv: Optional[Sequence[str]] = None) -> int:
    """Run the parser workflow; return process exit code."""
    args = parse_args(argv)
    files = list_input_files(args.inputs)

    rows: List[Dict[str, object]] = []
    for source_file in files:
        try:
            used_format, records = parse_file_records(source_file, args.format)
        except ParseFileError as exc:
            warn(str(exc))
            continue

        for record in records:
            summary = summarize_record(record, source_file, used_format)
            if int(summary["length"]) < args.min_length:
                continue
            rows.append(summary)

    rows.sort(key=lambda row: (str(row["source_file"]), str(row["id"])))
    write_summary_tsv(rows, Path(args.out))
    return 0


def main() -> None:
    """CLI entry point."""
    raise SystemExit(run())


if __name__ == "__main__":
    main()
