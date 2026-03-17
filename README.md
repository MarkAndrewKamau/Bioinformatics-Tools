# parse_fasta mini-project

## Usage

`parse_fasta.py` reads one or more FASTA/GenBank files (or directories containing them) and writes a summary TSV with per-record metadata:

- record id and description
- source file and parsed format
- length, inferred molecule type, GC%
- MD5 digest of the raw sequence

### Install dependencies

```bash
pip install -r requirements.txt
```

Or directly:

```bash
pip install biopython pytest
```

### Run the script

Parse a single file:

```bash
python parse_fasta.py input.fa
```

Parse multiple files from a directory:

```bash
python parse_fasta.py dir_with_fastas -o results.tsv
```

Specify the format:

```bash
python parse_fasta.py sample.gb -f genbank
```

Use minimum length filtering:

```bash
python parse_fasta.py *.fa --min-length 200
```

### Output format

Header:

```tsv
id	source_file	format	description	length	mol_type	gc_percent	md5
```

Example row (from `example_inputs/example.fa`):

```tsv
seq_dna	example_inputs/example.fa	fasta	seq_dna Example DNA record	10	DNA	50.00	45aff2fecf7615d56bc0567dffab9fa8
```

### Run tests

```bash
pytest -q
```

### Exit codes and error handling

- Exit code is `0` on successful completion.
- Parse failures, empty/unsupported files, and missing paths are reported as warnings to `stderr`.
- The script continues processing other files when one file fails.
