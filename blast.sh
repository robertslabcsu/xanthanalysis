#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
xanthanalysis blast
Usage:
  blast.sh --query Q.fa --db genome.fa --out out.tsv [--threads 8] [--task blastn]
Requires: blastn, makeblastdb
EOF
}

# defaults
THREADS=8
TASK="blastn"
OUT="blast.tsv"

# parse args
QUERY=""; DB=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --query) QUERY="$2"; shift 2 ;;
    --db)    DB="$2"; shift 2 ;;
    --out)   OUT="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --task) TASK="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1"; usage; exit 1 ;;
  esac
done

# checks
command -v makeblastdb >/dev/null || { echo "makeblastdb not found"; exit 127; }
command -v blastn >/dev/null || { echo "blastn not found"; exit 127; }
[[ -f "$QUERY" ]] || { echo "Missing --query file"; exit 2; }
[[ -f "$DB" ]] || { echo "Missing --db file"; exit 2; }

# build db (cached)
DBDIR="$(dirname "$DB")"
DBNAME="$(basename "$DB")"
if ! ls "${DB}".nin >/dev/null 2>&1; then
  makeblastdb -in "$DB" -dbtype nucl -parse_seqids -out "$DB"
fi

# run blast
blastn -query "$QUERY" -db "$DB" -task "$TASK" -num_threads "$THREADS" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  > "$OUT"

echo "Done. Results: $OUT"
