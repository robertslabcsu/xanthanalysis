import typer
import subprocess
from pathlib import Path

app = typer.Typer(help="Xanthomonas analysis CLI")

@app.command()
def setup_db(db: Path):
    typer.echo(f"Building BLAST DB for {db}")
    subprocess.run(["makeblastdb", "-in", str(db), "-dbtype", "nucl", "-parse_seqids", "-out", str(db)], check=True)

@app.command()
def blast(query: Path, db: Path, out: Path = Path("blast.tsv"), threads: int = 8, task: str = "blastn"):
    subprocess.run([
        "blastn","-query",str(query),"-db",str(db),
        "-task",task,"-num_threads",str(threads),
        "-outfmt","6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-out",str(out)
    ], check=True)
    typer.echo(f"Wrote {out}")

if __name__ == "__main__":
    app()
