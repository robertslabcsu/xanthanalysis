# xanthanalysis

xanthanalysis

Tools and reproducible steps to compare Xanthomonas translucens effector repertoires (with emphasis on TAL effectors) across genomes and isolates—exactly as used in the Microbial Genomics study on recent Colorado isolates CO236 (Xtt) and CO237 (Xtu). 
PMC

<p align="center"> <em>“Custom bash and Python scripts for the effector analysis are available at this repository.”</em> :contentReference[oaicite:1]{index=1} </p>
Table of contents

At a glance

Install

Quick start

What it does (methods overview)

Inputs & outputs

Reproduce the paper analysis

HPC notes (SLURM)

Project structure

Cite us

Data availability

Contributing

License

At a glance

🔎 Effector comparison: screens genomes against curated effector references to summarize presence/absence and similarity (TALs + other T3Es).

🧬 Paper-aligned: reproduces analyses supporting the finding that T3Es are largely conserved, while TAL effector presence/sequence varies among isolates and correlates with virulence patterns. 
PMC

🧰 Minimal deps: BLAST+ plus standard UNIX/Python tools.

🧪 Tiny example: included so you can test the pipeline without downloading large genomes (see Quick start).

Install
Option A — Conda (recommended)
# 1) Create environment
mamba env create -f environment.yml   # or: conda env create -f environment.yml
mamba activate xanthanalysis          # or: conda activate xanthanalysis

# 2) Verify BLAST+ is available
blastn -version
makeblastdb -version

Option B — System packages

Install NCBI BLAST+ via your package manager and ensure blastn/makeblastdb are on PATH.

Quick start

This toy run uses small example FASTA files to show the workflow end-to-end.

# 0) Activate env (see Install)
mamba activate xanthanalysis

# 1) Build a BLAST database from a (small) genome FASTA
makeblastdb -in examples/genome.fa -dbtype nucl -parse_seqids -out examples/genome.fa

# 2) Run BLAST of effector queries against the genome DB (tabular fmt 6)
blastn -query examples/effectors.fa -db examples/genome.fa -num_threads 4 \
  -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" \
  > results/example_blast.tsv

# 3) Summarize hits
python scripts/summarize_hits.py results/example_blast.tsv \
  --identity 70 --length 200 --out results/example_summary.tsv

# 4) Inspect results
column -t -s$'\t' results/example_summary.tsv | less -S


If you’re working with real genomes, see Reproduce the paper analysis below for data download pointers and full steps.

What it does (methods overview)

The repository implements the core steps used in the paper’s effector analysis:

Assemble effector query sets (TALs, T3Es) from prior literature/curated sources.

Index each target genome with makeblastdb.

Search effectors against genomes with blastn (or tblastn for protein queries), capturing hits in tabular format (fmt 6).

Filter and summarize by identity/length coverage to call presence/absence and nearest-match identity.

Export tables suitable for downstream visualization and comparison.

These steps mirror the study’s approach: previously characterized Xanthomonas effectors were compared against putative effectors in recent vs. older isolates to evaluate conservation and divergence, with TAL effectors showing notable differences while most T3Es were conserved. 
PMC

Inputs & outputs
Inputs

Genome FASTA files (assembled chromosomes/contigs) for each isolate of interest.

Effector FASTA queries (DNA or protein) for TALs and/or other T3Es.

Optional: isolate metadata (CSV) for labeling.

Outputs

*.blast.tsv — raw BLAST tabular hits (fmt 6).

*_summary.tsv — one row per effector × genome with best hit, identity, length, and presence/absence call.

presence_absence_matrix.tsv — wide table for plotting or clustering.

Column definitions are included at the top of each TSV.

Reproduce the paper analysis

The study generated high-quality long-read assemblies for Xtt CO236 and Xtu CO237, then compared effector repertoires against published isolates. As reported:

CO236 (barley) belongs to pathovar translucens (Xtt), and CO237 (wheat) belongs to pathovar undulosa (Xtu). 
PMC

Genomes are available under BioProject accessions PRJNA1017868 (CO236) and PRJNA1017870 (CO237). 
PMC

Raw ONT reads are deposited in Dryad (doi:10.5061/dryad.d51c5b06q). 
PMC

1) Get genomes

Download assembled genomes (FASTA) for CO236/CO237 and any comparator isolates from GenBank using the BioProject accessions above. 
PMC

2) Build databases
makeblastdb -in data/genomes/Xtu_CO237.fa -dbtype nucl -parse_seqids -out data/genomes/Xtu_CO237.fa
makeblastdb -in data/genomes/Xtt_CO236.fa -dbtype nucl -parse_seqids -out data/genomes/Xtt_CO236.fa

3) Prepare effector queries

Place curated TAL and T3E sequences in data/effectors/ (FASTA). The paper compared previously characterized effectors to putative effectors identified in the new isolates. 
PMC

4) Run searches
# TALs vs CO237 (Xtu)
blastn -query data/effectors/tal_effectors.fa -db data/genomes/Xtu_CO237.fa \
  -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" \
  -num_threads 8 > results/xu237_tal.tsv

# TALs vs CO236 (Xtt)
blastn -query data/effectors/tal_effectors.fa -db data/genomes/Xtt_CO236.fa \
  -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" \
  -num_threads 8 > results/xt236_tal.tsv

5) Summarize & compare
python scripts/summarize_hits.py results/xu237_tal.tsv --identity 70 --length 200 \
  --out results/xu237_tal_summary.tsv

python scripts/summarize_hits.py results/xt236_tal.tsv --identity 70 --length 200 \
  --out results/xt236_tal_summary.tsv

python scripts/make_presence_matrix.py results/*_tal_summary.tsv \
  --out results/tal_presence_absence.tsv


The paper reports that T3E repertoires were highly conserved across Xtu isolates, whereas TAL effectors differed in presence/sequence in ways consistent with virulence differences. Use the same steps for T3Es to verify conservation across your panel. 
PMC

HPC notes (SLURM)

Example SLURM array to BLAST many genomes:

#!/usr/bin/env bash
#SBATCH -J xa-blast
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH -A your_account
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
#SBATCH --array=0-$(($(wc -l < genome_list.txt)-1))

GENOME=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" genome_list.txt)
makeblastdb -in "$GENOME" -dbtype nucl -parse_seqids -out "$GENOME"
blastn -query data/effectors/tal_effectors.fa -db "$GENOME" -num_threads $SLURM_CPUS_PER_TASK \
  -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" \
  > "results/$(basename "${GENOME%.*}")_tal.tsv"

Project structure
xanthanalysis/
├─ README.md
├─ environment.yml
├─ scripts/
│  ├─ summarize_hits.py
│  └─ make_presence_matrix.py
├─ examples/
│  ├─ genome.fa
│  └─ effectors.fa
├─ data/               # (you add: genomes/ and effectors/)
└─ results/            # (created at runtime; gitignored)


scripts/summarize_hits.py filters BLAST fmt-6 to best hits and presence/absence.

scripts/make_presence_matrix.py pivots summaries to a wide matrix.

If you prefer a single command-line app (e.g., xanthanalysis blast ...), see the issue tracker—there’s a plan to expose these steps via a Typer-based CLI and an optional Snakemake workflow.

Cite us

If you use this code or reproduce the analysis, please cite:

Gutierrez-Castillo, D. E., Barrett, E., & Roberts, R. (2024).
A recently collected Xanthomonas translucens isolate encodes TAL effectors distinct from older, less virulent isolates. Microbial Genomics, 10(1), 001177.
DOI: 10.1099/mgen.0.001177. PMCID: PMC10868612. 
PubMed

Data availability

CO236 (Xtt) and CO237 (Xtu) genome accessions: PRJNA1017868 and PRJNA1017870. 
PMC

ONT reads: Dryad doi:10.5061/dryad.d51c5b06q. 
PMC

See the paper’s Data Summary section for a full list of isolates and software environment details used in the study. 
PMC

Contributing

Bug reports and PRs welcome!
Please open an issue describing:

Your dataset and goal (e.g., TAL vs T3E comparison),

Exact command(s) run,

Environment details (blastn -version, OS, conda env export).

Planned enhancements:

Packaged CLI (xanthanalysis command),

Snakemake workflow for multi-genome panels,

Plotting utilities for presence/absence heatmaps.

License

This repository is intended for open scientific reuse. See LICENSE for terms.
If no license is present yet, choose and add one via GitHub’s Add file → Add license (MIT/GPL-3.0/Apache-2.0 are common in bioinformatics). 
GitHub Docs
Gist

Acknowledgements

This repository accompanies the Microbial Genomics article and makes the scripts used therein publicly available for community reuse and extension. 
PMC
