#!/usr/bin/env bash

# Set the minimum identity threshold for BLAST hits, recommended 75 for bacterial effectors
min_identity=75

# Create the blastdb directory if it doesn't already exist
mkdir -p blastdb

# Create the blast_results directory if it doesn't already exist
mkdir -p blast_results

# Process each protein FASTA file in the Xtgenomes_proteins directory
for genome_fasta in Xtgenomes_proteins/*.faa; do
    # Extract the genome name from the FASTA filename
    genome=$(basename $genome_fasta | cut -d "." -f1)

    # Create a BLAST database for the current genome
    makeblastdb -in $genome_fasta -dbtype prot -out blastdb/$genome.db >> blastdb.log

    # Process each effector protein FASTA file in the allXanthomonasEffectors directory
    for effector_fasta in allXanthomonasEffectors/*; do
        # Extract the effector name from the FASTA filename
        effector_name=$(basename $effector_fasta | cut -d "." -f1)

        # Perform a BLASTP search of the effector against the genome database
        blastp -query $effector_fasta -db blastdb/$genome.db -outfmt "6 qseqid sseqid pident" > blast_tmp_result

        # Extract the identity percentage from the BLAST output
        identity=$(cut -f 3 blast_tmp_result | head -n1)

        # If the identity is greater than or equal to the minimum threshold, record the hit
        if [ $identity -ge $min_identity ]; then
            echo "$effector_name $identity" >> blast_results/$genome.results
        fi

        # Remove temporary files
        rm blast_tmp_result
    done
done
