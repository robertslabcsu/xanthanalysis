#!/usr/bin/env bash

for line in Xtgenomes_proteins/*.faa
do
    echo "$line" > tmp.txt
    tmpgenome=$(cut -d "/" -f 2 tmp.txt)
    makeblastdb -in $line -dbtype prot >> blastdb.log
    for effector in allXanthomonasEffectors/*
    do
        echo "$effector" > tmpname
        effect_name=$(cut -d "/" -f 2 tmpname)
        blastp -query "$effector" -db $line -outfmt "6 qseqid sseqid pident" > blast_tmp_result
        cut -f 3 blast_tmp_result | head -n1 > blast_tmp_ident
        tmp_pident=$(cat blast_tmp_ident)
        tmp_pident2=${tmp_pident%.*}
        if [[ $tmp_pident2 -ge 75 ]]
        then
            echo "$effect_name "$tmp_pident2"" >> blast_results/$tmpgenome.results
            rm blast_tmp_ident
            rm tmpname
        else
            echo "$effect_name FALSE" >> blast_results/$tmpgenome.results
            rm blast_tmp_ident
            rm tmpname
        fi
        rm blast_tmp_result
    done
    rm tmp.txt
done
