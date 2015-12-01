#! /bin/sh

src/generate_data.py setup
src/overlapping.py setup

cd reference-genomes
tar -xzvf ecoli_K12_MG1655_U00096.3.fasta.tar.gz
tar -xzvf saccharomyces_cerevisiae.fa.tar.gz
cd ..

src/generate_data.py run reference-genomes/ecoli_K12_MG1655_U00096.3.fasta reads-simulated/ 30
src/generate_data.py run reference-genomes/saccharomyces_cerevisiae.fa reads-simulated/ 30

./run-overlaps.sh
