#! /bin/sh

src/overlapping.py run reads-simulated/ecoli_K12_MG1655_U00096.3-e0.00/reads-m4.fa reads-simulated/ecoli_K12_MG1655_U00096.3-e0.00/reads.m4 results/ecoli_K12_MG1655_U00096.3-e0.00/
src/overlapping.py run reads-simulated/ecoli_K12_MG1655_U00096.3-e0.05/reads-m4.fa reads-simulated/ecoli_K12_MG1655_U00096.3-e0.05/reads.m4 results/ecoli_K12_MG1655_U00096.3-e0.05/
src/overlapping.py run reads-simulated/ecoli_K12_MG1655_U00096.3-e0.10/reads-m4.fa reads-simulated/ecoli_K12_MG1655_U00096.3-e0.10/reads.m4 results/ecoli_K12_MG1655_U00096.3-e0.10/
src/overlapping.py run reads-simulated/ecoli_K12_MG1655_U00096.3-e0.15/reads-m4.fa reads-simulated/ecoli_K12_MG1655_U00096.3-e0.15/reads.m4 results/ecoli_K12_MG1655_U00096.3-e0.15/
src/overlapping.py run reads-simulated/ecoli_K12_MG1655_U00096.3-e0.20/reads-m4.fa reads-simulated/ecoli_K12_MG1655_U00096.3-e0.20/reads.m4 results/ecoli_K12_MG1655_U00096.3-e0.20/

cat results/ecoli_K12_MG1655_U00096.3-e0.00/summary.csv > results/summary.csv
cat results/ecoli_K12_MG1655_U00096.3-e0.05/summary.csv >> results/summary.csv
cat results/ecoli_K12_MG1655_U00096.3-e0.10/summary.csv >> results/summary.csv
cat results/ecoli_K12_MG1655_U00096.3-e0.15/summary.csv >> results/summary.csv
cat results/ecoli_K12_MG1655_U00096.3-e0.20/summary.csv >> results/summary.csv
