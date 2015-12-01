#! /bin/sh

ref=ecoli_K12_MG1655_U00096.3
src/overlapping.py run reads-simulated/${ref}-e0.00/reads-m4.fa reads-simulated/${ref}-e0.00/reads.m4 results/${ref}-e0.00/
src/overlapping.py run reads-simulated/${ref}-e0.05/reads-m4.fa reads-simulated/${ref}-e0.05/reads.m4 results/${ref}-e0.05/
src/overlapping.py run reads-simulated/${ref}-e0.10/reads-m4.fa reads-simulated/${ref}-e0.10/reads.m4 results/${ref}-e0.10/
src/overlapping.py run reads-simulated/${ref}-e0.15/reads-m4.fa reads-simulated/${ref}-e0.15/reads.m4 results/${ref}-e0.15/
src/overlapping.py run reads-simulated/${ref}-e0.20/reads-m4.fa reads-simulated/${ref}-e0.20/reads.m4 results/${ref}-e0.20/
cat results/${ref}-e0.00/summary.csv > results/summary-${ref}.csv
cat results/${ref}-e0.05/summary.csv > results/summary-${ref}.csv
cat results/${ref}-e0.10/summary.csv > results/summary-${ref}.csv
cat results/${ref}-e0.15/summary.csv > results/summary-${ref}.csv
cat results/${ref}-e0.20/summary.csv > results/summary-${ref}.csv

ref=saccharomyces_cerevisiae
src/overlapping.py run reads-simulated/${ref}-e0.00/reads-m4.fa reads-simulated/${ref}-e0.00/reads.m4 results/${ref}-e0.00/
src/overlapping.py run reads-simulated/${ref}-e0.05/reads-m4.fa reads-simulated/${ref}-e0.05/reads.m4 results/${ref}-e0.05/
src/overlapping.py run reads-simulated/${ref}-e0.10/reads-m4.fa reads-simulated/${ref}-e0.10/reads.m4 results/${ref}-e0.10/
src/overlapping.py run reads-simulated/${ref}-e0.15/reads-m4.fa reads-simulated/${ref}-e0.15/reads.m4 results/${ref}-e0.15/
src/overlapping.py run reads-simulated/${ref}-e0.20/reads-m4.fa reads-simulated/${ref}-e0.20/reads.m4 results/${ref}-e0.20/
cat results/${ref}-e0.00/summary.csv > results/summary-${ref}.csv
cat results/${ref}-e0.05/summary.csv > results/summary-${ref}.csv
cat results/${ref}-e0.10/summary.csv > results/summary-${ref}.csv
cat results/${ref}-e0.15/summary.csv > results/summary-${ref}.csv
cat results/${ref}-e0.20/summary.csv > results/summary-${ref}.csv
