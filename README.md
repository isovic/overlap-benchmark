## simdata
Scripts to automate generating of simulated sequencing reads and running various tests on the data.  
Currently implemented tests only aim at benchmarking overlappers.  

This repo is in early stages of development.  

## Installation
Setup the tools required to simulate the data:  
```  
src/generate_data.py setup  
```  

Setup the tools that will be tested for overlapping:  
```  
src/overlapping.py setup  
```  

Unpack the sample reference provided in this repo:  
```  
cd reference-genomes  
tar -xzvf ecoli_K12_MG1655_U00096.3.fasta.tar.gz  
cd ..  
```  

## Usage  
Run the simulations for a given reference genome:  
```  
src/generate_data.py run reference-genomes/ecoli_K12_MG1655_U00096.3.fasta reads-simulated/ 30  
```  
This will simulate a ~30x coverage of the reference FASTA file, and output the results to reads-simulated/ folder.  

Run the overlapping process on a single dataset (this will run all overlappers):  
```  
src/overlapping.py run reads-simulated/ecoli_K12_MG1655_U00096.3-e0.10/reads-m4.fa reads-simulated/ecoli_K12_MG1655_U00096.3-e0.10/reads.m4 results/ecoli_K12_MG1655_U00096.3-e0.10/  
```  

Run a batch of overlapping tests pre-scripted in a Bash script:  
```  
./run-overlaps.sh  
```  
