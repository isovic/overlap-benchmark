## overlap-benchmark  
Scripts to automate generating of simulated sequencing reads and running various tests on the data.  
Currently implemented tests aim at benchmarking overlappers.  

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

Or, you can run the ```example.sh``` script which begins by setting up the required tools, simulates the data and runs all overlaps.  

## Info about simulations
Simulated datasets were generated using PBSim. The error profile used was "55:17:28" which are the same as the 2d parameters determined using LAST from the E. Coli R7.3 data from 2014. (See [GraphMap paper](http://biorxiv.org/content/early/2015/06/10/020719)).
To be verbose, simulations were performed with parameters specified as:
```  
length_mean = 5600  
length_sd = 3500  
length_min = 50  
length_max = 100000  
accuracy_mean = (1.0 - error_rate)  
accuracy_sd = 0.09  
accuracy_min = 0.40  
difference_ratio = '55:17:28' (mismatch:insertion:deletion)  
```  
