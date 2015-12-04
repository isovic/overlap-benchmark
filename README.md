### overlap-benchmark  
Scripts to automate generating of simulated sequencing reads and running various tests on the data.  
Currently implemented tests aim at benchmarking overlappers.  

This repo is in early stages of development.  

### Simulations and evaluation
Simulated datasets are generated using PBSim. The error profile used is "55:17:28" which is the same as the 2d parameters determined using LAST from the E. Coli R7.3 data from 2014. (See [GraphMap paper](http://biorxiv.org/content/early/2015/06/10/020719)).  
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
The alignments of simulated reads are converted from MAF to SAM format using LAST's ```maf-convert``` script.  

To evaluate overlaps, MHAP's ```edu.umd.marbl.mhap.main.EstimateROC``` module is used.  
Since ```EstimateROC``` requires either BLASR's M4 file (or an MHAP file) containing truth overlaps, we use BLASR's tool ```samtom4``` to generate the appropriate truth file. The truth overlaps are located in files named ```reads.m4```.  
Also, since ```EstimateROC``` requires an the input FASTA reads file to be provided, headers of the simulated reads are replaced with their ordinal number of appearance (so that they can be matched to the read ID's in the MHAP files). Final reads are given in files named ```reads-m4.fa```.  

Additionally, to convert Minimap's overlaps to MHAP format, a script ```paf2mhap.pl``` is used ([https://github.com/lh3/miniasm/blob/master/misc/paf2mhap.pl](https://github.com/lh3/miniasm/blob/master/misc/paf2mhap.pl)).

At this moment, only GraphMap, Minimap and MHAP are tested automatically.  
Thanks to ```da2paf.pl``` script provided on the Miniasm repo ([https://github.com/lh3/miniasm/blob/master/misc/da2paf.pl](https://github.com/lh3/miniasm/blob/master/misc/da2paf.pl)), DALIGNER should soon be added to the benchmark as well.  


### Installation
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

### Usage  
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

### Results  
Running the ```src/overlapping.py``` script will generate a CSV file in the specified output path.  
For example, should the script be run with:  
```  
src/overlapping.py run reads-simulated/ecoli_K12_MG1655_U00096.3-e0.10/reads-m4.fa reads-simulated/ecoli_K12_MG1655_U00096.3-e0.10/reads.m4 results/ecoli_K12_MG1655_U00096.3-e0.10/  
```  
A summary output file will be placed in:  
```
results/ecoli_K12_MG1655_U00096.3-e0.10/summary.csv  
```  

Description of the output lines:  
The first line of the CSV file contains the specified output path where the results can be found.  
The next line describes the columns of the results (a header line).  
The following lines contain results for each tested overlapper.  
