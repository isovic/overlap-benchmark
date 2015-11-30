## simdata
Scripts to automate generating of simulated sequencing reads and running various tests on the data.  
Examples of tests include testing:  
- overlappers  
- aligners  
- assemblers  

This repo is in early stages of development.  

## Usage  
Setup the tools required to simulate the data:  
```  
src/generate_data.py setup  
```  

Run the simulations for a given reference genome:  
```  
src/generate_data.py run reference.fa reads-simulated/ 30  
```  
This will simulate a ~30x coverage of the reference.fa FASTA file, and output the results to reads-simulated/ folder.  

