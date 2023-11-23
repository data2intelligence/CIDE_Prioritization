Example code of prioritize genes using risk scores computed from cancer immunotherapy datasets.    

# Pre-requisites:  
Python >= 3.9 (you may install [Anaconda](https://www.anaconda.com/download))  

# Usage:  
CD into folder "src", and type "./rank.py geneset_input". For example, run "./rank.py ../data/secreted_proteins.txt".  


#output:  
1, \*.rank.xlsx: Top genes passing the threshold  
2, \*.rank.stat.xlsx: Statistical test results for all genes in the input set.  
