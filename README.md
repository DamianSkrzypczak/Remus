# Remus  
Remus is an online tool which helps with identification of 
regulatory regions relevant to a given monogenic disease phenotype.  


Starting from a small set of genes implicated in the disease,
Remus allows iterative building of a tissue-specific set of regions that likely play
a role in regulating expression of the input genes. 
 
The growing inventory of data available in Remus at 
the moment includes coordinates of tissue-specific enhancers, 
transcription start sites, and regions of accessible chromatin 
from large scale public datasets (ENCODE, FANTOM5, GTEx). 
 
 
In upcoming releases Remus will also enable inclusion of relevant binding sites 
of transcription factors and miRNAs, as well as locations of 
TAD boundaries adjacent to the input genes.

Remus is in active development and will be made freely available in fall 2018.

### Installation
##### Dependencies:
    pip install -r requirements.txt

or if development mode:  
    
    pip install -r requirements-dev.txt
##### Data tree preparation:
In application root run:  
(ensure that at least "pandas" library is installed)  
    
    bash make_data_tree.sh
