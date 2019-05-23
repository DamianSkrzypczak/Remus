
# ![RemusLogo](remus/static/img/remus_logo_mini.png) Remus  

Remus is a web tool helping in identification of regulatory regions relevant to a given monogenic disease phenotype.  
Starting from a small set of genes implicated in the disease, Remus allows iterative building of a tissue-specific set of regions that likely play
a role in regulating expression of the input genes. 
After the list is finalized, it can be downloaded as a BED file with genomic coordinates or used directly to filter variants in a VCF file.

The growing inventory of regulatory data available in Remus at the moment includes coordinates of:

 - tissue-specific enhancers from [ENCODE](https://www.encodeproject.org) and [FANTOM5](http://fantom.gsc.riken.jp/5) repositories
 - transcription start sites ([SlideBase](http://slidebase.binf.ku.dk) / [FANTOM5](http://fantom.gsc.riken.jp/5)) 
 - regions of accessible chromatin ([ENCODE](https://www.encodeproject.org))
 - microRNA - mRNA interactions from [miRTarBase](http://mirtarbase.mbc.nctu.edu.tw/) and [miRWalk](http://mirwalk.umm.uni-heidelberg.de)
  
In upcoming releases Remus will also enable inclusion of other regulatory features.

Remus is developed at [BTM](https://biostat.umed.pl), Medical Univeristy of Lodz, Poland. 

### Usage

1. Select genome build (hg19, GRCh38). 
   Regulatory features available in primary sources only in one genome build (e.g. hg19 for FANTOM5 data) have been liftedOver to the other genome build.
   
2. Type-in symbols of genes relevant for the phenotype.

3. Choose tissues and/or cell-types relevant for the phenotype. 
   In parenthesis next to name of the tissue/cell-type, symbols for available datatypes are shown, e.g. 
   ENH_ENC - ENCODE enhancers, TSS_F5 - transcription-start sites from FANTOM5, CHRM - accessible chromatin.

4. Select types of regulatory features to include. 
   Set maximal distance, up- and downstream from transcript start (UCSC transcript coordinates are used).
   Choose if regulatory features present in any of selected tissues (permissive) or in all of them (strict) should be used.

   Please note that miRNA-gene interactions are filtered against accessible chromatin regions in selected tissues, i.e.
   only miRNAs encoded in accessible parts of the genome will be included. 

5. Download resulting BED file or filter your variants.
   Note that the VCF file is filtered in your browser - it is NOT sent or uploaded anywhere.

##### VCF filtering

Remus allows for in-browser filtering of a VCF file using the regions-of-interest BED file.
Variants falling into the regions-of-interest are selected and returned in plain text VCF file.
Input must be provided as sorted plain-text VCF, and filtering large files takes only a few seconds (~5s on 500M VCF).
Filtering BGZipped & Tabix'ed files was considerably slower in tests and although implemented, has been disabled for the time being.

In-browser filtering means that the variant file does not leave your computer, which is great for sensitive data.
The downside is that the functionality can be somewhat limited.
Currently the VCF file is read in one piece ([to be changed](https://github.com/seru71/Remus/issues/15)), and empirically tested size limitations were following:

 - plain text VCF in Chrome 68 can be upto 1GB,
 - plain text VCF limit in Firefox 64.0 was ~250MB


### Installation of local instance of Remus

##### Installing dependencies:

    pip install -r requirements.txt

or if development mode:  
    
    pip install -r requirements-dev.txt

##### Data tree preparation:

After reading `exterenal_resources/README.md`, download liftOver and chain files.
In application root run:

    bash external_resources/download.sh
      
Next download the data. This step can take long time because of large amount of data needed to be downloaded.

    bash make_data_tree.sh

##### Running

In application root run:

    python3 app.py
    
The website is available at `127.0.0.1:5000`



### Setup in Docker container

##### Building and running container

    docker build -t remus .
    docker run -p LOCAL_PORT_1:80 -p LOCAL_PORT_2:22 --name remus remus

##### Configuring container

    ssh root@localhost -p LOCAL_PORT_2
    # password=root
    
Change the password, or better, disable root password login and use key-based login

##### Data preparation
    
    cd /var/www/remus
    
After reading `exterenal_resources/README.md`, download liftOver and chains, run:

    bash external_resources/download.sh

    # This step can take long time because of large amount of data needed to be downloaded.    
    bash make_data_tree.sh
    
##### Starting Apache
    
    /usr/sbin/apache2ctl -D BACKGROUND

##### Accessing Remus website

If whole process of building and configuring the container succeded, exit the container. Remus website should be accessible under `127.0.0.1:LOCAL_PORT_1`
    
    
### Acknowledgements

Codebase for Remus was written by [Damian Skrzypczak](https://github.com/DamianSkrzypczak) as part of his MSc project. 
Since then, it has been extended by me.

Code for in-browser filtering of tabixed VCF files was adopted from [js-local-vcf](https://github.com/jsa-aerial/js-local-vcf) written by [Jon Anthony](https://github.com/jsa-aerial).

Liftover of genome coordinates is done using [liftOver](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/) tool developed by [Jim Kent](http://www.kentinformatics.com/about-us.html). 
Note that liftOver is [free only for academic use](external_resources/REAMDE.md).
 
Data used in Remus is downloaded from public databases and primary sources attributed in the description on top of the page.
