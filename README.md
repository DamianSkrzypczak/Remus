
# ![RemusLogo](remus/static/img/remus_logo_mini.png) Remus  

Remus is a web tool helping in identification of 
regulatory regions relevant to a given monogenic disease phenotype.  

Starting from a small set of genes implicated in the disease,
Remus allows iterative building of a tissue-specific set of regions that likely play
a role in regulating expression of the input genes. 
 
The growing inventory of data available in Remus at 
the moment includes coordinates of:

 - tissue-specific enhancers from [ENCODE](https://www.encodeproject.org) and [FANTOM5](http://fantom.gsc.riken.jp/5) repositories
 - transcription start sites ([SlideBase](http://slidebase.binf.ku.dk) / FANTOM5) 
 - regions of accessible chromatin (ENCODE)
 - microRNA - mRNA interactions from [miRTarBase](http://mirtarbase.mbc.nctu.edu.tw/) and [miRWalk](http://mirwalk.umm.uni-heidelberg.de)
  
In upcoming releases Remus will also enable inclusion of other regulatory features.

Remus is developed at [BTM](https://biostat.umed.pl), Medical Univeristy of Lodz, Poland. 

Code for in-browser filtering of tabixed VCF files was adopted from [js-local-vcf](https://github.com/jsa-aerial/js-local-vcf) written by [Jon Anthony](https://github.com/jsa-aerial).
### Usage

coming soon

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
    
    
