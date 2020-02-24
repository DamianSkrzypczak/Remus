
# ![RemusLogo](remus/static/img/remus_logo_mini.png) Remus  


Public Remus instance is available at [http://remus.btm.umed.pl](http://remus.btm.umed.pl)

-----

Remus is a web tool helping in identification of regulatory regions relevant to a given monogenic disease phenotype.
Starting from a small set of genes implicated in the disease, Remus allows iterative building of a tissue-specific set of regions that likely play
a role in regulating expression of the input genes. 
After the list is finalized, it can be downloaded as a BED file with genomic coordinates or used directly to filter variants in a VCF file.

The growing inventory of regulatory data available in Remus at the moment includes coordinates of:

 - tissue-specific enhancers from [FANTOM5](http://fantom.gsc.riken.jp/5) and [ENCODE](https://www.encodeproject.org) repositories, 
 including the [SCREEN](https://screen.encodeproject.org/) datasets
 - promoters ([SCREEN](https://screen.encodeproject.org/)) and transcription start sites ([SlideBase](http://slidebase.binf.ku.dk) / [FANTOM5](http://fantom.gsc.riken.jp/5)) 
 - regions of accessible chromatin ([ENCODE](https://www.encodeproject.org) and [SCREEN](https://screen.encodeproject.org/))
 - microRNA - mRNA interactions from [miRTarBase](http://mirtarbase.mbc.nctu.edu.tw/) and [miRWalk](http://mirwalk.umm.uni-heidelberg.de)
  
In upcoming releases Remus will also enable inclusion of other regulatory features.


### Usage

1. Select genome build (hg19, GRCh38). 
   Regulatory features available in the primary sources only in one genome build (e.g. hg19 for FANTOM5 data) have been liftedOver to the other genome build.
   
2. Type-in symbols of genes relevant for the phenotype.

3. Choose tissues and/or cell-types relevant for the phenotype. 
   In parenthesis next to name of the tissue/cell-type, symbols for available datatypes are shown, e.g. 
   ENH_E - enhancers from ENCODE, PR_F5 - FANTOM5 promoters, CHR_S - accessible chromatin from SCREEN.

4. Select types of regulatory features to include. 
   Set maximal distance, up- and downstream from the transcript start (RefSeq transcript coordinates are used).
   Choose if regulatory features present in any of selected tissues (permissive) or in all of them (strict) should be used.

   Please note that miRNA-gene interactions are filtered against accessible chromatin regions in selected tissues, i.e.
   only miRNAs encoded in accessible parts of the genome will be included. 

5. Download resulting BED file or filter your variants.
   Note that the VCF file is filtered in your browser - it is NOT sent or uploaded anywhere.

##### VCF filtering

Remus allows for in-browser filtering of a VCF file using the output BED file with regulatory regions.
Variants falling into the regions are selected and returned in a plain text VCF file.
Input must be provided as sorted plain-text VCF, and filtering large files takes only a few seconds (~5s on 500M VCF).
Filtering BGZipped & Tabix'ed files was considerably slower in tests, and although implemented, has been disabled for the time being.

In-browser filtering means that the variant file does not leave your computer - great feature if you are working with sensitive data.
The downside is that the functionality can be somewhat limited.
Currently the VCF file is read in one piece ([to be changed](https://github.com/seru71/Remus/issues/15)), and empirically tested size limitations were following:

 - plain text VCF in Chrome 68 can be upto 1GB,
 - plain text VCF limit in Firefox 64.0 was ~250MB


### Installation of a local instance of Remus

#### Building and running a Docker container

1. In the Remus repo (REMUS_DIR), build the docker image:

    `docker build -t remus .`
    
2. To prepare data for Remus, either: 
    
    (__shorter version__)
    
    Download archive with Remus data files from [here](http://remus.btm.umed.pl/data_download/remus_0.5_data.tar), and extract the archive in REMUS_DIR.

    or (__longer version__):

    Start Remus container interactively with write access to the repo directory (REMUS_DIR):
      
    `docker run --rm --name remus_databuild -v REMUS_DIR:/var/www/remus:rw -ti remus`

    After reading `exterenal_resources/README.md`, download liftOver and chains by:
      
    `cd external_resources && ./download.sh && cd ..`

    Next, launch
      
    `./make_data_tree.sh`
      
    This will download necessary files and fill REMUS_DIR/data with all Remus data.
    Now you can exit the container.

3. Start docker container with the app available at `http://localhost:LOCAL_PORT` 

    ```docker run --rm -d --name remus_app -v REMUS_DIR:/var/www/remus -p LOCAL_PORT:80 remus apachectl -D FOREGROUND```
        

#### Native or development install

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
    
The application is available at `127.0.0.1:5000`




### Credits

Remus has been developed at [BTM](https://biostat.umed.pl), Medical Univeristy of Lodz, Poland. 
Application's UI and initial work on its internals was done by [Damian Skrzypczak](https://github.com/DamianSkrzypczak) as part of his MSc project. 
Since then, it has been extended by [me](https://github.com/seru71).

Code for in-browser filtering of tabixed VCF files was adopted from [js-local-vcf](https://github.com/jsa-aerial/js-local-vcf) written by [Jon Anthony](https://github.com/jsa-aerial).

Liftover of genome coordinates is done using [liftOver](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/) tool developed by [Jim Kent](http://www.kentinformatics.com/about-us.html). 
Note that liftOver is [free only for academic use](external_resources/REAMDE.md).
 
Data used in Remus is downloaded from public databases and primary sources attributed in the description on top of the page.

This project is funded by NCN Polonez grant no 2016/23/P/NZ2/04251. This project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No 665778.

![eu](https://seru71.github.io/polonez-project/img/eu_logo.jpg)
![ncn](https://seru71.github.io/polonez-project/img/ncn_logo.png)

