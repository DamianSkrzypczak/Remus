
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

### Credits

Remus is developed at [BTM](https://biostat.umed.pl), Medical Univeristy of Lodz, Poland. Application's UI and initial work on its internals has been done by [Damian Skrzypczak](https://github.com/DamianSkrzypczak).

This project is funded by NCN Polonez grant no 2016/23/P/NZ2/04251. This project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No 665778.

![eu](https://seru71.github.io/polonez-project/img/eu_logo.jpg)
![ncn](https://seru71.github.io/polonez-project/img/ncn_logo.png)

    
    
