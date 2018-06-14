#
# Dockerfile for building docker container with the app
#

LABEL authors="Paweł Sztromwasser, Damian Skrzypczak"

FROM rastasheep/ubuntu-sshd

RUN apt update \
 && apt install -y \
    apache2 \
    libapache2-mod-wsgi \
	build-essential \
	python3-dev \
	python3-pip \
	zlib1g-dev \
	libbz2-dev \
	libssl-dev \
	liblzma-dev \
	samtools \
	vim \
	git \
 && pip3 install --timeout 36000 \
	pandas \
	pysam \
	pybedtools \
	flask \
 && wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz \
 && tar -xzf bedtools-2.27.1.tar.gz \
 && cd bedtools2 \
 && make \
 && ln -s /bedtools2/bin/bedtools /usr/local/bin
 && apt-get clean \
 && apt-get autoremove \
 && rm -rf /var/lib/apt/lists/*


COPY ./requirements.txt /var/www/remus/requirements.txt
RUN pip install -r /var/www/remus/requirements.txt

COPY apache-remus.conf /etc/apache2/sites-available/apache-remus.conf
RUN a2ensite remus
RUN a2enmod headers

COPY apache-remus.wsgi /var/www/remus/apache-remus.wsgi

COPY ./run.py /var/www/remus/app.py
COPY ./remus /var/www/remus/

RUN a2dissite 000-default.conf
RUN a2ensite apache-remus.conf

EXPOSE 80

WORKDIR /var/www/apache-remus

CMD  /usr/sbin/apache2ctl -D FOREGROUND