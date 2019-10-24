#
# Dockerfile for building docker container with the app
#

FROM ubuntu:18.04

LABEL authors="Paweł Sztromwasser, Damian Skrzypczak"

RUN apt update \
 && apt install -y \
    apache2 \
    libapache2-mod-wsgi-py3 \
	build-essential \
	python3-dev \
	python3-pip \
	zlib1g-dev \
	libbz2-dev \
	libssl-dev \
	liblzma-dev \
    p7zip-full \
	samtools \
    tabix \
	vim \
	git \
    locales \
 && locale-gen en_US.UTF-8 \
 && wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz \
 && tar -xzf bedtools-2.27.1.tar.gz \
 && cd bedtools2 \
 && make \
 && ln -s /bedtools2/bin/bedtools /usr/local/bin \
 && apt-get clean \
 && apt-get autoremove \
 && rm -rf /var/lib/apt/lists/*

ENV LANG en_US.UTF-8

WORKDIR /var/www/remus
COPY data_sources ./data_sources
COPY remus ./remus
COPY tests ./tests
COPY a* LICENSE make_data_tree.sh README.md requirements* ./


RUN pip3 install -r requirements.txt

COPY apache-remus.conf /etc/apache2/sites-available/apache-remus.conf
RUN a2enmod headers \
    && a2dissite 000-default.conf \
    && a2ensite apache-remus

EXPOSE 80

#
# Uncomment following lines to build Remus data as part of the image building process (bloats the image)
#
#RUN cd external_resources && bash download.sh && cd .. \
# && bash make_data_tree.sh

CMD systemctl reload apache2
CMD apachectl -D FOREGROUND