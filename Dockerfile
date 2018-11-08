#
# Dockerfile for building docker container with the app
#

FROM rastasheep/ubuntu-sshd

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

COPY ./ /var/www/remus
WORKDIR /var/www/remus

RUN pip3 install -r requirements.txt

COPY apache-remus.conf /etc/apache2/sites-available/apache-remus.conf
RUN a2ensite apache-remus
RUN a2enmod headers

RUN a2dissite 000-default.conf
RUN a2ensite apache-remus.conf

EXPOSE 80

#RUN bash external_resources/download.sh \
# && bash make_data_tree.sh

#CMD  /usr/sbin/apache2ctl -D BACKGROUND
