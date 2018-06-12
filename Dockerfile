#
# Dockerfile for building docker container with the app
#


FROM rastasheep/ubuntu-sshd

RUN apt update \
     && apt install -y \
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
