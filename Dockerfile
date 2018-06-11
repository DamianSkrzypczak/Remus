#
# Dockerfile for building docker container with the app
#


FROM rastasheep/ubuntu-sshd

RUN apt update && \
    apt install -y \
	build-essential \
	python3-dev \
	python3-pip \
	zlib1g-dev \
	liblzma-dev \
	bedtools \
	vim \
	git && \
    pip3 install \
	pandas \
	pybedtools

