#
# Dockerfile for building docker container with the app
#


FROM rastasheep/ubuntu-sshd

RUN apt update && \
    apt install -y \
	build-essentials \
	python-dev \
	python-pip \
	vim \
	git

