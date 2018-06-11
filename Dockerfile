#
# Dockerfile for building docker container with the app
#


FROM rastasheep/ubuntu-sshd

RUN apt update && \
    apt install -y \
	build-essential \
	python3-dev \
	python3-pip \
	vim \
	git

