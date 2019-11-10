FROM ubuntu:latest

RUN apt-get update \
	&& apt-get install -y \
		python3-pip \
		make \
		git \
		curl \
		zlib1g-dev \
		python3-all-dev \
		libhdf5-dev \
		libatlas-base-dev \
		libopenblas-base \
		libopenblas-dev \
		libbz2-dev \
		liblzma-dev \
		libffi-dev  \
		python-virtualenv \
		cmake \
		wget \
		bzip2 \
	&& curl -O https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh \
	&& apt-get install default-jre -y --no-install-recommends \
	&& wget https://github.com/broadinstitute/pilon/releases/download/v1.23/pilon-1.23.jar \
	&& cp pilon-1.23.jar /usr/local/bin \
	&& wget https://github.com/shenwei356/seqkit/releases/download/v0.11.0/seqkit_linux_amd64.tar.gz \
	&& tar -zxf seqkit_linux_amd64.tar.gz \
	&& cp seqkit /usr/local/bin \
	&& git clone https://github.com/lh3/minimap2 \
	&& pip3 install git+https://github.com/rrwick/Porechop \
	&& pip3 install pomoxis \
	&& virtualenv pomoxis --python=python3 --prompt "(pomoxis) "  
