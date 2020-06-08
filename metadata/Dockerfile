FROM gfzriesgos/shakyground-grid-file:20200605

# for not having interaction on installation process
ENV DEBIAN_FRONTEND noninteractive

WORKDIR /usr/share/git/shakyground

COPY requirements.txt .

RUN apt-get update && \
    apt-get install -y \
        libgirepository1.0-dev \
        python3-pip \
        gmt \
        gmt-dcw \
        gmt-gshhg \
        python3-tk \
        libxml2-utils \
        unzip && \
    pip3 install -r requirements.txt

COPY . .
