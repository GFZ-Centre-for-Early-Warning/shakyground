FROM ubuntu:18.04

# for not having interaction on installation process
ENV DEBIAN_FRONTEND noninteractive

WORKDIR /root

COPY requirements.txt .

RUN apt-get update && \
    apt-get install python3.7 python3-pip gmt gmt-dcw gmt-gshhg python3-tk unzip -y && \
    pip3 install -r requirements.txt

WORKDIR /usr/share/git/shakyground
RUN wget https://earthquake.usgs.gov/static/lfs/data/global_vs30_grd.zip && unzip global_vs30_grd.zip -d . && rm global_vs30_grd.zip
COPY . .
