FROM ubuntu:18.04

# for not having interaction on installation process
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get install python3.7 python3-pip git gmt gmt-dcw gmt-gshhg python3-tk wget unzip -y && \
    pip3 install pandas scipy lxml h5py shapely zmq mock matplotlib decorator scipy psutil pygmt

WORKDIR /usr/share/git/shakyground
RUN wget https://earthquake.usgs.gov/static/lfs/data/global_vs30_grd.zip && unzip global_vs30_grd.zip -d . && rm global_vs30_grd.zip
COPY . .
