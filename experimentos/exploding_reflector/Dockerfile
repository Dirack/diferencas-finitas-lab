FROM ghcr.io/dirack/madagascarenvironment:1.0.1
USER root
ENV PYTHONPATH=/home/tryitondocker/madagascar/lib/python2.7/dist-packages
ENV PATH=/home/tryitondocker/madagascar/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV RSFSRC=/home/tryitondocker/madagascar-3.0
ENV DATAPATH=/home/tryitondocker/rsfdata/
RUN echo "Installing madagascar..." && cd /home/tryitondocker/madagascar-3.0/user/fomels && scons && cd /home/tryitondocker/madagascar-3.0 && scons install && echo "...done"
WORKDIR /home/tryitondocker/experimentos
