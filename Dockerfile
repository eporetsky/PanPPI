# docker build --rm -f Dockerfile -t externelly/panppi:app .
# docker run --name=panppi -p 8080:8080 externelly/panppi:app 
# On the browser open: 0.0.0.0:8080

FROM mambaorg/micromamba:jammy
LABEL maintainer="externelly"

USER root

ARG MAMBA_DOCKERFILE_ACTIVATE=1


RUN apt-get update -y && \ 
    apt-get install -y wget curl git tar bzip2 unzip

RUN git clone https://github.com/eporetsky/PanPPI && \
    cd PanPPI && \
    micromamba update -y -n base --file environments.yml && \
    cd Dash && \
    wget -O panppi.zip "https://www.dropbox.com/scl/fi/h2dlndcilk2inh67oc9nv/panppi.zip?rlkey=qu3iqv7k0j416hqji7sswdzwt&dl=1" && \
    unzip panppi.zip && \
    rm panppi.zip && \ 
    mkdir tmp

WORKDIR /tmp/PanPPI/Dash

EXPOSE 8080

CMD ["python", "app.py"]