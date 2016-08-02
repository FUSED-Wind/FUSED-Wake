## To build
# $ docker build -t fwp .
#
## To run the Jupyter notebook:
# $ docker run -it -p 8898:8898 fwp
#
## To enter the running container
# $ docker exec -it fwp /bin/bash
#

FROM continuumio/anaconda

MAINTAINER Pierre-Elouan Rethore <pe@retho.re>

RUN apt-get update \
 && apt-get install -y \
    build-essential \
    gfortran

RUN conda install jupyter -y --quiet

ENV NOTEBOOKS    /opt/notebooks
ENV HOME         /opt/notebooks/home
ENV INSTALL      /opt/install

RUN mkdir $NOTEBOOKS
RUN mkdir $NOTEBOOKS/FUSED-Wake
RUN mkdir $HOME
RUN mkdir $INSTALL

# Installing the python dependencies -------------------------------------------
WORKDIR $INSTALL
COPY requirements_dev.txt $INSTALL

RUN pip install --upgrade pip && \
    pip install -r requirements_dev.txt

# Installing the FUSED-Wake code -------------------------------------------------
COPY . $NOTEBOOKS/FUSED-Wake

WORKDIR $NOTEBOOKS/FUSED-Wake
RUN python setup.py install

# Preparing to exit ------------------------------------------------------------
WORKDIR $NOTEBOOKS
EXPOSE 8898

CMD jupyter notebook \
    --notebook-dir=$NOTEBOOKS \
    --ip='*' \
    --port=8898 \
    --no-browser
