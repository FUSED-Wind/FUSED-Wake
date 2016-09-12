FROM continuumio/anaconda3
MAINTAINER Pierre Elouan Rethore <pe@retho.re>

RUN apt-get update \
 && apt-get install -y \
    build-essential \
    g++ \
    gfortran

RUN conda install -y mkl \
 && conda update conda \
 && conda install -y accelerate

RUN mkdir /install
WORKDIR /install
COPY setup.py /install/
COPY *.rst /install/
COPY fusedwake /install/fusedwake
COPY requirements_dev.txt /install/
RUN pip install -r /install/requirements_dev.txt
RUN python setup.py install

RUN pip install future
RUN pip install Theano
RUN pip install pydot-ng


CMD bash
