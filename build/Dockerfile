FROM python:3.6.8

RUN apt-get update && \
    apt-get install -y apt-utils && \
    apt-get install -y build-essential && \
    apt-get install -y vim && \
    apt-get install -y man && \
    apt-get install -y plink1.9 && \
    ln -s /usr/bin/plink1.9 /usr/bin/plink

WORKDIR /app

RUN mkdir /app/testData
RUN mkdir /app/build
RUN mkdir /app/scratch
RUN mkdir /app/src

ADD ./requirements.txt /app/build
ADD ./src/ /app/src
ADD ./download1kG.sh /app/testData
ADD ./compiler.py /app
ADD ./corr.pyx /app
ADD ./corr.cpp /app
ADD ./plink2_stats.cpp /app

RUN pip install --trusted-host pypi.python.org -r /app/build/requirements.txt
