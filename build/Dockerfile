FROM continuumio/miniconda3

RUN conda create -n env python=3.6.8
RUN echo "source activate env" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH
ENV SSH_PASSWD "root:Docker!"
RUN apt-get update && \
    apt-get install -y apt-utils && \
    apt-get install -y build-essential && \
    apt-get install -y vim && \
    apt-get install -y man && \
    apt-get install -y --no-install-recommends dialog && \
    apt-get install -y --no-install-recommends openssh-server && \
#    ln -s /usr/bin/plink1.9 /usr/bin/plink
    echo "$SSH_PASSWD" | chpasswd

WORKDIR /app

RUN mkdir /app/testData
RUN mkdir /app/build
RUN mkdir /app/scratch
RUN mkdir /app/src

ADD build/init.sh /app/build
ADD build/requirements.txt /app/build
ADD src/ /app/src
#ADD build/download1kG.sh /app/testData
ADD build/compiler.py /app
ADD build/corr.pyx /app
ADD build/corr.cpp /app
ADD build/plink2_stats.cpp /app
COPY build/sshd_config /etc/ssh/

RUN pip install --trusted-host pypi.python.org -r /app/build/requirements.txt
RUN conda install h5py==2.8.0 
ENV SSH_PASSWD "root:Docker!"
RUN chmod u+x /app/build/init.sh
ENV PORT 8080
ENV SSH_PORT 2222
EXPOSE 8080 2222 80
ENTRYPOINT ["/app/build/init.sh"]
