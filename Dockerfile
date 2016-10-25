FROM ubuntu:14.04

WORKDIR /opt

RUN apt-get update && apt-get install -y \
libcurl4-gnutls-dev \
libgnutls-dev \
python \
python-pip \ 
python-dev \ 
build-essential \
libncurses5-dev \
libncursesw5-dev \
pypy \
git \
wget \
r-base \
zlib1g-dev \
samtools

RUN pip install --upgrade pip
RUN pip install --upgrade virtualenv

RUN mkdir /opt/bin 
ENV PATH /opt/bin:$PATH

# Install R
RUN wget -O /tmp/ada_2.0-3.tar.gz https://cran.r-project.org/src/contrib/ada_2.0-5.tar.gz; \
R CMD INSTALL /tmp/ada_2.0-3.tar.gz

# Install BWA
RUN git clone https://github.com/lh3/bwa.git; \
cd bwa; \
make; \
ln -s /opt/bwa/bwa /opt/bin/

# Install seqtk
RUN git clone https://github.com/lh3/seqtk.git; \
cd seqtk; \
make; \
ln -s /opt/seqtk/seqtk /opt/bin/

# Install bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz; \
tar xvfz bedtools-2.26.0.tar.gz; \
cd bedtools2; \
make; \
ln -s /opt/bedtools2/bin/bedtools /opt/bin/

# Install BLAT
RUN wget -P /opt/bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v287/blat/blat; \
chmod +x /opt/bin/blat

# Install ericscript
RUN git clone https://github.com/cgrlab/EricScript.git; \
chmod +x EricScript/ericscript.pl
