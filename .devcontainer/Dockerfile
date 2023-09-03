FROM centos:8

LABEL maintainer="dimtsap@hotmail.com"

# Resolves no URLs in mirrorlist error (https://stackoverflow.com/questions/70926799/centos-through-a-vm-no-urls-in-mirrorlist)
RUN sed -i 's/mirrorlist/#mirrorlist/g' /etc/yum.repos.d/CentOS-Linux-* \
    && sed -i 's|#baseurl=http://mirror.centos.org|baseurl=http://vault.centos.org|g' /etc/yum.repos.d/CentOS-Linux-*

RUN dnf install centos-release-stream -y \
    && dnf swap centos-{linux,stream}-repos -y \
    && dnf distro-sync -y

# Install base utilities
RUN yum install wget -y

# Install miniconda
RUN mkdir -p ~/miniconda3 \
    && wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh \
    && bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3 \
    && rm -rf ~/miniconda3/miniconda.sh

RUN ~/miniconda3/bin/conda init bash \
    && ~/miniconda3/bin/conda init zsh

RUN dnf install libgfortran -y

RUN ~/miniconda3/bin/conda install -c conda-forge uqpy \
    && ~/miniconda3/bin/conda install -c conda-forge shapely \
    && ~/miniconda3/bin/conda install -c conda-forge seaborn
