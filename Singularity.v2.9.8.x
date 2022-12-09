Bootstrap: docker
From: ubuntu:20.04
IncludeCmd: yes

%labels
  MAINTAINER cschu (cschu1981@gmail.com)
  VERSION v.2.9.7

%environment
export LC_ALL=C
export PATH=$PATH:/opt/software/miniconda3/bin:/opt/software/mOTUs

%post
  apt-get update
  apt-get install -y apt-transport-https apt-utils software-properties-common
  apt-get install -y add-apt-key
  export DEBIAN_FRONTEND=noninteractive
  ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime
  apt-get install -y tzdata
  dpkg-reconfigure --frontend noninteractive tzdata

  apt-get install -y wget python3-pip git gawk

  mkdir -p /opt/software

  # all the conda packages won't work together.. ><;
  wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh 
  bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/software/miniconda3
  rm -f Miniconda3-latest-Linux-x86_64.sh

  /opt/software/miniconda3/bin/conda install -y -c bioconda 'samtools>=1.13' bwa minimap2 bwa-mem2
  git clone https://github.com/cschu/gff_quantifier.git
  cd gff_quantifier
  git checkout fix/deal_with_missing_ambig_counts_20221209
  /opt/software/miniconda3/bin/python setup.py bdist_wheel
  /opt/software/miniconda3/bin/pip install --force-reinstall dist/*whl
  # /opt/software/miniconda3/bin/pip install -e .
  cd ..
  rm -rf gff_quantifier/


# trigger trigger trigger trigger trigger trigger trigger
