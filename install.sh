#!/usr/bin/env bash

set -e
shopt -s nullglob

echo "Welcome to the microSALT installation script. Q to exit"
while true; do
    echo "What name would you like to give your microSALT environment?"
    read input
    if [[ $input = "q" ]] || [[ $input = "Q" ]]; then
        break
    else
        cname=$input
        break
    fi
done
while true; do
    echo "Would you like a 'release' or local 'source' environment?"
    read input
    if [[ $input = "q" ]] || [[ $input = "Q" ]]; then
        break
    elif [[ $input == "source" ]]  || [[ $input == "release" ]]; then
        type=$input
        validbranch=false
        while ! $validbranch; do
            echo "Name a branch to install (or just 'master')":
            read input
            branch=$input
            curl https://raw.githubusercontent.com/Clinical-Genomics/microSALT/$branch/LICENSE | grep -q 'License' && validbranch=true||echo "Invalid branch name"
        done
        break
    fi
done
echo "Thank you, setting up environment $cname!"

conda info| grep -q $rname && source deactivate || ""

#Accepts that environment doesnt exist
conda remove -y -n $cname --all || ""
conda create -y -n $cname python=3.6
#source deactivate
source activate $cname
conda config --add channels bioconda
conda install -y -c bioconda blast=2.9.0 bwa=0.7.17 picard=2.20.3 quast=5.0.2 samtools=1.9 spades=3.13.1 trimmomatic=0.39
if [[ $type == "release" ]]; then
    pip install -r https://raw.githubusercontent.com/Clinical-Genomics/microSALT/master/requirements.txt 
    pip install -U git+https://github.com/Clinical-Genomics/microSALT 
elif [[ $type == "source" ]]; then
  HERE=$PWD
  if [ -d ${HERE}/microSALT ]; then
    rm -rf microSALT
  fi
  git clone https://github.com/Clinical-Genomics/microSALT
  cd microSALT && git checkout $branch
  pip install -r requirements.txt && pip install -e . && cd ${HERE}
  echo "Source installed under ${HERE}/microSALT" 
fi 
echo "Installation Complete!"
while true; do
    echo "Configuration requires manual set-up as described in README.md ['ok']"
    read input
    if [[ $input = "ok" ]] || [[ $input = "OK" ]]; then
        break
    fi
done
