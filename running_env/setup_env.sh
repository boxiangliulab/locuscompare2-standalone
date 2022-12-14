#!/bin/bash

#set -xv

if [[ "$#" -ne 1 ]];then
  echo "usage: $0 path_to_conda_yml"
  exit 1
fi
if ! [[ -s "$1" ]];then
  echo "file $1 does not exist or is empty"
  exit 1
fi
if [[ -d "$1" ]];then
  echo "file $1 is a directory"
  exit 1
fi

echo "conda env yml path is: $1, creating env..."
conda env create -f $1
if [[ $? != 0 ]];then
	echo "error creating env from yml file"
	exit 1
fi

cond_env_name=$(cat $1|grep name -m 1| awk -F":" '{print $2}')
cond_env_name=$(echo $cond_env_name | tr -d '[:space:]')
echo "env name is $cond_env_name"

conda_base_env_path=$(conda info --base)
echo "base env path is $conda_base_env_path"
source ${conda_base_env_path}/etc/profile.d/conda.sh
conda activate ${cond_env_name}
if [[ $? != 0 ]];then
	echo "error activating env ${cond_env_name}"
	exit 1
fi

conda info

uname_out="$(uname -s)"
case "${uname_out}" in
    Linux*)     platform=Linux;;
    Darwin*)    platform=Mac;;
    *)          platform="Unknown"
esac

echo "Current platform is ${platform}"

if [[ "${platform}" == "Unknown" ]];then
	echo "Unsupported platform"
	exit 1
fi

#removing target directory if already exists
rm -rf ~/colotools_env
mkdir ~/colotools_env && cd ~/colotools_env
if [[ $? != 0 ]];then
	echo "error creating colotools_env directory"
	exit 1
fi

if [[ "${platform}" == "Mac" ]];then
	os_name=mac
else
	os_name=linux
fi

#retrieve conda bin diretory
current_conda_bin="${CONDA_PREFIX}/bin"
echo "current conda env bin path is ${current_conda_bin}"

# install JLIM
echo "Installing JLIM..."
jlim_file="jlimR_2.5.0.tar.gz"
curl --connect-timeout 10 --retry 3 -o ${jlim_file} https://biotech-coloc-hangzhou.oss-cn-hangzhou.aliyuncs.com/jlim/${jlim_file}
if [[ $? != 0 ]];then
	echo "error downloading JLIM"
	exit 1
fi
R CMD INSTALL ${jlim_file}
if [[ $? != 0 ]];then
	echo "error installing JLIM"
	exit 1
fi

# install SMR
echo "Installing SMR..."
smr_file="smr.zip"
smr_file_name=${smr_file%%.*}
curl --connect-timeout 10 --retry 3 -o ${smr_file} https://biotech-coloc-hangzhou.oss-cn-hangzhou.aliyuncs.com/smr/${os_name}/${smr_file}
if [[ $? != 0 ]];then
	echo "error downloading SMR"
	exit 1
fi
unzip -q ${smr_file}
cp ${smr_file_name} ${current_conda_bin}
if [[ $? != 0 ]];then
	echo "error installing SMR"
	exit 1
fi

#install predixcan
echo "Installing Predixcan..."
predixcan_file="predixcan.zip"
curl --connect-timeout 10 --retry 3 -o ${predixcan_file} https://codeload.github.com/hakyimlab/MetaXcan/zip/refs/heads/master
if [[ $? != 0 ]];then
	echo "error downloading Predixcan"
	exit 1
fi
unzip -q predixcan.zip
ln -s $(pwd)/MetaXcan-master/software/SPrediXcan.py ${current_conda_bin}/SPrediXcan.py
if [[ $? != 0 ]];then
	echo "error installing Predixcan"
	exit 1
fi
predixcan_path=$(which SPrediXcan.py)
echo "PrediXcan full path is ${predixcan_path}"

#install fastenloc
echo "Installing Fastenloc..."
torus_file="torus.zip"
torus_file_name=${torus_file%%.*}
fastenloc_file="fastenloc.zip"
fastenloc_file_name=${fastenloc_file%%.*}
curl --connect-timeout 10 --retry 3 -o ${torus_file} https://biotech-coloc-hangzhou.oss-cn-hangzhou.aliyuncs.com/fastenloc/${os_name}/${torus_file}
if [[ $? != 0 ]];then
	echo "error downloading Fastenloc"
	exit 1
fi
curl --connect-timeout 10 --retry 3 -o ${fastenloc_file} https://biotech-coloc-hangzhou.oss-cn-hangzhou.aliyuncs.com/fastenloc/${os_name}/${fastenloc_file}
if [[ $? != 0 ]];then
	echo "error downloading Fastenloc"
	exit 1
fi
unzip -q ${torus_file}
unzip -q ${fastenloc_file}
cp ${torus_file_name} ${fastenloc_file_name} ${current_conda_bin}
if [[ $? != 0 ]];then
	echo "error installing Fastenloc"
	exit 1
fi

#install finemap
echo "Installing finemap"
finemap_file="finemap"
curl --connect-timeout 10 --retry 3 -o ${finemap_file} https://biotech-coloc-hangzhou.oss-cn-hangzhou.aliyuncs.com/finemap/${os_name}/${finemap_file}
if [[ $? != 0 ]];then
    echo "error downloading finemap"
    exit 1
fi
chmod 777 ${finemap_file}
cp ${finemap_file} ${current_conda_bin}

#install intact
echo "Installing INTACT..."
intact_file="INTACT_0.99.0.tar.gz"
curl --connect-timeout 10 --retry 3 -o ${intact_file} https://biotech-coloc-hangzhou.oss-cn-hangzhou.aliyuncs.com/intact/${intact_file}
if [[ $? != 0 ]];then
	echo "error downloading INTACT"
	exit 1
fi
R CMD INSTALL ${intact_file}
if [[ $? != 0 ]];then
	echo "error installing INTACT"
	exit
fi

echo "All finished"
