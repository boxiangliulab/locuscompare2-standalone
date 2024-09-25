#!/bin/bash

#install fastenloc
if command -v fastenloc >/dev/null 2>&1; then
	echo "  "
else
	echo "Please refer to https://github.com/xqwen/fastenloc/tree/dev/src to install fastenloc"
	exit 1
fi

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
# rm -rf ~/colotools_env
# mkdir ~/colotools_env && cd ~/colotools_env
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

# install SMR
echo "Installing SMR..."

if [[ "${platform}" == "Mac" ]];then
	smr_file="software/smr_v1.3.1/Mac/smr"
else
	smr_file="software/smr_v1.3.1/Linux/smr"
fi

cp ${smr_file} ${current_conda_bin}
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
if [[ $? != 0 ]];then
	echo "error installing Predixcan"
	exit 1
fi
ln -s $(pwd)/MetaXcan-master/software/SPrediXcan.py ${current_conda_bin}/SPrediXcan.py
if [[ $? != 0 ]];then
	echo "error installing Predixcan"
	exit 1
fi
predixcan_path=$(which SPrediXcan.py)
echo "PrediXcan full path is ${predixcan_path}"

#install finemap
echo "Installing finemap"
if [[ "${platform}" == "Mac" ]];then
	finemap_file="software/finemap_v1.4.2/Mac/finemap"
else
	finemap_file="software/finemap_v1.4.2/Linux/finemap"
fi

chmod 777 ${finemap_file}
cp ${finemap_file} ${current_conda_bin}

#install intact
echo "Installing INTACT..."
intact_file="software/INTACT_0.99.0.tar.gz"
# curl --connect-timeout 10 --retry 3 -o ${intact_file} https://www.bioconductor.org/packages/release/bioc/src/contrib/INTACT_1.4.0.tar.gz
if [[ $? != 0 ]];then
	echo "error downloading INTACT"
	exit 1
fi
R CMD INSTALL ${intact_file}
if [[ $? != 0 ]];then
	echo "error installing INTACT"
	exit
fi

#install twas
echo "Installing TWAS..."
twas_file="twas.zip"
curl --connect-timeout 10 --retry 3 -o ${twas_file} https://codeload.github.com/gusevlab/fusion_twas/zip/refs/heads/master
if [[ $? != 0 ]];then
	echo "error downloading TWAS"
	exit 1
fi
unzip -q ${twas_file}
if [[ $? != 0 ]];then
	echo "error installing TWAS"
	exit
fi
chmod -R 755 fusion_twas-master
ln -s $(pwd)/fusion_twas-master/FUSION.compute_weights.R ${current_conda_bin}/FUSION.compute_weights.R
if [[ $? != 0 ]];then
	echo "error installing TWAS"
	exit
fi
ln -s $(pwd)/fusion_twas-master/FUSION.assoc_test.R ${current_conda_bin}/FUSION.assoc_test.R
if [[ $? != 0 ]];then
	echo "error installing TWAS"
	exit
fi
twas_assoc_path=$(which FUSION.assoc_test.R)
echo "TWAS assoc_test full path is ${twas_assoc_path}"

#install plink2R
echo "Installing plink2R..."
plink2R_file="software/plink2R_1.1.tar.gz"

R CMD INSTALL ${plink2R_file}
if [[ $? != 0 ]];then
	echo "error installing plink2R"
	exit
fi

echo "All finished"
