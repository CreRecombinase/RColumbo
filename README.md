#RColumbo

Installing RColumbo (without root privelages) (on RCC) 

## 0. Load required modules

`module load cmake`
`module load R`

## 1. Download and install HDF5

The module on RCC does not support C++

`wget https://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.0-patch1/src/hdf5-1.10.0-patch1.tar.gz`
`tar -xzf hdf5-1.10.0-patch1.tar.gz`
`cd hdf5-1.10.0-patch1`
`./configure --prefix=$HOME --enable-cxx --enable-build-mode=production`
`make install`

## 2. Download blosc


`git clone https://github.com/Blosc/c-blosc`
`cd c-blosc`
`mkdir build`
`cd build`
`cmake -DCMAKE_INSTALL_PREFIX=$HOME`

## 3. Install RColumbo

from within R
`#Install devools if necessary`
`install.packages(devtools)`
`library(devtools)`
`install_github("CreRecombinase/RColumbo",ref="release")`
`library(RColumbo)`

## 4. Read compressed HDF5 files really fast

cis eQTL