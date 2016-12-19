#RColumbo

Installing RColumbo (without root privelages) (on RCC) 

## 0. Load required modules and edit  `.bash_profile`

`module load cmake`

`module load R`

Add this line to the end of `~/.bash_profile`:

`export LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH`

then run 
`source ~/.bash_profile`


## 1. Download and install HDF5

The HDF5 module on RCC does not support C++

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

`cmake -DCMAKE_INSTALL_PREFIX=$HOME ..`

`make`

`make install`


## 3. Install RColumbo

from within R
`#Install devools if necessary`

`install.packages(devtools)`

`library(devtools)`

`#Install rhdf5 if necessary`

`source("https://bioconductor.org/biocLite.R")`

`biocLite("rhdf5")`

`devtools::install_github("CreRecombinase/RColumbo",ref="release")`

`library(RColumbo)`

## 4. Read compressed HDF5 files really fast

The eQTL results are stored as dataframes in HDF5 (.h5) files. There is one file for each chromosome. Each h5 file
has a `cis_eQTL` dataframe and a  `trans_eQTL` dataframe. The dataframes can be read using `RColumbo`'s `read_h5_df` function:

`library(RColumbo)`

`eqtl_file <- "/project/xinhe/eQTL/WB_Chr2_v6p_ortho_flip.h5"`

`cis_eqtl_df <- read_h5_df(eqtl_file,groupname="cis_eQTL")`

`trans_eqtl_df <- read_h5_df(eqtl_file,groupname="trans_eQTL")`

The columns of these dataframes are as follows
* `chrom`: chromosome
* `fgeneid`: a gene identifier (it's a basically the ensembl gene ID)
* `pos`: chromosome coordinate as given by GTEx
* `theta`:estimated eQTL regression coefficient (no intercept)
* `serr`: standard error of eQTL regression coefficient

