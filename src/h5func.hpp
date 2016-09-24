#ifndef H5FUNC_HPP
#define H5FUNC_HPP
#include "RcppArmadillo.h"
#include <vector>
#include <memory>
#include <H5Cpp.h>


using namespace H5;

bool file_exists (const std::string& name);

typedef std::shared_ptr<H5::H5File> H5FilePtr;

H5FilePtr create_or_open_file(const std::string &fname);


typedef std::shared_ptr<H5::Group> H5GroupPtr;

H5GroupPtr create_or_open_group(H5FilePtr &file,const std::string &groupname);


typedef std::shared_ptr<DataSet> H5DataSetPtr;

H5DataSetPtr create_or_open_dataset(H5GroupPtr &group,const std::string &dataname, const DataType &data_type,std::vector<hsize_t> &cdatadim,std::vector<hsize_t> &mdatadim,std::vector<hsize_t> &chunkdim,const unsigned int deflate_level);

ArrayType read_arraytype(const DataSet *dataset, const PredType pt);

std::vector<int> read_int_h5(const std::string h5file, const std::string groupname, const std::string dataname);

std::vector<unsigned int> read_uint_h5(const std::string h5file, const std::string groupname, const std::string dataname);

arma::uvec intersect_col(const std::string h5file1, const std::string h5groupname1, const std::string h5colname1, const std::string h5file2, const std::string h5groupname2, const std::string h5colname2);
arma::fmat read_fmat_h5(const std::string hap_h5file,const std::string groupname, const std::string dataname, size_t offset, size_t chunksize);
arma::mat read_dmat_h5(const std::string hap_h5file,const std::string groupname, const std::string dataname, const size_t offset, const size_t chunksize);
std::vector<float> read_float_h5(const std::string h5file, const std::string groupname, const std::string dataname,const size_t offset,const size_t chunksize);

arma::mat read_dmat_chunk_ind(const std::string h5file,const std::string groupname, const std::string dataname, const arma::uvec indvec);
arma::fmat read_fmat_chunk_ind(const std::string h5file,const std::string groupname, const std::string dataname, const arma::uvec indvec);
H5::DataType read_data_type(const std::string h5file, const std::string groupname, const std::string dataname);
std::vector<std::string> getObjects(const std::string h5file, const std::string groupname);



hsize_t get_arraysize(ArrayType &atype);

size_t get_arraydim(const std::string h5file, const std::string groupname,const std::string dataname);
size_t write_mat_h5(const std::string h5file, const std::string groupname, const std::string dataname,const hsize_t Nsnps, const hsize_t Nind, arma::mat &data,const unsigned int deflate_level);
size_t write_mat_h5(const std::string h5file, const std::string groupname, const std::string dataname,const hsize_t Nsnps, const hsize_t Nind, arma::fmat &data,const unsigned int deflate_level);

size_t write_umat_h5(const std::string h5file, const std::string groupname, const std::string dataname,const hsize_t Nsnps, const hsize_t Nind, arma::umat &data,const unsigned int deflate_level);

int write_dmatrix_h5(Rcpp::String h5file,Rcpp::String groupname, Rcpp::String dataname, Rcpp::IntegerVector Nsnps, Rcpp::IntegerVector Nind, Rcpp::NumericMatrix data,const unsigned int deflate_level);

int write_Rint_h5(const std::string h5file, const std::string groupname, const std::string dataname, Rcpp::IntegerVector &data,const unsigned int deflate_level);

int write_Rnumeric_h5(const std::string h5file, const std::string groupname, const std::string dataname, Rcpp::NumericVector &data,const unsigned int deflate_level);

size_t write_int_h5(const std::string h5file, const std::string groupname, const std::string dataname,const hsize_t Nsnps, arma::uvec &data,const unsigned int deflate_level);

size_t write_uint_h5(const std::string h5file, const std::string groupname, const std::string dataname,const hsize_t Nsnps, arma::uvec &data,const unsigned int deflate_level);

size_t get_rownum_h5(const std::string hap_h5file,const std::string groupname, const std::string dataname);

#endif
