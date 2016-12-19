#ifndef SNP_EXP_HPP
#define SNP_EXP_HPP
#include "RcppArmadillo.h"
#include <iterator>
#include "H5IO.hpp"
#include <unordered_map>
#include <string>
#include <H5Cpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>
#include "h5func.hpp"


class genepos_map{
private:
  size_t size;
  std::string gene_pos_file;
public:
  std::unordered_map<arma::uword,arma::uword> gp_map;
  arma::umat annomat;
  genepos_map(std::string tdbsnpfile);
};

class dbsnpmap{
private:
  size_t size;
  std::string dbsnpfile;
public:
  std::unordered_map<arma::uword,arma::uword> dbmap;
  arma::uvec posmap;
  dbsnpmap(std::string tdbsnpfile);
};

std::unordered_map<arma::uword,arma::uword> make_map(arma::uvec &chrom, arma::uvec &pos,arma::uvec &rsid,arma::uvec &posmap);


class Raw_Dataset{
private:
  const std::string h5filename;
  const std::string data_groupname;
  const std::string data_dataname;
  const std::string anno_groupname;
  const std::string anno_dataname;
  const bool isCont;

  size_t offset;
  arma::uvec data_index;
  DataType datatype;
public:
  Raw_Dataset(const std::string th5filename, const std::string tdata_groupname, const std::string tdata_dataname, const std::string tanno_groupname, const std::string tanno_dataname);
  Raw_Dataset(const std::string th5filename, const std::string tdata_groupname, const std::string tdata_dataname, const std::string tanno_groupname, const std::string tanno_dataname, const arma::uvec index);
  arma::fmat read_and_offset(const size_t num);
  size_t resize_chunk(const size_t rchunksize)const;
  arma::fmat read_chunk(const size_t num) const;
  arma::fvec read_ffeature(const std::string groupname,const std::string dataname)const;
  arma::fmat read_elem_ind(const arma::uvec ind) const;
  arma::fmat read_elem_direct_ind(const arma::uvec ind) const;
  size_t set_offset(const size_t off);
  size_t increment_offset(const size_t off);
  size_t reset_offset();
  size_t get_offset() const;
  arma::fmat read_elem_name(const arma::uvec ind) const;
  arma::uvec get_index()const;
  arma::uvec get_anno_index()const;
  size_t Nind;
  size_t P;
};



class LD_dataset:public Raw_Dataset{
private:
public:
  LD_dataset(const std::string h5filename);
  LD_dataset(const std::string h5filename, arma::uvec tindex);
  arma::fvec mapvec;
  arma::fvec get_map_chunk(const size_t chunksize);
};

class Gwas_Dataset:public Raw_Dataset{
private:
  LD_dataset LD_data;
public:
  Gwas_Dataset(const std::string h5filename);
  Gwas_Dataset(const std::string h5filename, arma::uvec tindex);
  arma::fvec betavec;
  arma::fvec get_map_chunk(const size_t chunksize);
};








#endif
