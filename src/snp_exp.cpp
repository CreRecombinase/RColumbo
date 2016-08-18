#include "RcppArmadillo.h"
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_match.hpp>
#include <boost/spirit/include/qi_uint.hpp>
#include <iterator>
#include <boost/unordered_map.hpp>
#include <numeric>
#include <fstream>
#include <string>
#include <H5Cpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>


using namespace H5;

bool file_exists (const std::string& name) {
  if (FILE *file = fopen(name.c_str(), "r")) {
    fclose(file);
    return true;
  } else {
    return false;
  }
}

typedef std::shared_ptr<H5File> H5FilePtr;

H5FilePtr create_or_open_file(const std::string &fname)
{
  H5::H5File* file;
  if(!file_exists(fname)){
    file = new H5::H5File(fname.c_str(), H5F_ACC_TRUNC);
  }else{
    file = new H5::H5File(fname.c_str(), H5F_ACC_RDWR);
  }
  return H5FilePtr(file);
}
typedef std::shared_ptr<Group> H5GroupPtr;

H5GroupPtr create_or_open_group(H5FilePtr &file,const std::string &groupname)
{
  Group* group;
  Group *rg = new Group(file->openGroup("/"));
  hsize_t objc= rg->getNumObjs();
  bool fgroup=false;
  if(objc!=0){
    for(hsize_t i=0; i<objc;i++){
      std::string tst=rg->getObjnameByIdx(i);

      if(tst==groupname){
        fgroup = true;
      }
    }
  }
  if(fgroup){
    group = new Group(file->openGroup(groupname));
  }else{
    group = new Group(file->createGroup("/"+groupname));
  }
  return H5GroupPtr(group);
}

typedef std::shared_ptr<DataSet> H5DataSetPtr;

H5DataSetPtr create_or_open_dataset(H5GroupPtr &group,const std::string &dataname, const DataType &data_type,std::vector<hsize_t> &cdatadim,std::vector<hsize_t> &mdatadim,std::vector<hsize_t> &chunkdim)
{
  H5::Exception::dontPrint();
  DataSet* dataset;
  DataSpace* fdataspace;
  hsize_t objc= group->getNumObjs();
  bool fdat=false;
  if(objc!=0){
    for(hsize_t i=0; i<objc;i++){
      std::string tst=group->getObjnameByIdx(i);
      if(tst==dataname){
        fdat = true;
      }
    }
  }
  if(fdat){
    dataset = new DataSet(group->openDataSet(dataname));
  }else{
    hsize_t *cumdima = new hsize_t[cdatadim.size()];
    hsize_t *chunkdima = new hsize_t[chunkdim.size()];
    hsize_t*mdima = new hsize_t[mdatadim.size()];

    for(int i=0; i<cdatadim.size();i++){
      cumdima[i]=cdatadim[i];
      mdima[i]=mdatadim[i];
      chunkdima[i]=chunkdim[i];
    }
    fdataspace= new DataSpace(cdatadim.size(),cumdima,mdima); //Create dataspace for dataset (on disk)
    DSetCreatPropList cparms; //Create chunksize file parameters
    cparms.setChunk(chunkdim.size(),chunkdima); //Set chunksize
    cparms.setDeflate(2);
    dataset = new DataSet(group->createDataSet(dataname,data_type,*fdataspace,cparms));

  }
  return H5DataSetPtr(dataset);
}


size_t write_mat_h5(const std::string h5file, const std::string groupname, const std::string dataname,const hsize_t Nsnps, const hsize_t Nind, arma::mat &data){

  hsize_t chunkstart=0;
  size_t chunksize = data.size()/Nind;
  std::vector<hsize_t> cumdim{0};
  std::vector<hsize_t> maxdim{Nsnps};
  std::vector<hsize_t> chunkdim{chunksize};
  hsize_t adim[1];//dimension of each array element (ncols)
  adim[0]=Nind; //Assign array size dimension
  ArrayType ahtypew(PredType::NATIVE_DOUBLE,1,adim); //Create array size type
  H5FilePtr file =create_or_open_file(h5file);
  H5GroupPtr group=create_or_open_group(file,groupname);
  H5DataSetPtr dataset = create_or_open_dataset(group,dataname,ahtypew,cumdim,maxdim,chunkdim);
  DataSpace* fdataspace = new DataSpace(dataset->getSpace());
  hsize_t datadim[1];
  fdataspace->getSimpleExtentDims(datadim,NULL);

  chunkstart = datadim[0];
  datadim[0]=datadim[0]+chunksize;
  dataset->extend(datadim);
  hsize_t memdim[]={chunksize};
  fdataspace->close();
  delete fdataspace;
  fdataspace = new DataSpace(dataset->getSpace());
  fdataspace->getSimpleExtentDims(datadim,NULL);

  DataSpace *mspace= new DataSpace(1,memdim); //Size of first dataset (in memory, can be bigger or smaller than size on disk, depending on how much you're writing)
  hsize_t odim[]={chunkstart};//dimension of each offset (current_chunk*chunksize)
  hsize_t stridea[]={1};
  hsize_t blocka[]={1};
  //  std::cout<<"Getting Data dimensions"<<std::endl;
  fdataspace->selectHyperslab( H5S_SELECT_SET, memdim, odim,stridea,blocka);
  std::cout<<"Starting to write data"<<std::endl;
  dataset->write(data.memptr(),ahtypew,*mspace,*fdataspace);

  dataset->close();
  fdataspace->close();
  mspace->close();
  group->close();
  file->close();
  delete fdataspace;
  delete mspace;
  std::cout<<"Data written and handles closed"<<std::endl;
  return(0);
}







size_t write_uint_h5(const std::string h5file, const std::string groupname, const std::string dataname,const size_t Nsnps, arma::uvec &data){

  typedef std::vector<unsigned int> uistd;
  std::vector<unsigned int>tdata = arma::conv_to<uistd>::from(data);
  hsize_t chunkstart=0;
  size_t chunksize = tdata.size();
  std::vector<hsize_t> cumdim{0};
  std::vector<hsize_t> maxdim{Nsnps};
  std::vector<hsize_t> chunkdim{chunksize};
  H5FilePtr file =create_or_open_file(h5file);
  H5GroupPtr group=create_or_open_group(file,groupname);
  H5DataSetPtr dataset = create_or_open_dataset(group,dataname,PredType::NATIVE_UINT,cumdim,maxdim,chunkdim);

  DataSpace* fdataspace = new DataSpace(dataset->getSpace());
  hsize_t datadim[1];
  fdataspace->getSimpleExtentDims(datadim,NULL);
  chunkstart = datadim[0];
  datadim[0]=datadim[0]+chunksize;
  dataset->extend(datadim);
  hsize_t memdim[]={chunksize};
  fdataspace->close();
  delete fdataspace;
  fdataspace = new DataSpace(dataset->getSpace());
  fdataspace->getSimpleExtentDims(datadim,NULL);

  DataSpace *mspace= new DataSpace(1,memdim); //Size of first dataset (in memory, can be bigger or smaller than size on disk, depending on how much you're writing)
  hsize_t odim[]={chunkstart};//dimension of each offset (current_chunk*chunksize)
  hsize_t stridea[]={1};
  hsize_t blocka[]={1};
  //  std::cout<<"Getting Data dimensions"<<std::endl;
  fdataspace->selectHyperslab( H5S_SELECT_SET, memdim, odim,stridea,blocka);
  std::cout<<"Starting to write data"<<std::endl;
  dataset->write(&tdata[0],PredType::NATIVE_UINT,*mspace,*fdataspace);
  return(0);
}


size_t write_int_h5(const std::string h5file, const std::string groupname, const std::string dataname,const size_t Nsnps, arma::uvec &data){


  typedef std::vector<int> istd;
  std::vector<int>tdata = arma::conv_to<istd>::from(data);
  hsize_t chunkstart=0;
  size_t chunksize = tdata.size();
  std::vector<hsize_t> cumdim{0};
  std::vector<hsize_t> maxdim{Nsnps};
  std::vector<hsize_t> chunkdim{chunksize};
  H5FilePtr file =create_or_open_file(h5file);
  H5GroupPtr group=create_or_open_group(file,groupname);
  H5DataSetPtr dataset = create_or_open_dataset(group,dataname,PredType::NATIVE_INT,cumdim,maxdim,chunkdim);

  DataSpace* fdataspace = new DataSpace(dataset->getSpace());
  hsize_t datadim[1];
  fdataspace->getSimpleExtentDims(datadim,NULL);
  chunkstart = datadim[0];
  datadim[0]=datadim[0]+chunksize;
  dataset->extend(datadim);
  hsize_t memdim[]={chunksize};
  fdataspace->close();
  delete fdataspace;
  fdataspace = new DataSpace(dataset->getSpace());
  fdataspace->getSimpleExtentDims(datadim,NULL);

  DataSpace *mspace= new DataSpace(1,memdim); //Size of first dataset (in memory, can be bigger or smaller than size on disk, depending on how much you're writing)
  hsize_t odim[]={chunkstart};//dimension of each offset (current_chunk*chunksize)
  hsize_t stridea[]={1};
  hsize_t blocka[]={1};
  //  std::cout<<"Getting Data dimensions"<<std::endl;
  fdataspace->selectHyperslab( H5S_SELECT_SET, memdim, odim,stridea,blocka);
  std::cout<<"Starting to write data"<<std::endl;
  dataset->write(&tdata[0],PredType::NATIVE_INT,*mspace,*fdataspace);
  return(0);
}




arma::uvec isFlip(std::vector<std::string> &ref,std::vector<std::string> &alt,arma::uvec indices){
  if(ref.size()!=alt.size()){
    Rcpp::Rcout<<"ref is of size:"<<ref.size()<<" alt is of size:"<<alt.size()<<std::endl;
    Rcpp::stop("size of ref not equal to size of alt!");
  }
  arma::uvec ret(indices.n_elem);
  for(size_t i=0;i<ret.size();i++){
    if(ref[indices[i]]<alt[indices[i]]){
      ret[i]=1;
    }else{
      ret[i]=0;
    }
  }
  return(ret);
}

size_t makeFlip(arma::uvec &isFlipv, arma::mat genomat){
  size_t ct;
  if(isFlipv.n_elem!=genomat.n_cols){
    Rcpp::Rcerr<<"isFlipv has:"<<isFlipv.n_elem<<" and genomat has colnum:"<<genomat.n_cols<<std::endl;
    Rcpp::stop("error in makeFlip!");
  }
  for(size_t i=0; i<isFlipv.n_elem;i++){
    if(isFlipv[i]==1){
      genomat.col(i)=arma::abs(2-genomat.col(i));
      ct++;
    }
  }
  return(ct);
}


size_t read_genotype_gz(boost::iostreams::filtering_istream &fs, const size_t Nsnps, const size_t Nind,const size_t chunksize,arma::uvec &mchroms,arma::uvec &mposs, std::vector<std::string> &refs, std::vector<std::string> &alts,arma::mat &mgenotypes){
  using namespace Rcpp;
  using namespace boost::spirit::qi;
  //  std::vector<double> genotypes;

  std::vector<int> chroms;
  std::vector<unsigned int> poss;
  std::vector<double> genotypes;

  chroms.reserve(chunksize);
  poss.reserve(chunksize);
  genotypes.reserve(chunksize*Nind);

  refs.clear();
  refs.reserve(chunksize);
  alts.clear();
  alts.reserve(chunksize);
  size_t ct=0;
  std::string line;
  while(getline(fs,line)){
    if(ct%100000==0){
      std::cout<<ct<<std::endl;
    }
    std::string::const_iterator sbeg = line.begin();
    std::string::const_iterator send = line.end();
    phrase_parse(sbeg,send,int_>>'_'>>uint_>>'_'>>as_string[+char_("ACTGN")]>>"_">>as_string[+char_("ACTGN")]>>"_">>"b37">>*double_,space,chroms,poss,refs,alts,genotypes);

    if(refs.size()!=alts.size()){
      Rcpp::Rcerr<<"refs and alts different sizes! at line "<<ct<<" ("<<refs.size()<<" "<<alts.size()<<")"<<std::endl;
      Rcpp::stop("error in read_genotype_gz");
    }
    ct++;
    if(ct==chunksize){
      break;
    }
  }

  mchroms=arma::conv_to<arma::uvec>::from(chroms);
  mposs = arma::conv_to<arma::uvec>::from(poss);
  mgenotypes = arma::mat(&genotypes[0],Nind,ct);
  arma::uvec sizes={mchroms.n_elem,mposs.n_elem,mgenotypes.n_cols,(arma::uword)refs.size(),(arma::uword)alts.size()};
  if(any(sizes!=mchroms.n_elem)){
    sizes.print();
    Rcpp::Rcerr<<"not all sizes are equal!:"<<std::endl;
    Rcpp::stop("error in read_genotype_gz");
  }
  std::cout<<"refs is of size:"<<refs.size()<<std::endl;
  std::cout<<"alts is of size:"<<alts.size()<<std::endl;
  return(ct);
}

// Rcpp::List filter_data(Rcpp::List dfl,Rcpp::DataFrame dbsnp,bool doFlip){
//
//   echrom = dfl.)
//
// }



std::vector<unsigned int> read_uint_h5(const std::string h5file, const std::string groupname, const std::string dataname){
  H5File* file;
  DataSet* dataset;
  try{
    file= new H5File( h5file.c_str(), H5F_ACC_RDONLY);
  }
  catch( FileIException error )
  {
    error.printError();
  }
  Group haplogroup = file->openGroup(groupname);
  try{
    dataset = new DataSet(haplogroup.openDataSet(dataname));
  }
  catch( DataSetIException error )
  {
    error.printError();
  }
  DataType dt = dataset->getDataType();
  DataSpace fspace =dataset->getSpace();
  hsize_t datadim[1];
  fspace.getSimpleExtentDims(datadim,NULL);
  hsize_t readrows=datadim[0];
  DataSpace memspace(1,datadim);
  std::vector<unsigned int> retdat(datadim[0]);

//  std::cout<<"Reading in int dataset dataname"<<std::endl;

  dataset->read(&retdat[0],dt,memspace,fspace);
//  std::cout<<"Last element of data is :"<<retdat[retdat.size()-1]<<std::endl;
  return(retdat);

}


std::vector<int> read_int_h5(const std::string h5file, const std::string groupname, const std::string dataname){
  H5File* file;
  DataSet* dataset;
  try{
    file= new H5File( h5file.c_str(), H5F_ACC_RDONLY);
  }
  catch( FileIException error )
  {
    error.printError();
  }
  Group haplogroup = file->openGroup(groupname);
  try{
    dataset = new DataSet(haplogroup.openDataSet(dataname));
  }
  catch( DataSetIException error )
  {
    error.printError();
  }
  DataType dt = dataset->getDataType();
  DataSpace fspace =dataset->getSpace();
  hsize_t datadim[1];
  fspace.getSimpleExtentDims(datadim,NULL);
  hsize_t readrows=datadim[0];
  DataSpace memspace(1,datadim);
  std::vector<int> retdat(datadim[0]);

  std::cout<<"Reading in int dataset dataname"<<std::endl;

  dataset->read(&retdat[0],dt,memspace,fspace);
  std::cout<<"Last element of data is :"<<retdat[retdat.size()-1]<<std::endl;
  return(retdat);

}


arma::uvec make_abs_pos(arma::uvec &chroms, arma::uvec&pos){
  arma::uvec posmap(23,arma::fill::zeros);

  std::cout<<"Creating posmap"<<std::endl;
  std::cout<<"Max of chroms is: "<<arma::max(chroms)<<std::endl;
  for(size_t i=0; i<chroms.size();i++){
    if(posmap[chroms[i]]<pos[i]){
      posmap[chroms[i]]= pos[i];
    }
  }
  arma::uvec cs=cumsum(posmap);
  std::cout<<"posmap:"<<std::endl;
  return(cs);
}

//[[Rcpp::export]]
arma::uvec make_long(arma::uvec &vchrom, arma::uvec &vpos,arma::uvec &posmap){
  if(vchrom.n_elem!=vpos.n_elem){
    Rcpp::stop("length of vectors not equal in make_long");
  }
  arma::uvec vlong =posmap.elem(vchrom-1)+vpos;
  std::cout<<"done with make_long, last element is "<<vlong[vlong.n_elem-1]<<std::endl;
  return(vlong);
}

boost::unordered_map<arma::uword,arma::uword> make_map(arma::uvec &chrom, arma::uvec &pos,arma::uvec &rsid,arma::uvec &posmap){
  arma::uvec hashes =make_long(chrom,pos,posmap);
  std::cout<<"initializing map of size:"<<hashes.size()<<std::endl;
  boost::unordered_map<arma::uword, arma::uword> retmap;
  std::cout<<"reserving space"<<std::endl;
  retmap.reserve(pos.size());
  std::cout<<"filling map space"<<std::endl;
  std::cout<<"hashes is of size:"<<hashes.size()<<" and rsid is of size"<<rsid.size()<<std::endl;
  if(hashes.size()!=rsid.size()){
    Rcpp::stop("rsid length different than hash length!");
  }
  for(size_t i=0; i<hashes.size(); i++){
    retmap[hashes[i]]=rsid[i];
  }
  std::cout<<"finished creating map"<<std::endl;
  return(retmap);
}

class dbsnpmap{
private:
  size_t size;
  std::string dbsnpfile;
public:
  boost::unordered_map<arma::uword,arma::uword> dbmap;
  arma::uvec posmap;
  dbsnpmap(std::string tdbsnpfile);
};

dbsnpmap::dbsnpmap(std::string tdbsnpfile){
  dbsnpfile=tdbsnpfile;
  std::cerr<<"this should be called first"<<std::endl;
  std::cout<<"Constructing map (really really)"<<std::endl;
  arma::uvec posv= arma::conv_to<arma::uvec>::from(read_uint_h5(dbsnpfile,"dbSNP","pos"));
  arma::uvec chromv= arma::conv_to<arma::uvec>::from(read_int_h5(dbsnpfile,"dbSNP","chrom"));
  arma::uvec rsidv = arma::conv_to<arma::uvec>::from(read_uint_h5(dbsnpfile,"dbSNP","rsid"));
  posmap =make_abs_pos(chromv,posv);
  dbmap = make_map(chromv,posv,rsidv,posmap);
}



arma::umat find_rsid(dbsnpmap &dbm, arma::uvec &chrv, arma::uvec &posv){
  std::vector<arma::uword> frsvec;
  std::vector<arma::uword> indexvec;
  frsvec.reserve(posv.size());
  indexvec.reserve(posv.size());
  arma::uvec hvec = make_long(chrv,posv,dbm.posmap);
  arma::uword vpos=0;
  boost::unordered_map<arma::uword,arma::uword>::iterator fit;
  std::cout<<"mapping rsid"<<std::endl;
  for(arma::uvec::iterator rit=hvec.begin(); rit!=hvec.end(); rit++){
    fit =dbm.dbmap.find(*rit);
    vpos = rit-hvec.begin();
    if(vpos%10000==0){
      std::cout<<"vpos:"<<vpos<<std::endl;
    }
    if(fit!=dbm.dbmap.end()){

      indexvec.push_back(vpos);
      frsvec.push_back(fit->second);
    }
  }
  arma::uvec rsv(&frsvec[0],frsvec.size());
  arma::uvec idv(&indexvec[0],indexvec.size());
  return(arma::join_horiz(rsv,idv));
}

std::vector<std::string> subset_string_vec(std::vector<std::string> &ostring, arma::uvec indices){
  std::vector<std::string> retvec(indices.n_elem);
  for(size_t i=0;i<indices.n_elem; i++){
    retvec[i]=ostring[indices[i]];
  }
  return(retvec);
}


// [[Rcpp::export]]
arma::uvec write_genotype_h5(const char* snpdatmat, size_t Nind, size_t Nsnps,size_t chunksize, const std::string h5file, bool doFlip,const std::string dbsnpfile){
  using namespace Rcpp;
  std::cout<<"Starting to map file"<<std::endl;
  boost::iostreams::mapped_file_source mapfile(snpdatmat);
  boost::iostreams::stream<boost::iostreams::mapped_file_source> textstream(mapfile);
  boost::iostreams::filtering_istream fs;
  fs.push(boost::iostreams::gzip_decompressor{});
  fs.push(textstream);
  std::string title;
  getline(fs,title);

  size_t sr=0;
  size_t scum=0;
  std::cout<<"Constructing map(this has changed)"<<std::endl;
  std::cout<<"Constructing map (really)"<<std::endl;
  dbsnpmap dbm(dbsnpfile);
  std::cout<<"Starting to read genotype data"<<std::endl;

  arma::uvec chroms;
  arma::uvec poss;
  arma::uvec retdoFlip;
  arma::mat genodat;
  std::vector<std::string> refs;
  std::vector<std::string> alts;
  arma::uvec rsidv;
  size_t count=0;

  while(scum<Nsnps){
    std::cout<<"Starting genotype data"<<std::endl;
    sr = read_genotype_gz(fs, Nsnps, Nind,chunksize,chroms,poss,refs,alts,genodat);

    scum=sr+scum;

    std::cout<<"Searching dbsnpmap"<<std::endl;

    arma::umat rsm =find_rsid(dbm,chroms,poss);
    arma::uvec indices = rsm.col(1);
    size_t retn =genodat.n_cols;
    genodat=genodat.cols(indices);
    if(doFlip){
      std::cout<<"Flipping alleles"<<std::endl;
      retdoFlip=isFlip(refs,alts,indices);
      makeFlip(retdoFlip,genodat);
    }
    std::cout<<"subsetting chroms"<<std::endl;
    chroms=chroms.elem(indices);
    poss=poss.elem(indices);

    rsidv=rsm.col(0);
    std::cout<<"Writing genotype matrix"<<std::endl;
    write_mat_h5(h5file, "SNPdata", "genotype",Nsnps, Nind,genodat);
    std::cout<<"Writing retchrom matrix"<<std::endl;
    write_int_h5(h5file,"SNPinfo","chrom",Nsnps,chroms);
    std::cout<<"Writing retpos matrix"<<std::endl;
    write_uint_h5(h5file,"SNPinfo","pos",Nsnps,poss);
    std::cout<<"Writing retrsid matrix"<<std::endl;
    write_uint_h5(h5file,"SNPinfo","rsid",Nsnps,rsidv);
    if(doFlip){
      std::cout<<"Writing doFlip"<<std::endl;
      write_int_h5(h5file,"SNPinfo","doFlip",Nsnps,retdoFlip);
    }

    count++;
  }

  // DataFrame retdf = DataFrame::create(_["chrom"]=IntegerVector(retchrom.begin(),retchrom.end()),
  //                 _["pos"]=IntegerVector(retpos.begin(),retpos.end()),
  //                 _["ref"]=StringVector(retref.begin(),retref.end()),
  //                 _["alt"]=StringVector(retalt.begin(),retalt.end()),
  //                 _["rsid"]=IntegerVector(retrsid.begin(),retrsid.end()));
  // arma::mat retgenodat(&genotypes[0],Nind,sr);
  // std::cout<<"Writing posmap"<<std::endl;
  // write_uint_h5(h5file,"SNPmap","posmap",22,dbm.posmap);
  // List retl= List::create(_["df"]=retdf,_["mat"]=retgenodat);
  return(dbm.posmap);
}

