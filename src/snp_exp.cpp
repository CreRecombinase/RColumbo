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
#include "H5IO.hpp"
#include <unordered_map>
#include <numeric>
#include <fstream>
#include <string>
#include <H5Cpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>
#include "snp_exp.hpp"
#include "h5func.hpp"

using namespace H5;


Raw_Dataset::Raw_Dataset(const std::string th5filename, const std::string tdata_groupname,
                         const std::string tdata_dataname, const std::string tanno_groupname,
                         const std::string tanno_dataname):h5filename(th5filename),data_groupname(tdata_groupname),
                         data_dataname(tdata_dataname),anno_groupname(tanno_groupname),anno_dataname(tanno_dataname),isCont(true){
  offset=0;

  H5::DataType pt=read_data_type(h5filename,data_groupname,data_dataname);
  datatype= pt.getSuper();
  P=get_rownum_h5(h5filename,data_groupname,data_dataname);
  Nind=get_arraydim(h5filename,data_groupname,data_dataname);
  data_index = arma::regspace<arma::uvec>(1,P);

}

Raw_Dataset::Raw_Dataset(const std::string th5filename, const std::string tdata_groupname,
                         const std::string tdata_dataname, const std::string tanno_groupname,
                         const std::string tanno_dataname,const arma::uvec index):h5filename(th5filename),data_groupname(tdata_groupname),
                         data_dataname(tdata_dataname),anno_groupname(tanno_groupname),anno_dataname(tanno_dataname),isCont(false),data_index(index){
  offset=0;
  P=data_index.n_elem;
  Nind=get_arraydim(h5filename,data_groupname,data_dataname);
}

arma::fmat Raw_Dataset::read_and_offset(const size_t num){
  arma::fmat retmat=read_chunk(num);
  increment_offset(num);
  return(retmat);
}

arma::fmat Raw_Dataset::read_chunk(const size_t num)const{

  using namespace H5;
  //H5DataSetPtr data_dataset = get_dataset(h5filename,data_groupname,data_dataname);
  size_t nnum = resize_chunk(num);
  // int doubtype = PredType::NATIVE_DOUBLE.getId();
  // int floattype=PredType::NATIVE_FLOAT.getId();
  if(isCont){
    // if(datatype==PredType::NATIVE_DOUBLE){
    //   return(arma::conv_to<arma::fmat>::from(read_dmat_h5(h5filename,data_groupname,data_dataname,offset,nnum)));
    // }
//    if(datatype==PredType::NATIVE_FLOAT){
      return(read_fmat_h5(h5filename,data_groupname,data_dataname,offset,nnum));
//    }
    // Rcpp::Rcerr<<"No method exists to read type: "<<datatype.fromClass()<<std::endl;
    // Rcpp::stop("Can't read from file "+h5filename+" dataset/"+data_groupname+"/"+data_dataname);
  }
  arma::size_t lind= std::min(offset+num-1,P-1);
  arma::uvec subind=data_index(arma::span(offset,lind));
  return(read_elem_direct_ind(subind));
}

arma::fmat Raw_Dataset::read_elem_ind(const arma::uvec ind) const{

  arma::uword maxind=arma::max(ind);
  if(maxind>P){
    Rcpp::stop("Can't read from file "+h5filename+" dataset/"+data_groupname+"/"+data_dataname+"Attempting to read outside range");
  }
  arma::uvec subvec=data_index.elem(ind-1);
  return(read_elem_direct_ind(subvec));
}

arma::fmat Raw_Dataset::read_elem_direct_ind(const arma::uvec ind) const{
  using namespace H5;
  // int doubtype = PredType::NATIVE_DOUBLE.getId();
  // int floattype=PredType::NATIVE_FLOAT.getId();
  // if(datatype==PredType::NATIVE_DOUBLE){
  //   return(arma::conv_to<arma::fmat>::from(read_dmat_chunk_ind(h5filename,data_groupname,data_dataname,ind)));
  // }
//  if(datatype==PredType::NATIVE_FLOAT){
    return(read_fmat_chunk_ind(h5filename,data_groupname,data_dataname,ind));
//  }
  // Rcpp::Rcerr<<"No method exists to read type: "<<datatype.fromClass()<<std::endl;
  // Rcpp::stop("Can't read from file "+h5filename+" dataset/"+data_groupname+"/"+data_dataname);
}

size_t Raw_Dataset::set_offset(const size_t off){
  if(off>P){
    Rcpp::stop("Can't set offset beyond P");
  }
  offset=off;
  return(off);
}
size_t Raw_Dataset::increment_offset(const size_t off){
  if(offset+off>P){
    Rcpp::warning("Setting offset to maximum(P-1)");
    offset=P;
  }else{
    offset=offset+off;
  }
  return(offset);
}
size_t Raw_Dataset::reset_offset(){
  offset=0;
  return(offset);
}
size_t Raw_Dataset::get_offset() const{
  return(offset);
}
arma::fmat Raw_Dataset::read_elem_name(const arma::uvec ind) const{
  // int doubtype = PredType::NATIVE_DOUBLE.getId();
  // int floattype=PredType::NATIVE_FLOAT.getId();
  // if(datatype==PredType::NATIVE_DOUBLE){
  //   return(arma::conv_to<arma::fmat>::from(read_dmat_rowname(h5filename,anno_groupname,anno_dataname,data_groupname,data_dataname,ind)));
  // }
//  if(datatype==PredType::NATIVE_FLOAT){
    return(read_fmat_rowname(h5filename,anno_groupname,anno_dataname,data_groupname,data_dataname,ind));
//  }
  // Rcpp::Rcerr<<"No method exists to read type: "<<datatype.fromClass()<<std::endl;
  // Rcpp::stop("Can't read from file "+h5filename+" dataset/"+data_groupname+"/"+data_dataname);
}

arma::uvec Raw_Dataset::get_index() const{
  return(data_index);
}

size_t Raw_Dataset::resize_chunk(const size_t rchunksize)const{
  size_t retchunk=std::min(rchunksize,P-offset);
  if(retchunk!=rchunksize){
      Rcpp::Rcerr<<"resizing chunk to "<<retchunk<<" from "<<rchunksize<<std::endl;
  }
  return(retchunk);
}

arma::uvec Raw_Dataset::get_anno_index() const{
  arma::uvec annocol = arma::conv_to<arma::uvec>::from(read_uint_h5(h5filename,anno_groupname,anno_dataname));
  if(isCont){
  return(annocol.elem(data_index-1));
  }
  return(annocol);
}

arma::fvec Raw_Dataset::read_ffeature(const std::string groupname,const std::string dataname)const {
  if(isCont){
    return(read_float_h5(h5filename,groupname,dataname,0,P));
  }
  arma::uword maxind=arma::max(data_index);
  arma::uword minind=arma::min(data_index);
  arma::uword tchunksize=maxind-minind+1;
  arma::fvec retvec= read_float_h5(h5filename,groupname,dataname,minind-1,tchunksize);
  return(retvec.elem(data_index-minind));
}


// Raw_Dataset(const std::string th5filename, const std::string tdata_groupname,
//             const std::string tdata_dataname, const std::string tanno_groupname,
//             const std::string tanno_dataname)

LD_dataset::LD_dataset(const std::string th5filename):Raw_Dataset(th5filename,"Haplotype","genotype","Legend","pos"){
  // chunksize=tchunksize;
  mapvec=read_ffeature("Legend","cummap");
}

LD_dataset::LD_dataset(const std::string th5filename,arma::uvec tindex):Raw_Dataset(th5filename,"Haplotype","genotype","Legend","pos",tindex){
  // chunksize=tchunksize;
  mapvec=read_ffeature("Legend","cummap");
}

arma::fvec LD_dataset::get_map_chunk(const size_t chunksize){
  size_t nchunk=resize_chunk(chunksize);
  arma::fvec retvec(nchunk);
  arma::uword toffset=this->get_offset();
  retvec=mapvec(arma::span(toffset,toffset+nchunk-1));
  return(retvec);
}



// class Gwas_Dataset:public Raw_Dataset{
// private:
//   Gwas_Dataset(const std::string h5filename);
//   Gwas_Dataset(const std::string h5filename, arma::uvec tindex);
//   arma::fvec mapvec;
//   arma::fvec get_map_chunk(const size_t chunksize);
// };

Gwas_Dataset::Gwas_Dataset(const std::string h5filename):Raw_Dataset(h5filename,"eQTL","genotype","Legend","pos"),LD_data(h5filename){
  betavec=read_ffeature("Legend","beta");

}
Gwas_Dataset::Gwas_Dataset(const std::string h5filename,arma::uvec tindex):Raw_Dataset(h5filename,"eQTL","genotype","Legend","pos",tindex),LD_data(h5filename,tindex){
  betavec=read_ffeature("Legend","beta");
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

  return(ct);
}


size_t read_genofile_gz(boost::iostreams::filtering_istream &fs,const size_t Nsnps,const size_t Nind, const size_t chunksize, arma::uvec &rsid, arma::uvec &chrom, arma::uvec &doFlip,arma::mat &dosagemat){

  using namespace Rcpp;
  using namespace boost::spirit::qi;
  std::vector<int> tchrom;
  std::vector<int> tnum;
  std::vector<std::string> tref;
  std::vector<std::string> talt;
  std::vector<unsigned int> trsid;
  std::vector<double> tprob;
  size_t ct=0;
  std::string line;

  tchrom.reserve(chunksize);
  tnum.reserve(chunksize);
  tref.reserve(chunksize);
  talt.reserve(chunksize);
  trsid.reserve(chunksize);
  tprob.reserve(chunksize*Nind*3);

  while(getline(fs,line)){

    std::string::const_iterator sbeg = line.begin();
    std::string::const_iterator send = line.end();
    phrase_parse(sbeg,send,int_>>"rs">>uint_>>int_>>as_string[char_("ACTGN")]>>as_string[char_("ACTGN")]>>*double_,space,tchrom,trsid,tnum,tref,talt,tprob);
    ct++;
    if(ct==chunksize){
      break;
    }
  }

  rsid=arma::conv_to<arma::uvec>::from(trsid);
  chrom = arma::conv_to<arma::uvec>::from(tchrom);
  arma::cube cprob(&tprob[0],3,Nind,ct);
  dosagemat.set_size(Nind,ct);
  arma::uvec indexv=arma::regspace<arma::uvec>(0,ct-1);
  doFlip=isFlip(tref,talt,indexv);
  arma::rowvec dvec({0,1,2});
  Rcout<<"computing dosage"<<std::endl;
  for(size_t i=0; i<ct; i++){

    dosagemat.col(i)=arma::trans(dvec*cprob.slice(i));
  }
  makeFlip(doFlip,dosagemat);

  arma::uvec sizes={rsid.n_elem,chrom.n_elem,dosagemat.n_cols};
  if(any(sizes!=rsid.n_elem)){
    sizes.print();
    Rcpp::Rcerr<<"not all sizes are equal!:"<<std::endl;
    Rcpp::stop("error in read_genofile");
  }
  return(ct);
}



size_t read_expression_gz(boost::iostreams::filtering_istream &fs, const size_t Ngenes, const size_t Nind,const size_t chunksize,arma::uvec &fgeneid,arma::uvec &sgeneid,arma::mat &expression){
  using namespace Rcpp;
  using namespace boost::spirit::qi;

  std::vector<int> tsgeneid;
  std::vector<unsigned int> tfgeneid;
  std::vector<double> texpression;


  tsgeneid.reserve(chunksize);
  tfgeneid.reserve(chunksize);
  texpression.reserve(chunksize*Nind);

  size_t ct=0;
  std::string line;
  while(getline(fs,line)){

    std::string::const_iterator sbeg = line.begin();
    std::string::const_iterator send = line.end();
    phrase_parse(sbeg,send,"ENSG">>uint_>>'.'>>int_>>*double_,space,tfgeneid,tsgeneid,texpression);
    ct++;
    if(ct==chunksize){
      break;
    }
  }

  fgeneid=arma::conv_to<arma::uvec>::from(tfgeneid);
  sgeneid = arma::conv_to<arma::uvec>::from(tsgeneid);
  expression = arma::mat(&texpression[0],Nind,ct);
  arma::uvec sizes={fgeneid.n_elem,sgeneid.n_elem,expression.n_cols};
  if(any(sizes!=fgeneid.n_elem)){
    sizes.print();
    Rcpp::Rcerr<<"not all sizes are equal!:"<<std::endl;
    Rcpp::stop("error in read_expression_gz");
  }
  return(ct);
}


size_t read_gene_position_gz(std::string gene_posfile,arma::umat &posmat){
  using namespace Rcpp;
  using namespace boost::spirit::qi;

  boost::iostreams::mapped_file_source mapfile(gene_posfile);
  boost::iostreams::stream<boost::iostreams::mapped_file_source> textstream(mapfile);
  boost::iostreams::filtering_istream fs;
  fs.push(boost::iostreams::gzip_decompressor{});
  fs.push(textstream);
  std::string title;
  getline(fs,title);


  std::vector<int> tsgeneid;
  std::vector<int> tchrom;
  std::vector<unsigned int> tfgeneid;
  std::vector<unsigned int> tstart;
  std::vector<unsigned int> tstop;

  size_t ct=0;
  std::string line;
  std::string chtest;
  symbols<char,int> sym;
  for(int i=1; i<23; i++){
    sym.add(std::to_string(i),i);
  }
  sym.add("X",23)
    ("Y",24)
    ("MT",25);

  while(getline(fs,line)){

    std::string::const_iterator sbeg = line.begin();
    std::string::const_iterator send = line.end();
    phrase_parse(sbeg,send,"ENSG">>uint_>>'.'>>int_>>sym>>uint_>>uint_,space,tfgeneid,tsgeneid,tchrom,tstart,tstop);
    ct++;
  }

  arma::uvec sizes={tfgeneid.size(),tsgeneid.size(),tchrom.size(),tstart.size(),tstop.size()};
  if(any(sizes!=tfgeneid.size())){
    sizes.print();
    Rcpp::Rcerr<<"not all sizes are equal!:"<<std::endl;
    Rcpp::stop("error in read_gene_position_gz");
  }
  posmat.set_size(tfgeneid.size(),5);
  posmat.col(0)=arma::conv_to<arma::uvec>::from(tfgeneid);
  posmat.col(1)=arma::conv_to<arma::uvec>::from(tsgeneid);
  posmat.col(2)=arma::conv_to<arma::uvec>::from(tchrom);
  posmat.col(3)=arma::conv_to<arma::uvec>::from(tstart);
  posmat.col(4) = arma::conv_to<arma::uvec>::from(tstop);

  return(ct);
}



genepos_map::genepos_map(std::string tdbsnpfile):annomat(){
  gene_pos_file=tdbsnpfile;
  size= read_gene_position_gz(gene_pos_file,annomat);
  gp_map.reserve(annomat.n_rows);
  for(arma::uword i=0; i<annomat.n_rows;i++){
    gp_map[annomat(i,0)]=i;
  }
}

arma::umat find_gene_pos(genepos_map &gpm,arma::uvec &fgeneid){
  arma::uvec indices(fgeneid.n_rows);
  std::unordered_map<arma::uword,arma::uword>::iterator fit;
  for(size_t i=0; i<fgeneid.n_elem;i++){
    fit =gpm.gp_map.find(fgeneid[i]);
    if(fit!=gpm.gp_map.end()){
      indices[i]=fit->second;
    }
    else{
      Rcpp::Rcerr<<"geneid:"<<fgeneid<<" not found!"<<std::endl;
      indices[i]=-1;
    }
  }
  arma::umat retmat= gpm.annomat.rows(indices).t();
  return(retmat);
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


//
// std::unordered_map<arma::uword, arma::uword> make_map(arma::uvec queryvec,arma::uvec targetvec){
//   std::unordered_map<arma::uword, arma::uword> retmap;
//   std::cout<<"reserving space"<<std::endl;
//   retmap.reserve(pos.size());
// }


std::unordered_map<arma::uword,arma::uword> make_map(arma::uvec &chrom, arma::uvec &pos,arma::uvec &rsid,arma::uvec &posmap){
  arma::uvec hashes =make_long(chrom,pos,posmap);
  std::cout<<"initializing map of size:"<<hashes.size()<<std::endl;
  std::unordered_map<arma::uword, arma::uword> retmap;
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



dbsnpmap::dbsnpmap(std::string tdbsnpfile){
  dbsnpfile=tdbsnpfile;
//  std::cerr<<"this should be called first"<<std::endl;
//  std::cout<<"Constructing map (really really)"<<std::endl;
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
  std::unordered_map<arma::uword,arma::uword>::iterator fit;
  std::cout<<"mapping rsid"<<std::endl;
  for(arma::uvec::iterator rit=hvec.begin(); rit!=hvec.end(); rit++){
    fit =dbm.dbmap.find(*rit);
    vpos = rit-hvec.begin();
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
//[[Rcpp::export]]
arma::rowvec colssd(const arma::mat &data){
  return(arma::stddev(data));
}



// [[Rcpp::export]]
arma::uvec write_genotype_h5(const char* snpdatmat, size_t Nind, size_t Nsnps,size_t chunksize, const std::string h5file, bool doFlip,const std::string dbsnpfile,const unsigned int deflate_level){
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
    write_mat_h5(h5file, "SNPdata", "genotype",Nsnps, Nind,genodat,deflate_level);
    std::cout<<"Writing retchrom matrix"<<std::endl;
    write_int_h5(h5file,"SNPinfo","chrom",Nsnps,chroms,deflate_level);
    std::cout<<"Writing retpos matrix"<<std::endl;
    write_uint_h5(h5file,"SNPinfo","pos",Nsnps,poss,deflate_level);
    std::cout<<"Writing retrsid matrix"<<std::endl;
    write_uint_h5(h5file,"SNPinfo","rsid",Nsnps,rsidv,deflate_level);
    if(doFlip){
      std::cout<<"Writing doFlip"<<std::endl;
      write_int_h5(h5file,"SNPinfo","doFlip",Nsnps,retdoFlip,deflate_level);
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





//[[Rcpp::export]]
size_t write_genofile_h5(const std::string genofile,const std::string dbsnpfile,size_t Nind,size_t Nsnps,size_t chunksize,const std::string h5file, const unsigned int deflate_level){

  using namespace Rcpp;
  std::cout<<"Starting to map file"<<std::endl;
  boost::iostreams::mapped_file_source mapfile(genofile);
  boost::iostreams::stream<boost::iostreams::mapped_file_source> textstream(mapfile);
  boost::iostreams::filtering_istream fs;
  fs.push(boost::iostreams::gzip_decompressor{});
  fs.push(textstream);
  size_t sr=0;
  size_t scum=0;
//  std::cout<<"Constructing map"<<std::endl;
//  dbsnpmap dbm(dbsnpfile);
  std::cout<<"Starting to read genotype data"<<std::endl;

  arma::uvec chrom;
  arma::uvec rsid;
  arma::uvec doFlip;
  arma::mat dosagemat;
  size_t count=0;

  while(scum<Nsnps){
    std::cout<<"Starting genofile data"<<std::endl;
    sr = read_genofile_gz(fs, Nsnps, Nind,chunksize,rsid,chrom,doFlip,dosagemat);
    scum=sr+scum;

    // std::cout<<"Writing genotype matrix"<<std::endl;
    write_mat_h5(h5file, "Genotype", "genotype",Nsnps, Nind,dosagemat,deflate_level);

    // // std::cout<<"Writing chrom vec"<<std::endl;
    // write_int_h5(h5file,"Legend","chrom",Nsnps,chrom,deflate_level);
    // // std::cout<<"Writing rsid vec"<<std::endl;
    // write_uint_h5(h5file,"Legend","rsid",Nsnps,rsid,deflate_level);
    // // std::cout<<"Writing doFlip vec"<<std::endl;
     write_uint_h5(h5file,"Legend","doFlip",Nsnps,doFlip,deflate_level);
    count++;
  }
}

// [[Rcpp::export]]
size_t write_expression_h5(const char* expdatmat, size_t Nind, size_t Ngenes,size_t chunksize, const std::string h5file,const std::string gene_posfile,const unsigned int deflate_level){
  using namespace Rcpp;
  std::cout<<"Starting to map file"<<std::endl;
  boost::iostreams::mapped_file_source mapfile(expdatmat);
  boost::iostreams::stream<boost::iostreams::mapped_file_source> textstream(mapfile);
  boost::iostreams::filtering_istream fs;
  fs.push(boost::iostreams::gzip_decompressor{});
  fs.push(textstream);
  std::string title;
  getline(fs,title);

  size_t sr=0;
  size_t scum=0;
  genepos_map gpm(gene_posfile);
  std::cout<<"Starting to read genotype data"<<std::endl;

  arma::uvec fgeneid;
  arma::uvec sgeneid;
  arma::mat expmat;
  size_t count=0;

  while(scum<Ngenes){
    std::cout<<"Starting genotype data"<<std::endl;
    sr = read_expression_gz(fs, Ngenes, Nind,chunksize,fgeneid,sgeneid,expmat);
    scum=sr+scum;

    std::cout<<"Searching dbsnpmap"<<std::endl;

    arma::umat annomat =find_gene_pos(gpm,fgeneid);

    std::cout<<"Writing expression matrix"<<std::endl;
    write_mat_h5(h5file, "EXPdata", "expression",Ngenes, Nind,expmat,deflate_level);
    std::cout<<"Writing gene_annotation matrix"<<std::endl;
    write_umat_h5(h5file,"EXPinfo","annomat",Ngenes,5,annomat,deflate_level);
    count++;
  }

  return(count);
}




















