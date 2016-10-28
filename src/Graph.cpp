#include "RcppArmadillo.h"
#include "Graph.hpp"
#include "h5func.hpp"
#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/sequential_vertex_coloring.hpp>
#include <utility>
#include <list>


//[[Rcpp::export]]
arma::umat find_adjmat(const arma::fmat &LDmat, const double LD_cutoff,const size_t row_offset,const size_t col_offset){
  arma::umat retmat=arma::ind2sub(arma::size(LDmat),arma::find(abs(LDmat)>LD_cutoff));
  retmat.row(0)+=row_offset;
  retmat.row(1)+=col_offset;

  return(retmat);
}

//[[Rcpp::export]]
Rcpp::DataFrame find_adjdf(const arma::fmat &LDmat,const double LD_cutoff,const size_t row_offset,const size_t col_offset){
  using namespace Rcpp;
  arma::uvec kelem =arma::find(abs(LDmat)>LD_cutoff);
  arma::umat retmat=arma::ind2sub(arma::size(LDmat),kelem);
  arma::fvec ldvec= LDmat.elem(kelem);
  retmat.row(0)+=row_offset;
  retmat.row(1)+=col_offset;
  arma::urowvec rowid= retmat.row(0);
  arma::urowvec colid=retmat.row(1);
  return(DataFrame::create(_["rowid"]= IntegerVector(rowid.begin(),rowid.end()),
                           _["colid"]= IntegerVector(colid.begin(),colid.end()),
                           _["LD"]=NumericVector(ldvec.begin(),ldvec.end())));
}

//[[Rcpp::export]]
arma::umat find_adjmat_chunk(const std::string h5file, const double LD_cutoff, const size_t chunksize){
  size_t dimension= get_rownum_h5(h5file,"LD_mat","LD");
  size_t numchunk= std::ceil((double)dimension/(double)chunksize);
  size_t rowoffset=0;
  size_t coloffset=0;
  arma::umat retmat;
  arma::fmat tLDmat;
  while(rowoffset<dimension){
    coloffset=0;
    while(coloffset<dimension){
      tLDmat=read_2dfmat_h5(h5file,"LD_mat","LD",rowoffset,coloffset,chunksize,chunksize);
      std::cout<<"Finding connected elements"<<std::endl;
      arma::umat tadjmat= find_adjmat(tLDmat,LD_cutoff,rowoffset,coloffset);
      coloffset+=tLDmat.n_cols-1;
      if(tadjmat.n_cols>0){
        if(retmat.n_cols==0){
          retmat=tadjmat;
        }else{
          std::cout<<"Joining matrices "<<std::endl;
          retmat=arma::join_horiz(retmat,tadjmat);
        }
      }else{
        break;
      }
    }
    rowoffset+=tLDmat.n_rows-1;
  }
  return(retmat);
}



//[[Rcpp::export]]
arma::umat sequential_coloring(const arma::umat &tedges){
  size_t n_edges=tedges.n_cols;
  arma::uvec vertices = arma::unique(arma::vectorise(tedges));
  arma::uword offset=arma::min(vertices);
  size_t n_vertices = vertices.n_elem;
  std::vector<Edge> edge_vec(n_edges);
  Rcpp::Rcout<<"Making pairs (number of vertices is "<<n_vertices<<")"<<std::endl;
  for(size_t i=0; i<n_edges;i++){
    edge_vec[i]=std::make_pair(tedges(0,i)-offset,tedges(1,i)-offset);
  }
  Rcpp::Rcout<<"Making graph"<<std::endl;
  boost_graph graph(edge_vec.begin(),edge_vec.end(),n_vertices,n_edges);
  Rcpp::Rcout<<"Coloring graph (sequentially)"<<std::endl;
  std::vector<vertices_size_type> color_vec(num_vertices(graph));
  iterator_property_map<vertices_size_type*,vertex_index_map> color(&color_vec.front(),get(vertex_index,graph));
  vertices_size_type num_colors= sequential_vertex_coloring(graph,color);
  Rcpp::Rcout<<"Number of colors is: "<<num_colors<<std::endl;
  arma::uvec coloring = arma::conv_to<arma::uvec>::from(color_vec);
  return(arma::join_horiz(vertices,coloring));
}

//can be easily parallelized up to the numbere of connected components
//(higher if you don't mind reading the same data twice)

//
// arma::fmat adj_lik(const arma::fvec rvec, const arma::uvec rowv, const arma::uvec colv,
//                    const arma::fvec betahat,const arma::fvec se, const float sa, const float pi){
//   arma::uvec usnps = arma::unique(arma::vectorise(join_horiz(rowv,colv)));
//   size_t n_vertices= usnps.n_elem;
//   size_t n_edges =rowv.n_elem;
//
//   if(rowv.n_elem!=colv.n_elem){
//     Rcpp::stop("Different number of elements in rowv and colv");
//   }
//
//   if(n_vertices!=betahat.n_elem){
//     Rcpp::stop("Number of unique SNPs in LD is not equal to extent of betahat!");
//   }
//
//   if(betahat.n_elem!=se.n_elem){
//     Rcpp::stop("Extent of betahat not equal to extent of se in adj_lik!");
//   }
//   arma::uword offset=arma::min(usnps);
//   std::vector<Edge> edge_vec(n_edges);
//   Rcpp::Rcout<<"Making pairs (number of vertices is "<<n_vertices<<")"<<std::endl;
//   for(size_t i=0; i<n_edges;i++){
//     edge_vec[i]=std::make_pair(rowv(i)-offset,colv(i)-offset);
//   }
//   std::vector<float> weights= arma::conv_to<std::vector<float>>::from(rvec);
//   Rcpp::Rcout<<"Making graph"<<std::endl;
//   wGraph graph(edge_vec.begin(),edge_vec.end(),&weights[0],n_vertices,n_edges);
//   wWeightMap weightindex= get(edge_weight,graph);
//   wIndexMap index = get(vertex_index,graph);
//   wVertex_iter vp;
//   wGraphTraits::out_edge_iterator edge_i, edge_end;
//   float likdat=0;
//   arma::fmat retmat(n_vertices,2,arma::fill::zeros);
//   std::cout<<"Computing Likelihood"<<std::endl;
// //  for(vp=vertices(graph).first;vp!=vertices(graph).second;++vp){
//   vp=vertices(graph).first;
// //    std::cout<<"*vp is:"<<*vp<<" index[*vp] is"<<index[*vp]<<std::endl;
//     float betah= betahat[index[*vp]];
//     float sih= se[index[*vp]];
//     int n=0;
//     int on=out_degree(*vp,graph);
//     likdat=0;
//     bool foundself=false;
//     for(boost::tie(edge_i,edge_end)=out_edges(*vp,graph);edge_i!=edge_end;++edge_i){
//       std::cout<<"Edge from"<<index[source(*edge_i,graph)]<<" to "<<index[target(*edge_i,graph)]<<std::endl;
//       if(index[source(*edge_i,graph)]==index[source(*edge_i,graph)])
//       if(!foundself&(index[source(*edge_i,graph)]))
//       float osih=se[index[target(*edge_i,graph)]];
//       float r =weightindex[*edge_i];
//       likdat+=R::dnorm4(betah,0,sqrt(sih*sih*(1+((r*r)/(osih*osih))*(sa*sa))),0);
//       n++;
//     }
//
//     retmat(index[*vp],0)=likdat/n;
//     retmat(index[*vp],1)=n;
// }
//   return(retmat);
// }

