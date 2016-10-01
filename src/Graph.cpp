#include "RcppArmadillo.h"
#include "Graph.hpp"
#include "h5func.hpp"
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/sequential_vertex_coloring.hpp>
#include <utility>
#include <list>


//[[Rcpp::export]]
arma::umat find_adjmat(const arma::fmat &LDmat, const double LD_cutoff,const size_t row_offset,const size_t col_offset){
  arma::umat retmat=arma::ind2sub(arma::size(LDmat),arma::find(abs(LDmat>LD_cutoff)));
  retmat.row(0)+=row_offset;
  retmat.row(1)+=col_offset;
  return(retmat);
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


