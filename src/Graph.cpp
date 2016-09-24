#include "RcppArmadillo.h"
#include "Graph.hpp"
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <utility>
#include <list>



//[[Rcpp::export]]
arma::umat find_subgraphs(const arma::umat &tedges,const size_t n_vertices,const arma::uword offset){
  size_t n_edges=tedges.n_rows;
  std::vector<Edge> edge_vec(n_edges);
  Rcpp::Rcout<<"Making pairs (number of vertices is "<<n_vertices<<")"<<std::endl;
  for(size_t i=0; i<n_edges;i++){
    edge_vec[i]=std::make_pair(tedges(i,0)-offset,tedges(i,1)-offset);
  }
  Rcpp::Rcout<<"Making graph"<<std::endl;
  boost_graph graph(edge_vec.begin(),edge_vec.end(),n_vertices,n_edges);
  Rcpp::Rcout<<"Finding connected components"<<std::endl;
  std::vector<arma::uword> components(n_vertices);
  int num_connected=connected_components(graph,&components[0]);
  Rcpp::Rcout<<"Number of connected components is: "<<num_connected<<std::endl;
  Rcpp::Rcout<<"Size of components is : "<<components.size()<<std::endl;
  arma::uvec vertices= arma::regspace<arma::uvec>(0,n_vertices-1)+offset;
  arma::uvec comp = arma::conv_to<arma::uvec>::from(components);
  return(arma::join_horiz(vertices,comp));
}


