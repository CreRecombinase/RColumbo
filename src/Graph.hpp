#ifndef GRAPH_HPP
#define GRAPH_HPP
#include"RcppArmadillo.h"
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <utility>
#include <list>

using namespace boost;
typedef adjacency_list<vecS,vecS,undirectedS> boost_graph;
typedef std::pair<arma::uword, arma::uword> Edge;

class Graph{
private:
  size_t n;
  size_t e;
  float LDcutoff;
  boost_graph graph;
  //  std::vector<bool> visited;
  //  std::vector<arma::uword> tsubgraph;
  //  std::list<arma::uvec> subgraphs;
  size_t subsize;
public:
  Graph(const arma::umat &tedges,const float LDcut, const size_t n_vertices);
  size_t NumberOfVertices()const{return(n);}
  size_t NumberOfEdges()const{return (e);}
  //  arma::uvec getAdj(arma::uword vertex);
  // void DFS(arma::uword vertex);
  // void push_subgraph();
  // std::list<arma::uvec> output_subgraphs();
  // void print_subgraphs();
};

Rcpp::List find_subgraphs(const arma::fmat &LDmat,const double ldcutoff);
#endif

