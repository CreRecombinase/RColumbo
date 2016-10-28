#ifndef GRAPH_HPP
#define GRAPH_HPP
#include"RcppArmadillo.h"
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <utility>
#include <list>

using namespace boost;
typedef adjacency_list<vecS,vecS,undirectedS> boost_graph;

typedef graph_traits<boost_graph>::vertex_descriptor vertex_descriptor;
typedef graph_traits<boost_graph>::vertices_size_type vertices_size_type;
typedef property_map<boost_graph,vertex_index_t>::const_type vertex_index_map;
typedef std::pair<arma::uword, arma::uword> Edge;

typedef adjacency_list<vecS,vecS,undirectedS,no_property,property<edge_weight_t,float>> wGraph;

typedef graph_traits<wGraph>::vertex_descriptor wVertex;
typedef property_map<wGraph,vertex_index_t>::type wIndexMap;
typedef property_map<wGraph,edge_weight_t>::type wWeightMap;
typedef graph_traits<wGraph>::vertex_iterator wVertex_iter;
typedef graph_traits<wGraph> wGraphTraits;



arma::umat find_subgraphs(const arma::fmat LDmat, const size_t n_vertices, const arma::uword offset);

#endif

