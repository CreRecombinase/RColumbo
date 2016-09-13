#include "RcppArmadillo.h"
#include "Graph.hpp"
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <utility>
#include <list>




//
// class Graph{
// private:
//   size_t n;
//   size_t e;
//   float LDcutoff;
//   boost_graph graph;
//   //  std::vector<bool> visited;
//   //  std::vector<arma::uword> tsubgraph;
//   //  std::list<arma::uvec> subgraphs;
//   size_t subsize;
// public:
//   Graph(const arma::umat &tedges,const float LDcut, const size_t n_vertices);
//   size_t NumberOfVertices()const{return(n);}
//   size_t NumberOfEdges()const{return (e);}
//   //  arma::uvec getAdj(arma::uword vertex);
//   // void DFS(arma::uword vertex);
//   // void push_subgraph();
//   // std::list<arma::uvec> output_subgraphs();
//   // void print_subgraphs();
// };


// Graph::Graph(const arma::fmat &LDmat, const float LDcut): LDcutoff(LDcut),adjmat(arma::size(LDmat),arma::fill::zeros){
//   adjmat.elem(arma::find(arma::abs(LDmat-eye(arma::size(LDmat)))>LDcutoff)).ones();
// //  std::cout<<"adjacency matrix is: "<<std::endl;
// //  adjmat.print();
//   n=adjmat.n_rows;
//   tsubgraph.reserve(n);
//   subsize=0;
//   visited.resize(n);
//   std::fill_n(visited.begin(),n,false);
// }


// arma::uvec Graph::getAdj(const arma::uword vertex){
//   return(arma::find(adjmat.col(vertex)!=0));
// }
//
// void Graph::DFS(arma::uword vertex){
//   arma::uvec neighbors=getAdj(vertex);
//   if(!visited[vertex]){
//     subsize++;
//   }
//   visited[vertex]=true;
//   for(arma::uword it=0;it!=neighbors.n_elem;it++){
//     if(!visited[neighbors[it]]){
//       tsubgraph.push_back(neighbors[it]);
//       DFS(neighbors[it]);
//     }
//   }
// }
//
// void Graph::push_subgraph(){
//   if(subsize>0){
// //    std::cout<<"Pushing subgraph:"<<std::endl;
// //    subgraph.print();
// //    std::fill_n(visited.begin(),n,false);
// //    adjmat.cols(subgraph).zeros();
// //    adjmat.rows(subgraph).zeros();
//     subgraphs.push_back(tsubgraph);
//     // std::cout<<"New adjacency matrix is:"<<std::endl;
//     // adjmat.print();
//     subsize=0;
//     tsubgraph.clear();
//   }
// }
//
// std::list<arma::uvec> Graph::output_subgraphs(){
//   for(size_t i=0; i<n;i++){
//     if(!visited[i]){
//       tsubgraph.push_back(i);
//       DFS(i);
//       push_subgraph();
//     }
//   }
//   print_subgraphs();
//   return(subgraphs);
// }
//
// void Graph::print_subgraphs(){
//   size_t iternum=0;
//   for(std::list<arma::uvec>::iterator it=subgraphs.begin(); it!=subgraphs.end();it++){
//     // std::cout<<"Subgraph "<<iternum<<std::endl;
//     iternum++;
//     for(auto jt =it->begin(); jt!= it->end(); jt++){
//       // std::cout<<*jt<<"\t";
//     }
//     // std::cout<<std::endl;
//   }
// }
//
// //[[Rcpp::export]]
// Rcpp::List find_subgraphs(const arma::fmat LDmat,const double ldcutoff){
//   arma::fmat tld=LDmat;
//    std::cout<<"Creating graph"<<std::endl;
//   Graph graph(tld,ldcutoff);
//    std::cout<<"finding subgraphs"<<std::endl;
//   return(Rcpp::wrap(graph.output_subgraphs()));
// }
//[[Rcpp::export]]
Rcpp::List find_subgraphs(const arma::umat &tedges,const float LDcut, const size_t n_vertices){

  e=tedges.n_rows;
  Edge* edge_array = new Edge[e];
  for(size_t i=0; i<e;i++){
    edge_array[i]=std::make_pair(tedges(i,0),tedges(i,1));
  }

}


