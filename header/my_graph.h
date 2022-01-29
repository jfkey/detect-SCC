#ifndef MY_GRAPH_H_
#define MY_GRAPH_H_
#include <iostream>
#include "boost/graph/adjacency_list.hpp"
using namespace boost;

typedef adjacency_list<vecS, vecS, bidirectionalS> CG;	// CG stands for citation graph.

#endif // !MY_GRAPH_H_

