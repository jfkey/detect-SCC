// "On-Line Graph Algorithms for Incremental Compilation", 
// Alberto Marchetti-Spaccamela, Umberto Nanni and Hans 
// Rohnert. Information Processing Letters, 1996.

// This is the implementation of the algorithm of MNR. It was first described in the  following:
//
// [1] Marchetti-Spaccamela, Alberto, Umberto Nanni, and Hans Rohnert. "On-line graph algorithms for 
//     incremental compilation." International Workshop on Graph-Theoretic Concepts in Computer Science. 
//     Springer, Berlin, Heidelberg, 1993.
// [2] Marchetti-Spaccamela, Alberto, Umberto Nanni, and Hans Rohnert. 1996. ¡°Maintaining a Topological 
//     Order under Edge Insertions.¡± Information Processing Letters 59 (1): 53¨C58.
// 
// The algorithm of the MNR implements $ord$ as a total and contiguous 
// ordering of vertices using an array of size |V|. Then, we need two 
// auxiliary arrays to implement the MNR. The first array n2i maps each 
// vertices to a unique interger in {1, 2, ..., |V|}. The second array 
// maps each index in the order to the corresponding vertex. 
// 

#ifndef DYNAMIC_MNR_SCC_H_
#define DYNAMIC_MNR_SCC_H_

#include "boost/property_map/property_map.hpp"
#include "disjoint_set.h"
#include "my_func.h"

template<class T>
class dynamic_mnr_scc{
private:
	typedef CG::vertex_descriptor vd_t;
	typedef CG::out_edge_iterator out_iterator;
	typedef CG::in_edge_iterator in_iterator;
public:
	CG & cg;														// citation graph.
	std::vector< T > & n2v;											// topo value of each node. 
	std::vector<vd_t> v2n; 
	std::vector<bool> visited, in_component;						
	disjoint_set &ds;
	std::vector<pair<unsigned int, unsigned int>> & inc_edges;
	std::string self_name;
	unsigned int np;

public:
	dynamic_mnr_scc();
	dynamic_mnr_scc(CG &_cg, disjoint_set &_ds, std::vector<T> &_n2v, std::vector<pair<unsigned int, unsigned int>> & _inc_edges, unsigned int _np);
	void add_edges();
	void add_edges(int start_index, int end_index);
	void add_edge(vd_t t, vd_t h);
	void mnr_dfs(vd_t h, unsigned int lb, unsigned int ub, std::vector<unsigned int>& reachable, bool &flag);
	void mnr_shift(unsigned int lb, unsigned int ub);
	void mnr_shift_cycle(unsigned int lb, unsigned int ub);

};

#endif



