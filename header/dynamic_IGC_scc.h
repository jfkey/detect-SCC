// This is the implementation of the algorithm of PK and IGC+. It was first described in the  following:
//
// [1] Pearce, D. J., & Kelly, P. H. (2007). A dynamic topological sort algorithm for 
//     directed acyclic graphs. Journal of Experimental Algorithmics (JEA), 11, 1-7.
// 
// [2] Fan, W., Hu, C., & Tian, C. (2017, May). Incremental graph computations: Doable and undoable. 
//     In Proceedings of the 2017 ACM International Conference on Management of Data (pp. 155-169).
//	 
// 

#ifndef DYNAMIC_IGC_SCC_H_
#define DYNAMIC_IGC_SCC_H_
#include "boost/graph/adjacency_list.hpp"
#include "disjoint_set.h"
#include "my_graph.h"
#include "disjoint_set.h"
#include "my_func.h"
#include <assert.h>

template <class T>			// dynamic_igc_scc<unsigned int>
class dynamic_igc_scc {
private:
	typedef CG::vertex_descriptor vd_t;
	typedef CG::out_edge_iterator out_iterator;
	typedef CG::in_edge_iterator in_iterator;
public:
	CG & cg;														// citation graph.
	std::vector< T > & n2v;											// topo value of each node. 
	std::vector<bool> visited, in_component;						// inAFF only use to indentify the affection region in scc
	disjoint_set &ds;
	std::vector<pair<unsigned int, unsigned int>> & inc_edges;
	std::string self_name;
	unsigned int np;

public:
	dynamic_igc_scc();
	dynamic_igc_scc(CG &cg, disjoint_set &_ds, std::vector<T> &_n2v, std::vector<pair<unsigned int, unsigned int>> & _inc_edges, unsigned int _np);
	void add_edges();
	void add_edges(int start_index, int end_index);
	void do_igc(vd_t t, vd_t h);
	void igc_scc_fwd_dfs(vd_t n, vd_t lb,  vd_t fo,  std::vector<unsigned int>&reachable, bool &flag);
	void igc_scc_back_dfs(vd_t n, vd_t ub, std::vector<unsigned int>&reaching);
	void igc_scc_reorder(std::vector<unsigned int> &reachable, std::vector<unsigned int> &reaching );
	~dynamic_igc_scc();
};

template <class T>
class dynamic_igc_scc_com
{
private:
	T & _n2v;

public:
	dynamic_igc_scc_com(T & v) : _n2v(v) {	}

	bool operator()(unsigned int x, unsigned int y) const
	{
		return _n2v[x] > _n2v[y];								// In descending order
	}
};


#endif // !dynamic_igc_scc