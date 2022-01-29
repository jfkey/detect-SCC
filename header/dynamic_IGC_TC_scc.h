// This is the implementation of the algorithm of Time-Coupled-PK(IGC). PK(IGC) was first described in the  following:
// Our algorithm TC-PK is base on the algorithm of PK and IGC.
//
// [1] Pearce, D. J., & Kelly, P. H. (2007). A dynamic topological sort algorithm for 
//     directed acyclic graphs. Journal of Experimental Algorithmics (JEA), 11, 1-7.
// 
// [2] Fan, W., Hu, C., & Tian, C. (2017, May). Incremental graph computations: Doable and undoable. 
//     In Proceedings of the 2017 ACM International Conference on Management of Data (pp. 155-169).

// We divide the citation graph into three types of subgraphs, i.e., Gm, Gsi, Gr
// The strongly connected components are only held in Gm and Gsi due to the time 
// property of the citation graph. Thus, we develop alg1, alg2, and alg3 to maintain 
// the Invariants of Gm, Gsi and Gr when doing incremental cycle detection.
//  

#ifndef DYNAMIC_TC_IGC_SCC_HPP
#define DYNAMIC_TC_IGC_SCC_HPP

#include <assert.h>
#include <iomanip>

#include "disjoint_set.h"
#include "scc_para.h"
#include "my_func.h"


template <class T >
class dynamic_tc_igc_scc;


template <typename T>
class dynamic_tc_igc_com
{
private:
	T & _n2v;

public:
	dynamic_tc_igc_com(T & v) : _n2v(v) {	}

	bool operator()(unsigned int x, unsigned int y) const
	{
		return _n2v[x] > _n2v[y];								// In descending order
	}
};


class dynamic_tc_igc_com_arr {
private:
	std::vector<unsigned int> & cc;
public:
	dynamic_tc_igc_com_arr(std::vector<unsigned int> & _cc) :cc(_cc) {}

	bool operator() (unsigned int i, unsigned int j) { return (cc[i] > cc[j]); }
};



template <class T >  
class dynamic_tc_igc_scc {
private:
	typedef CG::vertex_descriptor vd_t;
	typedef CG::out_edge_iterator out_iterator;
	typedef CG::in_edge_iterator in_iterator;

public:
	CG & cg;														// citation graph.
	std::vector< T > & n2v;											// topo value of each node. 
	std::vector<unsigned int> & n2y;
	unsigned int min_year, max_year;
	std::vector<bool> visited, in_component;						// inAFF only use to indentify the affection region in scc
	std::vector<bool> &is_in_unorder;
	std::vector<int> &is_in_same_year; 
	std::vector<unsigned int> &max_ccs; 
	disjoint_set &ds;
	std::vector<pair<unsigned int, unsigned int>> & inc_edges;
	std::string self_name;
	unsigned int np;

public:
	dynamic_tc_igc_scc();
	dynamic_tc_igc_scc(CG & _cg, disjoint_set &_ds, std::vector< T > & _n2v, std::vector<unsigned int> &_max_ccs, std::vector<bool> &_is_in_unorder,
		std::vector<int> &_is_in_same_year, std::vector<unsigned int> &_n2y, unsigned int _max_year, unsigned int _min_year, std::vector<pair<unsigned int, unsigned int>> & _inc_edges, unsigned int _np);
	void add_edges();
	void add_edges(int start_index, int end_index);
	void add_edge_Gm(vd_t t, vd_t h);
	void add_edge_alg1(vd_t t, vd_t h);
	void add_edge_alg2(vd_t t, vd_t h);
	void add_edge_alg3(vd_t t, vd_t h);
	void add_edge_Gsi(vd_t t, vd_t h);
	void tc_scc_reorder(std::vector<unsigned int> & reachable, std::vector<unsigned int> &reaching);
	void fwd_dfs_Gm(vd_t n, vd_t lb, std::vector<unsigned int> &reachable, bool &flag);
	/*fwd dfs in G from Gsi & Gr. Meanwhile the starts nodes have no topological, the cycle may cause in Gsi & Gr */
	void scan_GsGr_without_topo(vd_t n, vd_t lb, std::vector<unsigned int> &reachable, bool &flag);
	/*fwd dfs in G, from Gsi & Gr to Gm, NO cycle */
	void scan_GsGr(vd_t n, std::vector<unsigned int> &reachable);
	/*fwd dfs in G, from Gm to Gsi & Gr. If there exist the cycle, the added edge(s,t) must be in this cycle. */
	//void fwd_dfs_Gu2nGu(vd_t n, vd_t lb, std::vector<unsigned int> &reachable, bool &flag);
	// back dfs in Gm. 
	void back_dfs_Gm(vd_t n, vd_t ub, std::vector<unsigned int> &reaching);
	void fwd_dfs_Gsi(vd_t n, vd_t lb, std::vector<unsigned int> &reachable, bool &flag, unsigned int & year);
	void back_dfs_Gsi(vd_t n, vd_t ub, std::vector<unsigned int> &reaching, unsigned int & year);
	
	~dynamic_tc_igc_scc();

};

#endif


