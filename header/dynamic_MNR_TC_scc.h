// This is the implementation of the algorithm of  Time-Coupled-MNR. MNR was first described in the  following:
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
// maps each index in the order to the corresponding vertex. In fact, in
// order to maintain the scc, the $ord$ is not a contiguous ordering. 
// 
// We divide the citation graph into three types of subgraphs, i.e., Gm, Gsi, Gr
// The strongly connected components are only held in Gm and Gsi due to the time 
// property of the citation graph. Thus, we develop alg1, alg2, and alg3 to maintain 
// the Invariants of Gm, Gsi and Gr when doing incremental cycle detection.

#ifndef DYNAMIC_TC_MNR_SCC_H_
#define DYNAMIC_TC_MNR_SCC_H_

#include "boost/property_map/property_map.hpp"
#include "disjoint_set.h"
#include "my_func.h"

template<class T>
class dynamic_tc_mnr_scc {
private:
	typedef CG::vertex_descriptor vd_t;
	typedef CG::out_edge_iterator out_iterator;
	typedef CG::in_edge_iterator in_iterator;
public:
	CG & cg;														// citation graph.
	std::vector< T > & n2v;											// topo value of each node. 
	std::vector< std::vector<vd_t> > v2ns;
	std::vector<bool> visited, in_component;
	disjoint_set &ds;
	std::vector<bool> &is_in_unorder;
	std::vector<int> &is_in_same_year;
	std::vector<unsigned int> &max_ccs;
	std::vector<unsigned int> & n2y;
	std::vector<pair<unsigned int, unsigned int>> & inc_edges;
	std::string self_name;
	unsigned int max_year, min_year;
	unsigned int np;


public:
	dynamic_tc_mnr_scc();
	dynamic_tc_mnr_scc(CG & _cg, disjoint_set &_ds, std::vector< T > & _n2v, std::vector<unsigned int> &_max_ccs, std::vector<bool> &_is_in_unorder,
		std::vector<int> &_is_in_same_year, std::vector<unsigned int> &_n2y, unsigned int _max_year, unsigned int _min_year, std::vector<pair<unsigned int, unsigned int>> & _inc_edges, unsigned int _np);

	void add_edges();
	void add_edges(int start_index, int end_index);
	void add_edge_Gm(vd_t t, vd_t h);
	void add_edge_alg1(vd_t t, vd_t h);
	void add_edge_alg2(vd_t t, vd_t h);
	void add_edge_alg3(vd_t t, vd_t h);
	void add_edge_Gsi(vd_t t, vd_t h);

	void scan_GsGr(vd_t n, std::vector<unsigned int> &reachable);
	void scan_GsGr_without_topo(vd_t n, vd_t lb, std::vector<unsigned int> &reachable, bool &flag);
	void mnr_dfs_Gm(vd_t h, unsigned int lb, unsigned int ub, std::vector<unsigned int>& reachable, bool & flag);
	void mnr_dfs_Gsi(vd_t h, unsigned int lb, unsigned int ub, std::vector<unsigned int>& reachable, bool &flag, unsigned int year);
	void mnr_shift_Gm(unsigned int lb, unsigned int ub);
	void mnr_shift_Gsi(unsigned int lb, unsigned int ub, unsigned int year);
	void mnr_shift_cycle_Gm(unsigned int lb, unsigned int ub);
	void mnr_shift_cycle_Gsi(unsigned int lb, unsigned int ub, unsigned int year);
};
#endif
