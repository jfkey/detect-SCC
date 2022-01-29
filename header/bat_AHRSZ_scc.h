// This is the implementation of the algorithm of bat-AHRSZ. It was first described in the  following:
//
// [1] B. Alpern, R. Hoover, B. K. Rosen, P. F. Sweeney and F. K. 
//     Zadeck, "Incremental Evaluation of Computational Circuits", 
//     In Proc. of the First Annual ACM-SIAM Symposium on Discrete 
//     Algorithms, pages 32-42, 1990.
// 
// 

#ifndef DYNAMIC_BAT_AHRSZ_SCC_HPP
#define DYNAMIC_BAT_AHRSZ_SCC_HPP

#include <queue>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <random>
#include "boost/property_map/property_map.hpp"
#include "boost/graph/topological_sort.hpp"
#include "boost/lexical_cast.hpp"


#include "ordered_slist.h"
#include "ordered_slist2.h"
#include "disjoint_set.h"
#include "my_func.h"
#include "my_less.h"
#include "my_greater.h"

#if __GNUC__ >= 3
#include <ext/functional>
using __gnu_cxx::select2nd;
using __gnu_cxx::identity;
#else
#include "my_select.h"
#include "boost/compute/functional/identity.hpp"
#endif


template<typename T>
class batch_AHRSZ_scc {
private:
	typedef CG::vertex_descriptor vd_t;
	typedef CG::out_edge_iterator out_iterator;
	typedef CG::in_edge_iterator in_iterator;
	typedef std::priority_queue<vd_t, std::vector<vd_t>, ahrsz_priority_comp<std::vector< ahrsz_priority_value <T> >, std::less<ahrsz_priority_value<T> > > > max_priority_queue;
	typedef std::priority_queue<vd_t, std::vector<vd_t>, ahrsz_priority_comp<std::vector< ahrsz_priority_value <T> >, std::greater<ahrsz_priority_value<T> > > > min_priority_queue;

public:
	CG & cg;														// citation graph.
	std::vector< ahrsz_priority_value <T> > & n2v;					// topo value of each node. 
	T & _pspace;													// T: ordered_slist or ordered_slist2,  the priority space.
	std::vector<ahrsz_ext_priority_value<T> > _flooring;			// temporary storage needed by the reassignment stage
	std::vector<bool> _visited;
	std::vector<bool> _inK;
	std::vector<bool> _inB, _inF, in_component;						// use inB and inF to identify the cycle.
	std::vector<unsigned int> _indegree;
	//std::vector<bool> prev_visited, prev_inB, prev_inF;				// record the previously inserted edges which could arrive, to identify the cycle. 

	disjoint_set &ds;
	string self_name;
	std::vector<pair<unsigned int, unsigned int>> & inc_edges;
	unsigned int bs;

public:
	batch_AHRSZ_scc();
	batch_AHRSZ_scc(CG & cg, std::vector< ahrsz_priority_value <T> > & n2v,  T &pspace, disjoint_set & ds, std::vector<pair<unsigned int, unsigned int>> & inc_edges, const unsigned int np, unsigned int _bs);

	void add_edges();
	void add_edges(int start_index, int end_index);
	void add_edges_bat(std::vector<std::pair<unsigned int, unsigned int>> & invalid_edges);
	void add_edge_single(vd_t t, vd_t h, bool &has_cycle, std::vector<vd_t> &K, std::vector<unsigned int> & cycle_nodes);
	
	void discovery_bat(vd_t tail, vd_t head, std::vector<vd_t> &K, bool & has_cycle, std::vector<unsigned int> &cycle_nodes);
	void reassignment(std::vector<vd_t> &K);
	void reassignment_mini_priority(std::vector<vd_t> &K);
	void reassignment_tail(bool &has_cycle, std::vector<vd_t> &K, std::vector<unsigned int> & cycle_nodes);
	void refresh_invalid_edges(std::vector<std::pair<unsigned int, unsigned int>> &invalid_edges, int idx, std::vector<bool> &is_invalid);
	

	void grow_forw2(vd_t node, max_priority_queue & max_pq, bool & has_cycle, std::vector<vd_t> &K);
	void grow_back2(vd_t node, min_priority_queue & min_pq, bool & has_cycle, std::vector<vd_t> &K);
	void find_cycle(vd_t head);
	void compute_flooring(vd_t n, std::vector<vd_t> &rto);
	ahrsz_ext_priority_value<T> compute_ceiling(vd_t v);
	ahrsz_ext_priority_value<T> compute_priority(ahrsz_ext_priority_value<T> floor, ahrsz_ext_priority_value<T> ceiling);
	~batch_AHRSZ_scc();
	
	
	
};

#endif




