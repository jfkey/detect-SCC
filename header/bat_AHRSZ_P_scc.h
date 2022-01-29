// We propose an AHRSZ+ for incremental strongly connected component detection based on AHRSZ. 
//
// [1] B. Alpern, R. Hoover, B. K. Rosen, P. F. Sweeney and F. K. 
//     Zadeck, "Incremental Evaluation of Computational Circuits", 
//     In Proc. of the First Annual ACM-SIAM Symposium on Discrete 
//     Algorithms, pages 32-42, 1990.
// 
// The difference between our AHRSZ+ and AHRSZ is as follows:
// 1. we give the idea of how to detect the cycle and the scc maintance component is added in the AHRSZ 
// 2. the scc maintance could be done in O(C) complexity
// 3. the forward edges and backward edges are initialized to the indegree and outdegree of the node in the discovery phase, respectively
// 4. we simplify the reassignment phase by reducing the computation of the ceiling of the node but with higher efficiency 
// 5. batch-inc with a smaller space consumption.
// 6. we use topological sort while AHRSZ using pseudo topological sort 
// 
// 


#ifndef DYNAMIC_BAT_AHRSZ_P_SCC_HPP
#define DYNAMIC_BAT_AHRSZ_P_SCC_HPP

#include <queue>
#include <vector>
#include <stdexcept>
#include <cassert>
#include "boost/property_map/property_map.hpp"
#include "boost/graph/topological_sort.hpp"
#include "boost/lexical_cast.hpp"

#include "ordered_slist.h"
#include "ordered_slist2.h"
#include "disjoint_set.h"
#include "my_func.h"


#if __GNUC__ >= 3
#include <ext/functional>
using __gnu_cxx::select2nd;
using __gnu_cxx::identity;
#else
#include "my_select.h"
#include "boost/compute/functional/identity.hpp"
#endif


template<typename T>
class batch_AHRSZ_P_scc {
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
	std::vector<bool> prev_visited, prev_inB, prev_inF;				// record the previously inserted edges which could arrive, to identify the cycle. 

	disjoint_set &ds;
	string self_name;
	std::vector<pair<unsigned int, unsigned int>> & inc_edges;
	unsigned int bs;

	// the following four are for valid-edge aware bat alg.
	std::vector<std::vector<int>> start_invalid_index;				// store the index of the invalid edge's start node 
	std::vector<std::vector<int>> end_invalid_index;				// store the index of the invalid edge's end node 
	std::vector<bool> inv_rem_vis;									// invalid edges remain to be visited. 
	std::vector<bool> inv_has_vis;										// invalid edges have visited


public:
	batch_AHRSZ_P_scc();
	batch_AHRSZ_P_scc(CG & cg, std::vector< ahrsz_priority_value <T> > & n2v, T &pspace, disjoint_set & ds, std::vector<pair<unsigned int, unsigned int>> & inc_edges, const unsigned int np, unsigned int _bs);

	void add_edges();
	void add_edges(int start_index, int end_index);
	void add_edges_valid_edge_aware();
	void valid_edge_aware(std::vector<std::pair<unsigned int, unsigned int>>& invalid_edges, std::vector<std::pair<unsigned int, unsigned int>>& new_invalid_edges, std::vector<unsigned int>& unmarks, std::vector<unsigned int>unmarkt);
	void add_edges_bat(std::vector<std::pair<unsigned int, unsigned int>> & invalid_edges);

	void add_edge_single(vd_t t, vd_t h, bool &has_cycle, std::vector<vd_t> &K, std::vector<unsigned int> & cycle_nodes);

	void discovery_bat(vd_t tail, vd_t head, std::vector<vd_t> &K, bool & has_cycle, std::vector<unsigned int> &cycle_nodes);
	void reassignment(std::vector<vd_t> &K);
	void reassignment_tail(bool &has_cycle, std::vector<vd_t> &K, std::vector<unsigned int> & cycle_nodes);
	void refresh_invalid_edges(std::vector<std::pair<unsigned int, unsigned int>> &invalid_edges, int idx, std::vector<bool> &is_invalid);
	void refresh_invalid_edges(std::vector<std::pair<unsigned int, unsigned int>> &invalid_edges, int idx, std::vector<bool> &is_invalid, std::vector<std::pair<unsigned int, unsigned int>> &valid_edges);

	void grow_forw2(vd_t node, max_priority_queue & max_pq, bool & has_cycle, std::vector<vd_t> &K);
	void grow_back2(vd_t node, min_priority_queue & min_pq, bool & has_cycle, std::vector<vd_t> &K);
	void find_cycle(vd_t head);
	void compute_flooring(vd_t n, std::vector<vd_t> &rto);
	ahrsz_ext_priority_value<T> compute_ceiling(vd_t v);
	ahrsz_ext_priority_value<T> compute_priority(ahrsz_ext_priority_value<T> floor, ahrsz_ext_priority_value<T> ceiling);
	~batch_AHRSZ_P_scc();



};

#endif
