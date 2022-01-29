// (C) Copyright David James Pearce 2003. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
// Email: david.pearce@mcs.vuw.ac.nz
// 
// 
// Notes
// =====
// This is the implementation of the algorithm of AHRSZ+. AHRSZ was first described in the  following:
//
// [1] B. Alpern, R. Hoover, B. K. Rosen, P. F. Sweeney and F. K. 
//     Zadeck, "Incremental Evaluation of Computational Circuits", 
//     In Proc. of the First Annual ACM-SIAM Symposium on Discrete 
//     Algorithms, pages 32-42, 1990.
// 
// The strongly connected component maintenance is similar to the following paper:
// 
// [2] Haeupler, B., Kavitha, T., Mathew, R., Sen, S., & Tarjan, R. E. (2012).
//     Incremental cycle detection, topological ordering, and strong component 
//	   maintenance. ACM Transactions on Algorithms (TALG), 8(1), 1-33.
// 

#ifndef DYNAMIC_AHRSZ_SCC_H_
#define DYNAMIC_AHRSZ_SCC_H_

#include <queue>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <tuple>
#include "boost/property_map/property_map.hpp"

#include "ordered_slist.h"
#include "ordered_slist2.h"
#include "my_graph.h"
#include "disjoint_set.h"
#include "my_func.h"
#include "my_greater.h"
#include "my_less.h"


#if __GNUC__ >= 3
#include <ext/functional>
using __gnu_cxx::select2nd;
using __gnu_cxx::identity;
#else
#include "my_select.h"
#include "boost/compute/functional/identity.hpp"
#endif


template<typename T >
class dynamic_ahrsz_scc {
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
	disjoint_set &ds;
	string self_name;
	std::vector<pair<unsigned int, unsigned int>> & inc_edges; 
	

public: 
	dynamic_ahrsz_scc();
	dynamic_ahrsz_scc(CG & cg, std::vector< ahrsz_priority_value <T> > & n2v, T& pspace, disjoint_set & ds, std::vector<pair<unsigned int, unsigned int>> & inc_edges, const int np);
	
	void add_edges();
	void add_edges(int start_index, int end_index);
	void do_ahrsz(vd_t t, vd_t h);
	void find_cycle(vd_t head);
	void discovery(vd_t tail, vd_t head, std::vector<vd_t> &K, bool & has_cycle, std::vector<unsigned int> &cycle_nodes);
	void reassignment(std::vector<vd_t> &K);
	void reassignment_mini_priority(std::vector<vd_t> &K);
	void compute_flooring(vd_t n, std::vector<vd_t> &rto);
	ahrsz_ext_priority_value<T> compute_priority(ahrsz_ext_priority_value<T> floor, ahrsz_ext_priority_value<T> ceiling);
	ahrsz_ext_priority_value<T> compute_ceiling(vd_t v);



public:
	/*
	std::string a2str(ahrsz_ext_priority_value<T> x) {
		if (x.minus_inf()) {
			return "-oo";
		}
		else if (x.plus_inf()) {
			return "+oo";
		}
		uint64_t r(_pspace.order(x.base()));
		return boost::lexical_cast<std::string>(r);
	}

	std::string a2str(ahrsz_priority_value<T> x) {
		uint64_t r(_pspace.order(x.base()));
		return boost::lexical_cast<std::string>(r);
	}
	*/
};

#endif


