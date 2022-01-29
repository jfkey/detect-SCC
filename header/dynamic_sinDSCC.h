// This is the implementation of the algorithm of Time-Coupled-AHRSZ. AHRSZ was first described in the following:
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
// We divide the citation graph into three types of subgraphs, i.e., Gm, Gsi, Gr
// The strongly connected components are only held in Gm and Gsi due to the time 
// property of the citation graph. Thus, we develop alg1, alg2, and alg3 to maintain 
// the Invariants of Gm, Gsi and Gr when doing incremental cycle detection.
// 


#ifndef DYNAMIC_TC_AHRSZ_SCC_HPP
#define DYNAMIC_TC_AHRSZ_SCC_HPP

#include <queue>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <tuple>
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


template<class T>
class sinDSCC;

template<class T>
class sinDSCC {
private:
	typedef CG::vertex_descriptor vd_t;
	typedef CG::out_edge_iterator out_iterator;
	typedef CG::in_edge_iterator in_iterator;
	typedef std::priority_queue<vd_t, std::vector<vd_t>, ahrsz_priority_comp<std::vector< ahrsz_priority_value <T> >, std::less<ahrsz_priority_value<T> > > > max_priority_queue;
	typedef std::priority_queue<vd_t, std::vector<vd_t>, ahrsz_priority_comp<std::vector< ahrsz_priority_value <T> >, std::greater<ahrsz_priority_value<T> > > > min_priority_queue;

public:
	CG & cg;																	// citation graph.
	std::vector< ahrsz_priority_value <T> > & n2v;								// topo value of each node. 
	std::vector<unsigned int> &n2y;												// publish year of each paper.
	std::vector<T*> & pspaces; 													// T: ordered_slist or ordered_slist2,  the priority space.  sizeof pspaces: max_year - min_year + 2
	std::vector<ahrsz_ext_priority_value<T> > _flooring;						// temporary storage needed by the reassignment stage
	std::vector<bool> _visited;
	std::vector<bool> _inK;
	std::vector<bool> _inB, _inF, in_component;									// use inB and inF to identify the cycle.
	std::vector<bool> & is_in_unorder;
	std::vector<int> & is_in_same_year;
	disjoint_set &ds;
	string self_name;
	int min_year, max_year;
	std::vector<pair<unsigned int, unsigned int>> & inc_edges;

public:
	sinDSCC();
	sinDSCC(CG & cg, vector< ahrsz_priority_value <T> > & n2v, vector<T*> & _pspaces, disjoint_set & ds, 
		vector<unsigned int> & n2y, vector<int>& _is_in_same_year, vector<bool> &_is_in_unorder, int _min_year, 
		int _max_year, std::vector<pair<unsigned int, unsigned int>> & inc_edges, const int np);
	
	void find_cycle_Gm(vd_t head);
	void find_cycle_Gsi(vd_t head);
	void discover_Gm(vd_t tail, vd_t head, std::vector<vd_t> &K, bool & has_cycle, std::vector<unsigned int> &cycle_nodes);
	void discover_Gsi(vd_t tail, vd_t head, std::vector<vd_t> &K, bool & has_cycle, std::vector<unsigned int> &cycle_nodes);
	void maintain_Gm(std::vector<vd_t> &K);
	void maintain_Gsi(std::vector<vd_t> &K);
	void compute_flooring_Gm(vd_t n, std::vector<vd_t> &rto );
	void compute_flooring_Gsi(vd_t n, std::vector<vd_t> &rto);
	ahrsz_ext_priority_value<T> compute_priority(ahrsz_ext_priority_value<T> floor, ahrsz_ext_priority_value<T> ceiling, T & cur_pspace);
	ahrsz_ext_priority_value<T> compute_ceiling_Gm(vd_t v);
	ahrsz_ext_priority_value<T> compute_ceiling_Gsi(vd_t v);
	void add_edges();
	void add_edges(int start_index, int end_index);
	void add_edge_Gm(vd_t t, vd_t h);
	void add_edge_alg1(vd_t t, vd_t h);
	void add_edge_alg2(vd_t t, vd_t h);
	void add_edge_alg3(vd_t t, vd_t h);
	void scan_GsGr(vd_t n, std::vector<unsigned int> & reachable, ahrsz_ext_priority_value<T> &max_priority);
	void scan_GsGr_without_topo(vd_t n, vd_t lb, std::vector<unsigned int> &reachable, ahrsz_ext_priority_value<T> &max_priority, bool &flag);
	void scan_GsGr_without_topo2(vd_t n, vd_t lb, std::vector<unsigned int> &reachable, ahrsz_ext_priority_value<T> &max_priority, bool &flag);
	void add_edge_Gsi(vd_t t, vd_t h);

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
	}*/
};

#endif





