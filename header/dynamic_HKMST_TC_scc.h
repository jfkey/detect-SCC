// This is the implementation of the algorithm of Time-Coupled-HKMST. HKMST was first described in the  following:
//
// [1] Katriel, I., & Bodlaender, H. L. (2006). Online topological ordering. 
//     ACM Transactions on Algorithms (TALG), 2(3), 364-379.
// 
// The strongly connected component maintenance arises from the following paper:
// 
// [2] Haeupler, B., Kavitha, T., Mathew, R., Sen, S., & Tarjan, R. E. (2012). 
//     Incremental cycle detection, topological ordering, and strong component 
//     maintenance. ACM Transactions on Algorithms (TALG), 8(1), 1-33.
// 
// PK claims that: The arbitrary assignment algorithm is slightly simpler but does 
// not minimise the number of new priorities created. 
// This statement is true only when the size of the graph is very small, e.g., 2000 nodes and 2000 edges. 
// one could compare the efficiency of different ALGs, i.e., MNR, PK, PK2, AHRSZ, HKMST. 
// all algs are equipped with online topo sort and scc detection. 
// 
// We divide the citation graph into three types of subgraphs, i.e., Gm, Gsi, Gr
// The strongly connected components are only held in Gm and Gsi due to the time 
// property of the citation graph. Thus, we develop alg1, alg2, and alg3 to maintain 
// the Invariants of Gm, Gsi and Gr when doing incremental cycle detection.
//  

#ifndef	DYNAMIC_TC_HKMST_SCC_H_
#define DYNAMIC_TC_HKMST_SCC_H_

#include <queue>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <algorithm>
#include "boost/property_map/property_map.hpp"

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
class dynamic_tc_hkmst_scc {
private:
	typedef CG::vertex_descriptor vd_t;
	typedef CG::out_edge_iterator out_iterator;
	typedef CG::in_edge_iterator in_iterator;
public:
	class out_scc_iterator {
	public:
		std::vector<unsigned int>::iterator it_begin, it_end;
		out_iterator oei_begin, oei_end;
	};

	class in_scc_iterator {
	public:
		std::vector<unsigned int>::iterator it_begin, it_end;
		in_iterator iei_begin, iei_end;
	};

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

	//std::vector<std::pair<out_iterator, out_iterator>> out_its;
	std::vector<out_scc_iterator> out_its;
	std::vector<in_scc_iterator> in_its;
	std::vector<bool> is_out, is_in;

public:
	dynamic_tc_hkmst_scc();
	dynamic_tc_hkmst_scc(CG & cg, vector< ahrsz_priority_value <T> > & _n2v, vector<T*> & _pspaces, disjoint_set & _ds,
		vector<unsigned int> & _n2y, vector<int>& _is_in_same_year, vector<bool> &_is_in_unorder, int _min_year,
		int _max_year, std::vector<pair<unsigned int, unsigned int>> & _inc_edges, const int np);
	
	void add_edges();
	void add_edges(int start_index, int end_index);
	void add_edge_Gm(vd_t t, vd_t h);
	void add_edge_Gsi(vd_t t, vd_t h);
	void add_edge_alg1(vd_t t, vd_t h);
	void add_edge_alg2(vd_t t, vd_t h);
	void add_edge_alg3(vd_t t, vd_t h);
	void scan_GsGr(vd_t n, std::vector<unsigned int> &reachable, ahrsz_ext_priority_value<T> & max_priority);
	void scan_GsGr_without_topo(vd_t n, vd_t lb, std::vector<unsigned int> &reachable, ahrsz_ext_priority_value<T> & max_priority, bool & flag);
	

	void find_cycle_Gm(vd_t head);
	void find_cycle_Gsi(vd_t head);
			
	void search_step_Gm(vd_t & u, vd_t &z, std::list<vd_t>&FA, std::list<vd_t>&BA,
		std::list<vd_t>&forw_nodes, std::list<vd_t>&back_nodes, bool & has_cycle);
	void search_step_Gsi(vd_t & u, vd_t &z, std::list<vd_t>&FA, std::list<vd_t>&BA,
		std::list<vd_t>&forw_nodes, std::list<vd_t>&back_nodes, bool & has_cycle);

	void soft_threshold_search_Gm(vd_t & v, vd_t & w, std::list<vd_t> &forw_nodes, std::list<vd_t> &back_nodes,
		std::list<vd_t> & FA, std::list<vd_t> &FP, std::list<vd_t> &BA, std::list<vd_t> &BP, bool &has_cycle,
		std::vector<unsigned int> &cycle_nodes, vd_t & threshold);
	void soft_threshold_search_Gsi(vd_t & v, vd_t & w, std::list<vd_t> &forw_nodes, std::list<vd_t> &back_nodes,
		std::list<vd_t> & FA, std::list<vd_t> &FP, std::list<vd_t> &BA, std::list<vd_t> &BP, bool &has_cycle,
		std::vector<unsigned int> &cycle_nodes, vd_t & threshold);

	void maintain_Gm(std::list<vd_t> &forw_nodes, std::list<vd_t> &back_nodes, vd_t& v, vd_t& w);
	void maintain_Gsi(std::list<vd_t> &forw_nodes, std::list<vd_t> &back_nodes, vd_t& v, vd_t& w);
	void maintain_SCC_Gm(std::vector<vd_t> &K);
	void maintain_SCC_Gsi(std::vector<vd_t> &K);

	ahrsz_ext_priority_value<T> compute_priority(ahrsz_ext_priority_value<T> floor, ahrsz_ext_priority_value<T> ceiling, T & cur_pspace);
	void compute_flooring_Gm(vd_t n, std::vector<vd_t> &rto );
	void compute_flooring_Gsi(vd_t n, std::vector<vd_t> &rto);
	ahrsz_ext_priority_value<T> compute_ceiling_Gm(vd_t v);
	ahrsz_ext_priority_value<T> compute_ceiling_Gsi(vd_t v);

	//bool has_first_out(vd_t &node);
	bool has_first_out_Gm(vd_t &node);
	bool has_first_out_Gsi(vd_t &node);

	//bool has_first_in(vd_t &node);
	bool has_first_in_Gm(vd_t &node);
	bool has_first_in_Gsi(vd_t &node);

	//bool out_tail(out_scc_iterator &out);
	bool out_tail_Gm(out_scc_iterator &out);
	bool out_tail_Gsi(out_scc_iterator &out, vd_t & node);

	//bool in_tail(in_scc_iterator &in);
	bool in_tail_Gm(in_scc_iterator &in);
	bool in_tail_Gsi(in_scc_iterator &in,vd_t &node);

	//bool has_out_after_search(const vd_t node);
	bool has_out_after_search_Gm(vd_t node);
	bool has_out_after_search_Gsi(vd_t node);

	//auto first_out(vd_t &node);
	auto first_out_Gm(vd_t &node);
	auto first_out_Gsi(vd_t &node);

	//auto first_in(vd_t &node);
	auto first_in_Gm(vd_t &node);
	auto first_in_Gsi(vd_t &node);

	//bool has_next_out(out_scc_iterator &out);
	bool has_next_out_Gm(out_scc_iterator &out);
	bool has_next_out_Gsi(out_scc_iterator &out, vd_t node);

	//bool has_next_in(in_scc_iterator &in);
	bool has_next_in_Gm(in_scc_iterator &in);
	bool has_next_in_Gsi(in_scc_iterator  &in, vd_t node);

};

#endif





