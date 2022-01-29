// This is the implementation of the algorithm of HKMST. It was first described in the  following:
//
// [2] Haeupler, B., Kavitha, T., Mathew, R., Sen, S., & Tarjan, R. E. (2012). 
//     Incremental cycle detection, topological ordering, and strong component 
//     maintenance. ACM Transactions on Algorithms (TALG), 8(1), 1-33.
// 
// PK claims that: The arbitrary assignment algorithm is slightly simpler but 
// does not minimise the number of new priorities created. 

#ifndef	DYNAMIC_HKMST_SCC_H_
#define DYNAMIC_HKMST_SCC_H_

#include <queue>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <algorithm>
#include "boost/property_map/property_map.hpp"

#include "ordered_slist.h"
#include "ordered_slist2.h"
#include "disjoint_set.h"
#include "my_graph.h"
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
class dynamic_hkmst_scc{

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
	CG & cg;														// citation graph.
	std::vector< ahrsz_priority_value <T> > & n2v;					// topo value of each node. 
	T & _pspace;													// T: ordered_slist or ordered_slist2,  the priority space.
	std::vector<ahrsz_ext_priority_value<T> > _flooring;			// temporary storage needed by the reassignment stage
	std::vector<bool> _visited;
	std::vector<bool> _inK;
	std::vector<bool> _inB, _inF, in_component;						// use inB and inF to identify the cycle.
	disjoint_set &ds;
	string self_name;
	std::vector<pair<unsigned int, unsigned int>> & inc_edges;
	
	//std::vector<std::pair<out_iterator, out_iterator>> out_its;
	std::vector<out_scc_iterator> out_its;
	std::vector<in_scc_iterator> in_its;
	std::vector<bool> is_out, is_in;

public:
	dynamic_hkmst_scc();
	dynamic_hkmst_scc(CG & cg, std::vector< ahrsz_priority_value <T> > & n2v, T& pspace, disjoint_set & ds, std::vector<pair<unsigned int, unsigned int>> & inc_edges, const int np);
	void add_edges();
	void add_edges(int start_index, int end_index);
	void add_edge(vd_t t, vd_t h);
	bool has_out_after_search(const vd_t node);

	auto first_out(vd_t &node);
	bool out_tail(out_scc_iterator &out);
	bool has_first_out(vd_t &node);
	bool has_next_out(out_scc_iterator &out);

	
	auto first_in(vd_t &node);
	bool in_tail(in_scc_iterator &in);
	bool has_first_in(vd_t &node); 
	bool has_next_in(in_scc_iterator &in);


	void find_cycle(vd_t head);
	void search_step(vd_t & u, vd_t &z, std::list<vd_t>&FA, std::list<vd_t>&BA, 
		std::list<vd_t>&forw_nodes, std::list<vd_t>&back_nodes, bool & has_cycle);

	void soft_threshold_search(vd_t & v, vd_t & w, std::list<vd_t> &forw_nodes, std::list<vd_t> &back_nodes,
		std::list<vd_t> & FA, std::list<vd_t> &FP, std::list<vd_t> &BA, std::list<vd_t> &BP, bool &has_cycle,
		std::vector<unsigned int> &cycle_nodes, vd_t & threshold);

	void reassignment(std::list<vd_t> &forw_nodes, std::list<vd_t> &back_nodes, vd_t& v, vd_t& w);
	void reassignment_with_cycle(std::vector<vd_t> &K);
	void compute_flooring(vd_t n, std::vector<vd_t> &rto );
	ahrsz_ext_priority_value<T> compute_priority(ahrsz_ext_priority_value<T> floor, ahrsz_ext_priority_value<T> ceiling);
	ahrsz_ext_priority_value<T> compute_ceiling(vd_t& v );

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

/*
	public:
		unsigned int sn;		// scc node
	public:
		out_scc_iterator(unsigned int _sn):sn(_sn) {
			
			it_begin = ds.scc_nodes[sn].begin();
			it_end = ds.scc_nodes[sn].end();
			if (it_begin != it_end) {
				oei_begin = boost::out_edges(*it_begin, cg).first;
				oei_end = boost::out_edges(*it_begin, cg).second;
			}
		}

		bool has_first_out() {
			bool flag = false;
			while (it_begin != it_end) {
				if (oei_begin != oei_end) {
					flag = true;
					break;
				}
				it_begin++;
			}
			// reset state
			//it_begin = ds.scc_nodes[sn].begin();
			//it_end = ds.scc_nodes[sn].end();
			//if (it_begin != it_end) {
			//	oei_begin = boost::out_edges(*it_begin, cg).first;
			//	oei_end = boost::out_edges(*it_begin, cg).second;
			//}
			return flag;
		}
		

		
		std::pair<typename dynamic_hkmst_scc<T>::out_iterator, typename dynamic_hkmst_scc<T>::out_iterator> get_first_out(vd_t & _sn) {
			this->sn = _sn; 
			std::pair<out_iterator, out_iterator> scc_item_it; 
			it_begin = ds.scc_nodes[sn].begin();
			it_end = ds.scc_nodes[sn].end();
			if (it_begin != it_end) {
				oei_begin = boost::out_edges(*it_begin, cg).first;
				oei_end = boost::out_edges(*it_begin, cg).second;
			}

			if (is_out[sn] == false) {
				is_out[sn] = true;

				while (it_begin != it_end) {
					scc_item_it = boost::out_edges(*it_begin, cg);
					if (scc_item_it.first != scc_item_it.second) {
						out_its[sn] = scc_item_it; 
						break; 
					}
					it_begin++; 
				}
				return out_its[sn];
			}
			else {
				out_its[sn].first++;
				scc_item_it = out_its[sn];
				while (scc_item_it.first == scc_item_it.second) {
					it_begin++; 
					scc_item_it = boost::out_edges(*it_begin, cg);
				}
				out_its[sc] = scc_item_it; 
				return out_its[sc];
			}
		}
	
	
	
	
	};

	*/


};


#endif

