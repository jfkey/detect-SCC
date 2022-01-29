#ifndef STATIC_TOPO_HPP
#define STATIC_TOPO_HPP

#include <string>
#include <vector>
#include <fstream>
#include <iostream> 
#include <map> 
#include <algorithm>
#include <stack>
#include <set>
#include <time.h> 
#include "boost/graph/adjacency_list.hpp"

#include "my_func.h"
#include "disjoint_set.h"
#include "scc_para.h"
#include "utils.h"

template <typename T > class static_tarjan;

template <typename T>
class static_tarjan {
public:
	
	CG cg;															// citation graph.
	std::vector<T> n2v;												// topo value of each node. 
	std::vector<unsigned int> n2y;									// publish year of each paper.
	unsigned int np;												// number of paper.
	disjoint_set ds;												// disjoint set to represent the scc. 
	unsigned int max_cc_idx;										// the the size of cc. 

	std::string pyear;												// file name of pyear. 
	std::string pref;												// file name of pp.
	std::string self_name;											// the name of this alg. 
	std::vector<std::pair<unsigned int, unsigned int>> edgelist;	// all edges of the citation graph.
	std::vector<std::vector<unsigned int>> topo_scc;				// all of the sccs. 
	std::vector<unsigned int> cc;
	std::vector<bool> act_nodes, act_edges;									// active_nodes 
	std::vector<int> degrees;

	static_tarjan() {}

	static_tarjan(unsigned int _np, const std::string &_pyear, const std::string &_pref):cg(_np), np(_np), pyear(_pyear), pref(_pref), ds(_np) {
		n2y = std::vector<unsigned int>(np);
		n2v = std::vector<T>(np);
		degrees = std::vector <int>(np, 0);
		topo_scc.reserve(np);
		self_name = "Static-Tarjan-";
		myres.set_static_name(self_name);
	}


	// each scc (|scc| > 1) store in the vector.
	void get_scc_nodes(std::vector<std::vector<unsigned int>> &res) {
		std::vector<std::vector<unsigned int>>::iterator it1;
		std::vector<unsigned int>::iterator it2;
		for (it1 = topo_scc.begin(); it1 != topo_scc.end(); it1++) {
			if ((*it1).size() > 1) {
				res.push_back(*it1);
			}
		}
	}

	void eliminate_static_scc() {
		my_timer mt;
		std::vector<std::vector<unsigned int>> scc_nodes;
		get_scc_nodes(scc_nodes);
		//condense_sccs(scc_nodes, ds, cg);
		//update_property_map();
		myres.set_eliminate_scc_time(mt.elapsed());
	}


	void update_property_map() {
		assert((topo_scc.size() == (max_cc_idx - 1)));
		boost::graph_traits<CG>::vertex_iterator vi, vi_end;
		my_timer mt;
		for (std::tie(vi, vi_end) = vertices(cg); vi != vi_end; ++vi) {
			n2v[*vi] = cc[*vi];
		}
		
	}

	void topo_sort() {
		my_timer mt;
		scc_para para(np);
		boost::graph_traits<CG>::vertex_iterator vi, vi_end;
		for (std::tie(vi, vi_end) = boost::vertices(cg); vi != vi_end; ++vi) {
			if (para.num[*vi] == 0 && act_nodes[*vi] == true) { // && rm_nodes[*vi] == false
				strong_connected(*vi, para, topo_scc, cg);
			}
		}
		myres.set_scc_info(mt.elapsed());
		max_cc_idx = para.cc_idx;
		cc = para.cc;
		record_cur_per_info(myres);
	}



	void print_scc(const std::string & dir) {
		myres.gen_sccs(topo_scc);
	}

	void read_citation_graph() {
		unsigned int min_year = 9999, max_year = 0;
		iread_pyear(pyear, n2y, min_year, max_year);
		iread_edgelist(pref, np, edgelist, n2y, degrees);
		myres.set_read_data_info(np, edgelist.size(), max_year, min_year);
		act_nodes = std::vector<bool>(np, false);
		act_edges = std::vector<bool>(edgelist.size(), false);
		
		if (myres.scale_factor < 1 ) {
			gen_new_edges(edgelist, degrees, act_nodes, act_edges);
		}
		else {
			act_nodes = std::vector<bool>(np, true);
			act_edges = std::vector<bool>(edgelist.size(), true);
		}
	}

	void gen_citation_graph() {

		std::vector<std::pair<unsigned int, unsigned int>>::const_iterator it;
		int cnt = 0;
		for (it = edgelist.cbegin(); it != edgelist.cend(); it++) {
			if (act_edges[cnt] == true) {
				boost::add_edge(it->first, it->second, cg);
			}
			cnt++; 
		}
		edgelist.clear();
	}

};
#endif

