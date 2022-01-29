#ifndef STATIC_SCC_PEARCE2_H2
#define STATIC_SCC_PEARCE2_H2

#include <string>
#include <vector>
#include <fstream>
#include <iostream> 
#include <map> 
#include <time.h> 
#include "boost/graph/adjacency_list.hpp"


#include "my_func.h"
#include "disjoint_set.h"
#include "scc_para.h"
#include "utils.h"
#include "my_graph.h"


template <typename T > class static_pearce2;

template <typename T>
class static_pearce2 {

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
	ppara para;														// tmp parameter of Pearce2 

	std::vector<bool> act_nodes, act_edges;									// active_nodes 
	std::vector<int> degrees;


	static_pearce2() {}

	static_pearce2(unsigned int _np, const std::string &_pyear, const std::string &_pref) : cg(_np), np(_np), pyear(_pyear), pref(_pref), ds(_np), para(_np) {
		n2y = std::vector<unsigned int>(np);
		n2v = std::vector<T>(np);
		degrees = std::vector <int>(np, 0);
		topo_scc.reserve(np);
		self_name = "Static-Pearce2-";
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
		condense_sccs(scc_nodes, ds, cg);
		update_property_map();
		myres.set_eliminate_scc_time(mt.elapsed());
	}

	void update_property_map() {
		assert((topo_scc.size() == (max_cc_idx - 1)));
		my_timer mt;
		for (int i = 0; i < np; i++ ) {
			n2v[i] = np - para.rindex[i];
		}
	}
	
	void pearce_search2(CG::vertex_descriptor v, ppara & para) {
		myres.access_nodes++;
		boost::graph_traits<CG> ::out_edge_iterator ai, ai_end;
		int root = true;
		para.rindex[v] = para.index_p;
		para.index_p++;
		unsigned int w;

		for (std::tie(ai, ai_end) = out_edges(v, cg); ai != ai_end; ++ai) {
			w = boost::target(*ai, cg);
			myres.access_edges++;
			if (para.rindex[w] == 0) {
				pearce_search2(w, para);
			}

			if (para.rindex[w] < para.rindex[v]) {
				para.rindex[v] = para.rindex[w];
				root = false;
			}
		}
		if (root == true) {
			para.index_p--;
			std::vector<unsigned int> tmp;
			while (!para.S.empty() && (para.rindex[v] <= para.rindex[para.S.back()])) {
				int w = para.S.back();
				tmp.push_back(w);
				para.S.pop_back();
				para.rindex[w] = para.c;
				para.index_p--;
			}
			para.rindex[v] = para.c;
			tmp.push_back(v);
			para.c--;
			topo_scc.push_back(tmp);
		}
		else {
			para.S.push_back(v);
		}
	}


	void topo_sort() {
		my_timer mt;
		boost::graph_traits<CG>::vertex_iterator vi, vi_end;
		for (std::tie(vi, vi_end) = boost::vertices(cg); vi != vi_end; ++vi) {
			if (para.rindex[*vi] == 0 && act_nodes[*vi] == true) {
				pearce_search2(*vi, para);
			}
		}
		myres.set_scc_info(mt.elapsed());
		max_cc_idx = np - para.c;
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

		if (myres.scale_factor < 1) {
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
	

	~static_pearce2(){ 	}
};
#endif

