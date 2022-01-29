#ifndef STATIC_SCC_GABOW_H_
#define STATIC_SCC_GABOW_H_

#include <string>
#include <vector>
#include <fstream>
#include <iostream> 
#include <map> 
#include <time.h> 
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/topological_sort.hpp"
#include "boost/graph/strong_components.hpp"

#include "my_func.h"
#include "disjoint_set.h"
#include "scc_para.h"
#include "utils.h"



template <typename T > class static_gabow;

template <typename T>
class static_gabow {
public:
	CG cg;															// citation graph.
	std::vector<unsigned int> n2y;									// publish year of each paper.
	unsigned int np;												// number of paper.
	unsigned int max_cc_idx;										// the the size of cc. 

	std::string pyear;												// file name of pyear. 
	std::string pref;												// file name of pp.
	std::string self_name;											// the name of this alg. 
	std::vector<std::pair<unsigned int, unsigned int>> edgelist;	// all edges of the citation graph.
	std::vector<std::vector<unsigned int>> topo_scc;				// all of the sccs. 

	std::vector<bool> act_nodes, act_edges;									// active_nodes 
	std::vector<int> degrees;

	static_gabow() {}

	static_gabow(unsigned int _np, const std::string &_pyear, const std::string &_pref) :cg(_np), np(_np), pyear(_pyear), pref(_pref){
		n2y = std::vector<unsigned int>(np);
		degrees = std::vector <int>(np, 0);
		topo_scc.reserve(np);
		self_name = "Static-Gabow-";
		myres.static_name = self_name;
	}

	void eliminate_static_scc() {
	}

	void topo_sort() {
		my_timer mt;
		gabow_para para(np);
		boost::graph_traits<CG>::vertex_iterator vi, vi_end;
		for (std::tie(vi, vi_end) = boost::vertices(cg); vi != vi_end; ++vi) {
			if (para.num[*vi] == 0 && act_nodes[*vi] == true) {
				GabowDFS(*vi, para);
			}
		}
		max_cc_idx = para.cc_idx;
		myres.scc_time = mt.elapsed();
		record_cur_per_info(myres);
	}


	// void strong_connected(typename T::vertex_descriptor v, scc_para & para, std::vector<std::vector<unsigned int>> &topo_scc, T &graph) {
	void GabowDFS(CG::vertex_descriptor v, gabow_para& para) {
		boost::graph_traits<CG> ::out_edge_iterator ai, ai_end;
		myres.access_nodes++;
		para.num[v] = para.cc_idx;
		para.cc_idx++;
		para.P.push(v);
		para.S.push(v);
		unsigned int w;
		for (std::tie(ai, ai_end) = out_edges(v, cg); ai != ai_end; ++ai) {
			w = boost::target(*ai, cg);
			myres.access_edges++;
			if (para.num[w] == 0) {
				GabowDFS(w, para);
			}
			else if (!para.inscc[w]) {
				while (para.num[para.P.top()] > para.num[w]) {
					para.P.pop();
				}
			}
		}
		if (v == para.P.top()) {
			unsigned int x;
			std::vector<unsigned int> tmp;
			if (!para.inscc[v]) {
				tmp.push_back(v);
				para.inscc[v] = true;
			}
			while (v != para.S.top()) {
				x = para.S.top();
				if (!para.inscc[x]) {
					tmp.push_back(x);
					para.inscc[x] = true;
				}
				para.S.pop();
			}
			topo_scc.push_back(tmp);
			para.P.pop();
		}
	}


	void print_scc(const std::string & dir) {
		//print_scc_infile(dir, topo_scc, self_name, 0);
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


	~static_gabow() { 	}
};
#endif

