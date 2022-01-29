#ifndef STATIC_SCC_KOSARAJU_H_
#define STATIC_SCC_KOSARAJU_H_

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


template <typename T> class static_kosaraju;

template <typename T >
class static_kosaraju {
public:

	CG cg;															// citation graph.
	std::vector<unsigned int> n2y;									// publish year of each paper.
	unsigned int np;												// number of paper.

	std::string pyear;												// file name of pyear. 
	std::string pref;												// file name of pp.
	std::string self_name;											// the name of this alg. 
	std::vector<std::pair<unsigned int, unsigned int>> edgelist;	// all edges of the citation graph.
	std::vector<std::vector<unsigned int>> topo_scc;				// all of the sccs. 

	std::vector<bool> act_nodes, act_edges;									// active_nodes 
	std::vector<int> degrees;



	static_kosaraju() {}
	static_kosaraju(unsigned int _np, const std::string &_pyear, const std::string &_pref) :cg(_np), np(_np), pyear(_pyear), pref(_pref) {
		n2y = std::vector<unsigned int>(np);
		topo_scc.reserve(np);
		degrees = std::vector <int>(np, 0);
		self_name = "Static-Kosaraju-";
		myres.set_static_name(self_name);
	}

	void eliminate_static_scc() {}


	void topo_sort() {
		my_timer mt;
		kosaraju_para para(np);
		boost::graph_traits<CG>::vertex_iterator vi, vi_end;
		for (std::tie(vi, vi_end) = boost::vertices(cg); vi != vi_end; ++vi) {
			if (para.marked[*vi] == false && act_nodes[*vi] == true) {
				dfs_r(*vi, para);
			}
		}
		para.reset_marked();

		std::vector<unsigned int>::reverse_iterator ri(para.nodes.rbegin());
		unsigned int w;
		for (ri; ri != para.nodes.rend(); ri++) {
			if (para.marked[*ri] == false && act_nodes[*ri] == true) {
				std::vector<unsigned int> scc;
				dfs(*ri, para, scc);
				topo_scc.push_back(scc);
			}
		}
		myres.set_scc_info(mt.elapsed());
		record_cur_per_info(myres);
	}

	void dfs_r(CG::vertex_descriptor v, kosaraju_para & para) {
		myres.access_nodes++;
		para.marked[v] = true;
		unsigned int w;
		boost::graph_traits<CG> ::in_edge_iterator ai, ai_end;
		for (std::tie(ai, ai_end) = in_edges(v, cg); ai != ai_end; ++ai) {
			myres.access_edges++;
			w = boost::source(*ai, cg);
			if (para.marked[w] == false) {
				dfs_r(w, para);
			}
		}
		para.nodes.push_back(v);
	}

	void dfs(CG::vertex_descriptor v, kosaraju_para &para, std::vector<unsigned int>& cur_scc) {
		myres.access_nodes++;
		para.marked[v] = true;
		cur_scc.push_back(v);
		boost::graph_traits<CG> ::out_edge_iterator ai, ai_end;
		unsigned int w;
		for (std::tie(ai, ai_end) = out_edges(v, cg); ai != ai_end; ++ai) {
			myres.access_edges++;
			w = boost::target(*ai, cg);
			if (para.marked[w] == false) {
				dfs(w, para, cur_scc);
			}
		}
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

	~static_kosaraju(){	}
};
#endif

