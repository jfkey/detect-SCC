#ifndef STATIC_SCC_GABOW_TC_H_
#define STATIC_SCC_GABOW_TC_H_

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

template <typename T > class static_gabow_tc;

template <typename T >
class static_gabow_tc  {
public:
	CG cg;																		// citation graph.
	std::vector<unsigned int> n2y;												// publish year of each paper.
	unsigned int np;															// number of paper.
	unsigned int max_year = 0, min_year = 9999;									// max min year of the citation graph.
	std::vector<bool> is_in_unorder;											// identify if the node x is in the unorder_graph. 

	std::string pyear;															// file name of pyear. 
	std::string pref;															// file name of pp.
	std::string self_name;														// the name of this alg. 
	std::vector<std::pair<unsigned int, unsigned int>> normal_edgelist;			// all edges of the citation graph.
	std::vector<std::pair<unsigned int, unsigned int>> unorder_edgelist;
	std::vector<std::vector<std::pair<unsigned int, unsigned int> >>same_year_edgelist;
	std::vector<std::vector<unsigned int>> topo_scc;							// all of the sccs. 

	std::vector<bool> act_nodes, act_edges;									// active_nodes 
	std::vector<int> degrees;

	static_gabow_tc() {}
	static_gabow_tc(unsigned int _np, const std::string &_pyear, const std::string &_pref) :cg(_np), np(_np), pyear(_pyear), pref(_pref) {
		max_year = 0;
		min_year = 9999;
		n2y = std::vector<unsigned int>(np);
		is_in_unorder = std::vector<bool>(np, false);
		degrees = std::vector <int>(np, 0);
		topo_scc.reserve(np);
		self_name = "Static-Gabow-TC-";
		myres.set_static_name(self_name);
	}

	void eliminate_static_scc() {
	}

	void topo_sort() {
		gabow_para para(np);
		my_timer mt;
		for (std::vector<std::pair<unsigned int, unsigned int>>::const_iterator it = unorder_edgelist.cbegin(); it != unorder_edgelist.cend(); it++) {
			if (para.num[it->second] == 0 && act_nodes[it->second] == true) {
				GabowDFS_unorder(it->second, para);
			}
		}

		std::pair<unsigned int, unsigned int> edge;
		std::vector<std::pair<unsigned int, unsigned int>> cur_year_edges;
		unsigned int y = 0;
		for (int i = 0; i < same_year_edgelist.size(); i++) {
			cur_year_edges = same_year_edgelist[i];
			for (int j = 0; j < cur_year_edges.size(); j++) {
				edge = cur_year_edges[j];
				if (is_in_unorder[edge.second] == false && para.num[edge.second] == 0 && act_nodes[edge.second] == true) {
					GabowDFS_same_year(edge.second, para, n2y[edge.second]); // gabow in the original graph for same year subgraph 
				 
				}
			}
		}
		myres.set_scc_info(mt.elapsed());
		record_cur_per_info(myres);
	}

	void GabowDFS_same_year(CG::vertex_descriptor v, gabow_para& para, unsigned int &year) {
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
			if (para.num[w] == 0 && is_in_unorder[w] == false && n2y[w] == year) {
				GabowDFS_same_year(w, para, year);
			}
			else if (!para.inscc[w] && is_in_unorder[w] == false && n2y[w] == year) {
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


	void GabowDFS_unorder(CG::vertex_descriptor v, gabow_para& para) {
		is_in_unorder[v] = true;
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
				GabowDFS_unorder(w, para);
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
		iread_pyear(pyear, n2y, min_year, max_year);
		same_year_edgelist = std::vector<std::vector<std::pair<unsigned int, unsigned int> >>(max_year - min_year + 1);


		if (myres.scale_factor < 1.0) {
			std::vector<std::pair<unsigned int, unsigned int>> edgelist;
			iread_edgelist(pref, np, edgelist, n2y, degrees);
			act_nodes = std::vector<bool>(np, false);
			act_edges = std::vector<bool>(edgelist.size(), false);
			gen_new_edges(edgelist, degrees, act_nodes, act_edges);
			for (int i = 0; i < edgelist.size(); i++) {
				if (act_edges[i] == true) {
					if (n2y[edgelist[i].first] < n2y[edgelist[i].second]) {
						unorder_edgelist.push_back(edgelist[i]);
					}
					else if (n2y[edgelist[i].first] == n2y[edgelist[i].second]) {
						same_year_edgelist[n2y[edgelist[i].first] - min_year].push_back(edgelist[i]);
					}
					else {
						normal_edgelist.push_back(edgelist[i]);
					}
				}
			}

		}
		else {

			iread_edgelist(pref, np, normal_edgelist, unorder_edgelist, same_year_edgelist, n2y, min_year);
			int edgecnt = 0;
			for_each(same_year_edgelist.begin(), same_year_edgelist.end(), [&edgecnt](std::vector<std::pair<unsigned int, unsigned int>> a) {
				edgecnt += a.size();
			});
			act_nodes = std::vector<bool>(np, true);
			act_edges = std::vector<bool>(edgecnt, true);
		}

		int edgecnt = 0;
		for_each(same_year_edgelist.begin(), same_year_edgelist.end(), [&edgecnt](std::vector<std::pair<unsigned int, unsigned int>> a) {
			edgecnt += a.size();
		});

		myres.set_cites(unorder_edgelist.size(), edgecnt, normal_edgelist.size());
		myres.set_read_data_info(np, edgecnt + unorder_edgelist.size() + normal_edgelist.size(), max_year, min_year);
	}

	void gen_citation_graph() {
		std::vector<std::pair<unsigned int, unsigned int>>::const_iterator it;
		for (std::vector < std::vector<std::pair<unsigned int, unsigned int>>>::const_iterator it1 = same_year_edgelist.cbegin(); it1 != same_year_edgelist.cend(); it1++) {
			for (std::vector<std::pair<unsigned int, unsigned int>>::const_iterator it2 = it1->cbegin(); it2 != it1->cend(); it2++) {
				boost::add_edge((*it2).first, (*it2).second, cg);
			}
		}
		for (std::vector<std::pair<unsigned int, unsigned int>>::const_iterator it = unorder_edgelist.cbegin(); it != unorder_edgelist.cend(); ++it) {
			boost::add_edge((*it).first, (*it).second, cg);
		}
		for (std::vector<std::pair<unsigned int, unsigned int>>::const_iterator it = normal_edgelist.cbegin(); it != normal_edgelist.cend(); ++it) {
			boost::add_edge((*it).first, (*it).second, cg);
		}
		normal_edgelist.clear();
	}


	~static_gabow_tc(){ }
};
#endif

