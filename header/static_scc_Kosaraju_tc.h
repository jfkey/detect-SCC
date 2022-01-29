#ifndef STATIC_SCC_KOSARAJU_TC_HPP
#define STATIC_SCC_KOSARAJU_TC_HPP

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

template <typename T> class static_kosaraju_tc;

template <typename T>
class static_kosaraju_tc  {
public:
	CG cg;																		// citation graph.
	std::vector<unsigned int> n2y;												// publish year of each paper.
	unsigned int np;															// number of paper.
	unsigned int max_year = 0, min_year = 9999;									// max min year of the citation graph.
	std::vector<bool> is_in_unorder;											// identify if the node x is in the unorder_graph. 
	std::vector<bool> is_in_same_year;

	std::string pyear;															// file name of pyear. 
	std::string pref;															// file name of pp.
	std::string self_name;														// the name of this alg. 
	std::vector<std::pair<unsigned int, unsigned int>> normal_edgelist;			// all edges of the citation graph.
	std::vector<std::pair<unsigned int, unsigned int>> unorder_edgelist;
	std::vector<std::vector<std::pair<unsigned int, unsigned int> >>same_year_edgelist;
	std::vector<std::vector<unsigned int>> topo_scc;							// all of the sccs. 

	std::vector<bool> act_nodes, act_edges;									// active_nodes 
	std::vector<int> degrees;


	static_kosaraju_tc() {}
	static_kosaraju_tc( unsigned int _np, const std::string &_pyear, const std::string &_pref) :cg(_np), np(_np), pyear(_pyear), pref(_pref) {
		max_year = 0;
		min_year = 9999;
		n2y = std::vector<unsigned int>(np);
		is_in_unorder = std::vector<bool>(np, false);
		is_in_same_year = std::vector<bool>(np, false);
		degrees = std::vector <int>(np, 0);
		topo_scc.reserve(np);
		self_name = "Static-Kosaraju-TC-";
		myres.set_static_name(self_name);
	}



	void eliminate_static_scc() {
	}


	void topo_sort() {
		my_timer mt;
		kosaraju_para para(np);
		boost::graph_traits<CG>::vertex_iterator vi, vi_end;

		// get unoder subgraph scc. 
		for (std::vector<std::pair<unsigned int, unsigned int>>::const_iterator it = unorder_edgelist.cbegin(); it != unorder_edgelist.cend(); it++) {
			if (para.marked[it->second] == false && act_nodes[it->second] == true) {
				dfs_r(it->second, para);
			}
		}
		para.reset_marked();
		std::vector<unsigned int>::reverse_iterator ri(para.nodes.rbegin());
		unsigned int w;
		for (ri; ri != para.nodes.rend(); ri++) {
			if (para.marked[*ri] == false && is_in_unorder[*ri] == true) {
				std::vector<unsigned int> scc;
				dfs(*ri, para, scc);
				topo_scc.push_back(scc);
			}
		}
		para.reset_marked();
		para.nodes.clear();

		// get same year subgraph scc. 
		std::pair<unsigned int, unsigned int> edge;
		std::vector<std::pair<unsigned int, unsigned int>> cur_year_edges;
		unsigned int y = 0;
		for (int i = 0; i < same_year_edgelist.size(); i++) {
			cur_year_edges = same_year_edgelist[i];
			for (int j = 0; j < cur_year_edges.size(); j++) {
				edge = cur_year_edges[j];
				if (para.marked[edge.second] == 0 && is_in_unorder[edge.second] == false && act_nodes[edge.second] == true) {
					y = n2y[edge.second];
					dfs_r_same_year(edge.second, para, y);
				}
			}
		}
		para.reset_marked();
		ri = para.nodes.rbegin();
		for (ri; ri != para.nodes.rend(); ri++) {
			if (para.marked[*ri] == false) {
				std::vector<unsigned int> scc;
				y = n2y[*ri];
				dfs_same_year(*ri, para, scc, y);
				topo_scc.push_back(scc);
			}
		}

		myres.set_scc_info(mt.elapsed());
		record_cur_per_info(myres);
	}


	void dfs_r_same_year(CG::vertex_descriptor v, kosaraju_para & para, unsigned int &year) {
		myres.access_nodes++;
		para.marked[v] = true;
		is_in_same_year[v] = true;
		unsigned int w;
		boost::graph_traits<CG> ::in_edge_iterator ai, ai_end;
		for (std::tie(ai, ai_end) = in_edges(v, cg); ai != ai_end; ++ai) {
			myres.access_edges++;
			w = boost::source(*ai, cg);
			if (para.marked[w] == false && is_in_unorder[w] == false && n2y[w] == year) {
				dfs_r_same_year(w, para, year);
			}
		}
		para.nodes.push_back(v);
	}

	void dfs_same_year(CG::vertex_descriptor v, kosaraju_para &para, std::vector<unsigned int>& cur_scc, unsigned int & year) {
		myres.access_nodes++;
		para.marked[v] = true;
		cur_scc.push_back(v);
		boost::graph_traits<CG> ::out_edge_iterator ai, ai_end;
		unsigned int w;
		for (std::tie(ai, ai_end) = out_edges(v, cg); ai != ai_end; ++ai) {
			myres.access_edges++;
			w = boost::target(*ai, cg);
			if (para.marked[w] == false && is_in_unorder[w] == false && n2y[w] == year && is_in_same_year[w] == true) {
				dfs_same_year(w, para, cur_scc, year);
			}
		}
	}

	void dfs_r(CG::vertex_descriptor v, kosaraju_para & para) {
		is_in_unorder[v] = true;
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
			if (para.marked[w] == false && is_in_unorder[w] == true) {
				dfs(w, para, cur_scc);
			}
		}
	}

	void print_scc(const std::string & dir) {
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


	~static_kosaraju_tc(){	}
};
#endif

