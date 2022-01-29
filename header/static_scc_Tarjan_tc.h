#ifndef STATIC_TOPO_TC_H_
#define STATIC_TOPO_TC_H_

#include <string>
#include <vector>
#include <fstream>
#include <iostream> 
#include <algorithm>
#include <time.h> 
#include <iomanip>
#include <algorithm>
#include "boost/graph/adjacency_list.hpp"

#include "scc_para.h"
#include "disjoint_set.h"
#include "utils.h"

template <typename T > class static_tarjan_tc;

template <typename T>
class static_tarjan_tc {
public:
	CG cg;																		// citation graph.
	std::vector<T> n2v;															// topo value of each node. 
	std::vector<unsigned int> n2y;												// publish year of each paper.
	unsigned int np;															// number of paper.
	disjoint_set ds;															// disjoint set to represent the scc. 
	unsigned int max_year = 0, min_year = 9999;									// max min year of the citation graph.
	std::vector<bool> is_in_unorder;											// identify if the node x is in the unorder_graph. 
	std::vector<int> is_in_same_year;											// to check the node is in same_year_subgraph.  
	std::vector<unsigned int> max_ccs;											// the the size of cc of subgraphs. 	

	std::string pyear;															// file name of pyear. 
	std::string pref;															// file name of pp.
	std::string self_name;														// the name of this alg. 
	std::vector<std::pair<unsigned int, unsigned int>> normal_edgelist;			// all edges of the citation graph.
	std::vector<std::pair<unsigned int, unsigned int>> unorder_edgelist;
	std::vector<std::vector<std::pair<unsigned int, unsigned int> >>same_year_edgelist;
	std::vector<std::vector<unsigned int>> topo_scc;							// all of the sccs. 
	scc_para unorder_para, same_year_para;										// tmp parameter of Pearce2 

	std::vector<bool> act_nodes, act_edges;									// active_nodes 
	std::vector<int> degrees;

	

	static_tarjan_tc() {}
	static_tarjan_tc(unsigned int _np, const std::string &_pyear, const std::string &_pref):cg(_np), np(_np), pyear(_pyear), pref(_pref), ds(_np), unorder_para(_np), same_year_para(_np) {
		max_year = 0;
		min_year = 9999;
		n2y = std::vector<unsigned int>(np);
		n2v = std::vector<T>(np);
		degrees = std::vector <int>(np, 0);
		is_in_unorder = std::vector<bool>(np, false);
		is_in_same_year = std::vector<int>(np, 0);
		topo_scc.reserve(np);
		self_name = "Static-Topo-TC-";
		myres.set_static_name(self_name);
	}

	


	void topo_sort() {
		my_timer mt;
		for (std::vector<std::pair<unsigned int, unsigned int>>::const_iterator it = unorder_edgelist.cbegin(); it != unorder_edgelist.cend(); it++) {
			if (unorder_para.num[(*it).second] == 0 && act_nodes[it->second] == true) {
				strong_connected(it->second, is_in_unorder, unorder_para, topo_scc, cg);
			}
		}
		max_ccs.push_back(unorder_para.cc_idx - 1);

		std::pair<unsigned int, unsigned int> edge;
		std::vector<std::pair<unsigned int, unsigned int>> cur_year_edges;
		unsigned int y = 0;
		for (int i = 0; i < same_year_edgelist.size(); i++) {
			cur_year_edges = same_year_edgelist[i];
			for (int j = 0; j < cur_year_edges.size(); j++) {
				edge = cur_year_edges[j];
				if (unorder_para.num[edge.second] == 0 && same_year_para.num[edge.second] == 0 && act_nodes[edge.second] == true) {
					y = n2y[edge.second];
					strong_connected2(edge.second, is_in_same_year, same_year_para, unorder_para, topo_scc, cg, y, n2y);
				}
			}
			max_ccs.push_back(same_year_para.cc_idx - 1);
			same_year_para.cc_idx = 1;
		}
		myres.set_scc_info(mt.elapsed());
		record_cur_per_info(myres);
	}

	void eliminate_static_scc() {
		my_timer mt;
		std::vector<std::vector<unsigned int>> scc_nodes;
		get_scc_nodes(scc_nodes);
		// condense_sccs(scc_nodes, ds, cg);
		// update_property_map();
		myres.set_eliminate_scc_time(mt.elapsed());
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

	void update_property_map() {
		my_timer mt;
		boost::graph_traits<CG>::vertex_iterator vi, vi_end;
		for (std::tie(vi, vi_end) = vertices(cg); vi != vi_end; ++vi) {
			if (unorder_para.cc[*vi] != 0) {
				n2v[*vi] = unorder_para.cc[*vi];
			}
			else {
				n2v[*vi] = same_year_para.cc[*vi];
			}
		}
	}

	void print_scc(const std::string & dir) {
		myres.gen_sccs(topo_scc);
	}

	void read_citation_graph() {
		iread_pyear(pyear, n2y, min_year, max_year);
		same_year_edgelist = std::vector<std::vector<std::pair<unsigned int, unsigned int> >>(max_year - min_year + 1);
		

		if (myres.scale_factor < 1.0 ) {
			std::vector<std::pair<unsigned int, unsigned int>> edgelist;
			iread_edgelist(pref, np, edgelist, n2y, degrees);
			act_nodes = std::vector<bool>(np, false);
			act_edges = std::vector<bool>(edgelist.size(), false);
			gen_new_edges(edgelist, degrees, act_nodes, act_edges);
			for (int i = 0; i < edgelist.size(); i ++) {
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


	~static_tarjan_tc() {	}

};

#endif

