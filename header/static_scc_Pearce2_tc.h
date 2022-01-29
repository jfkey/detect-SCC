#ifndef STATIC_SCC_PEARCE2_TC_HPP
#define STATIC_SCC_PEARCE2_TC_HPP

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


template <typename T > class static_pearce2_tc;

template <typename T >
class static_pearce2_tc {

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
	ppara para;																	// tmp parameter of Pearce2 

	std::vector<bool> act_nodes, act_edges;									// active_nodes 
	std::vector<int> degrees;


	static_pearce2_tc() {}
	static_pearce2_tc(unsigned int _np, const std::string &_pyear, const std::string &_pref):cg(_np), np(_np), pyear(_pyear), pref(_pref), ds(_np), para(_np) {
		max_year = 0; 
		min_year = 9999;
		n2y = std::vector<unsigned int>(np);
		n2v = std::vector<T>(np);
		degrees = std::vector <int>(np, 0);
		is_in_unorder = std::vector<bool>(np, false);
		is_in_same_year = std::vector<int>(np, 0);
		topo_scc.reserve(np);
		self_name = "Static-Pearce2-TC-";
		myres.set_static_name(self_name);
	}

	void eliminate_static_scc() {
		// debug. generate same scc topo order. 
		//topo_scc.clear();
		//recalc_order();
		my_timer mt;
		std::vector<std::vector<unsigned int>> scc_nodes;
		get_scc_nodes(scc_nodes);
		condense_sccs(scc_nodes, ds, cg);
		update_property_map();
		myres.set_eliminate_scc_time(mt.elapsed());
	}

	void recalc_order() {
		para = ppara(np);
		boost::graph_traits<CG>::vertex_iterator vi, vi_end;
		for (std::tie(vi, vi_end) = boost::vertices(cg); vi != vi_end; ++vi) {
			if (para.rindex[*vi] == 0) {
				pearce_search2_t(*vi, para);
			}
		}
	}

	void pearce_search2_t(CG::vertex_descriptor v, ppara & para) {
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
				pearce_search2_t(w, para);
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
		for (int i = 0; i < np; i++) {
			if (para.rindex[i] != 0) {
				n2v[i] = np - para.rindex[i];
			}
			else {
				n2v[i] = 0;
			}
		}
	
	}
	
	// pearce2
	void topo_sort_pearce2() {
		// reset 
		topo_scc.clear();
		topo_scc.reserve(np);
		ppara para = ppara(np);
		my_timer mt;
		boost::graph_traits<CG>::vertex_iterator vi, vi_end;
		for (std::tie(vi, vi_end) = boost::vertices(cg); vi != vi_end; ++vi) {
			if (para.rindex[*vi] == 0) {
				pearce_search2_simple(*vi, para);
			}
		}
		myres.dn_eval_time = mt.elapsed();
		
	}

	void print_eo2n() {
		string path = "D:/datasets/Aminer/final/"; 
		string name = "Eo2n.txt"; 
		ofstream ofs (path + name);
		for (int i = 0; i < unorder_edgelist.size(); i++) {
			ofs << unorder_edgelist[i].first << "\t" << unorder_edgelist[i].second << std::endl; 
		}
		ofs.close();
	}
	// tc_pearce2
	void topo_sort() {
		para = ppara(np);
		int is = 0;
		max_ccs.clear();
		topo_scc.clear();
		topo_scc.reserve(np);
		my_timer mt;
		for (std::vector<std::pair<unsigned int, unsigned int>>::const_iterator it = unorder_edgelist.cbegin(); it != unorder_edgelist.cend(); it++) {
			is++;
			if (para.rindex[it->second] == 0 && act_nodes[it->second] == true) {
				pearce_search2(it->second, is_in_unorder, para, topo_scc, cg);
			}
		}
		max_ccs.push_back(np - para.c);
		// remain to be rewrite. 
		para.c = np - 1;
		para.index_p = 1;
		std::pair<unsigned int, unsigned int> edge;
		std::vector<std::pair<unsigned int, unsigned int>> cur_year_edges;
		unsigned int y = 0;
		for (int i = 0; i < same_year_edgelist.size(); i++) {
			cur_year_edges = same_year_edgelist[i];
			for (int j = 0; j < cur_year_edges.size(); j++) {
				edge = cur_year_edges[j];
				if(is_in_unorder[edge.second] == false && para.rindex[edge.second] == 0 && act_nodes[edge.second] == true) {
					pearce_search2(edge.second, is_in_unorder, is_in_same_year, para, topo_scc, cg, n2y, n2y[edge.second]);
				}
			}
			max_ccs.push_back(np - para.c);
			para.c = np - 1;
			para.index_p = 1;
		}
		myres.set_scc_info(mt.elapsed());
		
	}

	void print_scc(const std::string & dir) {
		myres.gen_sccs(topo_scc);
		//ofstream ofs(dir + "Eo2ng4.txt");
		//for (int i = 0; i < unorder_edgelist.size();i ++) {
		//	int s = unorder_edgelist[i].first; 
		//	int t = unorder_edgelist[i].second;
		//	if ((n2y[t] - n2y[s])>4) {
		//		ofs << s << "\t" << t << std::endl; 
		//	}
		//}
		//ofs.close();
		
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

	void add_inc_edges(std::vector<std::pair<unsigned int, unsigned int>> & inc_edges) {
		std::vector<std::pair<unsigned int, unsigned int>>::const_iterator it;
		
		unsigned int s, t; 
		for (it = inc_edges.cbegin(); it != inc_edges.cend(); it++) {
			s = it->first; 
			t = it->second; 
			if (n2y[s] < n2y[t]) {
				unorder_edgelist.push_back(std::make_pair(s, t));
			}
			else if (n2y[s] == n2y[t]) {
				same_year_edgelist[n2y[s] - min_year].push_back(std::make_pair(s, t));
			}
			
			boost::add_edge((*it).first, (*it).second, cg);
		}
	}

	void add_inc_edges(std::vector<std::pair<unsigned int, unsigned int>> & inc_edges, const int start_index, const int end_index) {
		unsigned int s, t;
		for (int i = start_index; i < end_index; i ++) {
			s = inc_edges[i].first;
			t = inc_edges[i].second;
			if (n2y[s] < n2y[t]) {
				unorder_edgelist.push_back(std::make_pair(s, t));
			}
			else if (n2y[s] == n2y[t]) {
				same_year_edgelist[n2y[s] - min_year].push_back(std::make_pair(s, t));
			}

			boost::add_edge(inc_edges[i].first, inc_edges[i].second, cg);
		}
	}

	//  used for pearce2
	void pearce_search2_simple(CG::vertex_descriptor v, ppara & para) {
		boost::graph_traits<CG> ::out_edge_iterator ai, ai_end;
		int root = true;
		para.rindex[v] = para.index_p;
		para.index_p++;
		unsigned int w;

		for (std::tie(ai, ai_end) = out_edges(v, cg); ai != ai_end; ++ai) {
			w = boost::target(*ai, cg);
			if (para.rindex[w] == 0) {
				pearce_search2_simple(w, para);
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

};
#endif

