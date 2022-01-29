#ifndef STATIC_SCC_PEARCE2_TC_OL_H_
#define STATIC_SCC_PEARCE2_TC_OL_H_

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



template <typename T > class static_pearce2_tc_ol;

template <typename T >
class static_pearce2_tc_ol {
public:
	CG cg;																		// citation graph.
	std::vector< ahrsz_priority_value <T> > n2v;								// topo value of each node. 
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
	std::vector<T *> pspaces;


	static_pearce2_tc_ol() {}
	static_pearce2_tc_ol(unsigned int _np, const std::string &_pyear, const std::string &_pref) :cg(_np), np(_np), pyear(_pyear), pref(_pref), ds(_np), para(_np) {
		max_year = 0;
		min_year = 9999;
		n2y = std::vector<unsigned int>(np);
		n2v = std::vector<ahrsz_priority_value<T> >(np);
		is_in_unorder = std::vector<bool>(np, false);
		is_in_same_year = std::vector<int>(np, 0);
		topo_scc.reserve(np);
		self_name = "Static-Pearce2-TC-OL-";
		myres.set_static_name(self_name);
	}


	void eliminate_static_scc() {
		// debug. generate same scc
		//topo_scc.clear();
		//recalc_order();
		my_timer mt;
		std::vector<std::vector<unsigned int>> scc_nodes;
		get_scc_nodes(scc_nodes);
		condense_sccs(scc_nodes, ds, cg);
		update_property_map();
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
		assert((max_year - min_year + 2) == pspaces.size());
		unsigned int v; int idx = 0;
		std::vector<std::vector<unsigned int>>::iterator it;

		std::vector<typename T::iterator> pits;
		std::for_each(pspaces.begin(), pspaces.end(), [&pits](T* a) {
			pits.push_back(a->begin());
		});
		for (it = topo_scc.begin(); it != topo_scc.end(); it++) {
			v = ds.find((*it).front());
			if (is_in_unorder[v] == true) {
				n2v[v] = ahrsz_priority_value<T>(pits[0], *(pspaces[0]));
				pits[0] ++;
			}
			else if (is_in_unorder[v] == false && is_in_same_year[v] != 0) {
				idx = n2y[v] - min_year + 1;
				n2v[v] = ahrsz_priority_value<T>(pits[idx], *(pspaces[idx]));
				pits[idx] ++;
			}
			else {
				//assert(0 == 1 && "Error Scc. in static_scc_Pearce2_tc_ol.h");
			}
		}

		// // test for test_datasets (V, E) = (28, 32)
		//std::cout << "(n2v[9] == n2v[3]) : " << (n2v[ds.find(9)] == n2v[ds.find(3)]) << std::endl;
		//std::cout << "(n2v[16] == n2v[22]) : " << (n2v[ds.find(16)] == n2v[ds.find(22)]) << std::endl;
		//std::cout << "(n2v[21] > n2v[23]) : " << (n2v[ds.find(21)] > n2v[ds.find(23)]) << std::endl;
		//std::cout << "(n2v[24] > n2v[25]) : " << (n2v[ds.find(24)] > n2v[ds.find(25)]) << std::endl;
	}


	void recalc_order() {
		ppara para_t(np);
		boost::graph_traits<CG>::vertex_iterator vi, vi_end;
		for (std::tie(vi, vi_end) = boost::vertices(cg); vi != vi_end; ++vi) {
			if (para_t.rindex[*vi] == 0) {
				pearce_search2_t(*vi, para_t);
			}
		}
	}

	void pearce_search2_t(CG::vertex_descriptor v, ppara & para) {
		//myres.access_nodes++;
		boost::graph_traits<CG> ::out_edge_iterator ai, ai_end;
		int root = true;
		para.rindex[v] = para.index_p;
		para.index_p++;
		unsigned int w;

		for (std::tie(ai, ai_end) = out_edges(v, cg); ai != ai_end; ++ai) {
			w = boost::target(*ai, cg);
			//myres.access_edges++;
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


	void topo_sort() {
		int is = 0;
		my_timer mt;
		for (std::vector<std::pair<unsigned int, unsigned int>>::const_iterator it = unorder_edgelist.cbegin(); it != unorder_edgelist.cend(); it++) {
			is++;
			if (para.rindex[it->second] == 0) {
				pearce_search2(it->second, is_in_unorder, para, topo_scc, cg);
			}
		}
		max_ccs.push_back(np - para.c);
		para.c = np - 1;
		para.index_p = 1;
		std::pair<unsigned int, unsigned int> edge;
		std::vector<std::pair<unsigned int, unsigned int>> cur_year_edges;
		unsigned int y = 0;
		for (int i = 0; i < same_year_edgelist.size(); i++) {
			cur_year_edges = same_year_edgelist[i];
			for (int j = 0; j < cur_year_edges.size(); j++) {
				edge = cur_year_edges[j];
				if (is_in_unorder[edge.second] == false && para.rindex[edge.second] == 0) { // para.rindex[v] == 0, i.e., v is not visited. 
					pearce_search2(edge.second, is_in_unorder, is_in_same_year, para, topo_scc, cg, n2y, n2y[edge.second]);
				}
			}
			max_ccs.push_back(np - para.c);
			para.c = np - 1;
			para.index_p = 1;
		}
		myres.set_scc_info(mt.elapsed());

		// create pspaces to store the topo value .
		std::for_each(max_ccs.begin(), max_ccs.end(), [this](int const &a) {
			this->pspaces.push_back(new T(a + 1));
		});
		record_cur_per_info(myres);
	}

	void print_scc(const std::string & dir) {
		myres.gen_sccs(topo_scc);
	}

	void read_citation_graph() {
		iread_pyear(pyear, n2y, min_year, max_year);
		same_year_edgelist = std::vector<std::vector<std::pair<unsigned int, unsigned int> >>(max_year - min_year + 1);
		iread_edgelist(pref, np, normal_edgelist, unorder_edgelist, same_year_edgelist, n2y, min_year);
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

};
#endif

