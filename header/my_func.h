#ifndef MY_FUNC_HPP
#define MY_FUNC_HPP
#include <iostream>
#include <vector>
#include <assert.h>
#include <fstream>
#include <time.h>
#include <string>
#include <map> 
#include <set>
#include <numeric>
#include <time.h> 
#include <iomanip>

#include "ordered_slist.h"
#include "ordered_slist2.h"
#include "boost/graph/adjacency_list.hpp"
#include "disjoint_set.h" 
#include "scc_para.h"
#include "my_graph.h"

#include "utils.h"
#define MAX_YEAR_GAP 10

using namespace std;

inline string get_inc_f(const std::string &insert_file, const  std::string & insert_type) {
	string inc_f;
	if (insert_file.length() == 0) {
		if (insert_type == SEQ_INSERT) {
			inc_f = "inc_s" + to_string(DELTAG_SIZE) + ".txt";
		}
		else if (insert_type == RAND_INSERT) {
			inc_f = "inc_r" + to_string(DELTAG_SIZE) + ".txt";
		} if (insert_type == COMBINE_INSERT) {
			inc_f = "inc_c" + to_string(DELTAG_SIZE) + ".txt";
		}
	}
	else {
		inc_f = insert_file;
	}
	return inc_f;
}

//inline string get_inc_f(const std::string &insert_file, int deltag_size, const  std::string & insert_type) {
//	string inc_f;
//	if (insert_file.length() == 0) {
//		if (insert_type == SEQ_INSERT) {
//			inc_f = "inc_s" + to_string(percentage) + ".txt";
//		}
//		else if (insert_type == RAND_INSERT) {
//			inc_f = "inc_r" + to_string(percentage) + ".txt";
//		}
//		else if () {
//
//		}
//	}
//	else {
//		inc_f = insert_file;
//	}
//	return inc_f;
//}


inline int get_alg_type(unsigned int s, unsigned int t, std::vector<int>& is_in_same_year, std::vector<bool> &is_in_unorder) {
	// sperate the original graph into three parts, e.g., Gm, Gsi, Gr
	// Gm, Gsi, Gr represents the unorder reference subgraph, the same year subgraph and not in G1 or Gsi, respectively.  
	if (is_in_unorder[s] == true && is_in_unorder[t] == true) {
		// s, t in Gm  
		return 1;
	}
	else if (is_in_unorder[s] == false && is_in_unorder[t] == true) {
		// s, t in Gsi & Gr, Gm  // algorithm 1.
		return 2;
	}
	else if (is_in_unorder[s] == true && is_in_unorder[t] == false) {
		// s, t in Gm, Gsi & Gr, //algorithm 2 
		return 3;
	}
	else if (is_in_unorder[s] == false && is_in_unorder[t] == false && is_in_same_year[s] != 0 && is_in_same_year[t] != 0) {
		// s, t in Gsi, Gsi. 
		return 4;	// algorithm 3
	}
	else if (is_in_unorder[t] == false && is_in_unorder[t] == false && (is_in_same_year[s] != 0 && is_in_same_year[t] == 0 || is_in_same_year[t] != 0 && is_in_same_year[s] == 0)) {
		// s, t in Gsi, Gr. OR  s, t in Gr, Gsi. // algorithm 4. 
		return 5;
	}
	else if (is_in_unorder[t] == false && is_in_unorder[t] == false && is_in_same_year[s] == 0 && is_in_same_year[t] == 0) {
		// s, t in Gr, Gr. // algorithm 5. 
		return 6;
	}
	else {
		assert(0 == 1 && "algorithm type error");
		return 7;
	}
}

inline void gen_new_edges(std::vector<std::pair<unsigned int, unsigned int>> &edgelist, std::vector<int> &degrees, std::vector<bool> &act_nodes, std::vector<bool> &act_edges) {
	vector<size_t> idx(degrees.size());
	iota(idx.begin(), idx.end(), 0);
	sort(idx.begin(), idx.end(), [&degrees](size_t i1, size_t i2) {return degrees[i1] > degrees[i2]; });

	std::vector<int> act_edge_idx = std::vector<int>();
	int act_node_size = degrees.size() * myres.scale_factor;
	int act_edge_size = edgelist.size() * myres.scale_factor;
	for (int i = 0; i < act_node_size; i ++) {
		act_nodes[ idx[i] ] = true; 
	}
	int cnt = 0; 
	for (int i = 0; i < edgelist.size(); i++ ) {
		if (act_nodes[edgelist[i].first] == true && act_nodes[edgelist[i].second] == true) {
			act_edges[i] = true; 
			act_edge_idx.push_back(i);
			cnt++; 
		}
		
	}
	if (cnt < (act_edge_size * 0.9 )) {
		std::cout << "The ratio of edges is less than " << myres.scale_factor << std::endl;
	} 
	if (cnt > act_edge_size) {
		std::srand(time(0));   
		std::random_shuffle(act_edge_idx.begin(), act_edge_idx.end());
		for (int i = 0; i < (cnt - act_edge_size); i ++ ) {
			act_edges[act_edge_idx[i]] = false; 
		}

	}
}


inline void assert_same_sccs(std::vector<std::pair<int, int>> &res, std::vector<std::pair<int, int>> &stand) {
	assert(res.size() == stand.size() && " error sccs sizes. ");
	for (int i = 0; i < res.size(); i++) {
		assert(res[i].first == stand[i].first && res[i].second == stand[i].second && "error sccs.");
	}
}

inline void print_scc_infile(const string & dir, std::vector<std::vector<unsigned int>> &scc_nodes, const string & alg_name, const int & inc_size) {
	std::vector<std::pair<int, int>> res_scc;

	assert(!scc_nodes.empty() && "The scc beign printed is EMPTY.");
	int np = scc_nodes.size();
	std::vector<std::vector<unsigned int>>::iterator it;
	std::vector<int>cc_size_cnt(np + 1, 0);
	for (it = scc_nodes.begin(); it != scc_nodes.end(); it++) {
		if ((*it).size() > 0) {
			cc_size_cnt[(*it).size()]++;
		}
	}
	std::ofstream ofs(dir + alg_name + "inc-" + to_string(inc_size) + ".txt");
	for (int i = 1; i < np + 1; i++) {
		if (cc_size_cnt[i] > 0) {
			ofs << i << "\t" << cc_size_cnt[i] << std::endl;
			res_scc.push_back(std::make_pair(i, cc_size_cnt[i]));
		}
	}
	ofs.close();
	std::vector<std::pair<int, int>> stand_scc;
	//stand_scc.push_back(std::make_pair(1, 11));
	//stand_scc.push_back(std::make_pair(2, 2));
	//stand_scc.push_back(std::make_pair(4, 1));
	//stand_scc.push_back(std::make_pair(9, 1));

	//stand_scc.push_back(std::make_pair(1, 10));
	//stand_scc.push_back(std::make_pair(2, 2));
	//stand_scc.push_back(std::make_pair(5, 1));
	//stand_scc.push_back(std::make_pair(9, 1));
	//if (inc_size != 0) assert_same_sccs(res_scc, stand_scc);
}


template <typename T, typename N2I>
inline void print_state_ol(T &g, typename boost::property_map<T, N2I>::type pmap) {
	std::vector<unsigned int> topval;
	for (int i = 0; i < boost::num_vertices(g); i++) {
		//std::cout << "i: " << i << "\t" << pmap[i]() << std::endl;
		topval.push_back(pmap[i]());
	}
	vector<size_t> idx(topval.size());
	iota(idx.begin(), idx.end(), 0);
	sort(idx.begin(), idx.end(), [&topval](size_t i1, size_t i2) {return topval[i1] < topval[i2]; });
	for (int i = 0; i < idx.size(); i++) {
		std::cout << std::setw(2) << idx[i] << ", ";
	}
	std::cout << std::endl << std::endl;
}



// std::vector<std::vector<unsigned int>> scc_nodes;
inline void cout_scc(std::vector<std::vector<unsigned int>> &scc_nodes) {
	std::vector<std::vector<unsigned int>>::iterator it1;
	std::vector<unsigned int>::iterator it2;
	unsigned int idx = 0;
	std::cout << " =============scc nodes==================" << std::endl;
	// 
	for (it1 = scc_nodes.begin(); it1 != scc_nodes.end(); it1++) {
		std::cout << idx++ << "\t";
		for (it2 = (*it1).begin(); it2 != (*it1).end(); it2++) {
			std::cout << *it2 << "\t";
		}
		std::cout << std::endl;
	}

	std::cout << " =============end print==================" << std::endl;

}

inline int out_degree_scc(unsigned int t, CG & graph, disjoint_set & ds) {
	int out_d = 0;
	std::vector<unsigned int>::iterator it; 
	for (it = ds.scc_nodes[t].begin(); it != ds.scc_nodes[t].end(); it++) {
		out_d += boost::out_degree(*it, graph);
	}
	return out_d;
}

inline int in_degree_scc(unsigned int t, CG & graph, disjoint_set & ds) {
	int in_d = 0;
	std::vector<unsigned int>::iterator it;
	for (it = ds.scc_nodes[t].begin(); it != ds.scc_nodes[t].end(); it++) {
		in_d += boost::in_degree(*it, graph);
	}
	return in_d;
}





//   inter-scc do NOT include any edges.
inline void condense_per_scc(unsigned int &s, std::vector<unsigned int>& s_scc, CG &graph, disjoint_set & ds) {// the node s and its scc. 
	boost::graph_traits<CG> ::out_edge_iterator oi, oi_end;
	boost::graph_traits<CG> ::in_edge_iterator  ii, ii_end;
	std::vector<unsigned int>::iterator scc_it;
	std::vector<std::pair<unsigned int, unsigned int>> rm_edges; 
	std::vector<unsigned int> in_nodes, out_nodes; 
	unsigned int w;
	for (scc_it = s_scc.begin(); scc_it != s_scc.end(); scc_it++) {
		for (std::tie(oi, oi_end) = boost::out_edges(*scc_it, graph); oi != oi_end; ++oi) {
			w = boost::target(*oi, graph);
			if (ds.find(s) == ds.find(w)) {
				rm_edges.push_back(std::make_pair(*scc_it, w));
			}
		}
	}

	for (int i = 0; i < rm_edges.size(); i++) {
		boost::remove_edge(rm_edges[i].first, rm_edges[i].second, graph);
	}
	rm_edges.clear();
	for (scc_it = s_scc.begin(); scc_it != s_scc.end(); scc_it++) {
		for (std::tie(ii, ii_end) = boost::in_edges(*scc_it, graph); ii != ii_end; ++ii) {
			w = boost::source(*ii, graph);
			if (ds.find(s) == ds.find(w)) {
				rm_edges.push_back(std::make_pair(w, *scc_it));
			}
		}
	}

	for (int i = 0; i < rm_edges.size(); i++) {
		boost::remove_edge(rm_edges[i].first, rm_edges[i].second, graph);
	}

}

inline void condense_sccs(std::vector<std::vector<unsigned int>> &scc_nodes, disjoint_set & ds, CG &graph) {
	std::vector<std::vector<unsigned int>>::iterator it1;
	std::vector<unsigned int>::iterator it2;
	unsigned int s, t;
	my_timer mt = my_timer();

	// we use a disjoint set to represent the vertex partition, e.g., one strong connected component in one group.   canonical node
	for (it1 = scc_nodes.begin(); it1 != scc_nodes.end(); it1++) {
		assert((*it1).size() > 1);
		s = (*it1).front();
		for (it2 = (*it1).begin(); it2 != (*it1).end(); it2++) {
			t = *it2;
			ds.join(s, t);
		}
		ds.collapse_scc(ds.find(s), *it1);
	}
	std::cout << " ds.collapse_scc(): " << mt.elapsed() << " sec." << std::endl;
	mt = my_timer();
	// collapse nodes and the canonical vertex of the component containing vertex v can be derive from find(v).
	for (it1 = scc_nodes.begin(); it1 != scc_nodes.end(); it1++) {
		//assert((*it1).size() > 1);		
		s = ds.find((*it1).front());
		condense_per_scc(s, *it1, graph, ds);
	}
	std::cout << " condense_per_scc(): " << mt.elapsed() << " sec." << std::endl;
}

// tarjan for original graph (utilized for Tarjan)
inline void strong_connected(CG::vertex_descriptor v, scc_para & para, std::vector<std::vector<unsigned int>> &topo_scc, CG &graph) {
	myres.access_nodes++;
	boost::graph_traits<CG> ::out_edge_iterator ai, ai_end;
	para.num[v] = para.lowlink[v] = para.time;
	para.time++;
	para.stack.push_back(v);
	para.in_stack[v] = true;
	int w;
	for (std::tie(ai, ai_end) = out_edges(v, graph); ai != ai_end; ++ai) {
		myres.access_edges++;
		w = boost::target(*ai, graph);
		if (para.num[w] == 0) {
			strong_connected(w, para, topo_scc, graph);
			if (para.lowlink[w] < para.lowlink[v]) para.lowlink[v] = para.lowlink[w];
		}
		else if (para.in_stack[w]) {
			if (para.num[w] < para.lowlink[v]) para.lowlink[v] = para.num[w];
		}
	}

	std::vector<unsigned int> tmp_component;
	if (para.lowlink[v] == para.num[v]) {
		unsigned int w = -1;
		do {
			w = para.stack.back(); para.stack.pop_back();
			para.cc[w] = para.cc_idx;
			para.in_stack[w] = false;
			tmp_component.push_back(w);
		} while (w != v);
		para.cc_idx++;
		topo_scc.push_back(tmp_component);
	}
}

// tarjan for unorder_subgraph graph on the original graph. (utilized for TC-Tarjan)
inline void strong_connected(CG::vertex_descriptor v, std::vector<bool> &is_in_unorder, scc_para & para,
	std::vector<std::vector<unsigned int>> &topo_scc, CG & graph) {
	myres.access_nodes++;
	is_in_unorder[v] = true;
	boost::graph_traits<CG> ::out_edge_iterator ai, ai_end;
	para.num[v] = para.lowlink[v] = para.time;
	para.time++;
	para.stack.push_back(v);
	para.in_stack[v] = true;
	int w;
	for (std::tie(ai, ai_end) = out_edges(v, graph); ai != ai_end; ++ai) {
		w = boost::target(*ai, graph);
		myres.access_edges++;
		if (para.num[w] == 0) {
			strong_connected(w, is_in_unorder, para, topo_scc, graph);
			if (para.lowlink[w] < para.lowlink[v]) para.lowlink[v] = para.lowlink[w];
		}
		else if (para.in_stack[w]) {
			if (para.num[w] < para.lowlink[v]) para.lowlink[v] = para.num[w];
		}
	}

	std::vector<unsigned int> tmp_component;
	if (para.lowlink[v] == para.num[v]) {
		unsigned int w = -1;
		do {
			w = para.stack.back(); para.stack.pop_back();
			para.cc[w] = para.cc_idx;
			para.in_stack[w] = false;
			tmp_component.push_back(w);
		} while (w != v);
		para.cc_idx++;
		topo_scc.push_back(tmp_component);
	}
}

// tarjan for same-year subgraph on the original graph. (utilized for TC-Tarjan)
inline void strong_connected2(CG::vertex_descriptor v, std::vector<int> &is_in_same_year, scc_para & para,
	scc_para &unorder_para, std::vector<std::vector<unsigned int>> &topo_scc, CG & graph, const int year, std::vector<unsigned int> &n2y ) {
	myres.access_nodes++;
	is_in_same_year[v] = year;
	boost::graph_traits<CG> ::out_edge_iterator ai, ai_end;
	para.num[v] = para.lowlink[v] = para.time;
	para.time++;
	para.stack.push_back(v);
	para.in_stack[v] = true;
	unsigned int w;
	for (std::tie(ai, ai_end) = boost::out_edges(v, graph); ai != ai_end; ++ai) {
		w = boost::target(*ai, graph);
		myres.access_edges++;
		if (para.num[w] == 0 && unorder_para.num[w] == 0 && n2y[w] == year) {
			strong_connected2(w, is_in_same_year, para, unorder_para, topo_scc, graph, year, n2y);
			if (para.lowlink[w] < para.lowlink[v]) para.lowlink[v] = para.lowlink[w];
		}
		else if (para.in_stack[w] && n2y[w] == year) {
			if (para.num[w] < para.lowlink[v]) para.lowlink[v] = para.num[w];
		}
	}

	std::vector<unsigned int> tmp_component;
	if (para.lowlink[v] == para.num[v]) {
		unsigned int w = -1;
		do {
			w = para.stack.back(); para.stack.pop_back();
			para.cc[w] = para.cc_idx;
			para.in_stack[w] = false;
			tmp_component.push_back(w);
		} while (w != v);
		para.cc_idx++;
		topo_scc.push_back(tmp_component);
	}
}


//  pearce_search2 for unorder_subgraph graph on the original graph 
inline void pearce_search2(CG::vertex_descriptor v, std::vector<bool> &is_in_unorder, ppara &para,
	std::vector<std::vector<unsigned int>> &topo_scc, CG &graph) {
	is_in_unorder[v] = true;
	myres.access_nodes++;
	boost::graph_traits<CG> ::out_edge_iterator ai, ai_end;
	int root = true;
	para.rindex[v] = para.index_p;
	para.index_p++;
	unsigned int w;

	for (std::tie(ai, ai_end) = out_edges(v, graph); ai != ai_end; ++ai) {
		w = boost::target(*ai, graph);
		myres.access_edges++;
		if (para.rindex[w] == 0) {
			pearce_search2(w, is_in_unorder, para, topo_scc, graph);
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


// pearce_search2 for same_year subgraph on the orginal graph 
inline void pearce_search2(CG::vertex_descriptor v, std::vector<bool> &is_in_unorder, std::vector<int> &is_in_same_year,
	ppara &para, std::vector<std::vector<unsigned int>> &topo_scc, CG &graph, std::vector<unsigned int> &n2y, unsigned int& year) {
	is_in_same_year[v] = year;
	myres.access_nodes++;
	boost::graph_traits<CG> ::out_edge_iterator ai, ai_end;
	int root = true;
	para.rindex[v] = para.index_p;
	para.index_p++;
	unsigned int w;

	for (std::tie(ai, ai_end) = out_edges(v, graph); ai != ai_end; ++ai) {
		w = boost::target(*ai, graph);
		myres.access_edges++;
		if (para.rindex[w] == 0 && is_in_unorder[w] == false && n2y[w] == year) {
			pearce_search2(w, is_in_unorder, is_in_same_year, para, topo_scc, graph, n2y, year);
		}

		if (para.rindex[w] < para.rindex[v] && is_in_unorder[w] == false && is_in_unorder[v] == false && n2y[w] == year) {
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

//  pearce_search2 for same_year subgraph graph.  generate same year subgraph.
template<typename T >
inline void pearce_search2(typename T::vertex_descriptor v, std::vector<bool>& is_in_unorder, std::vector<bool> &is_in_same_year,
	ppara &para, std::vector<std::vector<unsigned int>> &topo_scc, T &graph) {
	is_in_same_year[v] = true;
	myres.access_nodes++;
	typename boost::graph_traits<T> ::out_edge_iterator ai, ai_end;
	int root = true;
	para.rindex[v] = para.index_p;
	para.index_p++;
	unsigned int w;

	for (std::tie(ai, ai_end) = out_edges(v, graph); ai != ai_end; ++ai) {
		w = boost::target(*ai, graph);
		myres.access_edges++;
		if (para.rindex[w] == 0 && is_in_unorder[w] == false) {
			pearce_search2(w, is_in_unorder, is_in_same_year, para, topo_scc, graph);
		}
		// node v and node w are in the same_year_subgraph and not visited by unorder subgraphs.
		if (para.rindex[w] < para.rindex[v] && is_in_unorder[w] == false && is_in_unorder[v] == false) {
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



inline void iread_pyear(const std::string &pyear, std::vector<unsigned int> &n2y) {
	long timestamp = clock();
	std::string line;
	std::ifstream ifs(pyear);
	while (true) {
		getline(ifs, line);
		if (line.length() == 0)
			break;
		int pos = line.find('\t');
		n2y[stoi(line.substr(0, pos))] = stoi(line.substr(pos + 1));
	}
	ifs.close();
}

inline void iread_pyear(const std::string &pyear, std::vector<unsigned int> &n2y, unsigned int & min_year, unsigned int & max_year) {
	long timestamp = clock();
	std::string line;
	std::ifstream ifs(pyear);
	while (true)
	{
		getline(ifs, line);
		if (line.length() == 0)
			break;
		int pos = line.find('\t');
		int y = stoi(line.substr(pos + 1));
		if (y > max_year) max_year = y;
		if (y < min_year) min_year = y;
		n2y[stoi(line.substr(0, pos))] = y;
	}
	ifs.close();
}


inline void iread_edgelist(const std::string &pref, unsigned int np, std::vector<std::pair<unsigned int, unsigned int>> & edgelist, std::vector<unsigned int> &n2y) {
	edgelist.reserve(np << 2);
	std::ifstream ifs(pref);
	std::string line;
	while (true) {
		getline(ifs, line);
		if (line.length() == 0)
			break;
		int pos = line.find('\t');
		unsigned int s, t;
		s = stoi(line.substr(0, pos));
		t = stoi(line.substr(pos + 1));

		if (s != t) {
			if (n2y[s] < n2y[t]) {
				if (myres.simplify_data == false) {
					edgelist.push_back(std::make_pair(s, t));
				}
				else {
					if (n2y[t] - n2y[s] < MAX_YEAR_GAP) {
						edgelist.push_back(std::make_pair(s, t));
					}
				}
			}
			else {
				edgelist.push_back(std::make_pair(s, t));
			}
			
		}
	}
	ifs.close();
}



inline void iread_edgelist(const std::string &pref, unsigned int np, std::vector<std::pair<unsigned int, unsigned int>> & edgelist, std::vector<unsigned int> &n2y, std::vector<int> &degrees) {
	edgelist.reserve(np << 2);
	std::ifstream ifs(pref);
	std::string line;
	while (true) {
		getline(ifs, line);
		if (line.length() == 0)
			break;
		int pos = line.find('\t');
		unsigned int s, t;
		s = stoi(line.substr(0, pos));
		t = stoi(line.substr(pos + 1));

		if (s != t) {
			degrees[s] += 1;
			degrees[t] += 1;
			if (n2y[s] < n2y[t]) {
				if (myres.simplify_data == false) {
					edgelist.push_back(std::make_pair(s, t));
				}
				else {
					if (n2y[t] - n2y[s] < MAX_YEAR_GAP) {
						edgelist.push_back(std::make_pair(s, t));
					}
				}
			}
			else {
				edgelist.push_back(std::make_pair(s, t));
			}

		}
	}
	ifs.close();
}


inline void iread_edgelist(const std::string &pref, unsigned int np, std::vector<std::pair<unsigned int, unsigned int>> & normal_edgelist,
	std::vector<std::pair<unsigned int, unsigned int>> & unorder_edgelist, std::vector<std::vector<std::pair<unsigned int, unsigned int> >> & same_year_edgelist,
	std::vector<unsigned int> &n2y, const int &min_year) {
	normal_edgelist.reserve(np);
	unorder_edgelist.reserve(np);
	same_year_edgelist.reserve(np);

	std::ifstream ifs(pref);
	std::string line;
	while (true) {
		getline(ifs, line);
		if (line.length() == 0)
			break;
		int pos = line.find('\t');
		unsigned int s, t;
		s = stoi(line.substr(0, pos));
		t = stoi(line.substr(pos + 1));
		if (s != t) {
			if (n2y[s] < n2y[t]) {
				if (myres.simplify_data == false) {
					unorder_edgelist.push_back(std::make_pair(s, t));
				}
				else {
					if (n2y[t] - n2y[s] < MAX_YEAR_GAP) {
						unorder_edgelist.push_back(std::make_pair(s, t));
					}
				}
			}
			else if (n2y[s] == n2y[t]) {
				same_year_edgelist[n2y[s] - min_year].push_back(std::make_pair(s, t));
			}
			else {
				normal_edgelist.push_back(std::make_pair(s, t));
			}
		}
	}
	ifs.close();
}

//inline void iread_edgelist(const std::string &pref, unsigned int np, std::vector<std::pair<unsigned int, unsigned int>> & normal_edgelist,
//	std::vector<std::pair<unsigned int, unsigned int>> & unorder_edgelist, std::vector<std::vector<std::pair<unsigned int, unsigned int> >> & same_year_edgelist,
//	std::vector<unsigned int> &n2y, const int &min_year, const int & max_year, std::vector<std::pair<int, int>> &id2ids) {
//	normal_edgelist.reserve(np);
//	unorder_edgelist.reserve(np);
//	same_year_edgelist.reserve(np);
//
//
//	int i = 0, ui = 0, ni = 0;
//	std::vector<int> same_yearids = std::vector<int>(max_year - min_year + 1, 0);
//	std::ifstream ifs(pref);
//	std::string line;
//	while (true) {
//		getline(ifs, line);
//		if (line.length() == 0)
//			break;
//		int pos = line.find('\t');
//		unsigned int s, t;
//		s = stoi(line.substr(0, pos));
//		t = stoi(line.substr(pos + 1));
//		if (s != t) {
//			if (n2y[s] < n2y[t]) {
//				unorder_edgelist.push_back(std::make_pair(s, t));
//				id2ids[i] = std::make_pair((max_year - min_year) + 1, ui);
//				ui++;
//			}
//			else if (n2y[s] == n2y[t]) {
//				same_year_edgelist[n2y[s] - min_year].push_back(std::make_pair(s, t));
//				id2ids[i] = std::make_pair(n2y[s] - min_year, same_yearids[n2y[s] - min_year]);
//				same_yearids[n2y[s] - min_year] += 1;
//			}
//			else {
//				normal_edgelist.push_back(std::make_pair(s, t));
//				id2ids[i] = std::make_pair((max_year - min_year) + 2, ni);
//				ni++;
//			}
//		}
//		i++;
//	}
//	ifs.close();
//}


inline void iread_insert_edge(const std::string & fn, std::vector<std::pair<unsigned int, unsigned int>> &edgelist, std::vector<unsigned int> &n2y)
{
	std::ifstream ifs(fn);
	long start = clock();
	string line;
	while (true) {
		getline(ifs, line);
		if (line.length() == 0)
			break;
		int pos = line.find('\t');
		unsigned int s, t;
		s = stoi(line.substr(0, pos));
		t = stoi(line.substr(pos + 1));
		if (s != t) {
			if (n2y[s] < n2y[t]) {
				if (myres.simplify_data == false) {
					edgelist.push_back(std::make_pair(s, t));
				}
				else {
					if (n2y[t] - n2y[s] < MAX_YEAR_GAP) {
						edgelist.push_back(std::make_pair(s, t));
					}
				}
			}
			else {
				edgelist.push_back(std::make_pair(s, t));
			}
			
		}
	}
	ifs.close();
	//cout << "Load Insert Edges File Name: " << fn << " With Counts: " << edgelist.size() << " , in " << (clock() - start)* 0.001 << " sec. " << endl;
}

template<typename T, typename Comp>
class ahrsz_priority_comp {
private:
	//std::vector< ahrsz_priority_value <T> >  &n2v;
	T &n2v;
	Comp comp;

public:
	ahrsz_priority_comp(T & _n2v) : n2v(_n2v) {	}
	bool operator()(CG::vertex_descriptor l, CG::vertex_descriptor r) const {
		return comp(n2v[l], n2v[r]);
	}

};

//ahrsz_priority_comp<T, N2I, std::greater<ahrsz_priority_value<T> > >

//template<class T, class N2I, class Comp>
//class ahrsz_priority_comp {
//private:
//	typedef typename boost::property_map<T, N2I>::type priority_map;
//private:
//	Comp _comp;
//	priority_map _pmap; // efficiency  
//public:
//	ahrsz_priority_comp(priority_map pmap)
//		: _pmap(pmap) {
//	}
//
//	bool operator()(typename T::vertex_descriptor l,
//		typename T::vertex_descriptor r) const {
//		return _comp(_pmap[l], _pmap[r]);
//	}
//};

template<class PSPACE = ordered_slist<void> > //  
class ahrsz_priority_value {
private:
	typename PSPACE::iterator _iter;
	PSPACE  const *_ptr;
public:
	ahrsz_priority_value(void)
		: _iter(), _ptr(NULL) {
	}

	ahrsz_priority_value(typename PSPACE::iterator const &iter,
		PSPACE const &pspace)
		: _iter(iter), _ptr(&pspace) {
	}

	bool operator<(ahrsz_priority_value const &src) const {
		assert(_ptr != NULL);
		assert(_ptr == src._ptr);
		return _ptr->order_lt(_iter, src._iter);
	}

	//bool operator>(ahrsz_priority_value const &src) const {
	//	assert(_ptr != NULL);
	//	return _iter != src._iter && !(*this < src);
	//}

	// modify by Junfeng. 
	bool operator>(ahrsz_priority_value const &src) const {
		assert(_ptr != NULL);
		assert(_ptr == src._ptr);
		return _ptr->order_gt(_iter, src._iter);
	}

	bool operator!=(ahrsz_priority_value const &src) const {
		return _iter != src._iter;
	}

	uint32_t operator()() const {
		assert(_ptr != NULL);
		return _ptr->order(_iter);
	}

	bool operator==(ahrsz_priority_value const &src) const {
		// return _iter == src._iter; 
		assert(_ptr == src._ptr);
		return  _ptr->order_eq(_iter, src._iter); // modify by Junfeng. 
	}

	typename PSPACE::iterator &base(void) {
		return _iter;
	}

	ahrsz_priority_value const &operator++(int) {
		++_iter;
		return *this;
	}

	ahrsz_priority_value operator++(void) {
		ahrsz_priority_value r(*this);
		++_iter;
		return r;
	}
};

typedef enum { minus_infinity, plus_infinity } ahrsz_ext_priority_value_inf_t;

template<class PSPACE>
class ahrsz_ext_priority_value {
private:
	ahrsz_priority_value<PSPACE> _val;
	bool _minus_inf;
	bool _plus_inf;
public:
	ahrsz_ext_priority_value(ahrsz_ext_priority_value_inf_t f = minus_infinity)
		: _val(), _minus_inf(f == minus_infinity), _plus_inf(f == plus_infinity) {
	}

	ahrsz_ext_priority_value(ahrsz_priority_value<PSPACE> &e)
		: _val(e), _minus_inf(false), _plus_inf(false) {
	}

	ahrsz_ext_priority_value(typename PSPACE::iterator const &iter,
		PSPACE const &space)
		: _val(iter, space), _minus_inf(false), _plus_inf(false) {
	}

	bool minus_inf(void) const { return _minus_inf; }
	bool plus_inf(void) const { return _plus_inf; }

	bool operator<(ahrsz_ext_priority_value const &src) const {
		if (_minus_inf || src._minus_inf) {
			return _minus_inf && !src._minus_inf;
		}
		else if (_plus_inf || src._plus_inf) {
			return src._plus_inf && !_plus_inf;
		}
		return _val < src._val;
	}

	bool operator>=(ahrsz_ext_priority_value const &src) const {
		return !(*this < src);
	}

	bool operator==(ahrsz_ext_priority_value const &src) const {
		if (_minus_inf || src._minus_inf) {
			return _minus_inf == src._minus_inf;
		}
		else if (_plus_inf || src._plus_inf) {
			return _plus_inf == src._plus_inf;
		}
		return _val == src._val;
	}

	typename PSPACE::iterator &base(void) {
		return _val.base();
	}

	ahrsz_ext_priority_value const &operator++(int) {
		++_val;
		return *this;
	}

	ahrsz_ext_priority_value operator++(void) {
		ahrsz_ext_priority_value r(*this);
		++_val;
		return r;
	}
};


#endif