#include "dynamic_IGC_TC_scc.h"

template<class T>
dynamic_tc_igc_scc<T>::dynamic_tc_igc_scc() {}

template<class T>
dynamic_tc_igc_scc<T>::dynamic_tc_igc_scc(CG & _cg, disjoint_set & _ds, std::vector< T > & _n2v, std::vector<unsigned int> &_max_ccs,
	std::vector<bool> &_is_in_unorder, std::vector<int> &_is_in_same_year, std::vector<unsigned int> &_n2y, unsigned int _max_year, 
	unsigned int _min_year, std::vector<pair<unsigned int, unsigned int>> & _inc_edges, unsigned int _np)
	: cg(_cg), ds(_ds), n2v(_n2v), max_ccs(_max_ccs), is_in_same_year(_is_in_same_year), is_in_unorder(_is_in_unorder),
	n2y(_n2y), min_year(_min_year), max_year(_max_year), inc_edges(_inc_edges), np(_np)
{
	in_component = std::vector<bool>(np, false);
	visited = std::vector<bool>(np, false);
	std::cout << "Init Dynamic Time-Coupled IGC(PK) Algorithm" << std::endl;
	self_name = "Dynamic-TC-IGC-";
	myres.set_dynamic_name(self_name);
}

template<class T>
void dynamic_tc_igc_scc<T>::add_edges()
{
	my_timer mt; 
	int type = 0; 
	std::vector<pair<unsigned int, unsigned int>>::iterator it;
	for (it = inc_edges.begin(); it != inc_edges.end(); it++) {
		unsigned int s, t;
		s = ds.find(it->first);
		t = ds.find(it->second);
		if (s != t) {															// s == t, e.g., the edge in the same component. do not add edge. 
			my_timer mt = my_timer();
			boost::add_edge(it->first, it->second, cg);
			type = get_alg_type(s, t, is_in_same_year, is_in_unorder);

			if (type == 1) {
				add_edge_Gm(s, t);
			}
			else if (type == 2) {												// s, t in Gsi & Gr, Gm // algorithm1
				/*if (n2y[s] < n2y[t]) {
					add_edge_alg1(s, t);
				}*/
			}
			else if (type == 3) {												// s, t in Gm, Gsi & Gr, //algorithm 2 
				add_edge_alg2(s, t);
			}
			else if (type == 4 || type == 5 || type == 6) {						 // s in Gsi or Gr,  t in Gsi or Gr // algrithm 3
				add_edge_alg3(s, t);
			}
			else {
				assert(0 == 1 && "ERROR EDGE TYPE.");
			}
		}
	}

	myres.dn_eval_time = myres.search_time + myres.reassignment_time + myres.tc_comp_time;
	myres.dn_eval_time_AAN = mt.elapsed();
}

template<class T>
void dynamic_tc_igc_scc<T>::add_edges(int start_index, int end_index)
{
	myres.inc_edges += (end_index - start_index);
	my_timer mt;
	int type = 0;
	std::vector<pair<unsigned int, unsigned int>>::iterator it;
	for (int i = start_index; i < end_index; i++) {
		unsigned int s, t;
		s = ds.find(inc_edges[i].first);
		t = ds.find(inc_edges[i].second);
		if (s != t) {															// s == t, e.g., the edge in the same component. do not add edge. 
			my_timer mt = my_timer();
			boost::add_edge(inc_edges[i].first, inc_edges[i].second, cg);
			type = get_alg_type(s, t, is_in_same_year, is_in_unorder);

			if (type == 1) {
				add_edge_Gm(s, t);
			}
			else if (type == 2) {												// s, t in Gsi & Gr, Gm // algorithm1
				/*if (n2y[s] < n2y[t]) {
					add_edge_alg1(s, t);
				}*/
			}
			else if (type == 3) {												// s, t in Gm, Gsi & Gr, //algorithm 2 
				add_edge_alg2(s, t);
			}
			else if (type == 4 || type == 5 || type == 6) {						 // s in Gsi or Gr,  t in Gsi or Gr // algrithm 3
				add_edge_alg3(s, t);
			}
			else {
				assert(0 == 1 && "ERROR EDGE TYPE.");
			}
		}
	}
	myres.dn_eval_time = myres.search_time + myres.reassignment_time + myres.tc_comp_time;
	//myres.dn_eval_time_AAN = mt.elapsed();
	record_cur_per_info(myres);
}


template<class T>
void dynamic_tc_igc_scc<T>::add_edge_Gm(vd_t t, vd_t h)
{
	unsigned int tn2v(n2v[t]);
	unsigned int hn2v(n2v[h]);

	//tn2v == hn2v, add the edge(t,h) in the scc. does not affect the topological order.
	//tn2v > hn2v, the edge(t,h) is actual topological order.
	if (tn2v < hn2v) {
		myres.invalid_edges += 1;
		std::vector<unsigned int> reaching, reachable;
		bool flag = false;
		in_component[t] = true;
		my_timer mt;
		fwd_dfs_Gm(h, tn2v, reachable, flag);
		myres.search_time += mt.elapsed();

		if (flag == true) {		// there exists the cycle in the condensation graph. 
			for_each(reachable.begin(), reachable.end(), [this](const unsigned int & a) {
				visited[a] = false;
			});
			mt = my_timer(); 
			back_dfs_Gm(t, hn2v, reaching);
			myres.search_time += mt.elapsed();

			std::vector<unsigned int> cycle, reachable_tmp, reaching_tmp;
			std::vector<unsigned int>::iterator i, j;
			unsigned int s = 0;
			for (i = reachable.begin(); i != reachable.end(); ++i) {
				if (in_component[*i] == true) {
					cycle.push_back(*i);
				}
				else {
					reachable_tmp.push_back(*i);
				}
				visited[*i] = false;
			}
			for (i = reaching.begin(); i != reaching.end(); i++) {
				visited[*i] = false;
				if (in_component[*i] == false)
					reaching_tmp.push_back(*i);
			}
			mt = my_timer(); 
			assert(cycle.size() > 1 && "the number of vertices in the cycle must be greater than 1 Gm -> Gm.");
			s = cycle.front();
			for (i = cycle.begin(); i != cycle.end(); i++) {
				ds.join(s, *i);
				visited[*i] = false;			// unmark states. 
				in_component[*i] = false;
				in_component[s] = false;
			}
			s = ds.find(cycle.front());
			ds.collapse_scc(s, cycle);
			//condense_per_scc(s, cycle, g, g.ds);
			std::vector<unsigned int> cycle_nodes2;
			for (int i = 0; i < cycle.size(); i++) {
				cycle_nodes2.insert(cycle_nodes2.end(), ds.scc_nodes[cycle[i]].begin(), ds.scc_nodes[cycle[i]].end());
			}
			condense_per_scc(s, cycle_nodes2, cg, ds);
			myres.condensation_time += mt.elapsed();

			mt = my_timer();
			dynamic_tc_igc_com<std::vector<T>> pc(n2v);
			sort(reaching_tmp.begin(), reaching_tmp.end(), pc);
			sort(reachable_tmp.begin(), reachable_tmp.end(), pc);
			reaching_tmp.push_back(s);
			myres.aff_region += (reaching_tmp.size() + reachable_tmp.size());
			
			tc_scc_reorder(reachable_tmp, reaching_tmp);
			myres.reassignment_time += mt.elapsed();
			
		}
		else {
			my_timer mt; 
			back_dfs_Gm(t, hn2v, reaching);
			myres.search_time += mt.elapsed();
			mt = my_timer();
			dynamic_tc_igc_com<std::vector<T>> pc(n2v);
			sort(reaching.begin(), reaching.end(), pc);
			sort(reachable.begin(), reachable.end(), pc);
			myres.aff_region += (reachable.size() + reaching.size());
			tc_scc_reorder(reachable, reaching);
			myres.reassignment_time += mt.elapsed();

		}
	}
	in_component[t] = false;
}

template<class T>
void dynamic_tc_igc_scc<T>::add_edge_alg1(vd_t t, vd_t h)
{
	//my_timer mt; 
	//std::vector<unsigned int> reachable;
	//std::vector<unsigned int>::iterator it;
	//scan_GsGr(h, reachable);											// NOTE. tarverse the subgraph from h
	////scan_GsGr(t, reachable);											
	//// update notGu nodes info.
	//for (it = reachable.begin(); it != reachable.end(); it++) {
	//	is_in_unorder[*it] = true;
	//	n2v[*it] = max_ccs[0];
	//	max_ccs[0]++;
	//	visited[*it] = false;
	//}
	//myres.tc_comp_time += mt.elapsed();
}

// from Gm to Gsi & Gr, t in Gm and h in Gsi & Gr
template<class T>
void dynamic_tc_igc_scc<T>::add_edge_alg2(vd_t t, vd_t h)
{
	my_timer mt; 
	unsigned int tn2v(n2v[t]);
	std::vector<unsigned int> reachable2;
	std::vector<unsigned int>::iterator it, it_ds;
	scan_GsGr(h, reachable2);
	// update notGu nodes info.
	for (it = reachable2.begin(); it != reachable2.end(); it++) {
		is_in_unorder[*it] = true;
		n2v[*it] = max_ccs[0];
		max_ccs[0] ++;
		visited[*it] = false;
	}
	unsigned int hn2v(n2v[h]);
	myres.tc_comp_time += mt.elapsed();

	if (tn2v < hn2v) {
		mt = my_timer();
		std::vector<unsigned int> reaching, reachable;
		bool flag = false;
		in_component[t] = true;
		fwd_dfs_Gm(h, tn2v, reachable, flag);
		myres.tc_comp_time += mt.elapsed();
		if (flag == true) {		// there exists the cycle in the condensation graph. 
			mt = my_timer();
			for_each(reachable.begin(), reachable.end(), [this](const unsigned int & a) {
				visited[a] = false;
			});
			back_dfs_Gm(t, hn2v, reaching);
			myres.tc_comp_time += mt.elapsed();
			mt = my_timer();
			std::vector<unsigned int> cycle, reachable_tmp, reaching_tmp;
			std::vector<unsigned int>::iterator i, j;
			unsigned int s = 0;
			for (i = reachable.begin(); i != reachable.end(); ++i) {
				if (in_component[*i] == true) {
					cycle.push_back(*i);
				}
				else {
					reachable_tmp.push_back(*i);
				}
				visited[*i] = false;
			}
			for (i = reaching.begin(); i != reaching.end(); i++) {
				visited[*i] = false;
				if (in_component[*i] == false)
					reaching_tmp.push_back(*i);
			}
			
			assert(cycle.size() > 1 && "the number of vertices in the cycle must be greater than 1 Gm -> Gm.");
			s = cycle.front();
			for (i = cycle.begin(); i != cycle.end(); i++) {
				ds.join(s, *i);
				visited[*i] = false;			// unmark states. 
				in_component[*i] = false;
				in_component[s] = false;
			}
			s = ds.find(cycle.front());
			ds.collapse_scc(s, cycle);
			//condense_per_scc(s, cycle, g, g.ds);
			std::vector<unsigned int> cycle_nodes2;
			for (int i = 0; i < cycle.size(); i++) {
				cycle_nodes2.insert(cycle_nodes2.end(), ds.scc_nodes[cycle[i]].begin(), ds.scc_nodes[cycle[i]].end());
			}
			condense_per_scc(s, cycle_nodes2, cg, ds);
			myres.condensation_time += mt.elapsed();

			mt = my_timer();
			dynamic_tc_igc_com<std::vector<T>> pc(n2v);
			sort(reaching_tmp.begin(), reaching_tmp.end(), pc);
			sort(reachable_tmp.begin(), reachable_tmp.end(), pc);
			reaching_tmp.push_back(s);
			tc_scc_reorder(reachable_tmp, reaching_tmp );
			myres.tc_comp_time += mt.elapsed();
		}
		else {
			mt = my_timer();
			back_dfs_Gm (t, hn2v, reaching );
			dynamic_tc_igc_com<std::vector<T>> pc(n2v);
			sort(reaching.begin(), reaching.end(), pc);
			sort(reachable.begin(), reachable.end(), pc);
			tc_scc_reorder(reachable, reaching );
			myres.tc_comp_time += mt.elapsed();
		}
	}
	in_component[t] = false;

}

template<class T>
void dynamic_tc_igc_scc<T>::add_edge_alg3(vd_t t, vd_t h)
{
	my_timer mt; 
	if (n2y[t] == n2y[h]) {
		if (is_in_same_year[t] != 0 && is_in_same_year[h] == 0) {				// t in Gsi, h in Gr
			mt = my_timer();
			int yidx = n2y[h] - min_year + 1;
			unsigned int max_order = max_ccs[yidx];
			std::vector<unsigned int> reaching, reachable;
			n2v[h] = max_order;
			max_ccs[yidx] += 1;
			is_in_same_year[h] = n2y[t];
			myres.tc_comp_time += mt.elapsed();

			/*reachable.push_back(h);
			back_dfs_Gsi(t, max_order, reaching, n2y[h]);
			dynamic_tc_igc_com<std::vector<T>> pc(n2v);
			sort(reaching.begin(), reaching.end(), pc);
			sort(reachable.begin(), reachable.end(), pc);
			tc_scc_reorder(reachable, reaching);*/

			add_edge_Gsi(t, h);
		}
		else if (is_in_same_year[t] == 0 && is_in_same_year[h] != 0) {			// t in Gr, h in Gsi
		}
		else if (is_in_same_year[t] == 0 && is_in_same_year[h] == 0) {			// t in Gr, h in Gr
			mt = my_timer();
			is_in_same_year[h] = n2y[h];
			int yidx = n2y[h] - min_year + 1;
			n2v[h] = max_ccs[yidx];
			max_ccs[yidx] += 1;
			myres.tc_comp_time += mt.elapsed();
		}
		else if (is_in_same_year[t] != 0 && is_in_same_year[h] != 0) {// t, h in Gsi, Gsi
			// apply dynamic algorithm in Gsi. 
			add_edge_Gsi(t, h);
		}
		else {
			assert(0 == 1 && "Type Error. t, h NOT in Gsi OR Gr.");
		}

	}
	else if (n2y[t] < n2y[h]) {	 // e.g., t in 2017, h in 2010. The cycle must be in Gsi & Gr. 
		mt = my_timer();
		std::vector<unsigned int> reachable, cycle;
		std::vector<unsigned int>::iterator it, it_ds, i;
		bool flag = false; unsigned int in_scc_order = 0;
		/*in_component[t] = true;
		scan_GsGr_without_topo(t, t, reachable, flag);*/
		in_component[h] = true;
		scan_GsGr_without_topo(h, t, reachable, flag);
		myres.tc_comp_time += mt.elapsed();

		if (flag == true) {
			mt = my_timer();
			// update notGu nodes info.  
			for (it = reachable.begin(); it != reachable.end(); it++) {
				if (in_component[*it] == false) {
					is_in_unorder[*it] = true;
					n2v[*it] = max_ccs[0];
					max_ccs[0] += 1;
				}
				else {
					cycle.push_back(*it);
				}
				visited[*it] = false;
				in_component[*it] = false;
			}
			myres.tc_comp_time += mt.elapsed();
			mt = my_timer();
			unsigned int s = cycle.front();
			for (i = cycle.begin(); i != cycle.end(); i++) {
				ds.join(s, *i);
				visited[*i] = false;
				in_component[*i] = false;
				in_component[s] = false;
			}
			s = ds.find(cycle.front());
			ds.collapse_scc(s, cycle);
			//condense_per_scc(s, cycle, g, g.ds);
			// update cycle 's topological sort 
			std::vector<unsigned int> cycle_nodes2;
			for (int i = 0; i < cycle.size(); i++) {
				cycle_nodes2.insert(cycle_nodes2.end(), ds.scc_nodes[cycle[i]].begin(), ds.scc_nodes[cycle[i]].end());
			}
			condense_per_scc(s, cycle_nodes2, cg, ds);
			myres.tc_comp_time += mt.elapsed();

			is_in_unorder[s] = true;
			n2v[s] = max_ccs[0];
			max_ccs[0] += 1;
			in_component[t] = false;
		}
		else {
			// update notGu nodes info.  
			mt = my_timer();
			for (it = reachable.begin(); it != reachable.end(); it++) {
				is_in_unorder[*it] = true;
				n2v[*it] = max_ccs[0];
				max_ccs[0] += 1;
				visited[*it] = false;
				in_component[*it] = false;
			}
			myres.tc_comp_time += mt.elapsed();
		}
	}
	else {
		// only adding the edge(t, h) is OK!
	}

}



template<class T>
void dynamic_tc_igc_scc<T>::add_edge_Gsi(vd_t t, vd_t h)
{
	unsigned int tn2v(n2v[t]);
	unsigned int hn2v(n2v[h]);
	unsigned int year = n2y[t];
	my_timer mt;
	//tn2v == hn2v, add the edge(t,h) in the scc. does not affect the topological order.
	//tn2v > hn2v, the edge(t,h) is actual topological order.
	if (tn2v < hn2v) {
		myres.invalid_edges += 1; 
		mt = my_timer();
		std::vector<unsigned int> reaching, reachable;
		bool flag = false;
		in_component[t] = true;

		fwd_dfs_Gsi(h, tn2v, reachable, flag, year);
		myres.tc_comp_time += mt.elapsed();
		if (flag == true) {		// there exists the cycle in the condensation graph. 
			mt = my_timer();
			for_each(reachable.begin(), reachable.end(), [this](const unsigned int & a) {
				visited[a] = false;
			});
			back_dfs_Gsi(t, hn2v, reaching, year);
			myres.search_time += mt.elapsed();
			mt = my_timer();
			std::vector<unsigned int> cycle, reachable_tmp, reaching_tmp;
			std::vector<unsigned int>::iterator i, j;
			unsigned int s = 0;
			for (i = reachable.begin(); i != reachable.end(); ++i) {
				if (in_component[*i] == true) {
					cycle.push_back(*i);
				}
				else {
					reachable_tmp.push_back(*i);
				}
				visited[*i] = false;
			}
			for (i = reaching.begin(); i != reaching.end(); i++) {
				visited[*i] = false;
				if (in_component[*i] == false)
					reaching_tmp.push_back(*i);
			}

			assert(cycle.size() > 1 && "the number of vertices in the cycle must be greater than 1 Gm -> Gm.");
			s = cycle.front();
			for (i = cycle.begin(); i != cycle.end(); i++) {
				ds.join(s, *i);
				visited[*i] = false;
				in_component[*i] = false;
				in_component[s] = false;
			}
			s = ds.find(cycle.front());
			ds.collapse_scc(s, cycle);
			//condense_per_scc(s, cycle, g, g.ds);
			std::vector<unsigned int> cycle_nodes2;
			for (int i = 0; i < cycle.size(); i++) {
				cycle_nodes2.insert(cycle_nodes2.end(), ds.scc_nodes[cycle[i]].begin(), ds.scc_nodes[cycle[i]].end());
			}
			condense_per_scc(s, cycle_nodes2, cg, ds);
			myres.condensation_time += mt.elapsed();

			mt = my_timer();
			dynamic_tc_igc_com<std::vector<T>> pc(n2v);
			sort(reaching_tmp.begin(), reaching_tmp.end(), pc);
			sort(reachable_tmp.begin(), reachable_tmp.end(), pc);
			reaching_tmp.push_back(s);
			myres.aff_region += (reaching_tmp.size() + reachable_tmp.size());

			tc_scc_reorder(reachable_tmp, reaching_tmp );
			myres.search_time += mt.elapsed();
		}
		else {
			mt = my_timer();
			back_dfs_Gsi(t, hn2v, reaching, year);
			myres.search_time += mt.elapsed();
			mt = my_timer();
			dynamic_tc_igc_com<std::vector<T>> pc(n2v);
			sort(reaching.begin(), reaching.end(), pc);
			sort(reachable.begin(), reachable.end(), pc);
			myres.aff_region += (reaching.size() + reachable.size());
			tc_scc_reorder(reachable, reaching );
			myres.reassignment_time += mt.elapsed();
		}
	}
	in_component[t] = false;
}

template<class T>
void dynamic_tc_igc_scc<T>::tc_scc_reorder(std::vector<unsigned int>& reachable, std::vector<unsigned int>& reaching)
{
	std::vector<vd_t> tmp;
	tmp.reserve(reaching.size() + reachable.size());
	std::vector<unsigned int>::iterator iend(reaching.end());
	for (std::vector<unsigned int>::iterator i(reaching.begin()); i != iend; ++i) {
		tmp.push_back(*i);
		visited[*i] = false;
		*i = n2v[*i]; // dirty trick
	}
	iend = reachable.end();
	for (std::vector<unsigned int>::iterator i(reachable.begin()); i != iend; ++i) {
		tmp.push_back(*i);
		visited[*i] = false;
		*i = n2v[*i]; // dirty trick
	}
	sort(reachable.begin(), reachable.end(), [](const unsigned int &a, const unsigned int &b) {return a > b; });
	sort(reaching.begin(), reaching.end(), [](const unsigned int &a, const unsigned int &b) {return a > b; });

	std::vector<unsigned int>::iterator i(reachable.begin()), it_ds;
	std::vector<unsigned int>::iterator j(reaching.begin());
	std::vector<unsigned int>::iterator jend(reaching.end());
	unsigned int index(0);
	while (i != iend || j != jend) {
		unsigned int w;
		if (j == jend || (i != iend && *i > *j)) {
			w = *i; ++i;
		}
		else {
			w = *j; ++j;
		}
		vd_t n(tmp[index++]);
		// allocate n at w
		n2v[n] = w;
	}
}


template<class T>
void dynamic_tc_igc_scc<T>::fwd_dfs_Gm(vd_t n, vd_t lb, std::vector<unsigned int>& reachable, bool & flag)
{
	
	reachable.push_back(n);
	visited[n] = true;
	out_iterator i, iend;
	std::vector<unsigned int>::iterator it;
	for (it = ds.scc_nodes[n].begin(); it != ds.scc_nodes[n].end(); it++) {
		for (std::tie(i, iend) = boost::out_edges(*it, cg); i != iend; ++i) {
			//unsigned int w(target(*i, g));
			unsigned int w = ds.find(target(*i, cg));
			unsigned int wn2v(n2v[w]);
			if (wn2v == lb && !visited[w] && is_in_unorder[w] == true) {
				in_component[w] = true;
				fwd_dfs_Gm(w, lb, reachable, flag);
				flag = true;
			}
			else if (wn2v > lb && !visited[w] && is_in_unorder[w] == true) {
				fwd_dfs_Gm(w, lb, reachable, flag);
			}
			in_component[n] = in_component[n] || in_component[w];
		}
	}
}

template<class T>
void dynamic_tc_igc_scc<T>::scan_GsGr_without_topo(vd_t n, vd_t lb, std::vector<unsigned int>& reachable, bool & flag)
{
	out_iterator i, iend;
	visited[n] = true;
	std::vector<unsigned int>::iterator it;
	for (it = ds.scc_nodes[n].begin(); it != ds.scc_nodes[n].end(); it++) {
		for (std::tie(i, iend) = boost::out_edges(*it, cg); i != iend; ++i) {
			//unsigned int w(target(*i, g));
			unsigned int w = ds.find(target(*i, cg));
			if (w == lb) {
				in_component[w] = true;
				if (!visited[w] && is_in_unorder[w] == false) {
					scan_GsGr_without_topo(w, lb, reachable, flag);
				}
				flag = true;
			}
			else if (!visited[w] && is_in_unorder[w] == false) { // w not visited & w is not in Gm. Because w in Gm, the w reachable nodes also in Gm. 
				scan_GsGr_without_topo(w, lb, reachable, flag);					// Then the nodes that need to be reordered will not be traversed. 
			}
			in_component[n] = in_component[n] || in_component[w];
		}
	}
	reachable.push_back(n);
}

template<class T>
void dynamic_tc_igc_scc<T>::scan_GsGr(vd_t n, std::vector<unsigned int>& reachable)
{
	visited[n] = true;
	out_iterator i, iend;
	std::vector<unsigned int>::iterator it;
	for (it = ds.scc_nodes[n].begin(); it != ds.scc_nodes[n].end(); it++) {
		for (std::tie(i, iend) = boost::out_edges(*it, cg); i != iend; ++i) {
			//unsigned int w(target(*i, g));
			unsigned int w = ds.find(target(*i, cg));
			if (!visited[w] && is_in_unorder[w] == false) {	// if w not visited & w not in Gm
				scan_GsGr(w, reachable);
			}
		}
	}
	reachable.push_back(n);
}

//template<class T>
//void dynamic_tc_igc_scc<T>::fwd_dfs_Gu2nGu(vd_t n, vd_t lb, std::vector<unsigned int>& reachable, bool & flag)
//{
//	visited[n] = true;
//	out_iterator i, iend;
//	std::vector<unsigned int>::iterator it;
//	for (it = ds.scc_nodes[n].begin(); it != ds.scc_nodes[n].end(); it++) {
//		for (std::tie(i, iend) = boost::out_edges(*it, cg); i != iend; ++i) {
//			//unsigned int w(target(*i, g));
//			unsigned int w = ds.find(target(*i, cg));
//			unsigned int wn2v(unorder_para.cc[w]);
//			if (wn2v == lb && !visited[w]) {
//				in_component[w] = true;
//				fwd_dfs_Gu2nGu(w, lb, reachable, flag);
//				flag = true;
//			}
//			else if ((wn2v > lb || is_in_unorder[w] == false) && !visited[w]) {
//				fwd_dfs_Gu2nGu(w, lb, reachable, flag);
//			}
//			in_component[n] = in_component[n] || in_component[w];
//		}
//	}
//	reachable.push_back(n);
//}

template<class T>
void dynamic_tc_igc_scc<T>::back_dfs_Gm(vd_t n, vd_t ub, std::vector<unsigned int>& reaching)
{
	reaching.push_back(n);
	visited[n] = true;
	in_iterator i, iend;
	std::vector<unsigned int>::iterator it;
	for (it = ds.scc_nodes[n].begin(); it != ds.scc_nodes[n].end(); it++) {
		for (std::tie(i, iend) = in_edges(*it, cg); i != iend; ++i) {
			//unsigned int w(source(*i, g));
			unsigned int w = ds.find(source(*i, cg));
			unsigned int wn2v(n2v[w]);
			if (wn2v < ub && ! visited[w] && is_in_unorder[w] == true)		// NOTE, add g.unorder_para.cc[w] == true
			{
				back_dfs_Gm(w, ub, reaching);
			}
		}
	}
}

template<class T>
void dynamic_tc_igc_scc<T>::fwd_dfs_Gsi(vd_t n, vd_t lb, std::vector<unsigned int>& reachable, bool & flag, unsigned int & year)
{
	out_iterator i, iend;
	reachable.push_back(n);
	visited[n] = true;
	std::vector<unsigned int>::iterator it;
	for (it = ds.scc_nodes[n].begin(); it != ds.scc_nodes[n].end(); it++) {
		for (std::tie(i, iend) = boost::out_edges(*it,cg); i != iend; ++i) {
			//unsigned int w(target(*i, g));
			unsigned int w = ds.find(target(*i, cg));
			unsigned int wn2v(n2v[w]);
			if (wn2v == lb && !visited[w] && is_in_unorder[w] == false && is_in_same_year[w] == year && n2y[w] == year) {
				in_component[w] = true;
				fwd_dfs_Gsi(w, lb, reachable, flag, year);
				flag = true;
			}
			else if (wn2v > lb && !visited[w] && is_in_unorder[w] == false && is_in_same_year[w] == year && n2y[w] == year) {
				fwd_dfs_Gsi(w, lb, reachable, flag, year);
			}
			in_component[n] = in_component[n] || in_component[w];
		}
	}
}

template<class T>
void dynamic_tc_igc_scc<T>::back_dfs_Gsi(vd_t n, vd_t ub, std::vector<unsigned int>& reaching, unsigned int & year)
{
	in_iterator i, iend;
	reaching.push_back(n);
	visited[n] = true;
	std::vector<unsigned int>::iterator it;
	for (it = ds.scc_nodes[n].begin(); it != ds.scc_nodes[n].end(); it++) {
		for (std::tie(i, iend) = in_edges(*it, cg); i != iend; ++i) {
			//unsigned int w(source(*i, g));
			unsigned int w = ds.find(source(*i, cg));
			unsigned int wn2v(n2v[w]);
			if (wn2v < ub && ! visited[w] && is_in_unorder[w] == false && is_in_same_year[w] == year && n2y[w] == year) {
				back_dfs_Gsi(w, ub, reaching, year);
			}
		}
	}
}

template<class T>
dynamic_tc_igc_scc<T>::~dynamic_tc_igc_scc() { }
