#include "bat_DSCC_test.h"
template<class T>
batDSCC_test<T>::batDSCC_test() { }

template<class T>
batDSCC_test<T>::batDSCC_test(CG & cg, vector<ahrsz_priority_value<T>>& _n2v, vector<T*>& _pspaces, disjoint_set & _ds,
	vector<unsigned int>& _n2y, vector<int>& _is_in_same_year, vector<bool>& _is_in_unorder,
	int _min_year, int _max_year, std::vector<pair<unsigned int, unsigned int>>& _inc_edges, const int np, const int _bs) : cg(cg), n2v(_n2v), pspaces(_pspaces),
	ds(_ds), n2y(_n2y), is_in_same_year(_is_in_same_year), is_in_unorder(_is_in_unorder), min_year(_min_year), max_year(_max_year), inc_edges(_inc_edges), bs(_bs) {

	std::cout << "Init Dynamic BAT-DSCC-test- Algorithm" << std::endl;
	self_name = "BAT-DSCC-test-";
	myres.set_dynamic_name(self_name);
	_flooring = std::vector<ahrsz_ext_priority_value<T> >(np, minus_infinity);
	_visited = std::vector<bool>(np, false);
	_inK = std::vector<bool>(np, false);
	_inF = std::vector<bool>(np, false);
	_inB = std::vector<bool>(np, false);
	in_component = std::vector<bool>(np, false);
	inAffReg = std::vector<bool>(np, false);
}

template<class T>
void batDSCC_test<T>::add_edges()
{
	add_edges_bat();
}
template<class T>
void batDSCC_test<T>::add_edges_bat() {

	std::vector<pair<unsigned int, unsigned int>> invalid_edges;
	invalid_edges.reserve(bs);
	std::vector<pair<unsigned int, unsigned int>> bat_Gu, rem_bat;
	std::vector< std::vector<pair<unsigned int, unsigned int>>> bat_Gyi;		// 
	std::vector<int> bat_Gyi_idx;													// the index of bat_Gyi which remain edges to handle
	for (int i = min_year; i <= max_year; i++) {
		std::vector<pair<unsigned int, unsigned int>> tmp;
		bat_Gyi.push_back(tmp);
	}
	int tmp_b = 0, type = 0;

	my_timer mt;
	for (int i = 0; i < inc_edges.size(); ) {
		while (tmp_b < bs && i < inc_edges.size()) {
			unsigned int s, t;
			s = ds.find(inc_edges[i].first);
			t = ds.find(inc_edges[i].second);
			if (s != t) {
				type = get_alg_type(s, t, is_in_same_year, is_in_unorder);
				if (type == 1) {
					if (!(n2v[s] < n2v[t] || n2v[s] == n2v[t])) {
						boost::add_edge(inc_edges[i].first, inc_edges[i].second, cg);
					}
					else {
						bat_Gu.push_back(inc_edges[i]);
					}
				}
				else if (type == 4) {												// s, t are both in 		
					if (n2y[s] == n2y[t]) {											// push edge into bat_Gyi
						if (!(n2v[s] < n2v[t] || n2v[s] == n2v[t])) {
							boost::add_edge(inc_edges[i].first, inc_edges[i].second, cg);
						}
						else {
							if (bat_Gyi[n2y[t] - min_year].size() == 0) {
								bat_Gyi_idx.push_back(n2y[t] - min_year);
							}
							bat_Gyi[n2y[t] - min_year].push_back(inc_edges[i]);
						}
					}
					else if (n2y[s] < n2y[t]) {										// push into rem_bat	
						rem_bat.push_back(inc_edges[i]);
					}
					else {															// only add edge into Gsi
						boost::add_edge(inc_edges[i].first, inc_edges[i].second, cg);
					}
				}
				else {
					rem_bat.push_back(inc_edges[i]);
				}

			}
			i++;
			tmp_b++;
		}

		valid_aware_Gm(bat_Gu);
		valid_aware_Gsi(bat_Gyi, bat_Gyi_idx);
		add_edge_rem(rem_bat);
		bat_Gyi_idx.clear();
		tmp_b = 0;
	}
	myres.dn_eval_time = myres.search_time + myres.reassignment_time + myres.tc_comp_time;
	myres.dn_eval_time_AAN = mt.elapsed();
}

template<class T>
void batDSCC_test<T>::add_edges(int start_index, int end_index)
{
	myres.inc_edges += (end_index - start_index);
	std::vector<pair<unsigned int, unsigned int>> invalid_edges;
	invalid_edges.reserve(bs);
	std::vector<pair<unsigned int, unsigned int>> bat_Gu, rem_bat;
	std::vector< std::vector<pair<unsigned int, unsigned int>>> bat_Gyi;
	std::vector<pair<unsigned int, unsigned int>> eo2n, em2s, others; 
	
	std::vector<int> bat_Gyi_idx;													// the index of bat_Gyi which remain edges to handle
	for (int i = min_year; i <= max_year; i++) {
		std::vector<pair<unsigned int, unsigned int>> tmp;
		bat_Gyi.push_back(tmp);
	}
	int tmp_b = 0, type = 0;
	std::vector<bool> drop_oth = std::vector<bool>((end_index - start_index), false);

	my_timer mt;
	for (int i = start_index; i < end_index; ) {
		while (tmp_b < bs && i < end_index) {
			unsigned int s, t;
			s = ds.find(inc_edges[i].first);
			t = ds.find(inc_edges[i].second);
			if (s != t) {
				type = get_alg_type(s, t, is_in_same_year, is_in_unorder);
				if (type == 1) {
					if (!(n2v[s] < n2v[t] || n2v[s] == n2v[t])) {
						boost::add_edge(inc_edges[i].first, inc_edges[i].second, cg);
					}
					else {
						others.push_back(inc_edges[i]);
					}
				}
				else if (type ==  2) {
					boost::add_edge(inc_edges[i].first, inc_edges[i].second, cg);
				}
				else if (type == 3) {
					em2s.push_back(inc_edges[i]);
				} 
				else if (type == 4 || type == 5 || type == 6 ) {												
					if (n2y[s] == n2y[t]) {
						if(is_in_same_year[s] == 0 && is_in_same_year[t] != 0) {					// s in Gsi, t in Gr
							boost::add_edge(inc_edges[i].first, inc_edges[i].second, cg);
						}
						else if (is_in_same_year[s] != 0 && is_in_same_year[t] != 0) {				// s in Gsi, t in Gsi	
							if (!(n2v[s] < n2v[t] || n2v[s] == n2v[t])) {		
								boost::add_edge(inc_edges[i].first, inc_edges[i].second, cg);
							}
							else {
								others.push_back(inc_edges[i]);
							}
						}
						else {
							others.push_back(inc_edges[i]);
						}
					}
					else if (n2y[s] < n2y[t]) {
						eo2n.push_back(inc_edges[i]);
					}
					else {
						boost::add_edge(inc_edges[i].first, inc_edges[i].second, cg);
					}

				}
			}
			i++;
			tmp_b++;
		}
		
		add_eo2n(eo2n);
		
		for (int i = 0; i < others.size(); i++) {
			unsigned int s, t;
			s = ds.find(others[i].first);
			t = ds.find(others[i].second);
			if (s != t) {
				type = get_alg_type(s, t, is_in_same_year, is_in_unorder);
				if (type == 1) {
					if (!(n2v[s] < n2v[t] || n2v[s] == n2v[t])) {
						boost::add_edge(others[i].first, others[i].second, cg);
						drop_oth[i] = true;
					}
				}
				else if (type == 2) {
					boost::add_edge(others[i].first, others[i].second, cg);
					drop_oth[i] = true; 
				}
				else if (type == 3) {
					em2s.push_back(others[i]);
					drop_oth[i] = true;
				}
			}
		}
		for (int i = 0; i < em2s.size(); i++) {
			add_item(em2s[i]);
		}
		em2s.clear();
		for (int i = 0; i < others.size(); i++) {
			if (drop_oth[i] == false) {
				unsigned int s, t;
				s = ds.find(others[i].first);
				t = ds.find(others[i].second);
				if (s != t) {
					type = get_alg_type(s, t, is_in_same_year, is_in_unorder);
					if (type == 1) {
						if (!(n2v[s] < n2v[t] || n2v[s] == n2v[t])) {
							boost::add_edge(others[i].first, others[i].second, cg);
						}
						else {
							bat_Gu.push_back(others[i]);
						}
					}
					else if (type == 4) {												// s, t are both in 		
						if (n2y[s] == n2y[t]) {											// push edge into bat_Gyi
							if (!(n2v[s] < n2v[t] || n2v[s] == n2v[t])) {
								boost::add_edge(others[i].first, others[i].second, cg);
							}
							else {
								if (bat_Gyi[n2y[t] - min_year].size() == 0) {
									bat_Gyi_idx.push_back(n2y[t] - min_year);
								}
								bat_Gyi[n2y[t] - min_year].push_back(others[i]);
							}
						}
						else if (n2y[s] < n2y[t]) {										// push into rem_bat	
							rem_bat.push_back(others[i]);
						}
						else {															// only add edge into Gsi
							boost::add_edge(others[i].first, others[i].second, cg);
						}
					}
					else {
						rem_bat.push_back(others[i]);
					}
				}
				add_item(others[i]);
			}
			else {
				drop_oth[i] = false; 
			}
		}

		valid_aware_Gm(bat_Gu);
		valid_aware_Gsi(bat_Gyi, bat_Gyi_idx);
		add_edge_rem(rem_bat);
		
		bat_Gyi_idx.clear();
		eo2n.clear();
		em2s.clear();
		others.clear();
		tmp_b = 0;
	}
	myres.dn_eval_time = myres.search_time + myres.reassignment_time + myres.tc_comp_time;
	myres.dn_eval_time_AAN = mt.elapsed();
	record_cur_per_info(myres);
}

template<class T>
void batDSCC_test<T>::add_item(pair<unsigned int, unsigned int> & edge) { 
	unsigned int s, t;
	int type = -1;
	s = ds.find(edge.first);
	t = ds.find(edge.second);
	if (s == t) return;
	boost::add_edge(edge.first, edge.second, cg);										// NOTE: add_edge(it->first, it->second, cg); not add_edge(s, t, cg)
	
	type = get_alg_type(s, t, is_in_same_year, is_in_unorder);
	if (type == 1) {																// s, t in Gsi & Gr, Gsi & Gr 
		add_edge_Gm(s, t);
	}
	else if (type == 2) {															// s, t in Gsi & Gr, Gm // algorithm1
	}
	else if (type == 3) {															// s, t in Gm, Gsi & Gr, //algorithm 2 
		add_edge_alg2(s, t);
	}
	else if (type == 4 || type == 5 || type == 6) {									// s in Gsi or Gr,  t in Gsi or Gr // algrithm 3
		add_edge_alg3(s, t);
	}
	else {
		assert(0 == 1 && "ERROR EDGE TYPE.");
	}
}


template<class T>
void batDSCC_test<T>::add_eo2n(std::vector<pair<unsigned int, unsigned int>> & eo2n) {
	for (int i = 0; i < eo2n.size(); i ++) {
		boost::add_edge(eo2n[i].first, eo2n[i].second, cg);
	}
	std::vector<std::vector<unsigned int>> topo_scc;
	std::vector<std::vector<unsigned int>>::iterator it1;
	std::vector<unsigned int>::iterator it2;
	ppara para(boost::num_vertices(cg));
	ahrsz_ext_priority_value<T>  max_priority(minus_infinity);
	my_timer mt; 
	unsigned int t; 
	for (int i = 0; i < eo2n.size(); i++) {
		t = ds.find(eo2n[i].second);
		if (para.rindex[t] == 0) {
			pearce_gs(t, para, max_priority, topo_scc);
		}
	}
	
	ahrsz_ext_priority_value<T> floor, ceiling;
	T* u_space = (pspaces[0]);
	if (max_priority == ahrsz_ext_priority_value<T>(minus_infinity) || std::next(max_priority.base()) == u_space->end()) {
		floor = ahrsz_ext_priority_value<T>(u_space->previous(u_space->end()), *u_space);
		ceiling = ahrsz_ext_priority_value<T>(plus_infinity);
	}
	else {
		floor = ahrsz_ext_priority_value<T>(max_priority);
		ceiling = ahrsz_ext_priority_value<T>(std::next(max_priority.base()), *u_space);
	}
	unsigned int s; 
	for (it1 = topo_scc.begin(); it1 != topo_scc.end(); it1++) {
		if (it1->size() > 1) {
			s = it1->front();
			std::vector<unsigned int> cycle_nodes;
			for (it2 = it1->begin(); it2 != it1->end(); it2 ++) {
				ds.join(s, *it2);
				cycle_nodes.insert(cycle_nodes.end(), ds.scc_nodes[*it2].begin(), ds.scc_nodes[*it2].end() );
			}
			ds.collapse_scc(ds.find(s), *it1);
			s = ds.find(s);
			is_in_unorder[s] = true; 
			ahrsz_ext_priority_value<T> ep(compute_priority(floor, ceiling, *u_space));
			n2v[s] = ahrsz_priority_value<T>(ep.base(), *u_space);
			floor = n2v[s];
			condense_per_scc(s, cycle_nodes, cg, ds);
		}
		else {
			is_in_unorder[it1->front()] = true; 
			ahrsz_ext_priority_value<T> ep(compute_priority(floor, ceiling, *u_space));
			n2v[it1->front()] = ahrsz_priority_value<T>(ep.base(), *u_space);
			floor = n2v[it1->front()];
		}
	}
	myres.tc_comp_time += mt.elapsed();


}
template<class T>
void batDSCC_test<T>:: pearce_gs(vd_t v, ppara &para, ahrsz_ext_priority_value<T>& max_priority, std::vector<std::vector<unsigned int>> &topo_scc) {
	
	out_iterator ai, ai_end;
	int root = true;
	para.rindex[v] = para.index_p;
	para.index_p++;
	unsigned int w;

	for (std::tie(ai, ai_end) = out_edges(v, cg); ai != ai_end; ++ai) {
		w = ds.find(boost::target(*ai, cg));

		if (para.rindex[w] == 0 && is_in_unorder[w] == false ) { 
			pearce_gs(w, para, max_priority, topo_scc);
		}
		if (para.rindex[w] < para.rindex[v] && is_in_unorder[w] == false && is_in_unorder[v] == false) {
			para.rindex[v] = para.rindex[w];
			root = false;
		}
		if (is_in_unorder[w] == true) {
			if (max_priority < ahrsz_ext_priority_value<T>(n2v[w].base(), *(pspaces[0]))) {
				max_priority = n2v[w];
			}
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



// different from the implementation of bat_AHRSZ_P_scc. 
template<class T>
void batDSCC_test<T>::valid_aware_Gm(std::vector<pair<unsigned int, unsigned int>> &bat_Gu)
{
	std::vector<pair<unsigned int, unsigned int>> poss_valid_cand;
	unsigned int s, t;
	while (bat_Gu.size() > 0) {
		for (int i = 0; i < bat_Gu.size(); i++) {
			s = ds.find(bat_Gu[i].first);
			t = ds.find(bat_Gu[i].second);
			if (s == t) continue;
			if (inAffReg[s] || inAffReg[t]) {											// inAffReg only store the canonical vertex.
				poss_valid_cand.push_back(bat_Gu[i]);
			}
			else {
				boost::add_edge(bat_Gu[i].first, bat_Gu[i].second, cg);
				add_edge_Gm(s, t, inAffReg, allK);
			}
		}
		bat_Gu.clear();
		for (int i = 0; i < poss_valid_cand.size(); i++) {
			s = ds.find(poss_valid_cand[i].first);
			t = ds.find(poss_valid_cand[i].second);
			if (s == t) continue;
			if (!(n2v[s] < n2v[t] || n2v[s] == n2v[t])) {
				boost::add_edge(poss_valid_cand[i].first, poss_valid_cand[i].second, cg);
			}
			else {
				//bat_Gu.push_back(inc_edges[i]);
				bat_Gu.push_back(poss_valid_cand[i]);
			}
		}
		poss_valid_cand.clear();

		// clear inAffReg info. 
		for (int i = 0; i < allK.size(); i++) {
			inAffReg[allK[i]] = false;
		}
		allK.clear();
	}
}

template<class T>
void batDSCC_test<T>::valid_aware_Gsi(std::vector< std::vector<pair<unsigned int, unsigned int>>> &bat_Gyi, std::vector<int> & bat_Gyi_idx)
{

	std::vector<pair<unsigned int, unsigned int>> cur_Gyi;
	std::vector<pair<unsigned int, unsigned int>> poss_valid_cand;
	unsigned int s, t;

	for (int i = 0; i < bat_Gyi_idx.size(); i++) {
		cur_Gyi = bat_Gyi[bat_Gyi_idx[i]];
		while (cur_Gyi.size() > 0) {
			for (int j = 0; j < cur_Gyi.size(); j++) {
				s = ds.find(cur_Gyi[j].first);
				t = ds.find(cur_Gyi[j].second);
				if (s == t) continue;
				if (inAffReg[s] || inAffReg[t]) {											// inAffReg only store the canonical vertex.
					poss_valid_cand.push_back(cur_Gyi[j]);
				}
				else {
					boost::add_edge(cur_Gyi[j].first, cur_Gyi[j].second, cg);
					add_edge_Gsi(s, t, inAffReg, allK);
				}
			}
			cur_Gyi.clear();
			bat_Gyi[bat_Gyi_idx[i]].clear();
			for (int j = 0; j < poss_valid_cand.size(); j++) {
				s = ds.find(poss_valid_cand[j].first);
				t = ds.find(poss_valid_cand[j].second);
				if (s == t) continue;
				if (!(n2v[s] < n2v[t] || n2v[s] == n2v[t])) {
					boost::add_edge(poss_valid_cand[j].first, poss_valid_cand[j].second, cg);
				}
				else {
					cur_Gyi.push_back(poss_valid_cand[j]);
				}
			}
			poss_valid_cand.clear();
			// clear inAffReg info. 
			for (int j = 0; j < allK.size(); j++) {
				inAffReg[allK[j]] = false;
			}
			allK.clear();
		}
		// clear inAffReg info. 
		for (int j = 0; j < allK.size(); j++) {
			inAffReg[allK[j]] = false;
		}
		allK.clear();
		cur_Gyi.clear();
		bat_Gyi[bat_Gyi_idx[i]].clear();

	}
	bat_Gyi_idx.clear();

}

template<class T>
void batDSCC_test<T>::add_edge_rem(std::vector<pair<unsigned int, unsigned int>>& bat_rem)
{
	std::vector<pair<unsigned int, unsigned int>>::iterator it;
	for (it = bat_rem.begin(); it != bat_rem.end(); it++) {
		unsigned int s, t;
		int type = -1;
		s = ds.find(it->first);
		t = ds.find(it->second);
		if (s == t) continue;
		my_timer mt = my_timer();
		boost::add_edge(it->first, it->second, cg);										 // NOTE: add_edge(it->first, it->second, cg); not add_edge(s, t, cg)
		myres.add_edge_time += mt.elapsed();
		type = get_alg_type(s, t, is_in_same_year, is_in_unorder);
		if (type == 1) {																 //s, t are both in Gm, or Gsi
			add_edge_Gm(s, t);
		}
		else if (type == 2) {															 // s, t in Gsi & Gr, Gm // algorithm1
		}
		else if (type == 3) {															 // s, t in Gm, Gsi & Gr, //algorithm 2 	
			myres.invalid_part_cnt++;
			add_edge_alg2(s, t);
		}
		else if (type == 4 || type == 5 || type == 6) {									// s in Gsi or Gr,  t in Gsi or Gr // algrithm 3
			add_edge_alg3(s, t);
		}
		else {
			assert(0 == 1 && "ERROR EDGE TYPE.");
		}
	}
	bat_rem.clear();

}


template<class T>
void batDSCC_test<T>::add_edge_Gm(vd_t t, vd_t h) {
	// NOTE: the next line can be move to func of add_edges(); 
	if ((n2v[t] < n2v[h]) || (n2v[t] == n2v[h])) {				// t -> h, i.e,  n2v(t) > n2v(h). when not holds this priority, discovery, reassignment.
		std::vector<vd_t> K, remaining_K;
		std::vector<unsigned int> cycle_nodes;
		std::vector<vd_t>::iterator it;
		bool has_cycle = false;
		my_timer mt;
		discover_Gm(t, h, K, has_cycle, cycle_nodes);
		myres.set_search_time(myres.get_search_time() + mt.elapsed());
		myres.set_aff_region(myres.get_aff_region() + K.size());
		myres.set_invalid_edges(myres.get_invalid_edges() + 1);

		if (has_cycle) {
			assert(cycle_nodes.size() > 1 && "the number of vertices in the cycle must be greater than 1. dynamic_AHRSZ");
			unsigned int s = cycle_nodes.front();
			for (std::vector<unsigned int>::iterator it = cycle_nodes.begin(); it != cycle_nodes.end(); it++) {
				ds.join(s, *it);
			}
			s = ds.find(cycle_nodes.front());
			ds.collapse_scc(s, cycle_nodes);

			// NOTE: different with the fully condensation graph. i.e., V1. 
			mt = my_timer();
			std::vector<unsigned int> cycle_nodes2;
			for (int i = 0; i < cycle_nodes.size(); i++) {
				cycle_nodes2.insert(cycle_nodes2.end(), ds.scc_nodes[cycle_nodes[i]].begin(), ds.scc_nodes[cycle_nodes[i]].end());
			}
			condense_per_scc(s, cycle_nodes2, cg, ds);
			myres.condensation_time += (mt.elapsed());

			remaining_K.push_back(s);
			for (it = K.begin(); it != K.end(); it++) {
				in_component[*it] = false;												// unmark state info
				_visited[*it] = false;
				_inB[*it] = false;
				_inF[*it] = false;
				_inK[*it] = false;
				if (ds.find(*it) != ds.find(s)) {
					remaining_K.push_back(*it);
				}
			}
		}

		mt = my_timer();
		maintain_Gm(K);
		myres.set_reassign_time(myres.get_reassign_time() + mt.elapsed());
	}
}

template<class T>
void batDSCC_test<T>::add_edge_Gm(vd_t t, vd_t h, std::vector<bool> & inAffReg, std::vector<unsigned int> &allK) {
	// NOTE: the next line can be move to func of add_edges(); 
	if ((n2v[t] < n2v[h]) || (n2v[t] == n2v[h])) {				// t -> h, i.e,  n2v(t) > n2v(h). when not holds this priority, discovery, reassignment.
		myres.invalid_local_cnt++;
		std::vector<vd_t> K, remaining_K;
		std::vector<unsigned int> cycle_nodes;
		std::vector<vd_t>::iterator it;
		bool has_cycle = false;
		my_timer mt;
		discover_Gm(t, h, K, has_cycle, cycle_nodes);
		myres.set_search_time(myres.get_search_time() + mt.elapsed());
		myres.set_aff_region(myres.get_aff_region() + K.size());
		myres.set_invalid_edges(myres.get_invalid_edges() + 1);
		// record affreg info. 
		for (int i = 0; i < K.size(); i++) {
			inAffReg[K[i]] = true;
			allK.push_back(K[i]);
		}

		if (has_cycle) {
			assert(cycle_nodes.size() > 1 && "the number of vertices in the cycle must be greater than 1. dynamic_AHRSZ");
			unsigned int s = cycle_nodes.front();
			for (std::vector<unsigned int>::iterator it = cycle_nodes.begin(); it != cycle_nodes.end(); it++) {
				ds.join(s, *it);
			}
			s = ds.find(cycle_nodes.front());
			ds.collapse_scc(s, cycle_nodes);

			// NOTE: different with the fully condensation graph. i.e., V1. 
			mt = my_timer();
			std::vector<unsigned int> cycle_nodes2;
			for (int i = 0; i < cycle_nodes.size(); i++) {
				cycle_nodes2.insert(cycle_nodes2.end(), ds.scc_nodes[cycle_nodes[i]].begin(), ds.scc_nodes[cycle_nodes[i]].end());
			}
			condense_per_scc(s, cycle_nodes2, cg, ds);
			myres.condensation_time += (mt.elapsed());

			remaining_K.push_back(s);
			for (it = K.begin(); it != K.end(); it++) {
				in_component[*it] = false;												// unmark state info
				_visited[*it] = false;
				_inB[*it] = false;
				_inF[*it] = false;
				_inK[*it] = false;
				if (ds.find(*it) != ds.find(s)) {
					remaining_K.push_back(*it);
				}
			}
		}

		mt = my_timer();
		maintain_Gm(K);
		myres.set_reassign_time(myres.get_reassign_time() + mt.elapsed());
	}
}


template<class T>
void batDSCC_test<T>::add_edge_Gsi(vd_t t, vd_t h)
{
	if ((n2v[t] < n2v[h]) || (n2v[t] == n2v[h])) {									// t -> h, i.e,  pmap(t) > pmap(h). when not holds this priority, discovery, reassignment.
		myres.invalid_local_cnt++;
		std::vector<vd_t> K, remaining_K;
		std::vector<unsigned int> cycle_nodes;
		std::vector<vd_t>::iterator it;
		bool has_cycle = false;
		my_timer mt;
		discover_Gsi(t, h, K, has_cycle, cycle_nodes);
		myres.set_search_time(myres.get_search_time() + mt.elapsed());
		myres.set_aff_region(myres.get_aff_region() + K.size());
		myres.set_invalid_edges(myres.get_invalid_edges() + 1);

		if (has_cycle) {
			assert(cycle_nodes.size() > 1 && "the number of vertices in the cycle must be greater than 1. dynamic_AHRSZ");
			unsigned int s = cycle_nodes.front();
			for (std::vector<unsigned int>::iterator it = cycle_nodes.begin(); it != cycle_nodes.end(); it++) {
				ds.join(s, *it);
			}
			s = ds.find(cycle_nodes.front());
			ds.collapse_scc(s, cycle_nodes);

			mt = my_timer();
			// NOTE: different with the fully condensation graph. i.e., V1. 
			std::vector<unsigned int> cycle_nodes2;
			for (int i = 0; i < cycle_nodes.size(); i++) {
				cycle_nodes2.insert(cycle_nodes2.end(), ds.scc_nodes[cycle_nodes[i]].begin(), ds.scc_nodes[cycle_nodes[i]].end());
			}
			condense_per_scc(s, cycle_nodes2, cg, ds);
			myres.condensation_time += (mt.elapsed());

			remaining_K.push_back(s);
			for (it = K.begin(); it != K.end(); it++) {
				in_component[*it] = false;	// unmark state info
				_visited[*it] = false;
				_inB[*it] = false;
				_inF[*it] = false;
				_inK[*it] = false;
				if (ds.find(*it) != ds.find(s)) {
					remaining_K.push_back(*it);
				}
			}
		}

		mt = my_timer();
		maintain_Gsi(K);
		myres.set_reassign_time(myres.get_reassign_time() + mt.elapsed());

	}

}

template<class T>
void batDSCC_test<T>::add_edge_Gsi(vd_t t, vd_t h, std::vector<bool> & inAffReg, std::vector<unsigned int> &allK)
{
	if ((n2v[t] < n2v[h]) || (n2v[t] == n2v[h])) {									// t -> h, i.e,  pmap(t) > pmap(h). when not holds this priority, discovery, reassignment.
		std::vector<vd_t> K, remaining_K;
		std::vector<unsigned int> cycle_nodes;
		std::vector<vd_t>::iterator it;
		bool has_cycle = false;
		my_timer mt;
		discover_Gsi(t, h, K, has_cycle, cycle_nodes);
		myres.set_search_time(myres.get_search_time() + mt.elapsed());
		myres.set_aff_region(myres.get_aff_region() + K.size());
		myres.set_invalid_edges(myres.get_invalid_edges() + 1);
		// record affreg info. 
		for (int i = 0; i < K.size(); i++) {
			inAffReg[K[i]] = true;
			allK.push_back(K[i]);
		}

		if (has_cycle) {
			assert(cycle_nodes.size() > 1 && "the number of vertices in the cycle must be greater than 1. dynamic_AHRSZ");
			unsigned int s = cycle_nodes.front();
			for (std::vector<unsigned int>::iterator it = cycle_nodes.begin(); it != cycle_nodes.end(); it++) {
				ds.join(s, *it);
			}
			s = ds.find(cycle_nodes.front());
			ds.collapse_scc(s, cycle_nodes);

			mt = my_timer();
			// NOTE: different with the fully condensation graph. i.e., V1. 
			std::vector<unsigned int> cycle_nodes2;
			for (int i = 0; i < cycle_nodes.size(); i++) {
				cycle_nodes2.insert(cycle_nodes2.end(), ds.scc_nodes[cycle_nodes[i]].begin(), ds.scc_nodes[cycle_nodes[i]].end());
			}
			condense_per_scc(s, cycle_nodes2, cg, ds);
			myres.condensation_time += (mt.elapsed());

			remaining_K.push_back(s);
			for (it = K.begin(); it != K.end(); it++) {
				in_component[*it] = false;	// unmark state info
				_visited[*it] = false;
				_inB[*it] = false;
				_inF[*it] = false;
				_inK[*it] = false;
				if (ds.find(*it) != ds.find(s)) {
					remaining_K.push_back(*it);
				}
			}
		}

		mt = my_timer();
		maintain_Gsi(K);
		myres.set_reassign_time(myres.get_reassign_time() + mt.elapsed());

	}

}

template<class T>
void batDSCC_test<T>::add_edge_alg1(vd_t t, vd_t h)			// 
{
	my_timer mt;
	ahrsz_ext_priority_value<T>  max_priority(minus_infinity);
	std::vector<unsigned int> reachable;
	std::vector<unsigned int>::iterator it, it_ds;
	// NOTE: we find the affect region of |Gm| by traversing the G from h to facilitate to minimize the size of |Gm|

	//ahrsz_ext_priority_value<T> floor(n2v[h]), ceiling = minus_infinity;
	//scan_GsGr(t, reachable, max_priority);
	scan_GsGr(h, reachable, max_priority);

	ahrsz_ext_priority_value<T> floor(n2v[h]), ceiling = minus_infinity;
	if (floor < max_priority) {
		floor = max_priority;
	}
	if (std::next(floor.base()) == pspaces[0]->end()) {
		ceiling = plus_infinity;
	}
	else {
		ceiling = ahrsz_ext_priority_value<T>(std::next(floor.base()), *(pspaces[0]));
	}
	// update notGu nodes info.
	for (it = reachable.begin(); it != reachable.end(); it++) {
		_visited[*it] = false;
		is_in_unorder[*it] = true;
		ahrsz_ext_priority_value<T> ep(compute_priority(floor, ceiling, *(pspaces[0])));
		n2v[*it] = ahrsz_priority_value<T>(ep.base(), *(pspaces[0]));
		floor = n2v[*it];
	}
	myres.tc_comp_time += mt.elapsed();
}

template<class T>
void batDSCC_test<T>::add_edge_alg2(vd_t t, vd_t h)
{
	my_timer mt;
	ahrsz_ext_priority_value<T> max_priority(minus_infinity);
	std::vector<unsigned int> reachable;
	std::vector<unsigned int>::iterator it, it_ds;
	scan_GsGr(h, reachable, max_priority);
	myres.tc_comp_time += mt.elapsed();

	if (max_priority < ahrsz_ext_priority_value<T>(n2v[t])) {
		mt = my_timer();
		// no cycle. only update reachable order priority.
		ahrsz_ext_priority_value<T> floor, ceiling;
		if (max_priority == ahrsz_ext_priority_value<T>(minus_infinity)) {
			floor = ahrsz_ext_priority_value<T>(minus_infinity);
			ceiling = ahrsz_ext_priority_value<T>(pspaces[0]->begin(), *(pspaces[0]));
		}
		else {
			floor = max_priority;
			ceiling = ahrsz_ext_priority_value<T>(std::next(floor.base()), *(pspaces[0]));
		}
		for (it = reachable.begin(); it != reachable.end(); it++) {
			is_in_unorder[*it] = true;
			_visited[*it] = false;
			ahrsz_ext_priority_value<T> ep(compute_priority(floor, ceiling, *(pspaces[0])));
			n2v[*it] = ahrsz_priority_value<T>(ep.base(), *(pspaces[0]));
			floor = n2v[*it];
		}
		myres.tc_comp_time += mt.elapsed();
	}
	else {
		mt = my_timer();
		// the inserted edge( t -> h ) could cause the cycle. 
		ahrsz_ext_priority_value<T> floor, ceiling;
		floor = max_priority;
		if (std::next(max_priority.base()) == pspaces[0]->end()) {
			ceiling = plus_infinity;
		}
		else {
			ceiling = ahrsz_ext_priority_value<T>(std::next(max_priority.base()), *(pspaces[0]));
		}
		for (it = reachable.begin(); it != reachable.end(); it++) {
			_visited[*it] = false;
			is_in_unorder[*it] = true;
			ahrsz_ext_priority_value<T> ep(compute_priority(floor, ceiling, *(pspaces[0])));
			n2v[*it] = ahrsz_priority_value<T>(ep.base(), *(pspaces[0]));
			floor = n2v[*it];
		}
		myres.tc_comp_time += mt.elapsed();
		add_edge_Gm(t, h);		// t -> h is invalid edge. 
	}


}


template<class T>
void batDSCC_test<T>::add_edge_alg3(vd_t t, vd_t h)												// w.r.t. dfs start with t in Static-TC
{
	int yidx = n2y[t] - min_year + 1;
	T *y_pspace = (pspaces[yidx]);
	my_timer mt;
	if (n2y[t] == n2y[h]) {
		if (is_in_same_year[t] != 0 && is_in_same_year[h] == 0) {										// t in Gsi, h in Gr
			myres.invalid_part_cnt++;
			mt = my_timer();
			is_in_same_year[h] = n2y[t];
			ahrsz_ext_priority_value<T> max_priority(minus_infinity);
			ahrsz_ext_priority_value<T> min_priority(plus_infinity);
			//scan_GsGr_without_topo(h, t, reachable, max_priority, flag);
			out_iterator i, iend;
			std::vector<unsigned int>::iterator it;
			for (it = ds.scc_nodes[h].begin(); it != ds.scc_nodes[h].end(); it++) {
				for (std::tie(i, iend) = boost::out_edges(*it, cg); i != iend; ++i) {
					//unsigned int w(target(*i, g));
					unsigned int w = ds.find(target(*i, cg));
					if (is_in_unorder[w] == false && is_in_same_year[w] == is_in_same_year[h]) {
						if (max_priority < ahrsz_ext_priority_value<T>(n2v[w].base(), *y_pspace)) {
							max_priority = n2v[w];
						}
					}
				}
			}
			ahrsz_ext_priority_value<T> ep(compute_priority(max_priority, min_priority, *y_pspace));
			n2v[h] = ahrsz_priority_value<T>(ep.base(), *y_pspace);
			myres.tc_comp_time += mt.elapsed();
			if (n2v[t] < n2v[h]) {
				add_edge_Gsi(t, h);
			}
		}
		else if (is_in_same_year[t] == 0 && is_in_same_year[h] != 0) {									// t in Gr, h in Gsi
		}
		else if (is_in_same_year[t] == 0 && is_in_same_year[h] == 0) {									// t in Gr, h in Gr
			myres.invalid_part_cnt++;
			mt = my_timer();
			is_in_same_year[h] = n2y[t];
			ahrsz_ext_priority_value<T> min_priority(plus_infinity);
			ahrsz_ext_priority_value<T> max_priority(minus_infinity);
			out_iterator i, iend;
			std::vector<unsigned int>::iterator it;
			for (it = ds.scc_nodes[h].begin(); it != ds.scc_nodes[h].end(); it++) {
				for (std::tie(i, iend) = boost::out_edges(*it, cg); i != iend; ++i) {
					//unsigned int w(target(*i, g));
					unsigned int w = ds.find(target(*i, cg));
					if (is_in_unorder[w] == false && is_in_same_year[w] == is_in_same_year[h]) {
						if (max_priority < ahrsz_ext_priority_value<T>(n2v[w].base(), *y_pspace)) {
							max_priority = n2v[w];
						}
					}
				}
			}
			if (max_priority == ahrsz_ext_priority_value<T>(minus_infinity)) {
				ahrsz_ext_priority_value<T> cur_h(minus_infinity), cur_t(y_pspace->begin(), *y_pspace);
				ahrsz_ext_priority_value<T> ep(compute_priority(cur_h, cur_t, *y_pspace));
				n2v[h] = ahrsz_priority_value<T>(ep.base(), *y_pspace);
			}
			else {
				ahrsz_ext_priority_value<T> ep(compute_priority(max_priority, min_priority, *y_pspace));
				n2v[h] = ahrsz_priority_value<T>(ep.base(), *y_pspace);
			}
			myres.tc_comp_time += mt.elapsed();
		}
		else if (is_in_same_year[t] != 0 && is_in_same_year[h] != 0) {									// t, h in Gsi, Gsi
			// apply dynamic algorithm in Gsi. 
			add_edge_Gsi(t, h);
		}
		else {
			assert(0 == 1 && "Type Error. t, h NOT in Gsi OR Gr.");
		}

	}
	else if (n2y[t] < n2y[h]) {																		// e.g., t in 2017, h in 2010. The cycle must be in Gsi & Gr. 
		myres.invalid_part_cnt++;
		mt = my_timer();
		std::vector<unsigned int> reachable, cycle, reachable2;
		std::vector<unsigned int>::iterator it, it_ds, i;
		std::vector<unsigned int>::reverse_iterator rit;
		T* u_space = (pspaces[0]);
		ahrsz_ext_priority_value<T> max_priority(minus_infinity);
		bool flag = false, is_first = true;
		// NOTE debug. current version.
		in_component[h] = true;
		scan_GsGr_without_topo(h, t, reachable, max_priority, flag);
		//in_component[t] = true;
		//scan_GsGr_without_topo(t, t, reachable, max_priority, flag);

		myres.tc_comp_time += mt.elapsed();
		if (flag == true) {
			mt = my_timer();
			for (rit = reachable.rbegin(); rit != reachable.rend(); rit++) {
				if (in_component[*rit] == false) {
					reachable2.push_back(*rit);
				}
				else {
					cycle.push_back(*rit);
					if (is_first) {
						reachable2.push_back(*rit);
						is_first = false;
					}
				}
				in_component[*rit] = false;
				_visited[*rit] = false;
				is_in_unorder[*rit] = true;
			}
			unsigned int s = cycle.front();
			for (i = cycle.begin(); i != cycle.end(); i++) {
				ds.join(s, *i);
			}
			s = ds.find(cycle.front());
			ds.collapse_scc(s, cycle);
			std::vector<unsigned int> cycle_nodes2;
			for (int i = 0; i < cycle.size(); i++) {
				cycle_nodes2.insert(cycle_nodes2.end(), ds.scc_nodes[cycle[i]].begin(), ds.scc_nodes[cycle[i]].end());
			}
			condense_per_scc(s, cycle_nodes2, cg, ds);
			myres.condensation_time += mt.elapsed();

			mt = my_timer();
			ahrsz_ext_priority_value<T> floor, ceiling;
			if (max_priority == ahrsz_ext_priority_value<T>(minus_infinity) || std::next(max_priority.base()) == u_space->end()) {
				floor = ahrsz_ext_priority_value<T>(u_space->previous(u_space->end()), *u_space);
				ceiling = ahrsz_ext_priority_value<T>(plus_infinity);
			}
			else {
				floor = ahrsz_ext_priority_value<T>(max_priority);
				ceiling = ahrsz_ext_priority_value<T>(std::next(max_priority.base()), *u_space);
			}
			for (rit = reachable2.rbegin(); rit != reachable2.rend(); rit++) {
				ahrsz_ext_priority_value<T> ep(compute_priority(floor, ceiling, *u_space));
				s = ds.find(*rit);
				n2v[s] = ahrsz_priority_value<T>(ep.base(), *u_space);
				floor = n2v[s];
				is_in_unorder[s] = true;
			}
			myres.tc_comp_time += mt.elapsed();
		}
		else {
			mt = my_timer();
			// update notGu nodes info.  
			ahrsz_ext_priority_value<T> floor, ceiling;
			if (max_priority == ahrsz_ext_priority_value<T>(minus_infinity) || std::next(max_priority.base()) == u_space->end()) {
				floor = ahrsz_ext_priority_value<T>(u_space->previous(u_space->end()), *u_space);
				ceiling = ahrsz_ext_priority_value<T>(plus_infinity);
			}
			else {
				floor = ahrsz_ext_priority_value<T>(max_priority);
				ceiling = ahrsz_ext_priority_value<T>(std::next(max_priority.base()), *u_space);
			}
			for (it = reachable.begin(); it != reachable.end(); it++) {
				is_in_unorder[*it] = true;
				in_component[*it] = false;
				_visited[*it] = false;
				ahrsz_ext_priority_value<T> ep(compute_priority(floor, ceiling, *u_space));
				n2v[*it] = ahrsz_priority_value<T>(ep.base(), *u_space);
				floor = n2v[*it];
			}
			myres.tc_comp_time += mt.elapsed();
		}
	}
	else {
		// only adding the edge(t, h) is OK!
		// the edge (t, h) neither created the cycle nor affected the Invariant of the Three types of the subgraphs
	}
}

template<class T>
void batDSCC_test<T>::scan_GsGr(vd_t n, std::vector<unsigned int>& reachable, ahrsz_ext_priority_value<T>& max_priority)
{
	_visited[n] = true;
	std::vector<unsigned int>::iterator it;
	out_iterator i, iend;
	for (it = ds.scc_nodes[n].begin(); it != ds.scc_nodes[n].end(); it++) {
		for (std::tie(i, iend) = boost::out_edges(*it, cg); i != iend; ++i) {
			//unsigned int w(target(*i, g));
			unsigned int w = ds.find(target(*i, cg));
			if (is_in_unorder[w] == true) {
				if (max_priority < ahrsz_ext_priority_value<T>(n2v[w].base(), *(pspaces[0]))) {
					max_priority = n2v[w];
				}
			}
			if (!_visited[w] && is_in_unorder[w] == false) {	// if w not visited & w not in Gm
				scan_GsGr(w, reachable, max_priority);
			}
		}
	}
	reachable.push_back(n);

}

template<class T>
void batDSCC_test<T>::scan_GsGr_without_topo(vd_t n, vd_t lb, std::vector<unsigned int>& reachable, ahrsz_ext_priority_value<T>& max_priority, bool & flag)
{
	_visited[n] = true;
	out_iterator i, iend;
	std::vector<unsigned int>::iterator it;
	for (it = ds.scc_nodes[n].begin(); it != ds.scc_nodes[n].end(); it++) {
		for (std::tie(i, iend) = boost::out_edges(*it, cg); i != iend; ++i) {
			//unsigned int w(target(*i, g));
			unsigned int w = ds.find(target(*i, cg));
			if (is_in_unorder[w] == true) {
				if (max_priority < ahrsz_ext_priority_value<T>(n2v[w].base(), *(pspaces[0]))) {
					max_priority = n2v[w];
				}
			}
			if (w == lb) {
				in_component[w] = true;
				if (!_visited[w] && is_in_unorder[w] == false) {
					scan_GsGr_without_topo(w, lb, reachable, max_priority, flag);
				}
				flag = true;
			}
			else if (!_visited[w] && is_in_unorder[w] == false) {							// w not visited & w is not in Gm. Because w in Gm, the w reachable nodes also in Gm. 
				scan_GsGr_without_topo(w, lb, reachable, max_priority, flag);				// Then the nodes that need to be reordered will not be traversed. 
			}
			in_component[n] = in_component[n] || in_component[w];
		}
	}
	reachable.push_back(n);

}




template<class T>
void batDSCC_test<T>::find_cycle_Gm(vd_t head)
{
	out_iterator i, iend;
	_visited[head] = true;
	std::vector<unsigned int>::iterator it;
	for (it = ds.scc_nodes[head].begin(); it != ds.scc_nodes[head].end(); it++) {
		for (std::tie(i, iend) = boost::out_edges(*it, cg); i != iend; ++i) {
			//vd_t w = (target(*i, cg));
			unsigned int  w = ds.find(target(*i, cg));
			if (is_in_unorder[w] && _inK[w] && _visited[w] == false)
				find_cycle_Gm(w);
			in_component[head] = in_component[head] || in_component[w];
		}
	}
}

template<class T>
void batDSCC_test<T>::find_cycle_Gsi(vd_t head)
{
	out_iterator i, iend;
	std::vector<unsigned int>::iterator it;
	_visited[head] = true;
	for (it = ds.scc_nodes[head].begin(); it != ds.scc_nodes[head].end(); it++) {
		for (std::tie(i, iend) = boost::out_edges(*it, cg); i != iend; ++i) {
			//vd_t w = (target(*i, cg));
			unsigned int  w = ds.find(target(*i, cg));
			if (is_in_same_year[head] == is_in_same_year[w]) {
				if (_inK[w] && _visited[w] == false) find_cycle_Gsi(w);
			}
			in_component[head] = in_component[head] || in_component[w];
		}
	}
}

template<class T>
void batDSCC_test<T>::discover_Gm(vd_t tail, vd_t head, std::vector<vd_t>& K, bool & has_cycle, std::vector<unsigned int>& cycle_nodes)
{
	std::vector<unsigned int>::iterator it;
	vd_t f(head), b(tail), tmp;
	max_priority_queue ForwFron(n2v);											// the forward frontier
	min_priority_queue BackFron(n2v);											// the backward frontier
	unsigned int ForwEdges = out_degree_scc(head, cg, ds);
	unsigned int BackEdges = in_degree_scc(tail, cg, ds);
	ForwFron.push(head);
	BackFron.push(tail);

	_visited[head] = true;														// notice, I uses visited to indicate when  a node is on one of the queues.				
	_visited[tail] = true;
	_inB[tail] = true;															// indicate tail and head has been visited when detecting the cycle.
	_inF[head] = true;

	while ((n2v[b] < n2v[f] || n2v[b] == n2v[f]) && !ForwFron.empty() && !BackFron.empty()) {
		unsigned int u = std::min(ForwEdges, BackEdges);
		ForwEdges -= u;
		BackEdges -= u;

		if (ForwEdges == 0) {
			if (!_inK[f]) {
				K.push_back(f);
				_inK[f] = true;
			}
			ForwFron.pop();
			_visited[f] = false;

			out_iterator i, iend;
			for (it = ds.scc_nodes[f].begin(); it != ds.scc_nodes[f].end(); it++) {
				for (tie(i, iend) = out_edges(*it, cg); i != iend; ++i) {
					//vd_t w(target(*i, *this));
					unsigned int w = ds.find(target(*i, cg));
					if (is_in_unorder[w] == true) {
						if (_inB[w]) {	// cycle detected. 
							has_cycle = true;
							if (_visited[w] == false || (_visited[w] == true && _inF[w] == false)) {
								ForwFron.push(w);
								_visited[w] = true;
							}
							_inF[w] = true;
						}
						if (!_inF[w]) {
							if (_visited[w] == false || (_visited[w] == true && _inF[w] == false)) {
								ForwFron.push(w);
								_visited[w] = true;
							}
							_inF[w] = true;
						}
					}
				}
			}
			if (ForwFron.empty()) {
				f = tail;
			}
			else {
				f = ForwFron.top();
			}
			ForwEdges = out_degree_scc(f, cg, ds);
		}
		if (BackEdges == 0) {
			if (!_inK[b]) {
				K.push_back(b);
				_inK[b] = true;
			}
			BackFron.pop();
			_visited[b] = false;

			in_iterator i, iend;
			for (it = ds.scc_nodes[b].begin(); it != ds.scc_nodes[b].end(); it++) {
				for (tie(i, iend) = in_edges(*it, cg); i != iend; ++i) {
					//vd_t w(source(*i, *this));
					unsigned int w = ds.find(source(*i, cg));
					if (is_in_unorder[w] == true) {
						if (_inF[w]) {																		// cycle detected. 
							has_cycle = true;
							if (_visited[w] == false || (_visited[w] == true && _inB[w] == false)) {
								BackFron.push(w);
								_visited[w] = true;
							}
							_inB[w] = true;
						}
						if (!_inB[w]) {
							if (_visited[w] == false || (_visited[w] == true && _inB[w] == false)) {
								BackFron.push(w);
								_visited[w] = true;
							}
							_inB[w] = true;
						}
					}
				}
			}
			if (BackFron.empty()) {
				b = head;
			}
			else {
				b = BackFron.top();
			}
			//BackEdges = in_degree(b, *this);
			BackEdges = in_degree_scc(b, cg, ds);
		}
	}
	// Finally, unmark all remaining on queues
	while (!ForwFron.empty()) {
		tmp = ForwFron.top();
		_visited[tmp] = false;
		_inF[tmp] = false;
		ForwFron.pop();
	}
	while (!BackFron.empty()) {
		tmp = BackFron.top();
		_visited[tmp] = false;
		_inB[tmp] = false;
		BackFron.pop();
	}

	if (has_cycle) {
		in_component[tail] = true;
		find_cycle_Gm(head);
		for (typename std::vector<vd_t>::iterator i = K.begin(); i != K.end(); i++) {
			if (in_component[*i]) {
				cycle_nodes.push_back(*i);
				in_component[*i] = false;
			}
		}
	}
}

template<class T>
void batDSCC_test<T>::discover_Gsi(vd_t tail, vd_t head, std::vector<vd_t>& K, bool & has_cycle, std::vector<unsigned int>& cycle_nodes) {
	std::vector<unsigned int>::iterator it;
	vd_t f(head), b(tail), tmp;
	max_priority_queue ForwFron(n2v);														// the forward frontier
	min_priority_queue BackFron(n2v);														// the backward frontier
	unsigned int ForwEdges = out_degree_scc(head, cg, ds);
	unsigned int BackEdges = in_degree_scc(tail, cg, ds);
	ForwFron.push(head);
	BackFron.push(tail);

	_visited[head] = true;																	// notice, I uses visited to indicate when  a node is on one of the queues.
	_visited[tail] = true;
	_inB[tail] = true;																		// indicate tail and head has been visited when detecting the cycle.
	_inF[head] = true;

	while ((n2v[b] < n2v[f] || n2v[b] == n2v[f]) && !ForwFron.empty() && !BackFron.empty()) {
		unsigned int u = std::min(ForwEdges, BackEdges);
		ForwEdges -= u;
		BackEdges -= u;

		if (ForwEdges == 0) {
			if (!_inK[f]) {
				K.push_back(f);
				_inK[f] = true;
			}
			ForwFron.pop();
			_visited[f] = false;

			out_iterator i, iend;
			for (it = ds.scc_nodes[f].begin(); it != ds.scc_nodes[f].end(); it++) {
				for (tie(i, iend) = out_edges(*it, cg); i != iend; ++i) {
					//vd_t w(target(*i, *this));
					unsigned int w = ds.find(target(*i, cg));
					if (!is_in_unorder[w] && is_in_same_year[w] == is_in_same_year[tail]) {
						if (_inB[w]) {	// cycle detected. 
							has_cycle = true;
							if (_visited[w] == false || (_visited[w] == true && _inF[w] == false)) {
								ForwFron.push(w);
								_visited[w] = true;
							}
							_inF[w] = true;
						}
						if (!_inF[w]) {
							if (_visited[w] == false || (_visited[w] == true && _inF[w] == false)) {
								ForwFron.push(w);
								_visited[w] = true;
							}
							_inF[w] = true;
						}

					}
				}
			}
			if (ForwFron.empty()) {
				f = tail;
			}
			else {
				f = ForwFron.top();
			}
			//ForwEdges = out_degree(f, *this);
			ForwEdges = out_degree_scc(f, cg, ds);
		}
		if (BackEdges == 0) {
			if (!_inK[b]) {
				K.push_back(b);
				_inK[b] = true;
			}
			BackFron.pop();
			_visited[b] = false;

			in_iterator i, iend;
			for (it = ds.scc_nodes[b].begin(); it != ds.scc_nodes[b].end(); it++) {
				for (tie(i, iend) = in_edges(*it, cg); i != iend; ++i) {
					//vd_t w(source(*i, *this));
					unsigned int w = ds.find(source(*i, cg));
					if (!is_in_unorder[w] && is_in_same_year[w] == is_in_same_year[tail]) {
						if (_inF[w]) {	// cycle detected. 
							has_cycle = true;
							if (_visited[w] == false || (_visited[w] == true && _inB[w] == false)) {
								BackFron.push(w);
								_visited[w] = true;
							}
							_inB[w] = true;
						}
						if (!_inB[w]) {
							if (_visited[w] == false || (_visited[w] == true && _inB[w] == false)) {
								BackFron.push(w);
								_visited[w] = true;
							}
							_inB[w] = true;
						}
					}
				}
			}
			if (BackFron.empty()) {
				b = head;
			}
			else {
				b = BackFron.top();
			}
			//BackEdges = in_degree(b, *this);
			BackEdges = in_degree_scc(b, cg, ds);
		}
	}
	// Finally, unmark all remaining on queues
	while (!ForwFron.empty()) {
		tmp = ForwFron.top();
		_visited[tmp] = false;
		_inF[tmp] = false;
		ForwFron.pop();
	}
	while (!BackFron.empty()) {
		tmp = BackFron.top();
		_visited[tmp] = false;
		_inB[tmp] = false;
		BackFron.pop();
	}

	if (has_cycle) {
		in_component[tail] = true;
		find_cycle_Gsi(head);
		for (typename std::vector<vd_t>::iterator i = K.begin(); i != K.end(); i++) {
			if (in_component[*i]) {
				cycle_nodes.push_back(*i);
				in_component[*i] = false;
			}
		}
	}
}

template<class T>
void batDSCC_test<T>::maintain_Gm(std::vector<vd_t>& K)
{
	// Initialise temporary data 
	for (typename std::vector<vd_t>::iterator i(K.begin()); i != K.end(); ++i) {
		_flooring[*i] = minus_infinity;
		_inK[*i] = true;
		_inB[*i] = false;								// unmark states. 
		_inF[*i] = false;
	}
	// first pass - compute flooring Information
	std::vector<vd_t> rto; // reverse topological order
	for (typename std::vector<vd_t>::iterator i(K.begin()); i != K.end(); ++i) {
		if (!_visited[*i]) { compute_flooring_Gm(*i, rto); }
	}

	// second pass - perform the reassignment
	for (typename std::vector<vd_t>::reverse_iterator i(rto.rbegin()); i != rto.rend(); ++i) {
		ahrsz_ext_priority_value<T> ceiling;
		if (_flooring[*i].minus_inf()) {
			ceiling = ahrsz_ext_priority_value<T>(pspaces[0]->begin(), *pspaces[0]);
		}
		else {
			if (_flooring[*i].base() == pspaces[0]->end() || std::next(_flooring[*i].base()) == pspaces[0]->end()) {
				ceiling = ahrsz_ext_priority_value<T>(plus_infinity);
			}
			else {
				ceiling = ahrsz_ext_priority_value<T>(std::next(_flooring[*i].base()), *pspaces[0]);
			}
		}
		//ahrsz_ext_priority_value<T> ep(compute_priority(_flooring[*i], compute_ceiling_Gm(*i), *(pspaces[0])));
		ahrsz_ext_priority_value<T> ep(compute_priority(_flooring[*i], ceiling, *(pspaces[0])));
		n2v[*i] = ahrsz_priority_value<T>(ep.base(), *(pspaces[0]));
		_visited[*i] = false; // unmark state info. 
	}
	// reset visited information	unmark state info. 
	for (typename std::vector<vd_t>::iterator i(K.begin()); i != K.end(); ++i) {
		_visited[*i] = false;
		_inK[*i] = false;
	}
}

template<class T>
void batDSCC_test<T>::maintain_Gsi(std::vector<vd_t>& K)
{
	// Initialise temporary data 
	for (typename std::vector<vd_t>::iterator i(K.begin()); i != K.end(); ++i) {
		_flooring[*i] = minus_infinity;
		_inK[*i] = true;
		_inB[*i] = false;		// unmark states. 
		_inF[*i] = false;
	}
	// first pass - compute flooring Information
	std::vector<vd_t> rto; // reverse topological order
	for (typename std::vector<vd_t>::iterator i(K.begin()); i != K.end(); ++i) {
		if (!_visited[*i]) { compute_flooring_Gsi(*i, rto); }
	}
	int idx = n2y[K.front()] - min_year + 1;
	T* y_pspace = (pspaces[idx]);

	// second pass - perform the reassignment
	for (typename std::vector<vd_t>::reverse_iterator i(rto.rbegin()); i != rto.rend(); ++i) {
		ahrsz_ext_priority_value<T> ceiling;
		if (_flooring[*i].minus_inf()) {
			ceiling = ahrsz_ext_priority_value<T>(y_pspace->begin(), *y_pspace);
		}
		else {
			if (_flooring[*i].base() == y_pspace->end() || std::next(_flooring[*i].base()) == y_pspace->end()) {
				ceiling = ahrsz_ext_priority_value<T>(plus_infinity);
			}
			else {
				ceiling = ahrsz_ext_priority_value<T>(std::next(_flooring[*i].base()), *y_pspace);
			}
		}
		ahrsz_ext_priority_value<T> ep(compute_priority(_flooring[*i], ceiling, *y_pspace));
		//ahrsz_ext_priority_value<T> ep(compute_priority(_flooring[*i], compute_ceiling_Gsi(*i), *y_pspace));

		n2v[*i] = ahrsz_priority_value<T>(ep.base(), *y_pspace);
		_visited[*i] = false; // unmark state info. 
	}
	// reset visited information	unmark state info. 
	for (typename std::vector<vd_t>::iterator i(K.begin()); i != K.end(); ++i) {
		_visited[*i] = false;
		_inK[*i] = false;
	}
}

template<class T>
void batDSCC_test<T>::compute_flooring_Gm(vd_t n, std::vector<vd_t>& rto)
{
	_visited[n] = true;
	std::vector<unsigned int>::iterator it;
	out_iterator i, iend;
	for (it = ds.scc_nodes[n].begin(); it != ds.scc_nodes[n].end(); it++) {
		for (tie(i, iend) = out_edges(*it, cg); i != iend; ++i) {
			/*vd_t j(target(*i, *this));*/
			unsigned int j = ds.find(target(*i, cg));
			if (is_in_unorder[j]) {
				if (_inK[j]) {
					if (!_visited[j]) { compute_flooring_Gm(j, rto); }
					_flooring[n] = std::max(_flooring[j], _flooring[n]);
				}
				else {
					// successor is not in K
					_flooring[n] = std::max(_flooring[n], ahrsz_ext_priority_value<T>(n2v[j]));
				}
			}
		}
	}
	rto.push_back(n);

}

template<class T>
void batDSCC_test<T>::compute_flooring_Gsi(vd_t n, std::vector<vd_t>& rto)
{
	_visited[n] = true;
	std::vector<unsigned int>::iterator it;
	out_iterator i, iend;
	for (it = ds.scc_nodes[n].begin(); it != ds.scc_nodes[n].end(); it++) {
		for (tie(i, iend) = out_edges(*it, cg); i != iend; ++i) {
			//vd_t j(target(*i, *this));
			unsigned int j = ds.find(target(*i, cg));
			if (!is_in_unorder[j] && is_in_same_year[j] == is_in_same_year[n]) {
				if (_inK[j]) {
					if (!_visited[j]) { compute_flooring_Gsi(j, rto); }
					_flooring[n] = std::max(_flooring[j], _flooring[n]);
				}
				else {
					// successor is not in K
					_flooring[n] = std::max(_flooring[n], ahrsz_ext_priority_value<T>(n2v[j]));
				}
			}
		}
	}
	rto.push_back(n);
}

template<class T>
ahrsz_ext_priority_value<T> batDSCC_test<T>::compute_priority(ahrsz_ext_priority_value<T> floor, ahrsz_ext_priority_value<T> ceiling, T & cur_pspace)
{
	ceiling = ahrsz_ext_priority_value<T>(plus_infinity);
	assert(floor < ceiling);
	ahrsz_ext_priority_value<T> candidate;
	if (floor.minus_inf()) {
		cur_pspace.push_front();
		candidate = ahrsz_ext_priority_value<T>(cur_pspace.begin(), cur_pspace);
	}
	else {
		assert(!floor.plus_inf());
		candidate = ahrsz_ext_priority_value<T>(cur_pspace.insert_after(floor.base()), cur_pspace);
	}
	// std::cout << "compute priority" << a2str(floor) << "\t" << a2str(candidate) <<"\t"<< a2str(ceiling) << std::endl; 
	assert(floor < candidate && candidate < ceiling && "create new priority error");
	return candidate;
}

template<class T>
ahrsz_ext_priority_value<T> batDSCC_test<T>::compute_ceiling_Gm(vd_t v)
{
	ahrsz_ext_priority_value<T> ceiling(plus_infinity);
	std::vector<unsigned int>::iterator it;
	in_iterator j, jend;
	for (it = ds.scc_nodes[v].begin(); it != ds.scc_nodes[v].end(); it++) {
		for (tie(j, jend) = in_edges(*it, cg); j != jend; ++j) {
			unsigned int w = ds.find(source(*j, cg));
			if (is_in_unorder[w]) {
				ceiling = std::min(ceiling, ahrsz_ext_priority_value<T>(n2v[w]));
			}
		}
	}
	return ceiling;
}

template<class T>
ahrsz_ext_priority_value<T> batDSCC_test<T>::compute_ceiling_Gsi(vd_t v)
{
	ahrsz_ext_priority_value<T> ceiling(plus_infinity);
	std::vector<unsigned int>::iterator it;
	in_iterator j, jend;
	for (it = ds.scc_nodes[v].begin(); it != ds.scc_nodes[v].end(); it++) {
		for (tie(j, jend) = in_edges(*it, cg); j != jend; ++j) {
			unsigned int w = ds.find(source(*j, cg));
			if (is_in_same_year[w] == is_in_same_year[v]) {
				ceiling = std::min(ceiling, ahrsz_ext_priority_value<T>(n2v[w]));
			}
		}
	}
	return ceiling;
}


