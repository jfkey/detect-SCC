#include "dynamic_HKMST_TC_scc.h"

template<typename T>
dynamic_tc_hkmst_scc<T>::dynamic_tc_hkmst_scc()
{ }

template<typename T>
dynamic_tc_hkmst_scc<T>::dynamic_tc_hkmst_scc(CG & cg, vector<ahrsz_priority_value<T>>& _n2v, vector<T*>& _pspaces, disjoint_set & _ds, vector<unsigned int>& _n2y, vector<int>& _is_in_same_year, vector<bool>& _is_in_unorder, int _min_year, int _max_year, std::vector<pair<unsigned int, unsigned int>>& _inc_edges, const int np)
	: cg(cg), n2v(_n2v), pspaces(_pspaces), ds(_ds), n2y(_n2y), is_in_same_year(_is_in_same_year), is_in_unorder(_is_in_unorder), 
	min_year(_min_year), max_year(_max_year), inc_edges(_inc_edges), out_its(np), in_its( np ), is_out(np, false), is_in(np, false) {
	std::cout << "Init Time-Coupled Dynamic HKMST Algorithm" << std::endl;
	self_name = "Dynamic-TC-HKMST-";
	myres.set_dynamic_name(self_name);
	_flooring = std::vector<ahrsz_ext_priority_value<T> >(np, minus_infinity);
	_visited = std::vector<bool>(np, false);
	_inK = std::vector<bool>(np, false);
	_inF = std::vector<bool>(np, false);
	_inB = std::vector<bool>(np, false);
	in_component = std::vector<bool>(np, false);
}



template<typename T>
void dynamic_tc_hkmst_scc<T>::add_edges()
{
	my_timer mt = my_timer();
	std::vector<pair<unsigned int, unsigned int>>::iterator it;
	for (it = inc_edges.begin(); it != inc_edges.end(); it++) {
		unsigned int s, t;
		int type = -1;
		s = ds.find(it->first);
		t = ds.find(it->second);
		if (s == t) continue;
		my_timer mt = my_timer();
		boost::add_edge(it->first, it->second, cg);										// NOTE: add_edge(it->first, it->second, cg); not add_edge(s, t, cg)
		myres.add_edge_time += mt.elapsed();
		type = get_alg_type(s, t, is_in_same_year, is_in_unorder);
		if (type == 1) {																// s, t in Gsi & Gr, Gsi & Gr 
			add_edge_Gm(s, t);
		}
		else if (type == 2) {															// s, t in Gsi & Gr, Gm // algorithm1
			if (n2y[s] < n2y[t]) {
				add_edge_alg1(s, t);
			}
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
	myres.dn_eval_time += myres.search_time + myres.reassignment_time + myres.tc_comp_time;
	myres.dn_eval_time_AAN = mt.elapsed();
}

template<typename T>
void dynamic_tc_hkmst_scc<T>::add_edges(int start_index, int end_index)
{
	myres.inc_edges += (end_index - start_index);
	my_timer mt = my_timer();
	for (int i = start_index; i < end_index; i++) {
		unsigned int s, t;
		int type = -1;
		s = ds.find(inc_edges[i].first);
		t = ds.find(inc_edges[i].second);
		if (s == t) continue;
		my_timer mt = my_timer();
		boost::add_edge(inc_edges[i].first, inc_edges[i].second, cg);										// NOTE: add_edge(it->first, it->second, cg); not add_edge(s, t, cg)
		myres.add_edge_time += mt.elapsed();
		type = get_alg_type(s, t, is_in_same_year, is_in_unorder);
		if (type == 1) {																// s, t in Gsi & Gr, Gsi & Gr 
			add_edge_Gm(s, t);
		}
		else if (type == 2) {															// s, t in Gsi & Gr, Gm // algorithm1
			/*if (n2y[s] < n2y[t]) {
				add_edge_alg1(s, t);
			}*/
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
	myres.dn_eval_time += myres.search_time + myres.reassignment_time + myres.tc_comp_time;
	record_cur_per_info(myres);
}




template<typename T>
void dynamic_tc_hkmst_scc<T>::add_edge_Gm(vd_t t, vd_t h)
{
	if ((n2v[t] < n2v[h])) {				// the added edge: t -> h, i.e,  pmap(t) > pmap(h). when not holds this priority, perform discovery and reassignment pharses. 
		myres.set_invalid_edges(myres.get_invalid_edges() + 1);
		std::list<vd_t>::iterator it;
		vd_t threshold;

		std::list<vd_t> FA, FP, BA, BP, forw_nodes, back_nodes, F_greater, B_less;
		std::vector<vd_t>  remaining_K;
		std::vector<unsigned int> cycle_nodes;
		bool has_cycle = false;
		my_timer mt;
		soft_threshold_search_Gm(t, h, forw_nodes, back_nodes, FA, FP, BA, BP, has_cycle, cycle_nodes, threshold);
		myres.set_search_time(myres.get_search_time() + mt.elapsed());

		if (has_cycle) {
			mt = my_timer();
			assert(cycle_nodes.size() > 1);// && "the number of vertices in the cycle must be greater than 1. dynamic_hkmst_scc"
			unsigned int s = cycle_nodes.front();
			for (std::vector<unsigned int>::iterator it = cycle_nodes.begin(); it != cycle_nodes.end(); it++) {
				ds.join(s, *it);
			}
			s = ds.find(cycle_nodes.front());
			
			remaining_K.push_back(s);
			auto fun = [this, &remaining_K, s](vd_t t) {
				in_component[t] = false;	// unmark state info
				_visited[t] = false;
				_inB[t] = false;
				_inF[t] = false;
				is_in[t] = false;
				is_out[t] = false;
				if (ds.find(t) != ds.find(s)) {
					remaining_K.push_back(t);
				}
			};

			std::vector<vd_t> F_greater, B_less;
			vd_t threshold = t;

			for (it = forw_nodes.begin(); it != forw_nodes.end(); it++) {
				if (n2v[threshold] < n2v[*it] && has_out_after_search_Gm(*it)) { // is_in_unorder[*it] &&  
					threshold = *it;
				}
			}
			// F> 
			for (it = forw_nodes.begin(); it != forw_nodes.end(); it++) {
				if (n2v[threshold] < n2v[*it])			// is_in_unorder[*it] && 
					F_greater.push_back(*it);	// unmark states. 
				_inF[*it] = false;
				is_out[*it] = false;
			}
			// B< 
			for (it = back_nodes.begin(); it != back_nodes.end(); it++) {
				if (n2v[*it] < n2v[threshold])		//is_in_unorder[*it]  && 
					B_less.push_back(*it);
				_inB[*it] = false;
				is_in[*it] = false;
			}
			F_greater.push_back(t);
			std::for_each(F_greater.begin(), F_greater.end(), fun);
			std::for_each(B_less.begin(), B_less.end(), fun);

			// unmark state info
			auto fun2 = [this](vd_t t) {
				in_component[t] = false;
				_visited[t] = false;
				_inB[t] = false;
				_inF[t] = false;
				is_in[t] = false;
				is_out[t] = false;
			};
			std::for_each(forw_nodes.begin(), forw_nodes.end(), fun2);
			std::for_each(back_nodes.begin(), back_nodes.end(), fun2);

			// collapse scc operation could affect out_its and in_its.
			ds.collapse_scc(s, cycle_nodes);
			std::vector<unsigned int> cycle_nodes2;
			for (int i = 0; i < cycle_nodes.size(); i++) {
				cycle_nodes2.insert(cycle_nodes2.end(), ds.scc_nodes[cycle_nodes[i]].begin(), ds.scc_nodes[cycle_nodes[i]].end());
			}
			condense_per_scc(s, cycle_nodes2, cg, ds);
			myres.condensation_time += mt.elapsed();

			myres.set_aff_region(myres.get_aff_region() + remaining_K.size());
			my_timer mt;
			maintain_SCC_Gm(remaining_K);
			myres.set_reassign_time(myres.get_reassign_time() + mt.elapsed());
		}
		else {
			// reassignment for affected region without the cycle.  
			// threshold = max{v, F| x in F and out(x) is nuot null . 
			// threshold, F>  B< 
			//threshold = t; 
			my_timer mt;
			maintain_Gm(forw_nodes, back_nodes, t, h);
			myres.set_reassign_time(myres.get_reassign_time() + mt.elapsed());

		}
	}

}

template<typename T>
void dynamic_tc_hkmst_scc<T>::add_edge_Gsi(vd_t t, vd_t h)
{
	if ((n2v[t] < n2v[h])) {				// the added edge: t -> h, i.e,  pmap(t) > pmap(h). when not holds this priority, perform discovery and reassignment pharses. 
		myres.set_invalid_edges(myres.get_invalid_edges() + 1);
		std::list<vd_t>::iterator it;
		vd_t threshold;

		std::list<vd_t> FA, FP, BA, BP, forw_nodes, back_nodes, F_greater, B_less;
		std::vector<vd_t>  remaining_K;
		std::vector<unsigned int> cycle_nodes;
		bool has_cycle = false;
		my_timer mt;
		soft_threshold_search_Gsi(t, h, forw_nodes, back_nodes, FA, FP, BA, BP, has_cycle, cycle_nodes, threshold);
		myres.set_search_time(myres.get_search_time() + mt.elapsed());

		if (has_cycle) {
			mt = my_timer();
			assert(cycle_nodes.size() > 1);// && "the number of vertices in the cycle must be greater than 1. dynamic_hkmst_scc"
			unsigned int s = cycle_nodes.front();
			for (std::vector<unsigned int>::iterator it = cycle_nodes.begin(); it != cycle_nodes.end(); it++) {
				ds.join(s, *it);
			}
			s = ds.find(cycle_nodes.front());

			remaining_K.push_back(s);
			auto fun = [this, &remaining_K, s](vd_t t) {
				in_component[t] = false;	// unmark state info
				_visited[t] = false;
				_inB[t] = false;
				_inF[t] = false;
				is_in[t] = false;
				is_out[t] = false;
				if (ds.find(t) != ds.find(s)) {
					remaining_K.push_back(t);
				}
			};

			std::vector<vd_t> F_greater, B_less;
			vd_t threshold = t;

			for (it = forw_nodes.begin(); it != forw_nodes.end(); it++) {
				if (n2v[threshold] < n2v[*it] && has_out_after_search_Gsi(*it)) { // is_in_unorder[*it] &&  
					threshold = *it;
				}
			}
			// F> 
			for (it = forw_nodes.begin(); it != forw_nodes.end(); it++) {
				if (n2v[threshold] < n2v[*it])			// is_in_unorder[*it] && 
					F_greater.push_back(*it);	// unmark states. 
				_inF[*it] = false;
				is_out[*it] = false;
			}
			// B< 
			for (it = back_nodes.begin(); it != back_nodes.end(); it++) {
				if (n2v[*it] < n2v[threshold])		//is_in_unorder[*it]  && 
					B_less.push_back(*it);
				_inB[*it] = false;
				is_in[*it] = false;
			}
			F_greater.push_back(t);
			std::for_each(F_greater.begin(), F_greater.end(), fun);
			std::for_each(B_less.begin(), B_less.end(), fun);

			// unmark state info
			auto fun2 = [this](vd_t t) {
				in_component[t] = false;
				_visited[t] = false;
				_inB[t] = false;
				_inF[t] = false;
				is_in[t] = false;
				is_out[t] = false;
			};
			std::for_each(forw_nodes.begin(), forw_nodes.end(), fun2);
			std::for_each(back_nodes.begin(), back_nodes.end(), fun2);

			// collapse scc operation could affect out_its and in_its.
			ds.collapse_scc(s, cycle_nodes);
			std::vector<unsigned int> cycle_nodes2;
			for (int i = 0; i < cycle_nodes.size(); i++) {
				cycle_nodes2.insert(cycle_nodes2.end(), ds.scc_nodes[cycle_nodes[i]].begin(), ds.scc_nodes[cycle_nodes[i]].end());
			}
			condense_per_scc(s, cycle_nodes2, cg, ds);
			myres.condensation_time += mt.elapsed();

			myres.set_aff_region(myres.get_aff_region() + remaining_K.size());
			my_timer mt;
			maintain_SCC_Gsi(remaining_K);
			myres.set_reassign_time(myres.get_reassign_time() + mt.elapsed());
		}
		else {
			// reassignment for affected region without the cycle.  
			// threshold = max{v, F| x in F and out(x) is nuot null . 
			// threshold, F>  B< 
			//threshold = t; 
			my_timer mt;
			maintain_Gsi(forw_nodes, back_nodes, t, h);
			myres.set_reassign_time(myres.get_reassign_time() + mt.elapsed());

		}
	}
}

template<typename T>
void dynamic_tc_hkmst_scc<T>::add_edge_alg1(vd_t t, vd_t h)
{
	my_timer mt;
	ahrsz_ext_priority_value<T>  max_priority(minus_infinity);
	std::vector<unsigned int> reachable;
	std::vector<unsigned int>::iterator it, it_ds;
	//scan_GsGr(t, reachable, max_priority );
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

template<typename T>
void dynamic_tc_hkmst_scc<T>::add_edge_alg2(vd_t t, vd_t h)
{
	my_timer mt; 
	ahrsz_ext_priority_value<T> max_priority(minus_infinity);
	std::vector<unsigned int> reachable;
	std::vector<unsigned int>::iterator it, it_ds;
	scan_GsGr(h, reachable, max_priority );
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
			ahrsz_ext_priority_value<T> ep( compute_priority(floor, ceiling, *(pspaces[0])));
			n2v[*it] = ahrsz_priority_value<T>(ep.base(), *(pspaces[0]));
			floor = n2v[*it];
		}
		myres.tc_comp_time += mt.elapsed();
		// t -> h is the invalid edge. 
		add_edge_Gm(t, h );
	}
}
 

template<class T>
void dynamic_tc_hkmst_scc<T>::add_edge_alg3(vd_t t, vd_t h)												// w.r.t. dfs start with t in Static-TC
{
	int yidx = n2y[t] - min_year + 1;
	T *y_pspace = (pspaces[yidx]);
	my_timer mt;
	if (n2y[t] == n2y[h]) {
		if (is_in_same_year[t] != 0 && is_in_same_year[h] == 0) {										// t in Gsi, h in Gr
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
		mt= my_timer();
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

template<typename T>
void dynamic_tc_hkmst_scc<T>::scan_GsGr(vd_t n, std::vector<unsigned int>& reachable, ahrsz_ext_priority_value<T>& max_priority)
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

template<typename T>
void dynamic_tc_hkmst_scc<T>::scan_GsGr_without_topo(vd_t n, vd_t lb, std::vector<unsigned int>& reachable, ahrsz_ext_priority_value<T>& max_priority, bool & flag)
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



template<typename T>
void dynamic_tc_hkmst_scc<T>::find_cycle_Gm(vd_t head)
{
	out_iterator i, iend;
	_visited[head] = true;
	std::vector<unsigned int>::iterator it;
	for (it = ds.scc_nodes[head].begin(); it != ds.scc_nodes[head].end(); it++) {
		for (std::tie(i, iend) = boost::out_edges(*it, cg); i != iend; ++i) {
			//vd_t w = (target(*i, cg));
			unsigned int  w = ds.find(target(*i, cg));
			if (is_in_unorder[w] && (_inB[w] || _inF[w]) && _visited[w] == false)
				find_cycle_Gm(w);
			in_component[head] = in_component[head] || in_component[w];
		}
	}
}

template<typename T>
void dynamic_tc_hkmst_scc<T>::find_cycle_Gsi(vd_t head)
{
	out_iterator i, iend;
	std::vector<unsigned int>::iterator it;
	_visited[head] = true;
	for (it = ds.scc_nodes[head].begin(); it != ds.scc_nodes[head].end(); it++) {
		for (std::tie(i, iend) = boost::out_edges(*it, cg); i != iend; ++i) {
			//vd_t w = (target(*i, cg));
			unsigned int  w = ds.find(target(*i, cg));
			if ( !is_in_unorder[w] &&  (is_in_same_year[head] == is_in_same_year[w]) && (_inB[w] || _inF[w]) && _visited[w] == false) find_cycle_Gsi(w);
			in_component[head] = in_component[head] || in_component[w];
		}
	}
}

template<typename T>
void dynamic_tc_hkmst_scc<T>::search_step_Gm(vd_t & u, vd_t & z, std::list<vd_t>& FA, std::list<vd_t>& BA, std::list<vd_t>& forw_nodes, std::list<vd_t>& back_nodes, bool & has_cycle)
{
	out_scc_iterator out_u, out_x;
	in_scc_iterator in_z, in_y;
	vd_t out_u_x, in_z_y;
	bool has_no = false, has_ni = false;
	out_u = first_out_Gm(u);
	in_z = first_in_Gm(z);
	out_u_x = ds.find(boost::target(*(out_u.oei_begin), cg));
	in_z_y = ds.find(boost::source(*(in_z.iei_begin), cg));
	has_no = has_next_out_Gm(out_u);
	has_ni = has_next_in_Gm(in_z);

	if (has_no == false)
		FA.remove(u);
	if (has_ni == false)
		BA.remove(z);

	if (_inF[in_z_y] || _inB[out_u_x]) has_cycle = true;		// report the cycle . 

	if (!_inF[out_u_x]) {
		_inF[out_u_x] = true;
		forw_nodes.push_back(out_u_x);
		out_x = first_out_Gm(out_u_x);
		//if (out_x.first != out_x.second)
		if (out_tail_Gm(out_x)) {
			FA.push_back(out_u_x);
		}
		is_out[out_u_x] = false;
	}
	if (!_inB[in_z_y]) {
		_inB[in_z_y] = true;
		back_nodes.push_back(in_z_y);
		in_y = first_in_Gm(in_z_y);
		//if (in_y.first != in_y.second) {
		if (in_tail_Gm(in_y)) {
			BA.push_back(in_z_y);
		}
		is_in[in_z_y] = false;
	}

}

template<typename T>
void dynamic_tc_hkmst_scc<T>::search_step_Gsi(vd_t & u, vd_t & z, std::list<vd_t>& FA, std::list<vd_t>& BA, std::list<vd_t>& forw_nodes, std::list<vd_t>& back_nodes, bool & has_cycle)
{
	out_scc_iterator out_u, out_x;
	in_scc_iterator in_z, in_y;
	vd_t out_u_x, in_z_y;
	bool has_no = false, has_ni = false;
	out_u = first_out_Gsi(u);
	in_z = first_in_Gsi(z);
	out_u_x = ds.find(boost::target(*(out_u.oei_begin), cg));
	in_z_y = ds.find(boost::source(*(in_z.iei_begin), cg));
	has_no = has_next_out_Gsi(out_u, u);
	has_ni = has_next_in_Gsi(in_z, z);

	if (has_no == false)
		FA.remove(u);
	if (has_ni == false)
		BA.remove(z);

	if (_inF[in_z_y] || _inB[out_u_x]) has_cycle = true;		// report the cycle . 

	if (!_inF[out_u_x]) {
		_inF[out_u_x] = true;
		forw_nodes.push_back(out_u_x);
		out_x = first_out_Gsi(out_u_x);
		//if (out_x.first != out_x.second)
		if (out_tail_Gsi(out_x, u)) {
			FA.push_back(out_u_x);
		}
		is_out[out_u_x] = false;
	}
	if (!_inB[in_z_y]) {
		_inB[in_z_y] = true;
		back_nodes.push_back(in_z_y);
		in_y = first_in_Gsi(in_z_y);
		//if (in_y.first != in_y.second) {
		if (in_tail_Gsi(in_y, z)) {
			BA.push_back(in_z_y);
		}
		is_in[in_z_y] = false;
	}
}

template<typename T>
void dynamic_tc_hkmst_scc<T>::soft_threshold_search_Gm(vd_t & v, vd_t & w, std::list<vd_t>& forw_nodes, std::list<vd_t>& back_nodes, std::list<vd_t>& FA, std::list<vd_t>& FP, std::list<vd_t>& BA, std::list<vd_t>& BP, bool & has_cycle, std::vector<unsigned int>& cycle_nodes, vd_t & threshold)
{
		
	std::list<vd_t>::iterator it_list;
	threshold = v;

	forw_nodes.push_back(w);
	_inF[w] = true;
	back_nodes.push_back(v);
	_inB[v] = true;

	if (has_first_out_Gm(w)) {
		FA.push_back(w);
	}
	if (has_first_in_Gm(v)) {
		BA.push_back(v);
	}
	while (!FA.empty() && !BA.empty()) {
		vd_t u = FA.front();
		vd_t z = BA.front();
		if (n2v[z] < n2v[u] || n2v[z] == n2v[u]) { // n2v[z] == n2v[u] to identify the scc. 
			search_step_Gm(u, z, FA, BA, forw_nodes, back_nodes, has_cycle);
		}
		else {
			if (n2v[u] < n2v[threshold]) {
				FA.pop_front();
				FP.push_back(u);
			}
			if (n2v[threshold] < n2v[z]) {
				BA.pop_front();
				BP.push_back(z);
			}
			if (u == z && z == threshold) {
				break;
			}
		}
		if (FA.empty()) {
			BP.clear();
			BA.remove(threshold);							// NOTE: to be refine ?

			if (!FP.empty()) {
				threshold = FP.front();
				it_list = FP.begin();
				while (it_list != FP.end()) {
					if ((n2v[threshold] < n2v[*it_list]) || (n2v[*it_list] == n2v[threshold])) {
						FA.push_back(*it_list);
						FP.erase(it_list++);				// erase current and get next
					}
					else {
						++it_list;
					}
				}
			}
		}
		if (BA.empty()) {
			FP.clear();
			FA.remove(threshold);
			if (!BP.empty()) {
				threshold = BP.front();
				it_list = BP.begin();
				while (it_list != BP.end()) {
					if ((n2v[*it_list] < n2v[threshold]) || (n2v[*it_list] == n2v[threshold])) {
						BA.push_back(*it_list);
						BP.erase(it_list++);
					}
					else {
						++it_list;
					}
				}
			}
		}
	}

	if (has_cycle) {
		in_component[v] = true;
		find_cycle_Gm(w);
		auto fun = [this, &cycle_nodes](vd_t t) {
			if (this->in_component[t]) {
				cycle_nodes.push_back(t);
				this->in_component[t] = false;
			}
		};
		std::for_each(forw_nodes.begin(), forw_nodes.end(), fun);
		std::for_each(back_nodes.begin(), back_nodes.end(), fun);
	}
}

template<typename T>
void dynamic_tc_hkmst_scc<T>::soft_threshold_search_Gsi(vd_t & v, vd_t & w, std::list<vd_t>& forw_nodes, std::list<vd_t>& back_nodes, std::list<vd_t>& FA, std::list<vd_t>& FP, std::list<vd_t>& BA, std::list<vd_t>& BP, bool & has_cycle, std::vector<unsigned int>& cycle_nodes, vd_t & threshold)
{
	std::list<vd_t>::iterator it_list;
	threshold = v;

	forw_nodes.push_back(w);
	_inF[w] = true;
	back_nodes.push_back(v);
	_inB[v] = true;

	if (has_first_out_Gsi(w)) {
		FA.push_back(w);
	}
	if (has_first_in_Gsi(v)) {
		BA.push_back(v);
	}
	while (!FA.empty() && !BA.empty()) {
		vd_t u = FA.front();
		vd_t z = BA.front();
		if (n2v[z] < n2v[u] || n2v[z] == n2v[u]) { // n2v[z] == n2v[u] to identify the scc. 
			search_step_Gsi(u, z, FA, BA, forw_nodes, back_nodes, has_cycle);
		}
		else {
			if (n2v[u] < n2v[threshold]) {
				FA.pop_front();
				FP.push_back(u);
			}
			if (n2v[threshold] < n2v[z]) {
				BA.pop_front();
				BP.push_back(z);
			}
			if (u == z && z == threshold) {
				break;
			}
		}
		if (FA.empty()) {
			BP.clear();
			BA.remove(threshold);							// NOTE: to be refine ?

			if (!FP.empty()) {
				threshold = FP.front();
				it_list = FP.begin();
				while (it_list != FP.end()) {
					if ((n2v[threshold] < n2v[*it_list]) || (n2v[*it_list] == n2v[threshold])) {
						FA.push_back(*it_list);
						FP.erase(it_list++);				// erase current and get next
					}
					else {
						++it_list;
					}
				}
			}
		}
		if (BA.empty()) {
			FP.clear();
			FA.remove(threshold);
			if (!BP.empty()) {
				threshold = BP.front();
				it_list = BP.begin();
				while (it_list != BP.end()) {
					if ((n2v[*it_list] < n2v[threshold]) || (n2v[*it_list] == n2v[threshold])) {
						BA.push_back(*it_list);
						BP.erase(it_list++);
					}
					else {
						++it_list;
					}
				}
			}
		}
	}

	if (has_cycle) {
		in_component[v] = true;
		find_cycle_Gsi(w);
		auto fun = [this, &cycle_nodes](vd_t t) {
			if (this->in_component[t]) {
				cycle_nodes.push_back(t);
				this->in_component[t] = false;
			}
		};
		std::for_each(forw_nodes.begin(), forw_nodes.end(), fun);
		std::for_each(back_nodes.begin(), back_nodes.end(), fun);
	}

}

template<typename T>
void dynamic_tc_hkmst_scc<T>::maintain_Gm(std::list<vd_t>& forw_nodes, std::list<vd_t>& back_nodes, vd_t & v, vd_t & w)
{
	std::list<vd_t> F_greater, B_less;
	std::list<vd_t>::iterator it;
	std::list<vd_t>::reverse_iterator rit;
	ahrsz_ext_priority_value<T> val1, val2;
	vd_t threshold = v;

	for (it = forw_nodes.begin(); it != forw_nodes.end(); it++) {
		if (n2v[threshold] < n2v[*it] && has_out_after_search_Gm(*it)) {
			threshold = *it;
		}
	}
	// F> 
	for (it = forw_nodes.begin(); it != forw_nodes.end(); it++) {
		if (n2v[threshold] < n2v[*it])
			F_greater.push_back(*it);	// unmark states. 
		_inF[*it] = false;
		is_out[*it] = false;
	}
	// B< 
	for (it = back_nodes.begin(); it != back_nodes.end(); it++) {
		if (n2v[*it] < n2v[threshold])
			B_less.push_back(*it);
		_inB[*it] = false;
		is_in[*it] = false;
	}
	auto compare1 = [this](const vd_t &a, const vd_t &b) { return n2v[b] < n2v[a]; };
	auto compare2 = [this](const vd_t &a, const vd_t &b) { return n2v[a] < n2v[b]; };
	F_greater.sort(compare1);
	B_less.sort(compare2);
	T *u_pspace = pspaces[0];

	if (threshold != v) {
		if (n2v[threshold].base() == u_pspace->end() || std::next(n2v[threshold].base()) == u_pspace->end()) {
			val2 = plus_infinity;
		}
		else {
			val2 = ahrsz_ext_priority_value<T>(std::next(n2v[threshold].base()), *u_pspace);  // val_t = prev(S)
		}
		val1 = ahrsz_ext_priority_value<T>(n2v[threshold].base(), *u_pspace);
		for (rit = F_greater.rbegin(); rit != F_greater.rend(); rit++) {
			ahrsz_ext_priority_value<T> ep = compute_priority(val1, val2, *u_pspace);
			n2v[*rit] = ahrsz_priority_value<T>(ep.base(), *u_pspace);
			val1 = ahrsz_ext_priority_value<T>(n2v[*rit].base(), *u_pspace);
		}
		for (it = B_less.begin(); it != B_less.end(); it++) {
			ahrsz_ext_priority_value<T> ep = compute_priority(val1, val2, *u_pspace);
			n2v[*it] = ahrsz_priority_value<T>(ep.base(), *u_pspace);
			val1 = ahrsz_ext_priority_value<T>(n2v[*it].base(), *u_pspace);
		}
	}
	else {
		if (n2v[threshold].base() == u_pspace->begin()) {
			val1 = minus_infinity;
		}
		else {
			val1 = ahrsz_ext_priority_value<T>(u_pspace->previous(n2v[threshold].base()), *u_pspace);
		}
		val2 = ahrsz_ext_priority_value<T>(n2v[threshold].base(), *u_pspace);
		for (rit = F_greater.rbegin(); rit != F_greater.rend(); rit++) {
			ahrsz_ext_priority_value<T> ep = compute_priority(val1, val2, *u_pspace);
			n2v[*rit] = ahrsz_priority_value<T>(ep.base(), *u_pspace);
			val1 = ahrsz_ext_priority_value<T>(n2v[*rit].base(), *u_pspace);
		}
		for (it = B_less.begin(); it != B_less.end(); it++) {
			ahrsz_ext_priority_value<T> ep = compute_priority(val1, val2, *u_pspace);
			n2v[*it] = ahrsz_priority_value<T>(ep.base(), *u_pspace);
			val1 = ahrsz_ext_priority_value<T>(n2v[*it].base(), *u_pspace);
		}
	}
}

template<typename T>
void dynamic_tc_hkmst_scc<T>::maintain_Gsi(std::list<vd_t>& forw_nodes, std::list<vd_t>& back_nodes, vd_t & v, vd_t & w)
{
	std::list<vd_t> F_greater, B_less;
	typename std::list<vd_t>::iterator it;
	typename std::list<vd_t>::reverse_iterator rit;
	ahrsz_ext_priority_value<T> val1, val2;
	vd_t threshold = v;

	for (it = forw_nodes.begin(); it != forw_nodes.end(); it++) {
		if (n2v[threshold] < n2v[*it] && has_out_after_search_Gsi(*it)) {
			threshold = *it;
		}
	}
	// F> 
	for (it = forw_nodes.begin(); it != forw_nodes.end(); it++) {
		if (n2v[threshold] < n2v[*it])
			F_greater.push_back(*it);	// unmark states. 
		_inF[*it] = false;
		is_out[*it] = false;
	}
	// B< 
	for (it = back_nodes.begin(); it != back_nodes.end(); it++) {
		if (n2v[*it] < n2v[threshold])
			B_less.push_back(*it);
		_inB[*it] = false;
		is_in[*it] = false;
	}
	auto compare1 = [this](const vd_t &a, const vd_t &b) { return n2v[b] < n2v[a]; };
	auto compare2 = [this](const vd_t &a, const vd_t &b) { return n2v[a] < n2v[b]; };
	F_greater.sort(compare1);
	B_less.sort(compare2);

	int idx = 0;
	if (B_less.size() > 0) {
		idx = n2y[B_less.front()] - min_year + 1;
	}
	else {
		idx = n2y[F_greater.front()] - min_year + 1;
	}
	T* y_pspace = (pspaces[idx]);

	if (threshold != v) {
		if (n2v[threshold].base() == y_pspace->end() || std::next(n2v[threshold].base()) == y_pspace->end()) {
			val2 = plus_infinity;
		}
		else {
			val2 = ahrsz_ext_priority_value<T>(std::next(n2v[threshold].base()), *y_pspace);  // val_t = prev(S)
		}
		val1 = ahrsz_ext_priority_value<T>(n2v[threshold].base(), *y_pspace);
		for (rit = F_greater.rbegin(); rit != F_greater.rend(); rit++) {
			ahrsz_ext_priority_value<T> ep = compute_priority(val1, val2, *y_pspace);
			n2v[*rit] = ahrsz_priority_value<T>(ep.base(), *y_pspace);
			val1 = ahrsz_ext_priority_value<T>(n2v[*rit].base(), *y_pspace);
		}
		for (it = B_less.begin(); it != B_less.end(); it++) {
			ahrsz_ext_priority_value<T> ep = compute_priority(val1, val2, *y_pspace);
			n2v[*it] = ahrsz_priority_value<T>(ep.base(), *y_pspace);
			val1 = ahrsz_ext_priority_value<T>(n2v[*it].base(), *y_pspace);
		}
	}
	else {
		if (n2v[threshold].base() == y_pspace->begin()) {
			val1 = minus_infinity;
		}
		else {
			val1 = ahrsz_ext_priority_value<T>(y_pspace->previous(n2v[threshold].base()), *y_pspace);
		}
		val2 = ahrsz_ext_priority_value<T>(n2v[threshold].base(), *y_pspace);
		for (rit = F_greater.rbegin(); rit != F_greater.rend(); rit++) {
			ahrsz_ext_priority_value<T> ep = compute_priority(val1, val2, *y_pspace);
			n2v[*rit] = ahrsz_priority_value<T>(ep.base(), *y_pspace);
			val1 = ahrsz_ext_priority_value<T>(n2v[*rit].base(), *y_pspace);
		}
		for (it = B_less.begin(); it != B_less.end(); it++) {
			ahrsz_ext_priority_value<T> ep = compute_priority(val1, val2, *y_pspace);
			n2v[*it] = ahrsz_priority_value<T>(ep.base(), *y_pspace);
			val1 = ahrsz_ext_priority_value<T>(n2v[*it].base(), *y_pspace);
		}
	}

}

template<typename T>
void dynamic_tc_hkmst_scc<T>::maintain_SCC_Gm(std::vector<vd_t>& K)
{
	
	for (std::vector<vd_t>::iterator i(K.begin()); i != K.end(); ++i) {
		_flooring[*i] = minus_infinity;
		_inB[*i] = false;		// unmark states. 
		_inF[*i] = false;
		_inK[*i] = true;
	}

	// first pass - compute flooring Information
	std::vector<vd_t> rto; // reverse topological order
	for (std::vector<vd_t>::iterator i(K.begin()); i != K.end(); ++i) {
		if (!_visited[*i]) { compute_flooring_Gm(*i, rto); }
	}
	// second pass - perform the reassignment
	for (std::vector<vd_t>::reverse_iterator i(rto.rbegin()); i != rto.rend(); ++i) {
		ahrsz_ext_priority_value<T> ep(compute_priority(_flooring[*i], compute_ceiling_Gm(*i), *(pspaces[0])));
		n2v[*i] = ahrsz_priority_value<T>(ep.base(), *pspaces[0]);
		_visited[*i] = false; // unmark state info. 
	}
	// reset visited information	unmark state info. 
	for (std::vector<vd_t>::iterator i(K.begin()); i != K.end(); ++i) {
		_visited[*i] = false;
		_inK[*i] = false;
	}
}

template<typename T>
void dynamic_tc_hkmst_scc<T>::maintain_SCC_Gsi(std::vector<vd_t>& K)
{
	for (std::vector<vd_t>::iterator i(K.begin()); i != K.end(); ++i) {
		_flooring[*i] = minus_infinity;
		_inB[*i] = false;		// unmark states. 
		_inF[*i] = false;
		_inK[*i] = true;
	}
	int idx = n2y[K.front()] - min_year + 1;
	T* y_pspace = (pspaces[idx]);

	// first pass - compute flooring Information
	std::vector<vd_t> rto; // reverse topological order
	for (std::vector<vd_t>::iterator i(K.begin()); i != K.end(); ++i) {
		if (!_visited[*i]) { compute_flooring_Gsi(*i, rto ); }
	}
	// second pass - perform the reassignment
	for (std::vector<vd_t>::reverse_iterator i(rto.rbegin()); i != rto.rend(); ++i) {
		ahrsz_ext_priority_value<T> ep(compute_priority(_flooring[*i], compute_ceiling_Gsi(*i), *y_pspace));
		n2v[*i] = ahrsz_priority_value<T>(ep.base(), *y_pspace);
		_visited[*i] = false; // unmark state info. 
	}
	// reset visited information	unmark state info. 
	for (std::vector<vd_t>::iterator i(K.begin()); i != K.end(); ++i) {
		_visited[*i] = false;
		_inK[*i] = false;
	}

}

template<typename T>
ahrsz_ext_priority_value<T> dynamic_tc_hkmst_scc<T>::compute_priority(ahrsz_ext_priority_value<T> floor, ahrsz_ext_priority_value<T> ceiling, T & cur_pspace){
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

template<typename T>
void dynamic_tc_hkmst_scc<T>::compute_flooring_Gm(vd_t n, std::vector<vd_t>& rto)
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

template<typename T>
void dynamic_tc_hkmst_scc<T>::compute_flooring_Gsi(vd_t n, std::vector<vd_t>& rto)
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

template<typename T>
ahrsz_ext_priority_value<T> dynamic_tc_hkmst_scc<T>::compute_ceiling_Gm(vd_t v)
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

template<typename T>
ahrsz_ext_priority_value<T> dynamic_tc_hkmst_scc<T>::compute_ceiling_Gsi(vd_t v)
{
	ahrsz_ext_priority_value<T> ceiling(plus_infinity);
	std::vector<unsigned int>::iterator it;
	in_iterator j, jend;
	for (it = ds.scc_nodes[v].begin(); it != ds.scc_nodes[v].end(); it++) {
		for (tie(j, jend) = in_edges(*it, cg); j != jend; ++j) {
			unsigned int w = ds.find(source(*j, cg));
			if (!is_in_unorder[w] && is_in_same_year[w] == is_in_same_year[v]) {
				ceiling = std::min(ceiling, ahrsz_ext_priority_value<T>(n2v[w]));
			}
		}
	}
	return ceiling;
}

// to identify whether the adjacent vertices of 'node' are in Gm. called by soft_threshold_search_Gm
template<typename T>
bool dynamic_tc_hkmst_scc<T>::has_first_out_Gm(vd_t & node)
{
	std::vector<unsigned int>::iterator it_begin, it_end;
	out_iterator oei_begin, oei_end;
	it_begin = ds.scc_nodes[node].begin();
	it_end = ds.scc_nodes[node].end();
	if (it_begin != it_end) {
		oei_begin = boost::out_edges(*it_begin, cg).first;
		oei_end = boost::out_edges(*it_begin, cg).second;
	}
	bool flag = false;
	while (it_begin != it_end) {
		while (oei_begin != oei_end) {
			vd_t tmp = ds.find(boost::target(*oei_begin, cg));		
			if ( is_in_unorder[tmp] ) {
				flag = true;
				return flag;
			}
			oei_begin++; 
		}
		it_begin++;
		if (it_begin != it_end) {
			oei_begin = boost::out_edges(*it_begin, cg).first;
			oei_end = boost::out_edges(*it_begin, cg).second;
		}
	}
	return flag;
}

// to identify whether the adjacent vertices of 'node' are in Gsi. called by soft_threshold_search_Gsi
template<typename T>
bool dynamic_tc_hkmst_scc<T>::has_first_out_Gsi(vd_t & node)
{
	std::vector<unsigned int>::iterator it_begin, it_end;
	out_iterator oei_begin, oei_end;
	it_begin = ds.scc_nodes[node].begin();
	it_end = ds.scc_nodes[node].end();
	if (it_begin != it_end) {
		oei_begin = boost::out_edges(*it_begin, cg).first;
		oei_end = boost::out_edges(*it_begin, cg).second;
	}
	bool flag = false;
	while (it_begin != it_end) {
		while (oei_begin != oei_end) {
			vd_t tmp = ds.find(boost::target(*oei_begin, cg));	
			if (is_in_unorder[tmp] == false && is_in_same_year[tmp] == is_in_same_year[node] ) {
				flag = true;
				return flag; 
			}
			oei_begin ++; 
		}
		it_begin++;
		if (it_begin != it_end) {
			oei_begin = boost::out_edges(*it_begin, cg).first;
			oei_end = boost::out_edges(*it_begin, cg).second;
		}
	}
	return flag;
}

// to identify whether the adjacent vertices of 'node' are in Gm. called by soft_threshold_search_Gm
template<typename T>
bool dynamic_tc_hkmst_scc<T>::has_first_in_Gm(vd_t & node)
{
	std::vector<unsigned int>::iterator it_begin, it_end;
	in_iterator iei_begin, iei_end;
	it_begin = ds.scc_nodes[node].begin();
	it_end = ds.scc_nodes[node].end();
	if (it_begin != it_end) {
		iei_begin = boost::in_edges(*it_begin, cg).first;
		iei_end = boost::in_edges(*it_begin, cg).second;
	}
	bool flag = false;
	while (it_begin != it_end) {
		while (iei_begin != iei_end) {
			vd_t tmp = ds.find(boost::source(*iei_begin, cg));
			if (is_in_unorder[tmp]) {
				flag = true;
				return flag; 
			}
			iei_begin++; 
		}
		it_begin++;
		if (it_begin != it_end) {
			iei_begin = boost::in_edges(*it_begin, cg).first;
			iei_end = boost::in_edges(*it_begin, cg).second;
		}
	}
	return flag;
}

// to identify whether the adjacent vertices of 'node' are in Gsi. called by soft_threshold_search_Gsi
template<typename T>
bool dynamic_tc_hkmst_scc<T>::has_first_in_Gsi(vd_t & node)
{
	std::vector<unsigned int>::iterator it_begin, it_end;
	in_iterator iei_begin, iei_end;
	it_begin = ds.scc_nodes[node].begin();
	it_end = ds.scc_nodes[node].end();
	if (it_begin != it_end) {
		iei_begin = boost::in_edges(*it_begin, cg).first;
		iei_end = boost::in_edges(*it_begin, cg).second;
	}
	bool flag = false;
	while (it_begin != it_end) {
		while (iei_begin != iei_end) {
			vd_t tmp = ds.find(boost::source(*iei_begin, cg));		
			if (is_in_unorder[tmp] == false && is_in_same_year[tmp] == is_in_same_year[node]) {
				flag = true;
				return flag; 
			}
			iei_begin++; 
		}
		it_begin++;
		if (it_begin != it_end) {
			iei_begin = boost::in_edges(*it_begin, cg).first;
			iei_end = boost::in_edges(*it_begin, cg).second;
		}
	}
	return flag;
}

// whether the adjacent vertices of the out_scc_iterator still has a node in Gm. called by search_step_Gm
template<typename T>
bool dynamic_tc_hkmst_scc<T>::out_tail_Gm(out_scc_iterator & out)
{
	bool flag = false;
	out_iterator e1, e2;
	std::vector<unsigned int>::iterator i1, i2;
	e1 = out.oei_begin;
	e2 = out.oei_end;
	i1 = out.it_begin;
	i2 = out.it_end;
	while (i1 != i2) {
		while (e1 != e2) {
			vd_t tmp = ds.find(boost::target(*e1, cg));		
			if (is_in_unorder[tmp] ) {
				flag = true;
				return flag; 
			}
			e1++; 
		}
		i1 = std::next(i1, 1);
		if (i1 != i2) {
			e1 = boost::out_edges(*i1, cg).first;
			e2 = boost::out_edges(*i1, cg).second;
		}
	}
	return flag;
}

// whether the adjacent vertices of the out_scc_iterator still has a node in Gm. called by search_step_Gsi
template<typename T>
bool dynamic_tc_hkmst_scc<T>::out_tail_Gsi(out_scc_iterator & out, vd_t & node)
{
	bool flag = false;
	out_iterator e1, e2;
	std::vector<unsigned int>::iterator i1, i2;
	e1 = out.oei_begin;
	e2 = out.oei_end;
	i1 = out.it_begin;
	i2 = out.it_end;
	while (i1 != i2) {
		while (e1 != e2) {
			vd_t tmp = ds.find(boost::target(*e1, cg));				// note: do NOT need to call ds.find
			if (is_in_unorder[tmp]== false && is_in_same_year[tmp] == is_in_same_year[node]) {
				flag = true;
				return flag; 
			}
			e1++; 
		}
		i1 = std::next(i1, 1);
		if (i1 != i2) {
			e1 = boost::out_edges(*i1, cg).first;
			e2 = boost::out_edges(*i1, cg).second;
		}
	}
	return flag;
}

// whether the adjacent vertices of the in_scc_iterator still has a node in Gm. called by search_step_Gsi
template<typename T>
bool dynamic_tc_hkmst_scc<T>::in_tail_Gm(in_scc_iterator & in)
{
	bool flag = false;
	in_iterator e1, e2;
	std::vector<unsigned int>::iterator i1, i2;
	e1 = in.iei_begin;
	e2 = in.iei_end;
	i1 = in.it_begin;
	i2 = in.it_end;
	while (i1 != i2) {
		while (e1 != e2) {
			vd_t tmp = ds.find(boost::source(*e1, cg));
			if (is_in_unorder[tmp]) {
				flag = true;
				return flag; 
			}
			e1++; 
		}
		i1 = std::next(i1, 1);
		if (i1 != i2) {
			e1 = boost::in_edges(*i1, cg).first;
			e2 = boost::in_edges(*i1, cg).second;
		}
	}
	return flag;
}

template<typename T>
bool dynamic_tc_hkmst_scc<T>::in_tail_Gsi(in_scc_iterator & in, vd_t &node)
{
	bool flag = false;
	in_iterator e1, e2;
	std::vector<unsigned int>::iterator i1, i2;
	e1 = in.iei_begin;
	e2 = in.iei_end;
	i1 = in.it_begin;
	i2 = in.it_end;
	while (i1 != i2) {
		while (e1 != e2) {
			vd_t tmp = ds.find(boost::source(*e1, cg));
			if (is_in_unorder[tmp] == false && is_in_same_year[tmp] == is_in_same_year[node]) {
				flag = true;
				return flag; 
			}
			e1++;
		}
		i1 = std::next(i1, 1);
		if (i1 != i2) {
			e1 = boost::in_edges(*i1, cg).first;
			e2 = boost::in_edges(*i1, cg).second;
		}
	}
	return flag;
}

//called by maintain_Gm
template<typename T>
bool dynamic_tc_hkmst_scc<T>::has_out_after_search_Gm(vd_t node)
{
	out_scc_iterator node_it;
	if (is_out[node] == false) {
		return has_first_out_Gm(node);
	}
	else {
		// each next-out node updates by first-out function. Thus, we use std::next. 
		node_it = out_its[node];
		return has_next_out_Gm(node_it);
	}
}

template<typename T>
bool dynamic_tc_hkmst_scc<T>::has_out_after_search_Gsi(vd_t node)
{
	out_scc_iterator node_it;
	if (is_out[node] == false) {
		return has_first_out_Gsi(node);
	}
	else {
		// each next-out node updates by first-out function. Thus, we use std::next. 
		node_it = out_its[node];
		return has_next_out_Gm(node_it);
	}
}

// get the out_scc_iterator of the first_out node in Gm
template<typename T>
auto dynamic_tc_hkmst_scc<T>::first_out_Gm(vd_t & node)
{
	if (is_out[node] == false) {
		is_out[node] = true;
		std::vector<unsigned int>::iterator it_tmp;
		std::pair<out_iterator, out_iterator> scc_item_it;
		it_tmp = ds.scc_nodes[node].begin();
		out_its[node].it_begin = it_tmp;
		out_its[node].it_end = ds.scc_nodes[node].end();
		if (out_its[node].it_begin != out_its[node].it_end) {
			out_its[node].oei_begin = boost::out_edges(*it_tmp, cg).first;
			out_its[node].oei_end = boost::out_edges(*it_tmp, cg).second;
		}
		while (out_its[node].it_begin != out_its[node].it_end) {
			scc_item_it = boost::out_edges(*(out_its[node].it_begin), cg);
			while (scc_item_it.first != scc_item_it.second) {
				vd_t tmp = ds.find(boost::target(*(scc_item_it.first), cg));
				if (is_in_unorder[tmp]) {
					out_its[node].oei_begin = scc_item_it.first;
					out_its[node].oei_end = scc_item_it.second;
					return out_its[node];
				}
				scc_item_it.first++;
			}
			out_its[node].it_begin++;
		}
		return out_its[node];
	}
	else {
		out_its[node].oei_begin++;
		while (out_its[node].oei_begin != out_its[node].oei_end) {
			vd_t tmp = ds.find(boost::target(*(out_its[node].oei_begin), cg));
			if (is_in_unorder[tmp]) {
				return out_its[node];
			}
			out_its[node].oei_begin++;
		}
		if (out_its[node].oei_begin == out_its[node].oei_end) {
			while (out_its[node].it_begin != out_its[node].it_end) {
				out_its[node].it_begin++;
				if (out_its[node].it_begin != out_its[node].it_end) {
					out_its[node].oei_begin = boost::out_edges(*(out_its[node].it_begin), cg).first;
					out_its[node].oei_end = boost::out_edges(*(out_its[node].it_begin), cg).second;
					while (out_its[node].oei_begin != out_its[node].oei_end) {
						vd_t tmp = ds.find(boost::target(*(out_its[node].oei_begin), cg));		
						if (is_in_unorder[tmp]) {
							return out_its[node];
						}
						out_its[node].oei_begin++;
					}
				}
			}
		}
		return out_its[node];
	}

}

template<typename T>
auto dynamic_tc_hkmst_scc<T>::first_out_Gsi(vd_t & node)
{
	if (is_out[node] == false) {
		is_out[node] = true;
		std::vector<unsigned int>::iterator it_tmp;
		std::pair<out_iterator, out_iterator> scc_item_it;
		it_tmp = ds.scc_nodes[node].begin();
		out_its[node].it_begin = it_tmp;
		out_its[node].it_end = ds.scc_nodes[node].end();
		if (out_its[node].it_begin != out_its[node].it_end) {
			out_its[node].oei_begin = boost::out_edges(*it_tmp, cg).first;
			out_its[node].oei_end = boost::out_edges(*it_tmp, cg).second;
		}
		while (out_its[node].it_begin != out_its[node].it_end) {
			scc_item_it = boost::out_edges(*(out_its[node].it_begin), cg);
			while (scc_item_it.first != scc_item_it.second) {
				vd_t tmp = ds.find(boost::target(*(scc_item_it.first), cg));		
				if (is_in_unorder[tmp] == false && is_in_same_year[tmp] == is_in_same_year[node]) {
					out_its[node].oei_begin = scc_item_it.first;
					out_its[node].oei_end = scc_item_it.second;
					return out_its[node];
				}
				scc_item_it.first++;
			}
			out_its[node].it_begin++;
		}
		return out_its[node];
	}
	else {
		out_its[node].oei_begin++;
		while (out_its[node].oei_begin != out_its[node].oei_end) {
			vd_t tmp = ds.find(boost::target(*(out_its[node].oei_begin), cg));
			if (is_in_unorder[tmp] == false && is_in_same_year[tmp] == is_in_same_year[node]) {
				return out_its[node];
			}
			out_its[node].oei_begin++;
		}
		if (out_its[node].oei_begin == out_its[node].oei_end) {
			while (out_its[node].it_begin != out_its[node].it_end) {
				out_its[node].it_begin++;
				if (out_its[node].it_begin != out_its[node].it_end) {
					out_its[node].oei_begin = boost::out_edges(*(out_its[node].it_begin), cg).first;
					out_its[node].oei_end = boost::out_edges(*(out_its[node].it_begin), cg).second;
					while (out_its[node].oei_begin != out_its[node].oei_end) {
						vd_t tmp = ds.find( boost::target(*(out_its[node].oei_begin), cg));
						if (is_in_unorder[tmp] == false && is_in_same_year[tmp] == is_in_same_year[node]) {
							return out_its[node];
						}
						out_its[node].oei_begin++;
					}
				}
			}
		}
		return out_its[node];
	}
}

template<typename T>
auto dynamic_tc_hkmst_scc<T>::first_in_Gm(vd_t & node)
{
	if (is_in[node] == false) {
		is_in[node] = true;
		std::vector<unsigned int>::iterator it_tmp;
		std::pair<in_iterator, in_iterator> scc_item_it;
		it_tmp = ds.scc_nodes[node].begin();
		in_its[node].it_begin = it_tmp;
		in_its[node].it_end = ds.scc_nodes[node].end();
		if (in_its[node].it_begin != in_its[node].it_end) {
			in_its[node].iei_begin = boost::in_edges(*it_tmp, cg).first;
			in_its[node].iei_end = boost::in_edges(*it_tmp, cg).second;
		}
		while (in_its[node].it_begin != in_its[node].it_end) {
			scc_item_it = boost::in_edges(*(in_its[node].it_begin), cg);
			while (scc_item_it.first != scc_item_it.second) {
				vd_t tmp = ds.find(boost::source(*(scc_item_it.first), cg));
				if (is_in_unorder[tmp]) {
					in_its[node].iei_begin = scc_item_it.first;
					in_its[node].iei_end = scc_item_it.second;
					return in_its[node]; 
				}
				scc_item_it.first ++; 
			}
			in_its[node].it_begin++;
		}
		return in_its[node];
	}
	else {
		in_its[node].iei_begin++;
		// iterate the out_iterator of the current node to identify whether there is a node in Gm . 
		while (in_its[node].iei_begin != in_its[node].iei_end) {
			vd_t tmp = ds.find(boost::source(*(in_its[node].iei_begin), cg));
			if (is_in_unorder[tmp]) {
				return in_its[node];
			}
			in_its[node].iei_begin ++;
		}
		// iterate the out_iterator of the other nodes to identify whether there is a node in Gm . 
		if (in_its[node].iei_begin == in_its[node].iei_end) {
			while (in_its[node].it_begin != in_its[node].it_end) {
				in_its[node].it_begin++;
				if (in_its[node].it_begin != in_its[node].it_end) {
					in_its[node].iei_begin = boost::in_edges(*(in_its[node].it_begin), cg).first;
					in_its[node].iei_end = boost::in_edges(*(in_its[node].it_begin), cg).second;
					while (in_its[node].iei_begin != in_its[node].iei_end) {
						vd_t tmp = ds.find(boost::source(*(in_its[node].iei_begin), cg));
						if (is_in_unorder[tmp]) {
							return in_its[node];
						}
						in_its[node].iei_begin++; 
					}
				}
			}
		}
		return in_its[node];
	}

}

template<typename T>
auto dynamic_tc_hkmst_scc<T>::first_in_Gsi(vd_t & node)
{
	if (is_in[node] == false) {
		is_in[node] = true;
		std::vector<unsigned int>::iterator it_tmp;
		std::pair<in_iterator, in_iterator> scc_item_it;
		it_tmp = ds.scc_nodes[node].begin();
		in_its[node].it_begin = it_tmp;
		in_its[node].it_end = ds.scc_nodes[node].end();
		if (in_its[node].it_begin != in_its[node].it_end) {
			in_its[node].iei_begin = boost::in_edges(*it_tmp, cg).first;
			in_its[node].iei_end = boost::in_edges(*it_tmp, cg).second;
		}
		while (in_its[node].it_begin != in_its[node].it_end) {
			scc_item_it = boost::in_edges(*(in_its[node].it_begin), cg);
			while (scc_item_it.first != scc_item_it.second) {
				vd_t tmp = ds.find(boost::source(*(scc_item_it.first), cg));
				if (is_in_unorder[tmp] == false && is_in_same_year[tmp] == is_in_same_year[node]) {
					in_its[node].iei_begin = scc_item_it.first;
					in_its[node].iei_end = scc_item_it.second;
					return in_its[node];
				}
				scc_item_it.first++; 
			}
			in_its[node].it_begin++;
		}
		return in_its[node];
	}
	else {
		in_its[node].iei_begin++;
		while (in_its[node].iei_begin != in_its[node].iei_end) {
			vd_t tmp = ds.find(boost::source(*(in_its[node].iei_begin), cg));
			if (is_in_unorder[tmp] == false && is_in_same_year[tmp] == is_in_same_year[node]) {
				return in_its[node];
			}
			in_its[node].iei_begin++;
		}
		if (in_its[node].iei_begin == in_its[node].iei_end) {
			while (in_its[node].it_begin != in_its[node].it_end) {
				in_its[node].it_begin++;
				if (in_its[node].it_begin != in_its[node].it_end) {
					in_its[node].iei_begin = boost::in_edges(*(in_its[node].it_begin), cg).first;
					in_its[node].iei_end = boost::in_edges(*(in_its[node].it_begin), cg).second;
					while (in_its[node].iei_begin != in_its[node].iei_end) {
						vd_t tmp = ds.find(boost::source(*(in_its[node].iei_begin), cg));
						if (is_in_unorder[tmp]== false && is_in_same_year[tmp] == is_in_same_year[node]) {
							return in_its[node];
						}
						in_its[node].iei_begin++; 
					}
				}
			}
		}
		return in_its[node];
	}
}

template<typename T>
bool dynamic_tc_hkmst_scc<T>::has_next_out_Gm(out_scc_iterator & out)
{
	bool flag = false;
	out_iterator e1, e2;
	std::vector<unsigned int>::iterator i1, i2;
	e1 = out.oei_begin;
	e2 = out.oei_end;
	i1 = out.it_begin;
	i2 = out.it_end;
	if (e1 != e2 && std::next(e1, 1) != e2) {    // current node has next out ? 
		e1++;								// iterate next . 
		while (e1 != e2) {
			vd_t tmp = ds.find(boost::target(*e1, cg));
			if (is_in_unorder[tmp]) {
				flag = true;
				return flag;
			}
			e1 = std::next(e1, 1);
		}
	}
	if (i1 != i2 ) {
		i1 = std::next(i1, 1);			// iterate next . 
		if (i1 != i2) {
			e1 = boost::out_edges(*i1, cg).first;
			e2 = boost::out_edges(*i1, cg).second;
		}
		while (i1 != i2) {
			while (e1 != e2) {
				vd_t tmp = ds.find(boost::target(*(e1), cg));
				if (is_in_unorder[tmp]) {
					flag = true;
					return flag;
				}
				e1 = std::next(e1, 1);
			}
			i1 = std::next(i1, 1);
			if (i1 != i2) {
				e1 = boost::out_edges(*i1, cg).first;
				e2 = boost::out_edges(*i1, cg).second;
			}
		}
	}
	return flag;
}

template<typename T>
bool dynamic_tc_hkmst_scc<T>::has_next_out_Gsi(out_scc_iterator & out, vd_t node)
{
	bool flag = false;
	out_iterator e1, e2;
	std::vector<unsigned int>::iterator i1, i2;
	e1 = out.oei_begin;
	e2 = out.oei_end;
	i1 = out.it_begin;
	i2 = out.it_end;
	if (e1 != e2 && std::next(e1, 1) != e2) {    // current node has next out ? 
		e1++;
		while (e1 != e2) {
			vd_t tmp = ds.find( boost::target(*e1, cg));
			if (is_in_unorder[tmp] == false && is_in_same_year[tmp] == is_in_same_year[node]) {
				flag = true;
				return flag;
			}
			e1 = std::next(e1, 1);
		}
	}
	if (i1 != i2 ) {
		i1 = std::next(i1, 1);
		if (i1 != i2) {
			e1 = boost::out_edges(*i1, cg).first;
			e2 = boost::out_edges(*i1, cg).second;
		}
		while (i1 != i2) {
			while (e1 != e2) {
				vd_t tmp = ds.find(boost::target(*(e1), cg));
				if (is_in_unorder[tmp] == false && is_in_same_year[tmp] == is_in_same_year[node]) {
					flag = true;
					return flag;
				}
				e1 = std::next(e1, 1);
			}
			i1 = std::next(i1, 1);
			if (i1 != i2) {
				e1 = boost::out_edges(*i1, cg).first;
				e2 = boost::out_edges(*i1, cg).second;
			}
		}
	}
	return flag;
}

template<typename T>
bool dynamic_tc_hkmst_scc<T>::has_next_in_Gm(in_scc_iterator & in)
{
	bool flag = false;
	in_iterator e1, e2;
	std::vector<unsigned int>::iterator i1, i2;
	e1 = in.iei_begin;
	e2 = in.iei_end;
	i1 = in.it_begin;
	i2 = in.it_end;
	if (e1 != e2 && std::next(e1, 1) != e2) {			// current node has next in ? 
		e1++;
		while (e1 != e2) {
			vd_t tmp =ds.find( boost::source(*e1, cg));
			if (is_in_unorder[tmp]) {
				flag = true; 
				return flag; 
			}
			e1 = std::next(e1, 1);
		}
	}
	if (i1 != i2 ) {
		i1 = std::next(i1, 1);
		if (i1 != i2) {
			e1 = boost::in_edges(*i1, cg).first;
			e2 = boost::in_edges(*i1, cg).second;
		}
		while (i1 != i2) {
			while (e1 != e2) {
				vd_t tmp = ds.find(boost::source(*e1, cg));
				if (is_in_unorder[tmp]) {
					flag = true;
					return flag;
				}
				e1 = std::next(e1, 1);
			}
			i1 = std::next(i1, 1);
			if (i1 != i2) {
				e1 = boost::in_edges(*i1, cg).first;
				e2 = boost::in_edges(*i1, cg).second;
			}
		}
	}
	return flag;
}

template<typename T>
bool dynamic_tc_hkmst_scc<T>::has_next_in_Gsi(in_scc_iterator & in, vd_t node)
{
	bool flag = false;
	in_iterator e1, e2;
	std::vector<unsigned int>::iterator i1, i2;
	e1 = in.iei_begin;
	e2 = in.iei_end;
	i1 = in.it_begin;
	i2 = in.it_end;
	if (e1 != e2 && std::next(e1, 1) != e2) {			// current node has next in ? 
		e1++;
		while (e1 != e2) {
			vd_t tmp = ds.find( boost::source(*e1, cg));
			if (is_in_unorder[tmp] == false &&  is_in_same_year[tmp] == is_in_same_year[node]) {
				flag = true;
				return flag;
			}
			e1 = std::next(e1, 1);
		}
	}
	if (i1 != i2) {
		i1 = std::next(i1, 1);
		if (i1 != i2) {
			e1 = boost::in_edges(*i1, cg).first;
			e2 = boost::in_edges(*i1, cg).second;
		}
		while (i1 != i2) {
			while (e1 != e2) {
				vd_t tmp = ds.find(boost::source(*e1, cg));
				if (is_in_unorder[tmp] == false && is_in_same_year[tmp] == is_in_same_year[node]) {
					flag = true;
					return flag;
				}
				e1 = std::next(e1, 1);
			}
			i1 = std::next(i1, 1);
			if (i1 != i2) {
				e1 = boost::in_edges(*i1, cg).first;
				e2 = boost::in_edges(*i1, cg).second;
			}
		}
	}
	return flag;
}


