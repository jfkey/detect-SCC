#include "dynamic_IGC_scc.h"

template<class T>
dynamic_igc_scc<T>::dynamic_igc_scc()
{

}

template<class T>
dynamic_igc_scc<T>::dynamic_igc_scc(CG & _cg, disjoint_set & _ds, std::vector<T>& _n2v, std::vector<pair<unsigned int, unsigned int>> & _inc_edges, unsigned int _np):cg(_cg), ds(_ds), n2v(_n2v), inc_edges(_inc_edges), np(_np)
{
	in_component = std::vector<bool>(np, false);
	visited = std::vector<bool>(np, false);
	self_name = "Dynamic-IGC-";
	myres.set_dynamic_name(self_name);
}

template<class T>
void dynamic_igc_scc<T>::add_edges()
{
	my_timer mt; 
	std::vector<pair<unsigned int, unsigned int>>::iterator it;
	for (it = inc_edges.begin(); it != inc_edges.end(); it++) {
		//std::cout << it->first << ", " << it->second << std::endl;
		unsigned int s, t;
		s = ds.find(it->first);
		t = ds.find(it->second);
		if (s != t) {													// s == t, e.g., the edge in the same component. do not add edge. 
			my_timer mt = my_timer();
			boost::add_edge(it->first, it->second, cg);
			do_igc(s, t);
		}
	}
	myres.dn_eval_time = myres.search_time + myres.reassignment_time;
	myres.dn_eval_time_AAN = mt.elapsed();
}

template<class T>
void dynamic_igc_scc<T>::add_edges(int start_index, int end_index)
{
	myres.inc_edges += (end_index - start_index);
	my_timer mt;
	for (int i = start_index; i < end_index; i++) {
		unsigned int s, t;
		s = ds.find(inc_edges[i].first);
		t = ds.find(inc_edges[i].second);
		if (s != t) {													// s == t, e.g., the edge in the same component. do not add edge. 
			my_timer mt = my_timer();
			boost::add_edge(inc_edges[i].first, inc_edges[i].second, cg);
			do_igc(s, t);
		}
	}
	myres.dn_eval_time = myres.search_time + myres.reassignment_time;
	record_cur_per_info(myres);
}



template<class T>
void dynamic_igc_scc<T>::do_igc(vd_t t, vd_t h)
{
	unsigned int tn2v(n2v[t]);
	unsigned int hn2v(n2v[h]);
	my_timer mt;

	// tn2v == hn2v, add the edge(t,h) in the scc. does not affect the topological order.
	// tn2v > hn2v, the edge(t,h) is actual topological order.
	if (tn2v < hn2v) {
		myres.set_invalid_edges(myres.get_invalid_edges() + 1);

		std::vector<unsigned int> reaching, reachable;
		bool flag = false;
		in_component[t] = true;
		
		mt = my_timer();
		igc_scc_fwd_dfs(h, tn2v, hn2v, reachable, flag);
		myres.set_search_time(myres.get_search_time() + mt.elapsed());

		if (flag == true) {		// there exists the cycle in the condensation graph. 
			// reset search status. 
			for_each(reachable.begin(), reachable.end(), [this](const unsigned int & a) {
				this->visited[a] = false; 
			});
			mt = my_timer();
			igc_scc_back_dfs(t, hn2v, reaching);
			myres.set_search_time(myres.get_search_time() + mt.elapsed());
			myres.set_aff_region(myres.get_aff_region() + reaching.size() + reachable.size());

			std::vector<unsigned int> cycle_nodes, reachable_tmp, reaching_tmp;
			std::vector<unsigned int>::iterator i, j;
			unsigned int s = 0;
			for (i = reachable.begin(); i != reachable.end(); ++i) {
				if (in_component[*i] == true) {
					cycle_nodes.push_back(*i);
				}
				else {
					reachable_tmp.push_back(*i);
				}
			}
			for (i = reaching.begin(); i != reaching.end(); i++) {
				if (in_component[*i] == false)
					reaching_tmp.push_back(*i);
			}

			mt = my_timer();
			assert(cycle_nodes.size() > 1 && "the number of vertices in the cycle must be greater than 1");
			s = cycle_nodes.front();
			for (i = cycle_nodes.begin(); i != cycle_nodes.end(); i++) {
				ds.join(s, *i);	
				visited[*i] = false;															// reset cycle status. 	
				in_component[*i] = false;
				in_component[s] = false;
			}
			s = ds.find(cycle_nodes.front());
			ds.collapse_scc(s, cycle_nodes);
			
			//condense_per_scc(s, cycle_nodes, g, g.ds);
			std::vector<unsigned int> cycle_nodes2;
			for (int i = 0; i < cycle_nodes.size(); i++) {
				cycle_nodes2.insert(cycle_nodes2.end(), ds.scc_nodes[cycle_nodes[i]].begin(), ds.scc_nodes[cycle_nodes[i]].end());
			}
			condense_per_scc(s, cycle_nodes2, cg, ds);

			myres.condensation_time += mt.elapsed();
			
			dynamic_igc_scc_com<std::vector<T>> pc(n2v);
			sort(reaching_tmp.begin(), reaching_tmp.end(), pc);
			sort(reachable_tmp.begin(), reachable_tmp.end(), pc);
			reaching_tmp.push_back(s);
			
			// debug
			//std::vector<unsigned int> aff_reg;
			//for (i = reaching_tmp.begin(); i != reaching_tmp.end(); i++) {
			//	aff_reg.push_back(*i);
			//}
			//for (i = reachable_tmp.begin(); i != reachable_tmp.end(); i++) {
			//	aff_reg.push_back(*i);
			//}
			//std::tuple<bool, unsigned int, unsigned int, std::vector<unsigned int> > invedge_affreg(true, t, h, aff_reg);
			//myres.invedge_affreg.push_back(invedge_affreg);
			
			igc_scc_reorder(reachable_tmp, reaching_tmp );

		}
		else {
			mt = my_timer();
			igc_scc_back_dfs(t, hn2v, reaching );
			myres.set_search_time(myres.get_search_time() + mt.elapsed());
			mt = my_timer();
			dynamic_igc_scc_com<std::vector<T>> pc(n2v);
			myres.set_aff_region(myres.get_aff_region() + reaching.size() + reachable.size());
			
			sort(reaching.begin(), reaching.end(), pc);
			sort(reachable.begin(), reachable.end(), pc);
			// debug
			//std::vector<unsigned int>::iterator i;
			//std::vector<unsigned int> aff_reg;
			//for (i = reaching.begin(); i != reaching.end(); i++) {
			//	aff_reg.push_back(*i);
			//}
			//for (i = reachable.begin(); i != reachable.end(); i++) {
			//	aff_reg.push_back(*i);
			//}
			//std::tuple<bool, unsigned int, unsigned int, std::vector<unsigned int> > invedge_affreg(true, t, h, aff_reg);
			//myres.invedge_affreg.push_back(invedge_affreg);
			
			igc_scc_reorder(reachable, reaching);
			myres.set_reassign_time(myres.get_reassign_time() + mt.elapsed());

		}
		in_component[t] = false; // unmark states. 
	}
}

template<class T>
void dynamic_igc_scc<T>::igc_scc_fwd_dfs(vd_t n, vd_t lb, vd_t fo, std::vector<unsigned int>& reachable, bool & flag)
{
	reachable.push_back(n);
	visited[n] = true;
	out_iterator i, iend;
	std::vector<unsigned int>::iterator it;
	for (it = ds.scc_nodes[n].begin(); it != ds.scc_nodes[n].end(); it++) {
		for (std::tie(i, iend) = boost::out_edges(*it, cg); i != iend; ++i) {
			//unsigned int w(target(*i, cg));
			unsigned int w = ds.find(target(*i, cg));
			unsigned int wn2v(n2v[w]);
			unsigned int fot(n2v[n]);
			if (wn2v == lb && !visited[w]) {
				in_component[w] = true;
				igc_scc_fwd_dfs(w, lb, fot, reachable, flag);
				flag = true;
			}
			else if (fo == wn2v && !visited[w]) {
				in_component[w] = true;
				igc_scc_fwd_dfs(w, lb, fot, reachable, flag);
			}
			else if (wn2v > lb && !visited[w]) {
				igc_scc_fwd_dfs(w, lb, fot, reachable, flag);
			}
			in_component[n] = in_component[n] || in_component[w];
		}
	}
}

template<class T>
void dynamic_igc_scc<T>::igc_scc_back_dfs(vd_t n, vd_t ub, std::vector<unsigned int>& reaching)
{
	reaching.push_back(n);
	visited[n] = true;
	in_iterator i, iend;
	std::vector<unsigned int>::iterator it;
	for (it = ds.scc_nodes[n].begin(); it != ds.scc_nodes[n].end(); it++) {
		for (std::tie(i, iend) = boost::in_edges(*it, cg); i != iend; ++i) {
			//unsigned int w(source(*i, g));
			unsigned int w = ds.find(source(*i, cg));
			unsigned int wn2v(n2v[w]);
			if (wn2v < ub && ! visited[w]) {
				igc_scc_back_dfs(w, ub, reaching );
			}
		}
	}
}

template<class T>
void dynamic_igc_scc<T>::igc_scc_reorder(std::vector<unsigned int>& reachable, std::vector<unsigned int>& reaching )
{
	std::vector<vd_t> tmp;
	tmp.reserve(reaching.size() + reachable.size());
	std::vector<unsigned int>::iterator iend(reaching.end());
	for (std::vector<unsigned int>::iterator i(reaching.begin()); i != iend; ++i) {
		tmp.push_back(*i);
		visited[*i] = false;
		*i = n2v[*i];			// dirty trick
	}
	iend = reachable.end();
	for (std::vector<unsigned int>::iterator i(reachable.begin()); i != iend; ++i) {
		tmp.push_back(*i);
		visited[*i] = false;
		*i = n2v[*i]; // dirty trick
	}
	
	// the canonical node of the scc is random selected by disjoint set. thus, the index of the reachable and reaching nodes 
	// are NOT the descending.
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
		//n2v[n] = w;
		//for (it_ds = ds.scc_nodes[n].begin(); it_ds != ds.scc_nodes[n].end(); it_ds++) {
			//n2v[*it_ds] = w;
		//}
		n2v[n] = w; 
	}


}

template<class T>
dynamic_igc_scc<T>::~dynamic_igc_scc()
{
}
