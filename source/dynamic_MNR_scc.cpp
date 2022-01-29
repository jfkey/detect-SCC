#include "dynamic_MNR_scc.h"

template<class T>
dynamic_mnr_scc<T>::dynamic_mnr_scc()
{
}

template<class T>
dynamic_mnr_scc<T>::dynamic_mnr_scc(CG & _cg, disjoint_set & _ds, std::vector<T>& _n2v, std::vector<pair<unsigned int, unsigned int>>& _inc_edges, unsigned int _np):cg(_cg), ds(_ds), n2v(_n2v), inc_edges(_inc_edges), np(_np)
{
	self_name = "Dynamic-MNR-";
	in_component = std::vector<bool>(np, false);
	visited = std::vector<bool>(np, false);
	v2n = std::vector<vd_t>(np +1, -1);
	
	for (unsigned int i = 0; i < np; i++) {
		v2n[n2v[i]] = ds.find(i);
	}
	
	myres.set_dynamic_name(self_name);
}

template<class T>
void dynamic_mnr_scc<T>::add_edges()
{
	my_timer mt; 
	std::vector<pair<unsigned int, unsigned int>>::iterator it;
	for (it = inc_edges.begin(); it != inc_edges.end(); it++) {
		unsigned int s, t;
		s = ds.find(it->first);
		t = ds.find(it->second);
		if (s != t) {													// s == t, e.g., the edge in the same component. do not add edge. 
			my_timer mt = my_timer();
			boost::add_edge(it->first, it->second, cg);
			add_edge(s, t);
		}
	}
	myres.dn_eval_time = myres.search_time + myres.reassignment_time;
	myres.dn_eval_time_AAN = mt.elapsed();
}

template<class T>
void dynamic_mnr_scc<T>::add_edges(int start_index, int end_index)
{
	myres.inc_edges += (end_index - start_index);
	my_timer mt;
	//std::vector<pair<unsigned int, unsigned int>>::iterator it;
	for (int i = start_index; i < end_index; i++) {
		unsigned int s, t;
		s = ds.find(inc_edges[i].first);
		t = ds.find(inc_edges[i].second);
		if (s != t) {													// s == t, e.g., the edge in the same component. do not add edge. 
			my_timer mt = my_timer();
			boost::add_edge(inc_edges[i].first, inc_edges[i].second, cg);
			add_edge(s, t);
		}
	}
	myres.dn_eval_time = myres.search_time + myres.reassignment_time;
	record_cur_per_info(myres);
}



template<class T>
void dynamic_mnr_scc<T>::add_edge(vd_t t, vd_t h)
{
	
	unsigned int tn2v(n2v[t]); // lower bound
	unsigned int hn2v(n2v[h]); // upper bound
	bool has_cycle = false;
	std::vector<unsigned int> reachable;

	if (tn2v < hn2v) {
		myres.set_invalid_edges(myres.get_invalid_edges() + 1);

		in_component[t] = true;
		my_timer mt;
		mnr_dfs(h, tn2v, hn2v, reachable, has_cycle);
		myres.set_search_time(myres.get_search_time() + mt.elapsed());
		myres.set_aff_region(myres.get_aff_region() + (hn2v - tn2v));

		if (has_cycle) {
			std::vector<unsigned int> cycle_nodes;
			std::vector<unsigned int>::iterator it;
			for (it = reachable.begin(); it != reachable.end(); it++) {
				if (in_component[*it]) {
					cycle_nodes.push_back(*it);
				}
			}
			assert(cycle_nodes.size() > 1);
			mt = my_timer();
			unsigned int s = cycle_nodes.front();
			for (it = cycle_nodes.begin(); it != cycle_nodes.end(); it++) {
				ds.join(s, *it);
			}
			s = ds.find(cycle_nodes.front());
			ds.collapse_scc(s, cycle_nodes);
			//condense_per_scc(s, cycle_nodes, g, g.ds);
			std::vector<unsigned int> cycle_nodes2;
			for (int i = 0; i < cycle_nodes.size(); i++) {
				cycle_nodes2.insert(cycle_nodes2.end(), ds.scc_nodes[cycle_nodes[i]].begin(), ds.scc_nodes[cycle_nodes[i]].end());
			}
			condense_per_scc(s, cycle_nodes2, cg, ds);
			myres.set_reassign_time(myres.get_reassign_time() + mt.elapsed());

			mnr_shift_cycle(tn2v, hn2v);
		}
		else {
			mt = my_timer();
			mnr_shift(tn2v, hn2v);
			myres.set_reassign_time(myres.get_reassign_time() + mt.elapsed());
		}
		in_component[t] = false;
		visited[t] = false;
		visited[h] = false;
		for_each(reachable.begin(), reachable.end(), [this](unsigned int a) {
			in_component[a] = false;
			visited[a] = false;
		});
	}

}

template<class T>
void dynamic_mnr_scc<T>::mnr_dfs(vd_t h, unsigned int lb, unsigned int ub, std::vector<unsigned int>& reachable, bool & flag)
{
	out_iterator i, iend;
	visited[h] = true;
	reachable.push_back(h);
	std::vector<unsigned int>::iterator it;
	for (it = ds.scc_nodes[h].begin(); it != ds.scc_nodes[h].end(); it++) {
		for (std::tie(i, iend) = boost::out_edges(*it, cg); i != iend; ++i) {
			//vd_t w(target(*i, g));
			unsigned int w = ds.find(target(*i, cg));
			unsigned int wn2v(n2v[w]);
			if (wn2v == lb) {
				flag = true;
				in_component[w] = true;
				if (!visited[w]) {
					mnr_dfs(w, lb, ub, reachable, flag);
				}
			}
			else if (wn2v > lb && !visited[w]) {
				mnr_dfs(w, lb, ub, reachable, flag);
			}
			in_component[h] = in_component[h] || in_component[w];
		}
	}
}

template<class T>
void dynamic_mnr_scc<T>::mnr_shift(unsigned int lb, unsigned int ub)
{
	std::vector<vd_t> dfs_node, normal_node;
	std::vector<unsigned int> aff_orders;
	unsigned int i;
	for (i = ub; i >= lb; i--) {
		vd_t w(v2n[i]);
		if (w == -1)  continue;
		aff_orders.push_back(i);
		in_component[w] = false;
		if (visited[w]) {
			visited[w] = false;
			dfs_node.push_back(w);
		}
		else {
			normal_node.push_back(w);
		}
	}
	std::vector<vd_t>::iterator it1(normal_node.begin());
	std::vector<unsigned int>::iterator it2(aff_orders.begin());
	assert((dfs_node.size() + normal_node.size()) <= aff_orders.size());
	for (; it1 != normal_node.end(); it1++) {
		n2v[*it1] = *it2;
		v2n[*it2] = *it1;
		it2++;
	}
	for (it1 = dfs_node.begin(); it1 != dfs_node.end(); it1++) {
		n2v[*it1] = *it2;
		v2n[*it2] = *it1;
		it2++;
	}
	while (it2 != aff_orders.end()) {
		v2n[*it2] = -1;
		it2++;
	}
}

template<class T>
void dynamic_mnr_scc<T>::mnr_shift_cycle(unsigned int lb, unsigned int ub)
{
	std::vector<vd_t> dfs_node, normal_node;
	std::vector<unsigned int> aff_orders;
	unsigned int i;
	bool is_first = true;
	for (i = ub; i >= lb; i--) {
		vd_t w(v2n[i]);
		if (w == -1)  continue;
		aff_orders.push_back(i);
		if (in_component[w] && is_first) {
			dfs_node.push_back(ds.find(w));
			is_first = false;
			in_component[w] = false;
			v2n[i] = -1;
		}
		else if (in_component[w] && !is_first) {
			in_component[w] = false;
			v2n[i] = -1;
		}
		else {
			if (visited[w]) {
				visited[w] = false;
				dfs_node.push_back(w);
			}
			else {
				normal_node.push_back(w);
			}
		}
	}
	std::vector<vd_t>::iterator it1(normal_node.begin());
	std::vector<unsigned int>::iterator it2(aff_orders.begin());
	assert((dfs_node.size() + normal_node.size()) <= aff_orders.size());
	for (; it1 != normal_node.end(); it1++) {
		n2v[*it1] = *it2;
		v2n[*it2] = *it1;
		it2++;
	}
	for (it1 = dfs_node.begin(); it1 != dfs_node.end(); it1++) {
		n2v[*it1] = *it2;
		v2n[*it2] = *it1;
		it2++;
	}
	while (it2 != aff_orders.end()) {
		v2n[*it2] = -1;
		it2++;
	}
}

