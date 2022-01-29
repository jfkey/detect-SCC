#include "bat_PK2_scc.h"


template<class T>
batch_PK2_scc<T>::batch_PK2_scc()
{
}

template<class T>
batch_PK2_scc<T>::batch_PK2_scc(CG & _cg, disjoint_set & _ds, std::vector<T>& _n2v, std::vector<pair<unsigned int, unsigned int>>& _inc_edges,
	unsigned int _np, unsigned int _bs) :cg(_cg), ds(_ds), n2v(_n2v), inc_edges(_inc_edges), np(_np), bs(_bs), para(np)
{
	
	topo_scc.reserve(np);
	self_name = "Batch-PK2-";
	
	visited = std::vector<bool>(np, false);
	in_reachables = std::vector<bool>(np, false);

	v2n = std::vector<vd_t>(np + 1, -1);
	for (unsigned int i = 0; i < np; i++) {
		v2n[n2v[i]] = ds.find(i);
	}
	myres.set_dynamic_name(self_name);
}

template<class T>
void batch_PK2_scc<T>::add_edges()
{
	std::vector<my_tuple> batch;
	batch.reserve(bs);
	int tmp_b = 0;
	my_timer mt; 
	for (int i = 0; i < inc_edges.size(); ) {
		while (tmp_b < bs && i < inc_edges.size()) {
			unsigned int s, t;
			s = ds.find(inc_edges[i].first);
			t = ds.find(inc_edges[i].second);
			if (s != t) {
				// both valid edges and invalid edges are added. 
				mt = my_timer();
				boost::add_edge(inc_edges[i].first, inc_edges[i].second, cg);
				myres.add_edge_time += mt.elapsed();
				if (n2v[s] < n2v[t]) {
					myres.invalid_edges++;
					batch.push_back(my_tuple(n2v[s], s, t));
				}
			}
			i++;
			tmp_b++;
		}
		if (batch.size() > 0) {
			add_edge_bat(batch);
			batch.clear();
			tmp_b = 0;
		}
		tmp_b = 0;
	}
	myres.dn_eval_time = myres.search_time + myres.reassignment_time;
}


template<class T>
void batch_PK2_scc<T>::add_edges(int start_index, int end_index)
{
	myres.inc_edges += (end_index - start_index);
	std::vector<my_tuple> batch;
	batch.reserve(bs);
	int tmp_b = 0;
	my_timer mt;
	for (int i = start_index; i < end_index; ) {
		while (tmp_b < bs && i < end_index ) {
			unsigned int s, t;
			s = ds.find(inc_edges[i].first);
			t = ds.find(inc_edges[i].second);
			if (s != t) {
				// both valid edges and invalid edges are added. 
				mt = my_timer();
				boost::add_edge(inc_edges[i].first, inc_edges[i].second, cg);
				myres.add_edge_time += mt.elapsed();
				if (n2v[s] < n2v[t]) {
					myres.invalid_edges++; 
					batch.push_back(my_tuple(n2v[s], s, t));
				}
			}
			i++;
			tmp_b++;
		}
		if (batch.size() > 0) {
			mt = my_timer();
			add_edge_bat(batch);
			batch.clear();
			tmp_b = 0;
			myres.dn_eval_time_AAN += mt.elapsed();
		}
		tmp_b = 0;
	}
	myres.dn_eval_time = myres.search_time + myres.reassignment_time;
	record_cur_per_info(myres);
}


template<class T>
void batch_PK2_scc<T>::eval_single(int start_index, int end_index)
{
	
}


template<class T>
void batch_PK2_scc<T>::add_edge_bat(std::vector<my_tuple> &batch)
{
	// first:  max reaching topo value of t 
	// second: the tail of the edge (s -> t)
	std::vector<pair<unsigned int, unsigned int> > reachables;
	auto comp = [this](const my_tuple &a, const my_tuple &b) {
		return a.first < b.first;
	};
	sort(batch.begin(), batch.end(), comp);
	unsigned int ub = 0;
	my_timer mt; 

	for (int i = 0; i < batch.size(); i++) {
		if (n2v[batch[i].second] > ub && !reachables.empty()) {
			shift(ub, reachables);
			reachables.clear();
		}
		if (!visited[ ds.find(batch[i].third)]) {
			//myres.set_invalid_edges(myres.get_invalid_edges() + 1);
			mt = my_timer();
			dfs_f( ds.find(batch[i].third), n2v[ds.find(batch[i].second)], reachables);
			myres.set_search_time(myres.get_search_time() + mt.elapsed());
			mt = my_timer();
			for (int i = 0; i < topo_scc.size(); i++) {
				simple_sccs(topo_scc[i]);
			}
			topo_scc.clear();
			myres.condensation_time += mt.elapsed();
		}
		ub = std::max(ub, n2v[batch[i].third]);
		//ub = std::max(ub, n2v[ ds.find(batch[i].third) ]);
	}
	if (!reachables.empty()) {
		shift(ub, reachables);
	}
	batch.clear();
}

template<class T>
void batch_PK2_scc<T>::dfs_f(vd_t n, unsigned int lb, vector<pair<unsigned int, unsigned int> > &reachables)
{
	visited[n] = true;							// visited[] could be replaced by para.num[]
	out_iterator i, iend;
	std::vector<unsigned int>::iterator it;
	para.time++;
	para.num[n] = para.lowlink[n] = para.time;
	para.stack.push_back(n);
	para.in_stack[n] = true;

	for (it = ds.scc_nodes[n].begin(); it != ds.scc_nodes[n].end(); it++) {
		for (std::tie(i, iend) = boost::out_edges(*it, cg); i != iend; ++i) {
			unsigned int w = ds.find(target(*i, cg));
			if (n2v[w] >= lb) {
				if (para.num[w] == 0) {
					dfs_f(w, lb, reachables);
					if (para.lowlink[w] < para.lowlink[n]) para.lowlink[n] = para.lowlink[w];
				}
				else if (para.in_stack[w]) {
					if (para.num[w] < para.lowlink[n]) para.lowlink[n] = para.num[w];
				}
			}
		}
	}

	std::vector<unsigned int> tmp_component;
	if (para.lowlink[n] == para.num[n]) {
		unsigned int w = -1;
		do {
			w = para.stack.back(); para.stack.pop_back();
			para.in_stack[w] = false;
			tmp_component.push_back(w);
		} while (w != n);
		
		if (tmp_component.size() > 1) {
			topo_scc.push_back(tmp_component);
		}
	}
	reachables.push_back(make_pair(lb, n));
}



template<class T>
void batch_PK2_scc<T>::shift(unsigned int max_tv, std::vector<pair<unsigned int, unsigned int> > reachables)
{
	my_timer mt; 
	unsigned int aff_reg = 0;
	std::vector<pair<unsigned int, unsigned int> >::reverse_iterator i(reachables.rbegin()), iend(reachables.rend());
	unsigned int lindex(max_tv);
	while (i != iend) {
		unsigned int w = v2n[max_tv];
		if (w != -1) {
			v2n[max_tv] = -1;				// if has iterated current slot, v2n is set to  -1 . 
			if (visited[w]) {
				visited[w] = false;
			}
			else {
				v2n[lindex] = w;
				n2v[w] = lindex;
				--lindex;
			}
			while (i != iend && i->first == max_tv) {
				unsigned int i_s = ds.find(i->second);
				if (in_reachables[i_s] == false) {
					in_reachables[i_s] = true;
					v2n[lindex] = i_s;
					n2v[i_s] = lindex;
				}
				else {
					v2n[lindex] = -1; 
					n2v[i->second] = -1;
				}
				++i;
				--lindex;
				aff_reg++;
			}
		}
		--max_tv;
		aff_reg++; 
	}
	i = reachables.rbegin();
	iend = reachables.rend();
	for (; i != iend; i++) {
		// reset para states. 
		para.num[i->second] = 0;
		para.lowlink[i->second] = 0;
		para.time = 0;
		para.in_stack[i->second] = false;
		visited[i->second] = false;
		in_reachables[i->second] = false; 
	}
	myres.set_reassign_time(myres.get_reassign_time() + mt.elapsed());
	myres.set_aff_region(myres.get_aff_region() + aff_reg);
}

template<class T>
void batch_PK2_scc<T>::simple_sccs(std::vector<unsigned int>& cycle_nodes)
{
	std::vector<unsigned int>::iterator it; 
	assert(cycle_nodes.size() > 1);
	unsigned int s = cycle_nodes.front();
	for (it = cycle_nodes.begin(); it != cycle_nodes.end(); it++) {
		ds.join(s, *it);
	}
	s = ds.find(cycle_nodes.front());
	ds.collapse_scc(s, cycle_nodes);
	std::vector<unsigned int> cycle_nodes2;
	for (int i = 0; i < cycle_nodes.size(); i++) {
		cycle_nodes2.insert(cycle_nodes2.end(), ds.scc_nodes[cycle_nodes[i]].begin(), ds.scc_nodes[cycle_nodes[i]].end());
	}
	condense_per_scc(s, cycle_nodes2, cg, ds);
	cycle_nodes.clear();
}


