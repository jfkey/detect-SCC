#include "bat_AHRSZ_P_scc.h"
template<typename T>
inline batch_AHRSZ_P_scc<T>::batch_AHRSZ_P_scc()
{
}

template<typename T>
batch_AHRSZ_P_scc<T>::batch_AHRSZ_P_scc(CG & cg, std::vector<ahrsz_priority_value<T>>& n2v, T & pspace, disjoint_set & ds,
	std::vector<pair<unsigned int, unsigned int>>& inc_edges, const unsigned int np, unsigned int _bs) : cg(cg), n2v(n2v), _pspace(pspace), 
	ds(ds), inc_edges(inc_edges), bs(_bs)
{
	std::cout << "Init Dynamic Bat-AHRSZ+ Algorithm" << std::endl;
	self_name = "Bat-AHRSZ+";
	myres.set_dynamic_name(self_name);
	_flooring = std::vector<ahrsz_ext_priority_value<T> >(np, minus_infinity);
	_visited = std::vector<bool>(np, false);
	_inK = std::vector<bool>(np, false);
	_inF = std::vector<bool>(np, false);
	_inB = std::vector<bool>(np, false);
	in_component = std::vector<bool>(np, false);
	
	start_invalid_index = std::vector<std::vector<int>>(np);
	end_invalid_index = std::vector<std::vector<int>>(np);
	inv_rem_vis = std::vector<bool>(bs, false);
	inv_has_vis = std::vector<bool>(bs, false);
}


template<typename T>
void batch_AHRSZ_P_scc<T>::add_edges()
{
	//std::vector<pair<unsigned int, unsigned int>> invalid_edges;
	//invalid_edges.reserve(bs);
	//int tmp_b = 0;
	//my_timer mt;
	//for (int i = 0; i < inc_edges.size(); ) {
	//	while (tmp_b < bs && i < inc_edges.size()) {
	//		unsigned int s, t;
	//		s = ds.find(inc_edges[i].first);
	//		t = ds.find(inc_edges[i].second);
	//		if (s != t) {
	//			if (!(n2v[s] < n2v[t] || n2v[s] == n2v[t])) {
	//				boost::add_edge(inc_edges[i].first, inc_edges[i].second, cg);
	//			}
	//			else {
	//				//invalid_edges.push_back(std::make_pair(s, t));
	//				invalid_edges.push_back(std::make_pair(inc_edges[i].first, inc_edges[i].second));
	//			}
	//		}
	//		i++;
	//		tmp_b++;
	//	}
	//	if (invalid_edges.size() > 0) {
	//		if (invalid_edges[0].first == 6 && invalid_edges[0].second ==3) {
	//			add_edges_bat(invalid_edges);
	//		}
	//		else {
	//			add_edges_bat(invalid_edges);
	//		}
	//		invalid_edges.clear();
	//		tmp_b = 0;
	//	}
	//	tmp_b = 0;
	//}
	//myres.set_dn_total_time(mt.elapsed());
	//std::cout << "number of visited: " << myres.test_cnt << std::endl;
	add_edges_valid_edge_aware();
}





template<typename T>
void batch_AHRSZ_P_scc<T>::add_edges(int start_index, int end_index)
{
}


// valid_edges aware strategy
// 1. split the inc-edges into batchs {B1, ..., Bn}
// 2. add all of the valid edges of the batch Bi
// 3. add the invalid edges one by one. do reorder after doing discovery. the affect region is AFF. 
// 4. select the candidate invalid edges C which start node or end node in the AFF from the remain edges
// 5. add all of the candidate invalid edges C. 
// 6. repeat 3, 4, 5 while B is empty
// problem:
// 1. if two or three or more edges share the same affect region, but their start node and end node are not visited by the invalid edge. 
template<typename T>
void batch_AHRSZ_P_scc<T>::add_edges_valid_edge_aware()
{
	std::vector<pair<unsigned int, unsigned int>> invalid_edges;
	std::vector<pair<unsigned int, unsigned int>> new_invalid_edges;
	std::vector<pair<unsigned int, unsigned int>> tmp_v;
	std::vector<unsigned int>unmarks, unmarkt; 

	invalid_edges.reserve(bs);
	int tmp_b = 0;
	my_timer mt;
	int inv_idx = 0;		// invalid edge index 
	for (int i = 0; i < inc_edges.size(); ) {
		while (tmp_b < bs && i < inc_edges.size()) {
			unsigned int s, t;
			s = ds.find(inc_edges[i].first);
			t = ds.find(inc_edges[i].second);
			if (s != t) {
				if (!(n2v[s] < n2v[t] || n2v[s] == n2v[t])) {
					boost::add_edge(inc_edges[i].first, inc_edges[i].second, cg);
				}
				else {
					// construct invalid edges 
					start_invalid_index[s].push_back(inv_idx);
					end_invalid_index[t].push_back(inv_idx);
					unmarks.push_back(s);
					unmarkt.push_back(t);
					invalid_edges.push_back(inc_edges[i]);
					inv_idx++;
				}
			}
			i++;
			tmp_b++;
		}
		if (invalid_edges.size() > 0) {
			while (invalid_edges.size() > 0) {
				valid_edge_aware(invalid_edges, new_invalid_edges, unmarks, unmarkt);
				invalid_edges.clear();
				tmp_v = invalid_edges;
				invalid_edges = new_invalid_edges;
				new_invalid_edges = tmp_v;
			}

			invalid_edges.clear();
			tmp_b = 0;
		}
		tmp_b = 0;
		inv_idx = 0;
		// unmark start_invalid_index, end_invalid_index
		for (int i = 0; i < unmarks.size(); i++) {
			start_invalid_index[unmarks[i]].clear();
		}
		for (int i = 0; i < unmarkt.size(); i++) {
			end_invalid_index[unmarkt[i]].clear();
		}
		unmarks.clear();
		unmarkt.clear();
	}
	//myres.set_dn_total_time(mt.elapsed());
	//std::cout << "number of visited: " << myres.test_cnt << std::endl;
	myres.dn_eval_time = myres.search_time + myres.reassignment_time;

}

template<typename T>
void batch_AHRSZ_P_scc<T>::valid_edge_aware(std::vector<std::pair<unsigned int, unsigned int>>& invalid_edges, std::vector<std::pair<unsigned int, unsigned int>>& new_invalid_edges, std::vector<unsigned int>& unmarks, std::vector<unsigned int>unmarkt)
{
	std::list<std::pair<unsigned int, unsigned int>> poss_valid_cand;					//possible valid edge candidate
	std::list<std::pair<unsigned int, unsigned int>>::iterator lit; 

	std::vector<int> cand_idx; 
	std::vector<vd_t> K;
	std::vector<unsigned int> cycle_nodes;
	bool has_cycle = false;
	std::pair<unsigned int, unsigned int> tmp;
	unsigned int s, t;
	my_timer mt;
	for (int i = 0; i < invalid_edges.size(); i++) {
		if (inv_rem_vis[i] == true) {													// this invalid edges may change to valid edge. 
			poss_valid_cand.push_back(invalid_edges[i]);
			continue;
		}
		inv_has_vis[i] = true;
		tmp = invalid_edges[i];
		s = ds.find(tmp.first);
		t = ds.find(tmp.second);
		if (s != t) {
			myres.set_invalid_edges(myres.get_invalid_edges() + 1);
			boost::add_edge(tmp.first, tmp.second, cg);
			add_edge_single(s, t, has_cycle, K, cycle_nodes);
			
			// mark the index of the possible valid edges
			for (int j = 0; j < K.size(); j++) {
				cand_idx = start_invalid_index[K[j]];
				for (int l = 0; l < cand_idx.size(); l ++) {
					inv_rem_vis[cand_idx[l]] = true;
				}
				cand_idx = end_invalid_index[K[j]];
				for (int l = 0; l < cand_idx.size(); l ++) {
					inv_rem_vis[cand_idx[l]] = true;
				}
			}

			myres.set_aff_region(myres.get_aff_region() + K.size());
			if (K.size() != 0) {
				// debug.
				std::vector<unsigned int> aff_reg;
				std::vector<vd_t>::iterator it;
				for (it = K.begin(); it != K.end(); it++) {
					aff_reg.push_back(*it);
				}
				std::tuple<bool, unsigned int, unsigned int, std::vector<unsigned int> > invedge_affreg(true, s, t, aff_reg);
				myres.invedge_affreg.push_back(invedge_affreg);
				reassignment_tail(has_cycle, K, cycle_nodes);
			}

		}
	}
	// unmarks, start_invalid_index
	// unmarkt, end_invalid_index
	for (int i = 0; i < unmarks.size();i ++) {
		start_invalid_index[unmarks[i]].clear();
	}
	for (int i = 0; i < unmarkt.size(); i++) {
		end_invalid_index[unmarkt[i]].clear();
	}
	unmarks.clear();
	unmarkt.clear();

	if (poss_valid_cand.size() > 0) {
		// unmark inv_rem_vis, inv_has_vis, inv_idx, start_invalid_index, end_invalid_index
		// invalid edges remain visited, invalid edges have visited. 
		inv_rem_vis = std::vector<bool>(bs);
		inv_has_vis = std::vector<bool>(bs);
		int inv_idx = 0;
		// identify valid edges and invalid edges from candidate
		for (lit = poss_valid_cand.begin(); lit != poss_valid_cand.end(); lit++) {
			unsigned int cs, ct;
			cs = ds.find(lit->first);
			ct = ds.find(lit->second);
			if (cs != ct) {
				if (!(n2v[cs] < n2v[ct] || n2v[cs] == n2v[ct])) {
					boost::add_edge(lit->first, lit->second, cg);
				}
				else {
					// reconstruct invalid edges 
					start_invalid_index[cs].push_back(inv_idx);
					end_invalid_index[ct].push_back(inv_idx);
					unmarks.push_back(cs);
					unmarks.push_back(ct);
					new_invalid_edges.push_back(*lit);
					inv_idx++;
				}
			}
		}
	}
	

}



template<typename T>
void batch_AHRSZ_P_scc<T>::add_edges_bat(std::vector<std::pair<unsigned int, unsigned int>> & invalid_edges)
{
	std::vector<std::pair<unsigned int, unsigned int>> lazy_group;
	std::vector<vd_t> K;
	std::vector<unsigned int> cycle_nodes;
	std::vector<bool> is_invalid(invalid_edges.size(), true);
	bool has_cycle = false;
	std::pair<unsigned int, unsigned int> tmp;
	unsigned int s, t;

	my_timer mt;
	for (int i = 0; i < invalid_edges.size(); i++) {
		if (is_invalid[i] == false) continue;
		tmp = invalid_edges[i];
		s = ds.find(tmp.first);
		t = ds.find(tmp.second);
		//std::cout << tmp.first << " -> " << tmp.second << std::endl;

		if (s != t) {
			myres.set_invalid_edges(myres.get_invalid_edges() + 1);
			boost::add_edge(tmp.first, tmp.second, cg);
			add_edge_single(s, t, has_cycle, K, cycle_nodes);
			if (has_cycle) {
				myres.set_aff_region(myres.get_aff_region() + K.size());
				add_edge_single(s, t, has_cycle, K, cycle_nodes);
				refresh_invalid_edges(invalid_edges, i, is_invalid);
			}
		}
	}


	myres.set_aff_region(myres.get_aff_region() + K.size());
	if (K.size() != 0) {
		// debug.
		std::vector<unsigned int> aff_reg;
		std::vector<vd_t>::iterator it;
		for (it = K.begin(); it != K.end(); it++) {
			aff_reg.push_back(*it);
		}
		std::tuple<bool, unsigned int, unsigned int, std::vector<unsigned int> > invedge_affreg(true, s, t, aff_reg);
		myres.invedge_affreg.push_back(invedge_affreg);

		reassignment_tail(has_cycle, K, cycle_nodes);
	}
}


template<typename T>
void batch_AHRSZ_P_scc<T>::add_edge_single(vd_t t, vd_t h, bool & has_cycle, std::vector<vd_t>& K, std::vector<unsigned int>& cycle_nodes)
{
	if (has_cycle == false) {
		my_timer mt;
		long t2 = clock();
		if ((n2v[t] < n2v[h]) || (n2v[t] == n2v[h])) {			// t -> h, i.e,  n2v(t) > n2v(h). when not holds this priority, discovery, reassignment.
			discovery_bat(t, h, K, has_cycle, cycle_nodes);
			// debug.
			std::vector<vd_t>::iterator it;
			std::vector<unsigned int> aff_reg;
			for (it = K.begin(); it != K.end(); it++) {
				aff_reg.push_back(*it);
			}
			std::tuple<bool, unsigned int, unsigned int, std::vector<unsigned int> > invedge_affreg(true, t, h, aff_reg);
			myres.invedge_affreg.push_back(invedge_affreg);

		}
		else {
		}
		myres.set_search_time(myres.get_search_time() + mt.elapsed());
	}
	else {
		std::vector<vd_t> remaining_K;
		std::vector<vd_t>::iterator it;
		assert(cycle_nodes.size() > 1 && "the number of vertices in the cycle must be greater than 1. dynamic_AHRSZ");
		unsigned int s = cycle_nodes.front();
		for (std::vector<unsigned int>::iterator it = cycle_nodes.begin(); it != cycle_nodes.end(); it++) {
			ds.join(s, *it);
		}
		s = ds.find(cycle_nodes.front());
		ds.collapse_scc(s, cycle_nodes);
		std::vector<unsigned int> cycle_nodes2;
		for (int i = 0; i < cycle_nodes.size(); i++) {
			cycle_nodes2.insert(cycle_nodes2.end(), ds.scc_nodes[cycle_nodes[i]].begin(), ds.scc_nodes[cycle_nodes[i]].end());
		}
		// collapse scc in the original graph. 
		condense_per_scc(s, cycle_nodes2, cg, ds);

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
		my_timer mt = my_timer();
		reassignment(K);
		//myres.set_reassign_time(myres.get_reassign_time() + mt.elapsed());
		// debug.
		std::vector<unsigned int> aff_reg;
		for (it = K.begin(); it != K.end(); it++) {
			aff_reg.push_back(*it);
		}
		std::tuple<bool, unsigned int, unsigned int, std::vector<unsigned int> > invedge_affreg(true, t, h, aff_reg);
		myres.invedge_affreg.push_back(invedge_affreg);

		has_cycle = false;
		cycle_nodes.clear();
		K.clear();
	}
}



template<typename T>
void batch_AHRSZ_P_scc<T>::discovery_bat(vd_t tail, vd_t head, std::vector<vd_t>& K, bool & has_cycle, std::vector<unsigned int>& cycle_nodes)
{
	vd_t f(head), b(tail), tmp;
	max_priority_queue ForwFron(n2v);  // the forward frontier
	min_priority_queue BackFron(n2v);  // the backward frontier
	unsigned int ForwEdges = out_degree_scc(head, cg, ds);
	unsigned int BackEdges = in_degree_scc(tail, cg, ds);

	ForwFron.push(head);
	BackFron.push(tail);

	// notice, I use visited to indicate when a node is on one of the queues.
	_inB[tail] = true;
	_inF[head] = true;
	_visited[tail] = true;
	_visited[head] = true;

	while ((n2v[b] < n2v[f] || n2v[b] == n2v[f]) && !ForwFron.empty() && !BackFron.empty()) {
		unsigned int u = std::min(ForwEdges, BackEdges);
		ForwEdges -= u;
		BackEdges -= u;
		if (_inK[b]) BackEdges = 0;		// visited info ?
		if (_inK[f]) ForwEdges = 0;

		if (ForwEdges == 0) {
			if (!_inK[f]) {
				K.push_back(f);
				_inK[f] = true;
			}
			ForwFron.pop();
			_visited[f] = false;
			grow_forw2(f, ForwFron, has_cycle, K);
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
			grow_back2(b, BackFron, has_cycle, K);
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

	if (has_cycle) {
		// reset visited info, to find the cycle. 
		for (std::vector<vd_t>::iterator i = K.begin(); i != K.end(); i++) {
			_visited[*i] = false;
		}
		in_component[tail] = true;
		find_cycle(head);
		for (std::vector<vd_t>::iterator i = K.begin(); i != K.end(); i++) {
			if (in_component[*i]) {
				cycle_nodes.push_back(*i);
				in_component[*i] = false;
			}
		}
	}

	// unmark status.
	
	while (!ForwFron.empty()) {
		tmp = ForwFron.top();
		_visited[tmp] = false;
		_inF[tmp] = false;
		ForwFron.pop();
		myres.num_visited++;
	}
	while (!BackFron.empty()) {
		tmp = BackFron.top();
		_visited[tmp] = false;
		_inB[tmp] = false;
		BackFron.pop();
		myres.num_visited++; 
	}

	for (std::vector<vd_t>::iterator i = K.begin(); i != K.end(); i++) {
		_visited[*i] = false;
		_inB[*i] = false;
		_inF[*i] = false;
	}
	
}

template<typename T>
void batch_AHRSZ_P_scc<T>::reassignment(std::vector<vd_t>& K)
{
	my_timer mt;
	// Initialise temporary data 
	for (typename std::vector<vd_t>::iterator i(K.begin()); i != K.end(); ++i) {
		_flooring[*i] = minus_infinity;
		_inK[*i] = true;
		_inB[*i] = false;		// unmark states. 
		_inF[*i] = false;
		_visited[*i] = false;
	}

	// first pass - compute flooring Information
	std::vector<vd_t> rto;	// reverse topological order
	for (std::vector<vd_t>::iterator i(K.begin()); i != K.end(); ++i) {
		if (!_visited[*i]) { compute_flooring(*i, rto); }
	}
	
	int nv = 0; 
	// second pass - perform the reassignment
	for (std::vector<vd_t>::reverse_iterator i(rto.rbegin()); i != rto.rend(); ++i) {
		nv++;
		ahrsz_ext_priority_value<T> ceiling; 
		if (_flooring[*i].minus_inf()) {
			ceiling = ahrsz_ext_priority_value<T>(_pspace.begin(), _pspace);
		}
		else {
			if (_flooring[*i].base() == _pspace.end() || std::next(_flooring[*i].base()) == _pspace.end()) {
				ceiling = ahrsz_ext_priority_value<T>(plus_infinity);
			}
			else {
				ceiling = ahrsz_ext_priority_value<T>(std::next(_flooring[*i].base()), _pspace);
			}
		}
		ahrsz_ext_priority_value<T> ep(compute_priority(_flooring[*i], ceiling));
		n2v[*i] = ahrsz_priority_value<T>(ep.base(), _pspace);
		_visited[*i] = false; // unmark state info. 
	}
	// reset visited information	unmark state info. 
	for (std::vector<vd_t>::iterator i(K.begin()); i != K.end(); ++i) {
		_visited[*i] = false;
		_inK[*i] = false;
	}
	myres.set_reassign_time(myres.get_reassign_time() + mt.elapsed());
}

template<typename T>
void batch_AHRSZ_P_scc<T>::reassignment_tail(bool & has_cycle, std::vector<vd_t>& K, std::vector<unsigned int>& cycle_nodes)
{
	if (has_cycle) {
		std::vector<vd_t> remaining_K;
		std::vector<vd_t>::iterator it;
		assert(cycle_nodes.size() > 1 && "the number of vertices in the cycle must be greater than 1. dynamic_AHRSZ");
		unsigned int s = cycle_nodes.front();
		for (std::vector<unsigned int>::iterator it = cycle_nodes.begin(); it != cycle_nodes.end(); it++) {
			ds.join(s, *it);
		}
		s = ds.find(cycle_nodes.front());

		ds.collapse_scc(s, cycle_nodes);
		std::vector<unsigned int> cycle_nodes2;
		for (int i = 0; i < cycle_nodes.size(); i++) {
			cycle_nodes2.insert(cycle_nodes2.end(), ds.scc_nodes[cycle_nodes[i]].begin(), ds.scc_nodes[cycle_nodes[i]].end());
		}
		// collapse scc in the original graph. 
		condense_per_scc(s, cycle_nodes2, cg, ds);

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
		my_timer mt;
		reassignment(remaining_K);
		//myres.set_reassign_time(myres.get_reassign_time() + mt.elapsed());

	}
	else {
		my_timer mt;
		reassignment(K);
		//myres.set_reassign_time(myres.get_reassign_time() + mt.elapsed());
	}

	K.clear();
	cycle_nodes.clear();
	has_cycle = false;
}

template<typename T>
void batch_AHRSZ_P_scc<T>::refresh_invalid_edges(std::vector<std::pair<unsigned int, unsigned int>>& invalid_edges, int idx, std::vector<bool>& is_invalid)
{
	assert(idx < invalid_edges.size());
	unsigned int s, t;
	for (; idx < invalid_edges.size(); idx++) {
		s = ds.find(invalid_edges[idx].first);
		t = ds.find(invalid_edges[idx].second);
		if (!((n2v[s] < n2v[t]) || (n2v[s] == n2v[t]))) {
			// NOTE: to be modified. 
			//boost::add_edge(s, t, cg);
			boost::add_edge(invalid_edges[idx].first, invalid_edges[idx].second, cg);
			is_invalid[idx] = false;
		}
	}
}

template<typename T>
void batch_AHRSZ_P_scc<T>::refresh_invalid_edges(std::vector<std::pair<unsigned int, unsigned int>>& invalid_edges, int idx, std::vector<bool>& is_invalid, std::vector<std::pair<unsigned int, unsigned int>>& valid_edges)
{
	assert(idx < invalid_edges.size());
	unsigned int s, t;
	for (; idx < invalid_edges.size(); idx++) {
		s = ds.find(invalid_edges[idx].first);
		t = ds.find(invalid_edges[idx].second);
		if (!((n2v[s] < n2v[t]) || (n2v[s] == n2v[t]))) {
			// NOTE: to be modified. 
			//boost::add_edge(s, t, cg);
			//boost::add_edge(invalid_edges[idx].first, invalid_edges[idx].second, cg);
			//valid_edges.push_back(invalid_edges[idx]);
			is_invalid[idx] = false;
		}
	}
}



template<typename T>
void batch_AHRSZ_P_scc<T>::grow_forw2(vd_t node, max_priority_queue & max_pq, bool & has_cycle, std::vector<vd_t>& K)
{
	out_iterator i, iend;
	std::vector<unsigned int>::iterator it;
	for (it = ds.scc_nodes[node].begin(); it != ds.scc_nodes[node].end(); ++it) {
		for (std::tie(i, iend) = boost::out_edges(*it, cg); i != iend; ++i) {
			//vd_t w = (target(*i, *this));
			vd_t w = ds.find(target(*i, cg));
			if (_inB[w]) {
				has_cycle = true;
				// _visited[w] == false i.e., node w is neither in max_pq nor min_pq. 
				// (_visited[w] == true && _inF[w] == false)  i.e., node w is in min_pq but not in max_pq. 
				if (_visited[w] == false || (_visited[w] == true && _inF[w] == false)) {
					max_pq.push(w);
					_visited[w] = true;
				}
				_inF[w] = true;
			}

			if (!_inF[w]) {
				/*if (_inK[w]) {
					_inF[w] = true;
					grow_forw2(w, max_pq, has_cycle, K);	// bat-ahrsz-1-11-7-6  if( ! has_cycle )   -> if( ! _inF[w] )  
				}
				else {
					if (_visited[w] == false || (_visited[w] == true && _inF[w] == false)) {
						max_pq.push(w);
						_visited[w] = true;
					}
					_inF[w] = true;
				}*/
				if (_visited[w] == false || (_visited[w] == true && _inF[w] == false)) {
					max_pq.push(w);
					_visited[w] = true;
				}
				_inF[w] = true;
			}
		}
	}
}

template<typename T>
void batch_AHRSZ_P_scc<T>::grow_back2(vd_t node, min_priority_queue & min_pq, bool & has_cycle, std::vector<vd_t>& K)
{
	std::vector<unsigned int>::iterator it;
	in_iterator i, iend;
	for (it = ds.scc_nodes[node].begin(); it != ds.scc_nodes[node].end(); it++) {
		for (std::tie(i, iend) = boost::in_edges(*it, cg); i != iend; ++i) {
			//vd_t w = (boost::source(*i, *this));
			vd_t w = ds.find(boost::source(*i, cg));
			if (_inF[w]) {
				has_cycle = true;
				if (_visited[w] == false || (_visited[w] == true && _inB[w] == false)) {
					_visited[w] = true;
					_inB[w] = true;
					min_pq.push(w);
				}

			}
			if (!_inB[w]) {
				/*if (_inK[w]) {
					_inB[w] = true;
					grow_back2(w, min_pq, has_cycle, K);
				}
				else {
					if (_visited[w] == false || (_visited[w] == true && _inB[w] == false)) {
						_visited[w] = true;
						min_pq.push(w);
					}
					_inB[w] = true;
				}*/
				if (_visited[w] == false || (_visited[w] == true && _inB[w] == false)) {
					min_pq.push(w);
					_visited[w] = true;
				}
				_inB[w] = true;
			}
		}
	}
}

template<typename T>
void batch_AHRSZ_P_scc<T>::find_cycle(vd_t head)
{
	std::vector<unsigned int>::iterator it;
	out_iterator i, iend;
	_visited[head] = true;
	for (it = ds.scc_nodes[head].begin(); it != ds.scc_nodes[head].end(); it++) {
		for (std::tie(i, iend) = boost::out_edges(*it, cg); i != iend; ++i) {
			//vd_t w = (target(*i, *this));
			unsigned int w = ds.find(boost::target(*i, cg));
			if ((_inK[w]) && _visited[w] == false) find_cycle(w);
			in_component[head] = in_component[head] || in_component[w];
		}
	}
}

template<typename T>
void batch_AHRSZ_P_scc<T>::compute_flooring(vd_t n, std::vector<vd_t>& rto)
{
	_visited[n] = true;
	std::vector<unsigned int>::iterator it;
	out_iterator i, iend;
	for (it = ds.scc_nodes[n].begin(); it != ds.scc_nodes[n].end(); it++) {
		for (tie(i, iend) = out_edges(*it, cg); i != iend; ++i) {
			//vd_t j(target(*i, cg));
			unsigned int j = ds.find(target(*i, cg));
			if (_inK[j]) {
				// successor is in K
				if (!_visited[j]) { compute_flooring(j, rto); }
				_flooring[n] = std::max(_flooring[j], _flooring[n]);
			}
			else {
				// successor is not in K
				_flooring[n] = std::max(_flooring[n], ahrsz_ext_priority_value<T>(n2v[j]));
			}
		}
	}
	rto.push_back(n);
}

template<typename T>
ahrsz_ext_priority_value<T> batch_AHRSZ_P_scc<T>::compute_ceiling(vd_t v)
{
}

// the ceiling may cause bug in debug model.  
// due to the case of computing ceiling: std::next(_flooring[*i].base()) == _pspace.end() not considered
template<typename T>
ahrsz_ext_priority_value<T> batch_AHRSZ_P_scc<T>::compute_priority(ahrsz_ext_priority_value<T> floor, ahrsz_ext_priority_value<T> ceiling)
{
	assert(floor < ceiling);
	ahrsz_ext_priority_value<T> candidate;
	if (floor.minus_inf()) {
		_pspace.push_front();
		candidate = ahrsz_ext_priority_value<T>(_pspace.begin(), _pspace);
	}
	else {
		assert(!floor.plus_inf());
		candidate = ahrsz_ext_priority_value<T>(_pspace.insert_after(floor.base()), _pspace);
	}
	// std::cout << "compute priority" << a2str(floor) << "\t" << a2str(candidate) <<"\t"<< a2str(ceiling) << std::endl; 
	assert(floor < candidate && candidate < ceiling && "create new priority error");
	return candidate;
}

template<typename T>
batch_AHRSZ_P_scc<T>::~batch_AHRSZ_P_scc()
{
}

