#include "bat_AHRSZ_scc.h"
template<typename T>
inline batch_AHRSZ_scc<T>::batch_AHRSZ_scc()
{
}

template<typename T>
batch_AHRSZ_scc<T>::batch_AHRSZ_scc(CG & cg, std::vector<ahrsz_priority_value<T>>& n2v, T & pspace, disjoint_set & ds, 
	std::vector<pair<unsigned int, unsigned int>>& inc_edges, const unsigned int np, unsigned int _bs) : cg(cg), n2v(n2v), _pspace(pspace), ds(ds), inc_edges(inc_edges), bs(_bs)
{
	std::cout << "Init Dynamic Bat-AHRSZ Algorithm" << std::endl;
	self_name = "Bat-AHRSZ-";
	myres.set_dynamic_name(self_name);
	_flooring = std::vector<ahrsz_ext_priority_value<T> >(np, minus_infinity);
	_visited = std::vector<bool>(np, false);
	_inK = std::vector<bool>(np, false);
	_inF = std::vector<bool>(np, false);
	_inB = std::vector<bool>(np, false);
	in_component = std::vector<bool>(np, false);
	_indegree = std::vector<unsigned int >(np, 0);

}

template<typename T>
void batch_AHRSZ_scc<T>::add_edges()
{
	std::vector<pair<unsigned int, unsigned int>> invalid_edges;
	invalid_edges.reserve(bs);
	int tmp_b = 0;
	my_timer mt;
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
					invalid_edges.push_back(std::make_pair(inc_edges[i].first, inc_edges[i].second));
				}
			}
			i++;
			tmp_b++;
		}
		if (invalid_edges.size() > 0) {
			add_edges_bat(invalid_edges);
			invalid_edges.clear();
			tmp_b = 0;
		}
		tmp_b = 0;
	}
	myres.dn_eval_time = myres.search_time + myres.reassignment_time;
	myres.dn_eval_time_AAN = mt.elapsed();

}

template<typename T>
void batch_AHRSZ_scc<T>::add_edges(int start_index, int end_index)
{
	myres.inc_edges += (end_index - start_index);
	std::vector<pair<unsigned int, unsigned int>> invalid_edges;
	invalid_edges.reserve(bs);
	int tmp_b = 0;
	my_timer mt; 
	for (int i = start_index; i < end_index;  ) {
		while (tmp_b < bs && i < end_index ) {
			unsigned int s, t;
			s = ds.find(inc_edges[i].first);
			t = ds.find(inc_edges[i].second);
			if (s != t) {
				if (!(n2v[s] < n2v[t] || n2v[s] == n2v[t])) {
					boost::add_edge(inc_edges[i].first, inc_edges[i].second, cg);
				}
				else {
					invalid_edges.push_back(std::make_pair(inc_edges[i].first, inc_edges[i].second));
				}
			}
			i++;
			tmp_b++;
		}
		if (invalid_edges.size() > 0) {
			add_edges_bat(invalid_edges);
			invalid_edges.clear();
			tmp_b = 0;
		}
		tmp_b = 0;
	}
	myres.dn_eval_time = myres.search_time + myres.reassignment_time;
	//myres.dn_eval_time_AAN = mt.elapsed();
	record_cur_per_info(myres);
}



template<typename T>
void batch_AHRSZ_scc<T>::add_edges_bat(std::vector<std::pair<unsigned int, unsigned int>> & invalid_edges)
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
		reassignment_tail(has_cycle, K, cycle_nodes);
	}
}



template<typename T>
void batch_AHRSZ_scc<T>::add_edge_single(vd_t t, vd_t h, bool & has_cycle, std::vector<vd_t>& K, std::vector<unsigned int>& cycle_nodes)
{
	if (has_cycle == false) {
		
		long t2 = clock();
		if ((n2v[t] < n2v[h]) || (n2v[t] == n2v[h])) {			// t -> h, i.e,  n2v(t) > n2v(h). when not holds this priority, discovery, reassignment.
			my_timer mt;
			discovery_bat(t, h, K, has_cycle, cycle_nodes);
			myres.set_search_time(myres.get_search_time() + mt.elapsed());
		}
		else {
		}
	}
	else {
		my_timer mt;
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
		myres.condensation_time += mt.elapsed();

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


		mt = my_timer();
		//reassignment(K);
		reassignment_mini_priority(K);
		myres.set_reassign_time(myres.get_reassign_time() + mt.elapsed());

		has_cycle = false;
		cycle_nodes.clear();
		K.clear();
	}
}



template<typename T>
void batch_AHRSZ_scc<T>::discovery_bat(vd_t tail, vd_t head, std::vector<vd_t>& K, bool & has_cycle, std::vector<unsigned int>& cycle_nodes)
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
void batch_AHRSZ_scc<T>::reassignment(std::vector<vd_t>& K)
{
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
	// second pass - perform the reassignment
	for (std::vector<vd_t>::reverse_iterator i(rto.rbegin()); i != rto.rend(); ++i) {
		ahrsz_ext_priority_value<T> ep(compute_priority(_flooring[*i], compute_ceiling(*i)));
		n2v[*i] = ahrsz_priority_value<T>(ep.base(), _pspace);
		_visited[*i] = false; // unmark state info. 
	}
	
	// reset visited information	unmark state info. 
	for (std::vector<vd_t>::iterator i(K.begin()); i != K.end(); ++i) {
		_visited[*i] = false;
		_inK[*i] = false;
	}
}

template<typename T>
void batch_AHRSZ_scc<T>::reassignment_tail(bool & has_cycle, std::vector<vd_t>& K, std::vector<unsigned int>& cycle_nodes)
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
		//reassignment(remaining_K);
		reassignment_mini_priority(remaining_K);
		myres.set_reassign_time(myres.get_reassign_time() + mt.elapsed());
		
	}
	else {
		my_timer mt;
		//reassignment(K);
		reassignment_mini_priority(K);
		myres.set_reassign_time(myres.get_reassign_time() + mt.elapsed());
	}

	K.clear();
	cycle_nodes.clear();
	has_cycle = false;
}

// Minimize the priority of creating 
template<typename T>
void batch_AHRSZ_scc<T>::reassignment_mini_priority(std::vector<vd_t>& K)
{

	typedef std::pair<vd_t, ahrsz_ext_priority_value<T> > Q_e;
	typedef my_less<Q_e, select2nd<Q_e> > Q_c;
	typedef std::priority_queue<Q_e, std::vector<Q_e>, Q_c> Q_t;
	Q_t Q;

	// Initialise temporary data 
	for (std::vector<vd_t>::iterator i(K.begin()); i != K.end(); ++i) {
		_flooring[*i] = minus_infinity;
		_inK[*i] = true;
		_inB[*i] = false;																						// unmark states. 
		_inF[*i] = false;
	}

	// first pass - compute flooring Information
	std::vector<vd_t> rto;																						// reverse topological order
	for (std::vector<vd_t>::iterator i(K.begin()); i != K.end(); ++i) {
		if (!_visited[*i]) { compute_flooring(*i, rto); }
	}

	// init Q. 
	std::vector<unsigned int>::iterator it;
	for (std::vector<vd_t>::iterator i(K.begin()); i != K.end(); ++i) {
		in_iterator j, jend;
		ahrsz_ext_priority_value<T> ceiling(plus_infinity);
		unsigned int kid = 0;
		for (it = ds.scc_nodes[*i].begin(); it != ds.scc_nodes[*i].end(); it++) {
			for (tie(j, jend) = in_edges(*it, cg); j != jend; ++j) {
				unsigned int w = ds.find(source(*j, cg));
				ceiling = std::min(ceiling, ahrsz_ext_priority_value<T>(n2v[w]));
				if (_inK[w]) { kid++; }
			}
		}
		_indegree[*i] = kid;
		if (kid == 0) {
			Q.push(Q_e(*i, ceiling));
		}
	}



	std::vector<vd_t> Z;
	while (!Q.empty()) {
		Z.clear();
		ahrsz_ext_priority_value<T> Z_ceil = Q.top().second;
		ahrsz_ext_priority_value<T> Z_floor = minus_infinity;
		while (!Q.empty() && Q.top().second == Z_ceil) {
			vd_t x = Q.top().first;
			Z_floor = std::max(Z_floor, _flooring[x]);
			Z.push_back(x);
			Q.pop();
		}

		assert(Q.empty() || (Q.top().second < Z_ceil));

		// all of Z must get same priority
		ahrsz_ext_priority_value<T> Z_p = compute_priority(Z_floor, Z_ceil);

		for (std::vector<vd_t>::const_iterator i(Z.begin()); i != Z.end(); ++i) {
			n2v[*i] = ahrsz_priority_value<T>(Z_p.base(), _pspace);
			out_iterator oi, oiend;
			for (it = ds.scc_nodes[*i].begin(); it != ds.scc_nodes[*i].end(); it++) {
				for (tie(oi, oiend) = out_edges(*it, cg); oi != oiend; ++oi) {
					//vd_t j(target(*i, cg));
					unsigned int j = ds.find(target(*oi, cg));
					if (_inK[j] && --_indegree[j] == 0) {
						Q.push(std::make_pair(j, compute_ceiling(j)));
					}
				}
			}
		}
	}

	// reset visited information	unmark state info. 
	for (std::vector<vd_t>::iterator i(K.begin()); i != K.end(); ++i) {
		_visited[*i] = false;
		_inK[*i] = false;
		_indegree[*i] = 0;
	}
}



template<typename T>
void batch_AHRSZ_scc<T>::refresh_invalid_edges(std::vector<std::pair<unsigned int, unsigned int>>& invalid_edges, int idx, std::vector<bool>& is_invalid)
{
	assert(idx < invalid_edges.size());
	unsigned int s, t;
		for (; idx < invalid_edges.size(); idx++) {
		s = ds.find(invalid_edges[idx].first);
		t = ds.find(invalid_edges[idx].second);
		if (!((n2v[s] < n2v[t]) || (n2v[s] == n2v[t]))) {
			if (is_invalid[idx] == true) {
				boost::add_edge(invalid_edges[idx].first, invalid_edges[idx].second, cg);
				is_invalid[idx] = false;
			}
		}
	}
}



template<typename T>
void batch_AHRSZ_scc<T>::grow_forw2(vd_t node, max_priority_queue & max_pq, bool & has_cycle, std::vector<vd_t>& K)
{
	out_iterator i, iend;
	std::vector<unsigned int>::iterator it;
	for (it = ds.scc_nodes[node].begin(); it != ds.scc_nodes[node].end(); ++it) {
		for (std::tie(i, iend) = boost::out_edges(*it, cg); i != iend; ++i) {
			//vd_t w = (target(*i, *this));
			vd_t w = ds.find(target(*i, cg));
			if (_inB[w] ) {
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
				if (_inK[w]) {
					_inF[w] = true;
					grow_forw2(w, max_pq, has_cycle, K);	
				}
				else {
					if (_visited[w] == false || (_visited[w] == true && _inF[w] == false)) {
						max_pq.push(w);
						_visited[w] = true;
					}
					_inF[w] = true;
				}
			}
		}
	}
}

template<typename T>
void batch_AHRSZ_scc<T>::grow_back2(vd_t node, min_priority_queue & min_pq, bool & has_cycle, std::vector<vd_t>& K)
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
				if (_inK[w]) {
					_inB[w] = true;
					grow_back2(w, min_pq, has_cycle, K);
				}
				else {
					if (_visited[w] == false || (_visited[w] == true && _inB[w] == false)) {
						_visited[w] = true;
						min_pq.push(w);
					}
					_inB[w] = true;
				}
			}
			
			
		}
	}
}

template<typename T>
void batch_AHRSZ_scc<T>::find_cycle(vd_t head)
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
void batch_AHRSZ_scc<T>::compute_flooring(vd_t n, std::vector<vd_t>& rto)
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
ahrsz_ext_priority_value<T> batch_AHRSZ_scc<T>::compute_ceiling(vd_t v)
{
	ahrsz_ext_priority_value<T> ceiling(plus_infinity);
	std::vector<unsigned int>::iterator it;
	in_iterator j, jend;
	for (it = ds.scc_nodes[v].begin(); it != ds.scc_nodes[v].end(); it++) {
		for (tie(j, jend) = in_edges(*it, cg); j != jend; ++j) {
			unsigned int w = ds.find(source(*j, cg));
			ceiling = std::min(ceiling, ahrsz_ext_priority_value<T>(n2v[w]));
		}
	}
	return ceiling;
	
}

template<typename T>
ahrsz_ext_priority_value<T> batch_AHRSZ_scc<T>::compute_priority(ahrsz_ext_priority_value<T> floor, ahrsz_ext_priority_value<T> ceiling)
{
	assert(floor < ceiling);
	ahrsz_ext_priority_value<T> candidate;
	if (floor.minus_inf()) {
		candidate = ahrsz_ext_priority_value<T>(_pspace.begin(), _pspace);
		if (candidate >= ceiling) {
			_pspace.push_front();
			candidate = ahrsz_ext_priority_value<T>(_pspace.begin(), _pspace);
		}
	}
	else {
		assert(!floor.plus_inf());
		candidate = floor;
		++candidate;
		if (candidate.base() == _pspace.end() || candidate >= ceiling) {
			// we need a new priority
			candidate = ahrsz_ext_priority_value<T>(_pspace.insert_after(floor.base()), _pspace);
		}
	}
	assert(floor < candidate && candidate < ceiling && "create new priority error");
	return candidate;
}

template<typename T>
batch_AHRSZ_scc<T>::~batch_AHRSZ_scc()
{
}

