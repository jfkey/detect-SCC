#include "dynamic_AHRSZ_scc.h"

template<typename T>
dynamic_ahrsz_scc<T>::dynamic_ahrsz_scc() { }

template<typename T>
dynamic_ahrsz_scc<T>::dynamic_ahrsz_scc(CG & cg, std::vector<ahrsz_priority_value<T>>& n2v, T & pspace, disjoint_set & ds, 
	std::vector<pair<unsigned int, unsigned int>>& inc_edges, const int np): cg(cg), n2v(n2v), _pspace(pspace), ds (ds), inc_edges(inc_edges)
{
	std::cout << "Init Dynamic AHRSZ Algorithm" << std::endl;
	self_name = "Dynamic-AHRSZ-";
	myres.set_dynamic_name(self_name);
	_flooring = std::vector<ahrsz_ext_priority_value<T> >(np, minus_infinity);
	_visited = std::vector<bool>(np, false);
	_inK = std::vector<bool>(np, false);
	_inF = std::vector<bool>(np, false);
	_inB = std::vector<bool>(np, false);
	in_component = std::vector<bool>(np, false);
	_indegree = std::vector<unsigned int>(np, 0);
}


template<typename T>
void dynamic_ahrsz_scc<T>::add_edges() 
{
	my_timer mtf; 
	int inv = 0; 
	std::vector<pair<unsigned int, unsigned int>>::iterator it; 
	for (it = inc_edges.begin(); it != inc_edges.end(); it++) {
		//std::cout << it->first << ", " << it->second << std::endl;
		unsigned int s, t;
		s = ds.find(it->first);
		t = ds.find(it->second);
		if (s != t) {													// s == t, e.g., the edge in the same component. do not add edge. 
			my_timer mt = my_timer();
			boost::add_edge(it->first, it->second, cg);
			myres.add_edge_time += (mt.elapsed());
			if ((n2v[s] < n2v[t]) || (n2v[s] == n2v[t])) {			// t -> h, i.e,  pmap(t) > pmap(h). when not holds this priority, discovery, reassignment.
				my_timer mtt = my_timer();
				inv++;
				if (it->first ==2 && it->second ==0) {
					do_ahrsz(s, t);
				}
				else {
					do_ahrsz(s, t);
				}
				myres.invalid_all_cnt++; 
			}
		}
	}
	//myres.set_dn_total_time(mtf.elapsed());
	myres.dn_eval_time = myres.search_time + myres.reassignment_time;
	std::cout << "inv: " << inv << std::endl; 
}


// start index is included, end index is not include 
template<typename T>
void dynamic_ahrsz_scc<T>::add_edges(int start_index, int end_index) {
	myres.inc_edges += (end_index - start_index);
	for (int i = start_index; i < end_index; i++) {
		unsigned int s, t;
		s = ds.find(inc_edges[i].first);
		t = ds.find(inc_edges[i].second);
		if (s != t) {													// s == t, e.g., the edge in the same component. do not add edge. 
			my_timer mt = my_timer();
			boost::add_edge(inc_edges[i].first, inc_edges[i].second, cg);
			myres.add_edge_time += (mt.elapsed());
			if ((n2v[s] < n2v[t]) || (n2v[s] == n2v[t])) {			// t -> h, i.e,  pmap(t) > pmap(h). when not holds this priority, discovery, reassignment.
				do_ahrsz(s, t);
			}
		}
	}
	myres.dn_eval_time = myres.search_time + myres.reassignment_time;
	record_cur_per_info(myres);

}


template<typename T>
void dynamic_ahrsz_scc<T>::do_ahrsz(vd_t t, vd_t h)
{
	std::vector<vd_t> K, remaining_K;
	std::vector<unsigned int> cycle_nodes;
	std::vector<vd_t>::iterator it;
	bool has_cycle = false;
	my_timer mt;
	discovery(t, h, K, has_cycle, cycle_nodes);
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
		ds.collapse_scc(s, cycle_nodes);											// collapse scc in the original graph. 
		
		my_timer mt;
		// NOTE: different with the fully condensation graph. i.e., V1. 
		std::vector<unsigned int> cycle_nodes2;
		for (int i = 0; i < cycle_nodes.size(); i++) {
			cycle_nodes2.insert(cycle_nodes2.end(), ds.scc_nodes[cycle_nodes[i]].begin(), ds.scc_nodes[cycle_nodes[i]].end() );
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
	//reassignment(K);
	reassignment_mini_priority(K);
	myres.set_reassign_time(myres.get_reassign_time() + mt.elapsed());

}

template<typename T>
void dynamic_ahrsz_scc<T>::find_cycle(vd_t head)
{
	out_iterator i, iend;
	std::vector<unsigned int>::iterator it; 
	_visited[head] = true; 

	for (it = ds.scc_nodes[head].begin(); it != ds.scc_nodes[head].end(); it++) {
		for (std::tie(i, iend) = boost::out_edges(*it, cg); i != iend; ++i) {
			//vd_t w = (target(*i, cg));
			unsigned int w = ds.find(target(*i, cg));
			if (_inK[w] && _visited[w] == false) find_cycle(w);
			in_component[head] = in_component[head] || in_component[w];
		}
	}
	
}

template<typename T>
void dynamic_ahrsz_scc<T>::discovery(vd_t tail, vd_t head, std::vector<vd_t>& K, bool & has_cycle, std::vector<unsigned int>& cycle_nodes)
{
	vd_t f(head), b(tail), tmp;
	std::vector<unsigned int>::iterator it; 
	max_priority_queue ForwFron(n2v );											// the forward frontier										
	min_priority_queue BackFron(n2v );											// the backward frontie
	unsigned int ForwEdges = out_degree_scc(head, cg, ds);
	unsigned int BackEdges = in_degree_scc(tail, cg, ds);
	ForwFron.push(head);
	BackFron.push(tail);

	// notice, I uses visited to indicate when  a node is on one of the queues.
	_visited[head] = true;
	_visited[tail] = true;
	_inB[tail] = true;	// indicate tail and head has been visited when detecting the cycle.
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
					// vd_t w(target(*i, cg));
					unsigned int w = ds.find(target(*i, cg));
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

			if (ForwFron.empty()) {
				f = tail;
			}
			else {
				f = ForwFron.top();
			}
			//ForwEdges = out_degree(f, cg);
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
					//vd_t w(source(*i, cg));
					unsigned int w = ds.find(source(*i, cg));
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
			if (BackFron.empty()) {
				b = head;
			}
			else {
				b = BackFron.top();
			}
			//BackEdges = in_degree(b, cg);
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
		find_cycle(head);
		for (typename std::vector<vd_t>::iterator i = K.begin(); i != K.end(); i++) {
			if (in_component[*i]) {
				cycle_nodes.push_back(*i);
				in_component[*i] = false;
			}
		}
	}

}



template<typename T>
void dynamic_ahrsz_scc<T>::reassignment(std::vector<vd_t>& K)
{
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

	// second pass - perform the reassignment
	for (std::vector<vd_t>::reverse_iterator i(rto.rbegin()); i != rto.rend(); ++i) {
		ahrsz_ext_priority_value<T> ep(compute_priority(_flooring[*i], compute_ceiling(*i )));
		n2v[*i] = ahrsz_priority_value<T>(ep.base(), _pspace);
		_visited[*i] = false;																					// unmark state info. 
	}

	// reset visited information	unmark state info. 
	for (std::vector<vd_t>::iterator i(K.begin()); i != K.end(); ++i) {
		_visited[*i] = false;
		_inK[*i] = false;
	}
}



// Minimize the priority of creating 
template<typename T>
void dynamic_ahrsz_scc<T>::reassignment_mini_priority(std::vector<vd_t>& K)
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
				if (_inK[w]) {kid++; }
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

		assert(Q.empty() || (Q.top().second < Z_ceil) );

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
						Q.push(std::make_pair(j, compute_ceiling(j )));
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
void dynamic_ahrsz_scc<T>::compute_flooring(vd_t n, std::vector<vd_t>& rto)
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
ahrsz_ext_priority_value<T> dynamic_ahrsz_scc<T>::compute_priority(ahrsz_ext_priority_value<T> floor, ahrsz_ext_priority_value<T> ceiling)
{
	my_timer mt;
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
ahrsz_ext_priority_value<T> dynamic_ahrsz_scc<T>::compute_ceiling(vd_t v)
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
