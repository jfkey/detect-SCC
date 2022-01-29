#include "dynamic_HKMST_scc.h"
#include <algorithm>

template<typename T>
dynamic_hkmst_scc<T>::dynamic_hkmst_scc()
{
}

template<typename T>
dynamic_hkmst_scc<T>::dynamic_hkmst_scc(CG & cg, std::vector<ahrsz_priority_value<T>>& n2v, T & pspace, disjoint_set & ds, std::vector<pair<unsigned int, unsigned int>>& inc_edges, const int np)
	: cg(cg), n2v(n2v), _pspace(pspace), ds(ds), inc_edges(inc_edges), out_its(np), in_its(np)
{

	std::cout << "Init Dynamic HKMST Algorithm" << std::endl;
	self_name = "Dynamic-HKMST-";
	myres.set_dynamic_name(self_name);
	_flooring = std::vector<ahrsz_ext_priority_value<T> >(np, minus_infinity);
	_visited = std::vector<bool>(np, false);
	_inK = std::vector<bool>(np, false);
	_inF = std::vector<bool>(np, false);
	_inB = std::vector<bool>(np, false);
	is_in = std::vector<bool>(np, false);
	is_out = std::vector<bool>(np, false);
	in_component = std::vector<bool>(np, false);

}



template<typename T>
void dynamic_hkmst_scc<T>::add_edges()
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
			myres.add_edge_time += (mt.elapsed());
			if ((n2v[s] < n2v[t]) || (n2v[s] == n2v[t])) {				// t -> h, i.e,  pmap(t) > pmap(h). when not holds this priority, discovery, reassignment.
				add_edge(s, t);
			}
		}
	}
	myres.dn_eval_time = myres.search_time + myres.reassignment_time;
	myres.dn_eval_time_AAN = mt.elapsed();
}


template<typename T>
void dynamic_hkmst_scc<T>::add_edges(int start_index, int end_index)
{
	
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
				add_edge(s, t);
			}
		}
	}
	myres.dn_eval_time = myres.search_time + myres.reassignment_time;
	record_cur_per_info(myres);
}



template<typename T>
void dynamic_hkmst_scc<T>::add_edge(vd_t t, vd_t h)
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
		soft_threshold_search(t, h, forw_nodes, back_nodes, FA, FP, BA, BP, has_cycle, cycle_nodes, threshold);
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
				if (n2v[threshold] < n2v[*it] && has_out_after_search(*it)) {
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
			F_greater.push_back(t);
			std::for_each(F_greater.begin(), F_greater.end(), fun);
			std::for_each(B_less.begin(), B_less.end(), fun);
			
			// unmark state info
			auto fun2 = [this ](vd_t t) {
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

			myres.set_aff_region( myres.get_aff_region() + remaining_K.size());
			my_timer mt;
			reassignment_with_cycle(remaining_K);
			myres.set_reassign_time(myres.get_reassign_time() + mt.elapsed());
		}
		else {
			// reassignment for affected region without the cycle.  
			// threshold = max{v, F| x in F and out(x) is nuot null . 
			// threshold, F>  B< 
			//threshold = t; 
			my_timer mt;
			reassignment(forw_nodes, back_nodes, t, h);
			myres.set_reassign_time(myres.get_reassign_time() + mt.elapsed());

		}
	}
}


template<typename T>
bool dynamic_hkmst_scc<T>::has_out_after_search(vd_t node)
{
	out_scc_iterator node_it;
	if (is_out[node] == false) {
		return has_first_out(node);
	}
	else {
		return has_next_out( out_its[node]);
	}
}

template<typename T>
auto dynamic_hkmst_scc<T>::first_out(vd_t & node)
{
	//return out_scc_iterator(); 
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
			if (scc_item_it.first != scc_item_it.second) {
				out_its[node].oei_begin = scc_item_it.first;
				out_its[node].oei_end = scc_item_it.second;
				break;
			}
			out_its[node].it_begin++;
		}
		return out_its[node];
	}
	else {
		/*out_its[node].oei_begin++;
		while (out_its[node].oei_begin == out_its[node].oei_end && out_its[node].it_begin != out_its[node].it_end) {
			out_its[node].it_begin++;
			if (out_its[node].it_begin != out_its[node].it_end) {
				out_its[node].oei_begin = boost::out_edges(*(out_its[node].it_begin), cg).first;
				out_its[node].oei_end = boost::out_edges(*(out_its[node].it_begin), cg).second;
			}
		}
		return out_its[node];*/
		out_its[node].oei_begin++;
		if (out_its[node].oei_begin == out_its[node].oei_end) {
			while (out_its[node].it_begin != out_its[node].it_end) {
				out_its[node].it_begin++;
				if (out_its[node].it_begin != out_its[node].it_end) {
					out_its[node].oei_begin = boost::out_edges(*(out_its[node].it_begin), cg).first;
					out_its[node].oei_end = boost::out_edges(*(out_its[node].it_begin), cg).second;
					if (out_its[node].oei_begin != out_its[node].oei_end) {
						break;
					}
				}
			}
		}
		return out_its[node];
	}

}

template<typename T>
bool dynamic_hkmst_scc<T>::out_tail(out_scc_iterator & out)
{
	bool flag = false;
	out_iterator e1, e2;
	std::vector<unsigned int>::iterator i1, i2;

	e1 = out.oei_begin;
	e2 = out.oei_end;
	i1 = out.it_begin;
	i2 = out.it_end;

	while (i1 != i2) {
		if (e1 != e2) {
			flag = true;
			break;
		}
		i1 = std::next(i1, 1);
		if (i1 != i2) {
			e1 = boost::out_edges(*i1, cg).first;
			e2 = boost::out_edges(*i1, cg).second;
		}
	}
	return flag;
}



template<typename T>
bool dynamic_hkmst_scc<T>::has_first_out(vd_t & node)
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
		if (oei_begin != oei_end) {
			flag = true;
			break;
		}
		it_begin++;
		if (it_begin != it_end) {
			oei_begin = boost::out_edges(*it_begin, cg).first;
			oei_end = boost::out_edges(*it_begin, cg).second;
		}
	}
	return flag;
}

template<typename T>
bool dynamic_hkmst_scc<T>::has_next_out(out_scc_iterator & out)
{
	bool flag = false;
	out_iterator e1, e2, e3, e4;
	std::vector<unsigned int>::iterator i1, i2;

	e1 = out.oei_begin;
	e2 = out.oei_end;
	i1 = out.it_begin;
	i2 = out.it_end;
	// current node has next out ? 
	if (e1 != e2 && std::next(e1, 1) != e2) {
		flag = true; 
		return flag; 
	}
	// other nodes have next out ? 
	i1 = std::next(i1, 1);
	if (i1 != i2) {
		e1 = boost::out_edges(*i1, cg).first;
		e2 = boost::out_edges(*i1, cg).second;
	}
	while (i1 != i2) {
		if (e1 != e2) {
			flag = true;
			break;
		}
		i1 = std::next(i1, 1);
		if (i1 != i2) {
			e1 = boost::out_edges(*i1, cg).first;
			e2 = boost::out_edges(*i1, cg).second;
		}
	}
	return flag;

}

template<typename T>
auto dynamic_hkmst_scc<T>::first_in(vd_t & node)
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
			if (scc_item_it.first != scc_item_it.second) {
				in_its[node].iei_begin = scc_item_it.first;
				in_its[node].iei_end = scc_item_it.second;
				break;
			}
			in_its[node].it_begin++;
		}
		return in_its[node];
	}
	else {
		/*in_its[node].iei_begin++;
		while (in_its[node].iei_begin == in_its[node].iei_end && in_its[node].it_begin != in_its[node].it_end) {
			in_its[node].it_begin++;
			if (in_its[node].it_begin != in_its[node].it_end) {
				in_its[node].iei_begin = boost::in_edges(*(in_its[node].it_begin), cg).first;
				in_its[node].iei_end = boost::in_edges(*(in_its[node].it_begin), cg).second;
			}
		}
		return in_its[node];*/
		in_its[node].iei_begin++;
		if (in_its[node].iei_begin == in_its[node].iei_end ) {
			while (in_its[node].it_begin != in_its[node].it_end) {
				in_its[node].it_begin++;
				if (in_its[node].it_begin != in_its[node].it_end) {
					in_its[node].iei_begin = boost::in_edges(*(in_its[node].it_begin), cg).first;
					in_its[node].iei_end = boost::in_edges(*(in_its[node].it_begin), cg).second;
					if (in_its[node].iei_begin != in_its[node].iei_end) {
						break;
					}
				}
			}
		}
		return in_its[node];
	}
}

template<typename T>
bool dynamic_hkmst_scc<T>::in_tail(in_scc_iterator & in)
{
	bool flag = false;
	in_iterator e1, e2;
	std::vector<unsigned int>::iterator i1, i2;

	e1 = in.iei_begin;
	e2 = in.iei_end;
	i1 = in.it_begin;
	i2 = in.it_end;

	while (i1 != i2) {
		if (e1 != e2) {
			flag = true;
			break;
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
bool dynamic_hkmst_scc<T>::has_first_in(vd_t & node)
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
		if (iei_begin != iei_end) {
			flag = true;
			break;
		}
		it_begin++;
		if (it_begin != it_end) {
			iei_begin = boost::in_edges(*it_begin, cg).first;
			iei_end = boost::in_edges(*it_begin, cg).second;
		}
	}
	return flag;
}

template<typename T>
bool dynamic_hkmst_scc<T>::has_next_in(in_scc_iterator & in)
{
	bool flag = false;
	in_iterator e1, e2;
	std::vector<unsigned int>::iterator i1, i2;

	e1 = in.iei_begin;
	e2 = in.iei_end;
	i1 = in.it_begin;
	i2 = in.it_end;
	// current node has next in ? 
	if (e1 != e2 && std::next(e1, 1) != e2) {			
		flag = true;
		return flag;
	}
	// other nodes has next in ? 		
	i1 = std::next(i1, 1);
	if (i1 != i2) {
		e1 = boost::in_edges(*i1, cg).first;
		e2 = boost::in_edges(*i1, cg).second;
	}
	while (i1 != i2) {
		if (e1 != e2) {
			flag = true;
			break;
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
void dynamic_hkmst_scc<T>::find_cycle(vd_t head)
{
	out_iterator i, iend;
	_visited[head] = true;
	std::vector<unsigned int>::iterator it;
	for (it = ds.scc_nodes[head].begin(); it != ds.scc_nodes[head].end(); it++) {
		for (std::tie(i, iend) = boost::out_edges(*it, cg); i != iend; ++i) {
			//vd_t w = (target(*i, *this));
			unsigned int w = ds.find(target(*i, cg));
			if ((_inB[w] || _inF[w]) && _visited[w] == false) find_cycle(w);
			in_component[head] = in_component[head] || in_component[w];
		}
	}
}

template<typename T>
void dynamic_hkmst_scc<T>::search_step(vd_t & u, vd_t & z, std::list<vd_t>& FA, std::list<vd_t>& BA, std::list<vd_t>& forw_nodes, std::list<vd_t>& back_nodes, bool & has_cycle)
{
	out_scc_iterator out_u, out_x;
	in_scc_iterator in_z, in_y;
	vd_t out_u_x, in_z_y;
	bool has_no = false, has_ni = false;
	out_u = first_out(u);
	in_z = first_in(z);
	out_u_x = ds.find(boost::target(*(out_u.oei_begin), cg));
	in_z_y = ds.find(boost::source(*(in_z.iei_begin), cg));
	has_no = has_next_out(out_u);
	has_ni = has_next_in(in_z);
	
	if (has_no == false)
		FA.remove(u);
	if (has_ni == false)
		BA.remove(z);

	if (_inF[in_z_y] || _inB[out_u_x]) has_cycle = true;		// report the cycle . 

	if (!_inF[out_u_x]) {
		_inF[out_u_x] = true;
		forw_nodes.push_back(out_u_x);
		out_x = first_out(out_u_x);
		//if (out_x.first != out_x.second)
		if (out_tail(out_x)) {
			FA.push_back(out_u_x);
		}
		is_out[out_u_x] = false;
	}
	if (!_inB[in_z_y]) {
		_inB[in_z_y] = true;
		back_nodes.push_back(in_z_y);
		in_y = first_in(in_z_y);
		//if (in_y.first != in_y.second) {
		if( in_tail(in_y) ){
			BA.push_back(in_z_y);
		}
		is_in[in_z_y] = false;
	}

}

template<typename T>
void dynamic_hkmst_scc<T>::soft_threshold_search(vd_t & v, vd_t & w, std::list<vd_t>& forw_nodes, std::list<vd_t>& back_nodes, std::list<vd_t>& FA, std::list<vd_t>& FP, std::list<vd_t>& BA, std::list<vd_t>& BP, bool & has_cycle, std::vector<unsigned int>& cycle_nodes, vd_t & threshold)
{
	std::list<vd_t>::iterator it_list;
	threshold = v;

	forw_nodes.push_back(w);
	_inF[w] = true;
	back_nodes.push_back(v);
	_inB[v] = true;

	if ( has_first_out(w)) {
		FA.push_back(w);
	}
	if (has_first_in(v)) {
		BA.push_back(v);
	}
	while (!FA.empty() && !BA.empty()) {
		vd_t u = FA.front();
		vd_t z = BA.front();
		if (n2v[z] < n2v[u] || n2v[z] == n2v[u]) { // n2v[z] == n2v[u] to identify the scc. 
			search_step(u, z, FA, BA, forw_nodes, back_nodes, has_cycle);
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
		find_cycle(w);
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
void dynamic_hkmst_scc<T>::reassignment(std::list<vd_t>& forw_nodes, std::list<vd_t>& back_nodes, vd_t & v, vd_t & w)
{
	std::list<vd_t> F_greater, B_less;
	std::list<vd_t>::iterator it;
	std::list<vd_t>::reverse_iterator rit;
	ahrsz_ext_priority_value<T> val1, val2;
	vd_t threshold = v;

	for (it = forw_nodes.begin(); it != forw_nodes.end(); it++) {
		if (n2v[threshold] < n2v[*it] && has_out_after_search(*it)) {
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
	myres.set_aff_region(myres.get_aff_region() + F_greater.size() + B_less.size());

	if (threshold != v) {
		if (n2v[threshold].base() == _pspace.end() || std::next(n2v[threshold].base()) == _pspace.end()) {
			val2 = plus_infinity;
		}
		else {
			val2 = ahrsz_ext_priority_value<T>(std::next(n2v[threshold].base()), _pspace);  // val_t = prev(S)
		}
		val1 = ahrsz_ext_priority_value<T>(n2v[threshold].base(), _pspace);
		for (rit = F_greater.rbegin(); rit != F_greater.rend(); rit++) {
			ahrsz_ext_priority_value<T> ep = compute_priority(val1, val2);
			n2v[*rit] = ahrsz_priority_value<T>(ep.base(), _pspace);
			val1 = ahrsz_ext_priority_value<T>(n2v[*rit].base(), _pspace);
		}
		for (it = B_less.begin(); it != B_less.end(); it++) {
			ahrsz_ext_priority_value<T> ep = compute_priority(val1, val2);
			n2v[*it] = ahrsz_priority_value<T>(ep.base(), _pspace);
			val1 = ahrsz_ext_priority_value<T>(n2v[*it].base(), _pspace);
		}
	}
	else {
		if (n2v[threshold].base() == _pspace.begin()) {
			val1 = minus_infinity;
		}
		else {
			val1 = ahrsz_ext_priority_value<T>(_pspace.previous(n2v[threshold].base()), _pspace);
		}
		val2 = ahrsz_ext_priority_value<T>(n2v[threshold].base(), _pspace);
		for (rit = F_greater.rbegin(); rit != F_greater.rend(); rit++) {
			ahrsz_ext_priority_value<T> ep = compute_priority(val1, val2);
			n2v[*rit] = ahrsz_priority_value<T>(ep.base(), _pspace);
			val1 = ahrsz_ext_priority_value<T>(n2v[*rit].base(), _pspace);
		}
		for (it = B_less.begin(); it != B_less.end(); it++) {
			ahrsz_ext_priority_value<T> ep = compute_priority(val1, val2);
			n2v[*it] = ahrsz_priority_value<T>(ep.base(), _pspace);
			val1 = ahrsz_ext_priority_value<T>(n2v[*it].base(), _pspace);
		}
	}
}

template<typename T>
void dynamic_hkmst_scc<T>::reassignment_with_cycle(std::vector<vd_t>& K)
{
	for (std::vector<vd_t>::iterator i(K.begin()); i != K.end(); ++i) {
		_flooring[*i] = minus_infinity;
		_inB[*i] = false;					// unmark states. 
		_inF[*i] = false;
		_inK[*i] = true;
	}

	// first pass - compute flooring Information
	std::vector<vd_t> rto; // reverse topological order
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
void dynamic_hkmst_scc<T>::compute_flooring(vd_t n, std::vector<vd_t>& rto)
{
	_visited[n] = true;
	out_iterator i, iend;
	std::vector<unsigned int>::iterator it;
	for (it = ds.scc_nodes[n].begin(); it != ds.scc_nodes[n].end(); it++) {
		for (tie(i, iend) = out_edges(*it, cg); i != iend; ++i) {
			//vd_t j(target(*i, *this));
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
ahrsz_ext_priority_value<T> dynamic_hkmst_scc<T>::compute_priority(ahrsz_ext_priority_value<T> floor, ahrsz_ext_priority_value<T> ceiling)
{
	my_timer mt2; 
	ceiling = ahrsz_ext_priority_value<T>(plus_infinity);
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
ahrsz_ext_priority_value<T> dynamic_hkmst_scc<T>::compute_ceiling(vd_t & v)
{
	ahrsz_ext_priority_value<T> ceiling(plus_infinity);
	in_iterator j, jend;
	std::vector<unsigned int>::iterator it;
	for (it = ds.scc_nodes[v].begin(); it != ds.scc_nodes[v].end(); it++) {
		for (tie(j, jend) = in_edges(*it, cg); j != jend; ++j) {
			unsigned int w = ds.find(source(*j, cg));
			ceiling = std::min(ceiling, ahrsz_ext_priority_value<T>(n2v[w]));
		}
	}
	
	return ceiling;
}


