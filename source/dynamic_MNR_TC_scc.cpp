#include "dynamic_MNR_TC_scc.h"

template<class T>
 dynamic_tc_mnr_scc<T>::dynamic_tc_mnr_scc()
{
}

 template<class T>
 dynamic_tc_mnr_scc<T>::dynamic_tc_mnr_scc(CG & _cg, disjoint_set & _ds, std::vector<T>& _n2v, std::vector<unsigned int>& _max_ccs, 
	 std::vector<bool>& _is_in_unorder, std::vector<int>& _is_in_same_year, std::vector<unsigned int>& _n2y, unsigned int _max_year, 
	 unsigned int _min_year, std::vector<pair<unsigned int, unsigned int>>& _inc_edges, unsigned int _np) 
	 : cg(_cg), ds(_ds), n2v(_n2v), max_ccs(_max_ccs), is_in_same_year(_is_in_same_year), is_in_unorder(_is_in_unorder),
	 n2y(_n2y), min_year(_min_year), max_year(_max_year), inc_edges(_inc_edges), np(_np)
 {
	 in_component = std::vector<bool>(np, false);
	 visited = std::vector<bool>(np, false);
	 std::cout << "Init Dynamic Time-Coupled MNR Algorithm" << std::endl;
	 self_name = "Dynamic-TC-MNR-";
	 myres.set_dynamic_name(self_name);
	 
	 std::vector<vd_t> v2n_u = std::vector<vd_t>(np, -1);
	 v2ns.push_back(v2n_u);
	 for (int i = 1; i < max_ccs.size(); i++) {
		 max_ccs[i] += 1;														// NOTE THAT, the value of the array n2v is from 1
		 std::vector<vd_t> v2n_y = std::vector<vd_t>(max_ccs[i], -1);
		 v2ns.push_back(v2n_y);
	 }

	 for (unsigned int i = 0; i < np; i++) {
		 if (is_in_unorder[i]) {
			 v2ns[0][n2v[i]] = ds.find(i);
		 }
		 else if (is_in_same_year[i] != 0) {
			 int idx = is_in_same_year[i] - min_year + 1;
			 v2ns[idx][n2v[i]] = ds.find(i);									// NOTE THAT, the value of the array n2v is from 1
		 }
	 }

 }


 template<class T>
 void dynamic_tc_mnr_scc<T>::add_edges()
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
				 if (n2y[s] < n2y[t]) {
					 add_edge_alg1(s, t);
				 }
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
 void dynamic_tc_mnr_scc<T>::add_edges(int start_index, int end_index)
 {
	 myres.inc_edges += (end_index - start_index);
	 my_timer mt;
	 int type = 0;
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
				 if (n2y[s] < n2y[t]) {
					 add_edge_alg1(s, t);
				 }
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
 void dynamic_tc_mnr_scc<T>::add_edge_Gm(vd_t t, vd_t h)
 {
	 unsigned int tn2v(n2v[t]); // lower bound
	 unsigned int hn2v(n2v[h]); // upper bound
	 bool has_cycle = false;
	 std::vector<unsigned int> reachable;
	 my_timer mt; 
	 if (tn2v < hn2v) {
		 in_component[t] = true;
		 mnr_dfs_Gm(h, tn2v, hn2v, reachable, has_cycle);
		 myres.invalid_edges += 1; 
		 myres.aff_region += (hn2v-tn2v);
		 myres.search_time += mt.elapsed();

		 if (has_cycle) {
			 // identify the scc and condense scc. 
			 mt = my_timer();
			 std::vector<unsigned int> cycle_nodes;
			 std::vector<unsigned int>::iterator it;
			 for (it = reachable.begin(); it != reachable.end(); it++) {
				 if (in_component[*it]) {
					 cycle_nodes.push_back(*it);
				 }
			 }
			 assert(cycle_nodes.size() > 1);
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
			 myres.condensation_time += mt.elapsed();
			 mt = my_timer();
			 mnr_shift_cycle_Gm(tn2v, hn2v);
			 myres.reassignment_time += mt.elapsed();
		 }
		 else {
			 mt = my_timer();
			 mnr_shift_Gm(tn2v, hn2v );
			 myres.reassignment_time += mt.elapsed();
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
 void dynamic_tc_mnr_scc<T>::add_edge_alg1(vd_t t, vd_t h)
 {
	 /*
	 my_timer mt; 
	 std::vector<unsigned int> reachable;
	 std::vector<unsigned int>::iterator it, it_ds;
	 //scan_GsGr(t, reachable);
	 scan_GsGr(h, reachable);
	 // update notGu nodes info.
	 for (it = reachable.begin(); it != reachable.end(); it++) {
		 is_in_unorder[*it] = true;
		 n2v[*it] = max_ccs[0];
		 v2ns[0][max_ccs[0]] = *it;
		 max_ccs[0] += 1;
		 visited[*it] = false;
		 //  if node in Gsi, then set the i2n[*it] of Gsi as  -1
	 }
	 myres.tc_comp_time += mt.elapsed();
	 */
 }

 template<class T>
 void dynamic_tc_mnr_scc<T>::add_edge_alg2(vd_t t, vd_t h)
 {
	 my_timer mt; 
	 unsigned int tn2v(n2v[t]);
	 std::vector<unsigned int> reachable2;
	 std::vector<unsigned int>::iterator it, it_ds;
	 scan_GsGr(h, reachable2 );
	 
	 // update notGu nodes info.
	 for (it = reachable2.begin(); it != reachable2.end(); it++) {
		 is_in_unorder[*it] = true;
		 n2v[*it] = max_ccs[0];
		 v2ns[0][max_ccs[0]] = *it;
		 max_ccs[0] += 1;
		 visited[*it] = false;
	 }
	 unsigned int hn2v(n2v[h]);
	 myres.tc_comp_time += mt.elapsed();
	 if (tn2v < hn2v) {
		 mt = my_timer();
		 std::vector<unsigned int> reaching, reachable;
		 bool has_cycle = false;
		 in_component[t] = true;
		 mnr_dfs_Gm(h, tn2v, hn2v, reachable, has_cycle);
		 myres.tc_comp_time += mt.elapsed();
		 if (has_cycle) {
			 // identify the scc and condense scc. 
			 mt = my_timer();
			 std::vector<unsigned int> cycle_nodes;
			 std::vector<unsigned int>::iterator it;
			 for (it = reachable.begin(); it != reachable.end(); it++) {
				 if (in_component[*it]) {
					 cycle_nodes.push_back(*it);
				 }
			 }
			 assert(cycle_nodes.size() > 1);
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
			 myres.condensation_time += mt.elapsed();
			 mt = my_timer();

			 mnr_shift_cycle_Gm(tn2v, hn2v );
			 myres.tc_comp_time += mt.elapsed();
		 }
		 else {
			 mt = my_timer();
			 mnr_shift_Gm(tn2v, hn2v);
			 myres.tc_comp_time += mt.elapsed();
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
 void dynamic_tc_mnr_scc<T>::add_edge_alg3(vd_t t, vd_t h)
 {
	 my_timer mt; 
	 std::vector<vd_t> v2n_y;
	 if (n2y[t] == n2y[h]) {
		 if (is_in_same_year[t] != 0 && is_in_same_year[h] == 0) { // t in Gsi, h in Gr
			 mt = my_timer();
			 is_in_same_year[h] = n2y[h];
			 int yidx = n2y[h] - min_year + 1;
			 n2v[h] = max_ccs[yidx]; // v2n_y[ g.max_ccs[yidx] ] = n2v[h]; due to current max_ccs[yidx] == v2n_y.size()
			 v2n_y = v2ns[yidx];
			 assert(max_ccs[yidx] == v2n_y.size());
			 // v2n_y.push_back( h ); 
			 v2ns[yidx].push_back(h);
			 max_ccs[yidx] += 1;

			 add_edge_Gsi(t, h);
		 }
		 else if (is_in_same_year[t] == 0 && is_in_same_year[h] != 0) { // t in Gr, h in Gsi
		 }
		 else if (is_in_same_year[t] == 0 && is_in_same_year[h] == 0) {// t in Gr, h in Gr
			 mt = my_timer();
			 is_in_same_year[h] = n2y[h];
			 int yidx = n2y[h] - min_year + 1;
			 n2v[h] = max_ccs[yidx];
			 v2n_y = v2ns[yidx];
			 assert(max_ccs[yidx] == v2n_y.size());
			 v2ns[yidx].push_back(h);
			 max_ccs[yidx] += 1;
			 myres.tc_comp_time += mt.elapsed();
		 }
		 else if (is_in_same_year[t] != 0 && is_in_same_year[h] != 0) {// t, h in Gsi, Gsi
			 // apply dynamic algorithm in Gsi. 
			 add_edge_Gsi(t, h );
		 }
		 else {
			 assert(0 == 1 && "Type Error. t, h NOT in Gsi OR Gr.");
		 }

	 }
	 else if (n2y[t] < n2y[h]) { // e.g., t in 2017, h in 2010. The cycle must be in Gsi & Gr. 
		 mt = my_timer();
		 std::vector<unsigned int> reachable, cycle;
		 std::vector<unsigned int>::iterator it, it_ds, i;
		 bool flag = false; unsigned int in_scc_order = 0;
		 //in_component[t] = true;
		 //scan_GsGr_without_topo(t, t, reachable, flag);
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
					 v2ns[0][max_ccs[0]] = *it;
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
			 myres.condensation_time += mt.elapsed();

			 is_in_unorder[s] = true;
			 n2v[s] = max_ccs[0];
			 v2ns[0][max_ccs[0]] = s;
			 max_ccs[0] += 1;

			 in_component[t] = false;
		 }

		 else {
			 mt = my_timer();
			 // update notGu nodes info.  
			 for (it = reachable.begin(); it != reachable.end(); it++) {
				 is_in_unorder[*it] = true;
				 n2v[*it] = max_ccs[0];
				 v2ns[0][max_ccs[0]] = *it;
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
 void dynamic_tc_mnr_scc<T>::add_edge_Gsi(vd_t t, vd_t h)
 {
	 unsigned int tn2v(n2v[t]); // lower bound
	 unsigned int hn2v(n2v[h]); // upper bound
	 bool has_cycle = false;
	 std::vector<unsigned int> reachable;
	 my_timer mt; 
	 if (tn2v < hn2v) {
		 in_component[t] = true;
		 mnr_dfs_Gsi(h, tn2v, hn2v, reachable, has_cycle, n2y[t]);

		 myres.invalid_edges += 1;
		 myres.aff_region += (hn2v - tn2v);
		 myres.search_time += mt.elapsed();

		 if (has_cycle) {
			 mt = my_timer();
			 // identify the scc and condense scc. 
			 std::vector<unsigned int> cycle_nodes;
			 std::vector<unsigned int>::iterator it;
			 for (it = reachable.begin(); it != reachable.end(); it++) {
				 if (in_component[*it]) {
					 cycle_nodes.push_back(*it);
				 }
			 }
			 assert(cycle_nodes.size() > 1);
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
			 myres.condensation_time += mt.elapsed();
			 mt = my_timer();

			 mnr_shift_cycle_Gsi(tn2v, hn2v, n2y[t]);
			 myres.reassignment_time += mt.elapsed();

		 }
		 else {
			 mt = my_timer();
			 mnr_shift_Gsi(tn2v, hn2v, n2y[t]);
			 myres.reassignment_time += mt.elapsed();
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
 void dynamic_tc_mnr_scc<T>::scan_GsGr(vd_t n, std::vector<unsigned int>& reachable)
 {
	 out_iterator i, iend;
	 visited[n] = true;
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

 template<class T>
 void dynamic_tc_mnr_scc<T>::scan_GsGr_without_topo(vd_t n, vd_t lb, std::vector<unsigned int>& reachable, bool & flag)
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
			 else if (!visited[w] && is_in_unorder[w] == false) {					// w not visited & w is not in Gm. Because w in Gm, the w reachable nodes also in Gm. 
				 scan_GsGr_without_topo(w, lb, reachable, flag);					// Then the nodes that need to be reordered will not be traversed. 
			 }
			 in_component[n] = in_component[n] || in_component[w];
		 }
	 }
	 reachable.push_back(n);

 }

 template<class T>
 void dynamic_tc_mnr_scc<T>::mnr_dfs_Gm(vd_t h, unsigned int lb, unsigned int ub, std::vector<unsigned int>& reachable, bool & flag)
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
			 if (wn2v == lb && ! visited[w] && is_in_unorder[w] == true) {
				 flag = true;
				 in_component[w] = true;
				 mnr_dfs_Gm(w, lb, ub, reachable, flag);
			 }
			 else if (wn2v > lb && !visited[w] && is_in_unorder[w] == true) {
				 mnr_dfs_Gm(w, lb, ub, reachable, flag);
			 }
			 in_component[h] = in_component[h] || in_component[w];
		 }
	 }

 }

 template<class T>
 void dynamic_tc_mnr_scc<T>::mnr_dfs_Gsi(vd_t h, unsigned int lb, unsigned int ub, std::vector<unsigned int>& reachable, bool & flag, unsigned int year)
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
			 if (wn2v == lb && !visited[w] && is_in_unorder[w] == false && is_in_same_year[w] == year && n2y[w] == year) {
				 flag = true;
				 in_component[w] = true;
				 mnr_dfs_Gsi(w, lb, ub, reachable, flag, year);
			 }
			 else if (wn2v > lb && ! visited[w] && is_in_unorder[w] == false &&  is_in_same_year[w] == year && n2y[w] == year) {
				 mnr_dfs_Gsi(w, lb, ub, reachable, flag, year);
			 }
			 in_component[h] = in_component[h] || in_component[w];
		 }
	 }
 }

 template<class T>
 void dynamic_tc_mnr_scc<T>::mnr_shift_Gm(unsigned int lb, unsigned int ub)
 {
	 std::vector<vd_t> & v2n_u = v2ns[0];
	 std::vector<vd_t> dfs_node, normal_node;
	 std::vector<unsigned int> aff_orders;
	 unsigned int i;
	 for (i = ub; i >= lb; i--) { // unsigned int -1 is 2^32 -1 , i.e.,  4294967295
		 vd_t w(v2n_u[i]);
		 if (w == -1) {
			 if (i == 0) {
				 break;
			 }
			 else {
				 continue;
			 }
		 }
		 if (! is_in_unorder[w]) {
			 continue;
		 }
		 aff_orders.push_back(i);
		 in_component[w] = false;
		 if (visited[w]) {
			 visited[w] = false;
			 dfs_node.push_back(w);
		 }
		 else {
			 normal_node.push_back(w);
		 }
		 if (i == 0) break;	 // unsigned int -1 is 2^32 -1 , i.e.,  4294967295
	 }
	 std::vector<vd_t>::iterator it1(normal_node.begin());
	 std::vector<unsigned int>::iterator it2(aff_orders.begin());
	 assert((dfs_node.size() + normal_node.size()) <= aff_orders.size());
	 for (; it1 != normal_node.end(); it1++) {
		 n2v[*it1] = *it2;
		 v2n_u[*it2] = *it1;
		 it2++;
	 }
	 for (it1 = dfs_node.begin(); it1 != dfs_node.end(); it1++) {
		 n2v[*it1] = *it2;
		 v2n_u[*it2] = *it1;
		 it2++;
	 }
	 while (it2 != aff_orders.end()) {
		 v2n_u[*it2] = -1;
		 it2++;
	 }

 }

 template<class T>
 void dynamic_tc_mnr_scc<T>::mnr_shift_Gsi(unsigned int lb, unsigned int ub, unsigned int year)
 {
	 std::vector<vd_t> & v2n_y = v2ns[year - min_year + 1];
	 std::vector<vd_t> dfs_node, normal_node;
	 std::vector<unsigned int> aff_orders;
	 unsigned int i;
	 for (i = ub; i >= lb; i--) {
		 vd_t w(v2n_y[i]);
		 if (w == -1) {
			 if (i == 0) {
				 break;
			 }
			 else {
				 continue;
			 }
		 }
		 if (is_in_unorder[w]) {
			 continue;
		 }
		 aff_orders.push_back(i);
		 in_component[w] = false;
		 if (visited[w]) {
			 visited[w] = false;
			 dfs_node.push_back(w);
		 }
		 else {
			 normal_node.push_back(w);
		 }
		 if (i == 0) break;
	 }
	 std::vector<vd_t>::iterator it1(normal_node.begin());
	 std::vector<unsigned int>::iterator it2(aff_orders.begin());
	 assert((dfs_node.size() + normal_node.size()) <= aff_orders.size());
	 for (; it1 != normal_node.end(); it1++) {
		 n2v[*it1] = *it2;
		 v2n_y[*it2] = *it1;
		 it2++;
	 }
	 for (it1 = dfs_node.begin(); it1 != dfs_node.end(); it1++) {
		 n2v[*it1] = *it2;
		 v2n_y[*it2] = *it1;
		 it2++;
	 }
	 while (it2 != aff_orders.end()) {
		 v2n_y[*it2] = -1;
		 it2++;
	 }
 }

 template<class T>
 void dynamic_tc_mnr_scc<T>::mnr_shift_cycle_Gm(unsigned int lb, unsigned int ub)
 {
	 std::vector<vd_t> dfs_node, normal_node;
	 std::vector<unsigned int> aff_orders;
	 std::vector<vd_t> & v2n_u = v2ns[0];
	 unsigned int i;
	 bool is_first = true;
	 for (i = ub; i >= lb; i--) {
		 vd_t w(v2n_u[i]);
		 if (w == -1) {
			 if (i == 0) {
				 break;
			 }
			 else {
				 continue;
			 }
		 }
		 if (!is_in_unorder[w]) {
			 continue;
		 }
		 aff_orders.push_back(i);
		 if (in_component[w] && is_first) {
			 dfs_node.push_back(ds.find(w));
			 is_first = false;
			 in_component[w] = false;
			 v2n_u[i] = -1;
		 }
		 else if (in_component[w] && !is_first) {
			 in_component[w] = false;
			 v2n_u[i] = -1;
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
		 if (i == 0) break;
	 }
	 std::vector<vd_t>::iterator it1(normal_node.begin());
	 std::vector<unsigned int>::iterator it2(aff_orders.begin());
	 assert((dfs_node.size() + normal_node.size()) <= aff_orders.size());
	 for (; it1 != normal_node.end(); it1++) {
		 n2v[*it1] = *it2;
		 v2n_u[*it2] = *it1;
		 it2++;
	 }
	 for (it1 = dfs_node.begin(); it1 != dfs_node.end(); it1++) {
		 n2v[*it1] = *it2;
		 v2n_u[*it2] = *it1;
		 it2++;
	 }
	 while (it2 != aff_orders.end()) {
		 v2n_u[*it2] = -1;
		 it2++;
	 }

 }

 template<class T>
 void dynamic_tc_mnr_scc<T>::mnr_shift_cycle_Gsi(unsigned int lb, unsigned int ub, unsigned int year)
 {
	 std::vector<vd_t> dfs_node, normal_node;
	 std::vector<unsigned int> aff_orders;
	 std::vector<vd_t> & v2n_y = v2ns[year - min_year + 1];
	 unsigned int i;
	 bool is_first = true;
	 for (i = ub; i >= lb; i--) {
		 vd_t w(v2n_y[i]);
		 if (w == -1) {
			 if (i == 0) {
				 break;
			 }
			 else {
				 continue;
			 }
		 }
		 if (is_in_unorder[w]) {
			 continue;
		 }
		 aff_orders.push_back(i);
		 if (in_component[w] && is_first) {
			 dfs_node.push_back(ds.find(w));
			 is_first = false;
			 in_component[w] = false;
			 v2n_y[i] = -1;
		 }
		 else if (in_component[w] && !is_first) {
			 in_component[w] = false;
			 v2n_y[i] = -1;
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
		 if (i == 0) break;
	 }
	 std::vector<vd_t>::iterator it1(normal_node.begin());
	 std::vector<unsigned int>::iterator it2(aff_orders.begin());
	 assert((dfs_node.size() + normal_node.size()) <= aff_orders.size());
	 for (; it1 != normal_node.end(); it1++) {
		 n2v[*it1] = *it2;
		 v2n_y[*it2] = *it1;
		 it2++;
	 }
	 for (it1 = dfs_node.begin(); it1 != dfs_node.end(); it1++) {
		 n2v[*it1] = *it2;
		 v2n_y[*it2] = *it1;
		 it2++;
	 }
	 while (it2 != aff_orders.end()) {
		 v2n_y[*it2] = -1;
		 it2++;
	 }

 }
