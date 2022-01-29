// This is the implementation of the algorithm of PK2. It was first described in the  following:
//
// [1] Pearce, David J., and Paul HJ Kelly. "A batch algorithm for maintaining a topological order." 
//     Proceedings of the Thirty-Third Australasian Conferenc on Computer Science-Volume 102. 2010. 79-88.
// 
// They claimed it is (PK2) the first batch version for maintaining online topological order.
// PK2 introduces a shift queue to handle batch inc-edges based on MNR (1996), which fails to 
// describe how to do cycle detection. Thus, We expanded PK2 to support batch-incremental 
// citation graph cycle detetion. 
// 

#ifndef BATCH_PK2_SCC_H_
#define BATCH_PK2_SCC_H_

#include "boost/property_map/property_map.hpp"
#include "disjoint_set.h"
#include "my_func.h"
#include "scc_para.h"

class my_tuple {
public:
	unsigned int first;
	unsigned int second;
	unsigned int third;
public:
	my_tuple(unsigned int _f, unsigned int _s, unsigned int _t): first(_f), second(_s), third(_t) {
	}
};


template<class T>
class batch_PK2_scc {
private:
	typedef CG::vertex_descriptor vd_t;
	typedef CG::out_edge_iterator out_iterator;
	typedef CG::in_edge_iterator in_iterator;
public:
	CG & cg;														// citation graph.
	std::vector< T > & n2v;											// topo value of each node. 
	std::vector<vd_t> v2n;
	std::vector<bool> visited, in_reachables;

	disjoint_set &ds;
	std::vector<pair<unsigned int, unsigned int>> & inc_edges;
	std::string self_name;
	unsigned int np;
	unsigned int bs;												// batch size.
	std::vector<std::vector<unsigned int>> topo_scc;
	scc_para para; 
	

public:
	batch_PK2_scc();
	batch_PK2_scc(CG &_cg, disjoint_set &_ds, std::vector<T> &_n2v, std::vector<pair<unsigned int, unsigned int>> & _inc_edges, unsigned int _np, unsigned int _bs);
	void add_edges();
	void add_edges(int start_index, int end_index);
	void eval_single(int start_index, int end_index);
	void add_edge_bat(std::vector<my_tuple> &batch);
	void dfs_f(vd_t n, unsigned int lb, vector<pair<unsigned int, unsigned int> > &reachables);
	void shift(unsigned int max_tv, std::vector<pair<unsigned int, unsigned int> > reachables);
	void simple_sccs(std::vector<unsigned int> &cycle_nodes);


};



#endif

