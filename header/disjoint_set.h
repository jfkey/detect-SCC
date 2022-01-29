#ifndef DISJOINT_SET_H_
#define DISJOINT_SET_H_

#include <iostream>
#include <vector>

class disjoint_set {
private:
	unsigned int * id;
	unsigned int cnt;

public:
	//  the inside vector reps. the strongly connected component. 
	//  the outsize vector reps. the node index of the Graph.
	// the node index of the Graph starts from 0. 
	std::vector<std::vector<unsigned int>> scc_nodes;

private:
	disjoint_set() {}
public:
	disjoint_set(const unsigned int _cnt) : cnt(_cnt) {
		id = new unsigned int[_cnt];
		for (unsigned int i = 0; i < cnt; i++) {
			id[i] = i;
		}
		for (unsigned int i = 0; i < cnt; i++) {
			std::vector<unsigned int> tmp(1, i);
			scc_nodes.push_back(tmp);
		}
	}

	~disjoint_set() {
		delete[]id;
		//std::cout << "~disjoint_set with cnt " << cnt << std::endl; 
	}

public:
	unsigned int find(unsigned int x) {
		unsigned int r = x;
		while (id[r] != r)
			r = id[r];
		unsigned int i = x, j;
		while (i != r) {			//path compression
			j = id[i];
			id[i] = r;
			i = j;
		}
		return r;
	}

	void join(unsigned int x, unsigned int y) {
		if (x == y) return;
		unsigned int fx = find(x), fy = find(y);
		if (fx != fy)
			id[fx] = fy;
	}

	// select x as canonical node in the scc. 
	void collapse_scc(unsigned int x, std::vector<unsigned int> scc) {
		std::vector<unsigned int>::iterator it;
		for (it = scc.begin(); it != scc.end(); it++) {
			if (x != *it) {
				if (scc_nodes[*it].size() == 0) {
					std::cout << "error collapse scc. " << std::endl;
				}
				//assert(scc_nodes[*it].size() == 0);
				scc_nodes[x].insert(scc_nodes[x].end(), scc_nodes[*it].begin(), scc_nodes[*it].end());
				scc_nodes[*it].resize(0);
			}

		}

	}
};

#endif // !DISJOINT_SET
