#ifndef PARA_HPP
#define PARA_HPP
#include <iostream>
#include <vector>
#include <stack>

// tarjan.
class scc_para {
public:
	std::vector <unsigned int> num, lowlink, cc, stack;
	std::vector<bool> in_stack;
	unsigned int n;
	unsigned int time;
	unsigned int cc_idx;

	scc_para(int np) : n(np), time(1), cc_idx(1) {
		num = std::vector<unsigned int>(np, 0);
		lowlink = std::vector<unsigned int>(np, 0);
		cc = std::vector<unsigned int>(np, 0);
		stack = std::vector<unsigned int>(0, 0);
		stack.reserve(np);
		in_stack = std::vector<bool>(np, false);
	}

	~scc_para() {

	}
};

// pearce2. 
class ppara {
public:
	std::vector<unsigned int> rindex;
	std::vector<unsigned int>S;
	int c;
	int index_p;
	int np;

	ppara(int _np) :np(_np), index_p(1), c(_np - 1) {
		rindex = std::vector<unsigned int>(np, 0);
		S.clear();
		S.reserve(np);
	}

};

// gabow 
class gabow_para {
public:
	std::vector <unsigned int> num;
	std::vector<bool> inscc;
	std::stack<unsigned int> S, P;
	int np;
	int cc_idx;
	gabow_para(int _np) :np(_np), num(_np, 0), inscc(_np, false) {
		cc_idx = 1;
	}
};

// kosaraju. 
class kosaraju_para {
public:
	std::vector<bool> marked;
	std::vector<unsigned int> nodes;
	int np;
	kosaraju_para(unsigned int _np) : np(_np), marked(_np, false) {
		nodes.reserve(np);
	}

	void reset_marked() {
		marked = std::vector<bool>(np, false);
	}
};

#endif
