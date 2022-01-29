#include <iostream> 
#include <string>
#include <time.h> 
#include <algorithm>
#include <assert.h>
#include <random>
#include <exception>

#include "boost/graph/adjacency_list.hpp"
#include "boost/program_options.hpp"
#include "boost/graph/graphviz.hpp"

#include "my_graph.h"
#include "my_func.h"
#include "utils.h"
#include "disjoint_set.h"
#include "memory.h"

#include "static_scc_Pearce2.h"
#include "static_scc_Pearce2_ol.h"
#include "static_scc_Pearce2_tc.h"
#include "static_scc_Pearce2_tc_ol.h"
#include "static_scc_Tarjan.h"
#include "static_scc_Tarjan_tc.h"
#include "static_scc_Gabow.h"
#include "static_scc_Gabow_tc.h"
#include "static_scc_Kosaraju.h"
#include "static_scc_Kosaraju_tc.h"

#include "dynamic_AHRSZ_scc.cpp"
#include "dynamic_sinDSCC.cpp"
#include "dynamic_HKMST_scc.cpp"
#include "dynamic_HKMST_TC_scc.cpp"
#include "dynamic_IGC_scc.cpp"
#include "dynamic_IGC_TC_scc.cpp"
#include "dynamic_MNR_scc.cpp"
#include "dynamic_MNR_TC_scc.cpp"

#include "bat_PK2_scc.cpp"
#include "bat_AHRSZ_scc.cpp"
#include "bat_AHRSZ_P_scc.cpp"
#include "bat_DSCC.cpp"
#include "bat_DSCC_test.cpp"


#define HKMST_GENERATE_STATS
#define IGC_GENERATE_STATS

#define NODE_NAME "indpyear.txt"
#define EDGE_NAME "indpp_refine.txt"


#define NEW_STACK_SIZE 5120			//MB. 
#define isShown true 
my_result myres(isShown);
#define EVAL_SINGLE_CNT 15000

std::vector<double> inc_percentage_vec = { 0.5, 1, 5, 10, 15, 20, 25, 30};



std::vector<bool> cond_scc_in;
std::vector<bool> cond_scc_out;

using namespace std;
using namespace boost;
namespace po = boost::program_options;

namespace {
	const size_t ERROR_IN_COMMAND_LINE = 1;
	const size_t SUCCESS = 0;
}

typedef static_pearce2<unsigned int> static_scc_pearce2;									// 2016 IPL PEARCE2.
typedef static_pearce2_tc<unsigned int> static_scc_pearce2_tc;
typedef static_pearce2_ol<ordered_slist<void>> static_scc_pearce2_ol;
typedef static_pearce2_tc_ol<ordered_slist<void> > static_scc_pearce2_tc_ol;

typedef static_tarjan<unsigned int> static_scc_tarjan;										// 1972 Tarjan.	
typedef static_tarjan_tc<unsigned int> static_scc_tarjan_tc;								// TC-Tarjan
typedef static_gabow<unsigned int> static_scc_gabow;										// Gabow
typedef static_gabow_tc<unsigned int> static_scc_gabow_tc;									// TC-Gabow
typedef static_kosaraju<unsigned int> static_scc_kosaraju;									// Kosaraju.
typedef static_kosaraju_tc<unsigned int> static_scc_kosaraju_tc;							// TC-Kosaraju.

typedef dynamic_ahrsz_scc<ordered_slist<void>> dn_ahrsz1_scc;								// 1990-SODA-AHRSZ
typedef sinDSCC<ordered_slist<void>> sin_dscc;												// our single 
typedef dynamic_hkmst_scc<ordered_slist<void>> dn_hkmst_scc;								// 2012-TALG-HKMST. 
typedef dynamic_tc_hkmst_scc<ordered_slist<void>> dn_tc_hkmst_scc;							// Time-Coupled + HKMST
typedef dynamic_igc_scc<unsigned int> dn_igc_scc;											// 2017-SIGMOD, 2006-PK
typedef dynamic_tc_igc_scc<unsigned int> dn_tc_igc_scc;										// Time-Coupled + PK
typedef dynamic_mnr_scc<unsigned int> dn_mnr_scc;											// MNR. 
typedef dynamic_tc_mnr_scc<unsigned int> dn_tc_mnr_scc;										// TC-MNR .
typedef batch_PK2_scc<unsigned int> bat_pk2_scc;											// 2010-PK2 

typedef batch_AHRSZ_scc<ordered_slist<void>> bat_ahrsz_scc;									// AHRSZ. i.e, one of the baselines 
typedef batch_AHRSZ_P_scc<ordered_slist<void>> bat_ahrszp_scc;								// bat-AHRSZ 
typedef batDSCC<ordered_slist<void>> bat_dscc;												//  our batch incremental 
typedef batDSCC_test<ordered_slist<void>> bat_dscc_test;							

Config parseParams(int argc, char** argv) {
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce help message")
		("data-folder,f", po::value<string>()->required(), "graph data folder")
		("algo-type,t", po::value<string>()->required(), "algorithm type")
		("algo-name,a", po::value<string>()->required(), "algorithm name")
		("graph-node,n", po::value<int>()->required(), "citation graph node")
		("insert-type,p", po::value<string>()->default_value(""), "insert edge type")
		("insert-file,l", po::value<string>()->default_value(""), "insert edge type")
		("node-name,d", po::value<string>()->default_value(""), "node file name")
		("edge-name,g", po::value<string>()->default_value(""), "edge file name")
		("batch-size,b", po::value<int>()->default_value(1), "batch size for batch algorithm")
		("inc-percentage,c", po::value<double>()->default_value(0), "increase the percentage of edges")
		("simplify-dataset,s", po::value<bool>()->default_value(false), "remove the edge: n2y[s] < n2y[t]")
		("scale-factor,r", po::value<double>()->default_value(1), "scale factor for static SCCs detection")
		;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm); // can throw 
	po::notify(vm);
	Config config;
	if (vm.count("help")) {
		cout << desc << '\n';
		exit(0);
	}
	if (vm.count("data-folder")) {
		config.folder = vm["data-folder"].as<string>();
	}
	if (vm.count("algo-type")) {
		config.alg_type = vm["algo-type"].as<string>();
	}
	if (vm.count("algo-name")) {
		config.alg_name = vm["algo-name"].as<string>();
	}
	if (vm.count("graph-node")) {
		config.n = vm["graph-node"].as<int>();
	}
	if (vm.count("insert-type")) {
		config.insert_type = vm["insert-type"].as<string>();
	}
	if (vm.count("insert-file")) {
		config.insert_file = vm["insert-file"].as<string>();
	}
	if (vm.count("node-name")) {
		config.node_name = vm["node-name"].as<string>();
	}
	if (vm.count("edge-name")) {
		config.edge_name = vm["edge-name"].as<string>();
	}
	if (vm.count("batch-size")) {
		config.bat_size = vm["batch-size"].as<int>();
	}
	if (vm.count("inc-percentage")) {
		config.inc_percentage = vm["inc-percentage"].as<double>();
	}
	if (vm.count("simplify-dataset")) {
		config.simple_data = vm["simplify-dataset"].as<bool>();
	}
	if (vm.count("inc-percentage")) {
		config.scale_factor = vm["scale-factor"].as<double>();
	}

	return config;
}


template <typename T>
void do_static_cd(T& graph, const std::string dir) {
	
	graph.read_citation_graph();
	graph.gen_citation_graph();
	size_t meg = Memory::toKB(Memory::getCurrentRSS());
	
	graph.topo_sort();
	size_t mescc = Memory::toKB(Memory::getCurrentRSS()) - meg;
	//graph.eliminate_static_scc();
	graph.print_scc(dir);
	
	myres.graph_memory = meg;
	myres.static_memory = mescc;
}


// MNR: dn_mnr_scc, static_scc_pearce2, subgraph1
// PK(IGC): dn_igc_scc, static_scc_pearce2, subgraph1
template <typename T, typename S>
void do_dynamic_cd(S & sg, const std::string& dir, const string & pyear, const string& pref,
	const string &insert_type, string & insert_file, const double percentage, const int B) {
	sg.read_citation_graph();
	sg.gen_citation_graph();
	size_t meg = Memory::toKB(Memory::getCurrentRSS());
	sg.topo_sort();
	size_t mescc = Memory::toKB(Memory::getCurrentRSS()) - meg;
	sg.eliminate_static_scc();
	sg.print_scc(dir);

	std::vector<std::pair<unsigned int, unsigned int>> edges;
	//iread_insert_edge(dir + get_inc_f(insert_file, percentage, insert_type), edges, sg.n2y);
	iread_insert_edge(dir + get_inc_f(insert_file, insert_type), edges, sg.n2y);
	int start_idx = 0, end_idx = 0; 
	T dynamic_graph(sg.cg, sg.ds, sg.n2v, edges, sg.np);
	size_t med = Memory::toKB(Memory::getCurrentRSS());

	//dynamic_graph.add_edges();
	if (percentage < 0) {
		for (int i = 0; i < inc_percentage_vec.size(); i ++) {
			end_idx = (int)(edges.size() * inc_percentage_vec[i] / DELTAG_SIZE);
			dynamic_graph.add_edges(start_idx, end_idx);
			start_idx = end_idx;
		}
	}
	else if (percentage == 0 ) {
		std::vector<std::pair<unsigned int, unsigned int>> single_cand;
		int step = (edges.size() - 1) / EVAL_SINGLE_CNT;
		for (int i = 0; i < EVAL_SINGLE_CNT; i++) {
			single_cand.push_back(edges[step * i]);
		}
		
		dynamic_graph.inc_edges = single_cand;
		dynamic_graph.add_edges(0, EVAL_SINGLE_CNT);
	}
	else {
		end_idx = (int)(edges.size() * percentage / DELTAG_SIZE);
		dynamic_graph.add_edges(start_idx, end_idx);
	}
	
	size_t meinc = Memory::toKB(Memory::getCurrentRSS()) - med;

	myres.gen_dynamic_scc(dynamic_graph.ds.scc_nodes);
	/*print_scc_infile(dir, dynamic_graph.ds.scc_nodes, dynamic_graph.self_name, edges.size());*/
	myres.dynamic_info(edges.size(), B );
	myres.graph_memory = meg; 
	myres.static_memory = mescc; 
	myres.dynamic_memory = meinc;
	myres.inc_fn = get_inc_f(insert_file, insert_type);

}

// AHRSZ: dn_ahrsz1_scc, static_scc_pearce2_ol 
// HKMST: dn_hkmst_scc, static_scc_pearce2_ol, subgraph2
template <typename T, typename S >
void do_dynamic_cd2(S & sg, const std::string& dir, const string & pyear, const string& pref,
	const string &insert_type, string & insert_file, const double percentage, const int B) {						// sg: stands for static graph.
	sg.read_citation_graph();
	sg.gen_citation_graph();
	size_t meg = Memory::toKB(Memory::getCurrentRSS());
	sg.topo_sort();
	size_t mescc = Memory::toKB(Memory::getCurrentRSS()) - meg;
	sg.eliminate_static_scc();
	sg.print_scc(dir);
	std::vector<std::pair<unsigned int, unsigned int>> edges;
	iread_insert_edge(dir + get_inc_f(insert_file, insert_type), edges, sg.n2y);
	int start_idx = 0, end_idx = 0;

	T dynamic_graph(sg.cg, sg.n2v, sg.p_space, sg.ds, edges, sg.np);
	size_t med = Memory::toKB(Memory::getCurrentRSS());
	//dynamic_graph.add_edges();
	if (percentage < 0) {
		for (int i = 0; i < inc_percentage_vec.size(); i++) {
			end_idx = (int)(edges.size() * inc_percentage_vec[i] / DELTAG_SIZE);
			dynamic_graph.add_edges(start_idx, end_idx);
			start_idx = end_idx;
		}
	}
	else if (percentage == 0) {
		std::vector<std::pair<unsigned int, unsigned int>> single_cand;
		int step = (edges.size() - 1) / EVAL_SINGLE_CNT;
		for (int i = 0; i < EVAL_SINGLE_CNT; i++) {
			single_cand.push_back(edges[step * i]);
		}
		
		dynamic_graph.inc_edges = single_cand;
		dynamic_graph.add_edges(0, EVAL_SINGLE_CNT);
	}
	else {
		end_idx = (int)(edges.size() * percentage / DELTAG_SIZE);
		dynamic_graph.add_edges(start_idx, end_idx);
	}

	size_t meinc = Memory::toKB(Memory::getCurrentRSS()) - med;
	myres.gen_dynamic_scc(dynamic_graph.ds.scc_nodes);
	myres.dynamic_info(edges.size(), B);
	myres.graph_memory = meg;
	myres.static_memory = mescc;
	myres.dynamic_memory = meinc;
	myres.inc_fn = get_inc_f(insert_file,  insert_type);
}

//TC-MNR:
//TC-PK: 
template <typename T, typename S>
void do_dynamic_cd3(S & sg, const std::string& dir, const string & pyear, const string& pref,
	const string &insert_type, string & insert_file, const double percentage, const int B) {
	sg.read_citation_graph();
	sg.gen_citation_graph();
	size_t meg = Memory::toKB(Memory::getCurrentRSS());
	sg.topo_sort();
	size_t mescc = Memory::toKB(Memory::getCurrentRSS()) - meg;
	sg.eliminate_static_scc();
	sg.print_scc(dir);

	std::vector<std::pair<unsigned int, unsigned int>> edges;
	iread_insert_edge(dir + get_inc_f(insert_file, insert_type), edges, sg.n2y);
	int start_idx = 0, end_idx = 0;

	T dynamic_graph(sg.cg, sg.ds, sg.n2v, sg.max_ccs, sg.is_in_unorder, sg.is_in_same_year, sg.n2y, sg.max_year, sg.min_year, edges, sg.np);
	size_t med = Memory::toKB(Memory::getCurrentRSS());

	if (percentage < 0) {
		for (int i = 0; i < inc_percentage_vec.size(); i++) {
			end_idx = (int)(edges.size() * inc_percentage_vec[i] / DELTAG_SIZE);
			dynamic_graph.add_edges(start_idx, end_idx);
			start_idx = end_idx;
		}
	}
	else if (percentage == 0) {
		std::vector<std::pair<unsigned int, unsigned int>> single_cand;
		int step = (edges.size() - 1) / EVAL_SINGLE_CNT;
		for (int i = 0; i < EVAL_SINGLE_CNT; i++) {
			single_cand.push_back(edges[step * i]);
		}
		
		dynamic_graph.inc_edges = single_cand;
		dynamic_graph.add_edges(0, EVAL_SINGLE_CNT);
	}
	else {
		end_idx = (int)(edges.size() * percentage / DELTAG_SIZE);
		dynamic_graph.add_edges(start_idx, end_idx);
	}

	size_t meinc = Memory::toKB(Memory::getCurrentRSS()) - med;
	myres.gen_dynamic_scc(dynamic_graph.ds.scc_nodes);
	myres.dynamic_info(edges.size(), B );
	myres.graph_memory = meg;
	myres.static_memory = mescc;
	myres.dynamic_memory = meinc;
	myres.inc_fn = get_inc_f(insert_file, insert_type);
}

// TC-HKMST.
// TC-AHRSZ
template <typename T, typename S>
void do_dynamic_cd4(S & sg, const std::string& dir, const string & pyear, const string& pref,
	const string &insert_type, string & insert_file, const double percentage, const int B) {
	sg.read_citation_graph();
	sg.gen_citation_graph();
	size_t meg = Memory::toKB(Memory::getCurrentRSS());
	sg.topo_sort();
	size_t mescc = Memory::toKB(Memory::getCurrentRSS()) - meg;
	sg.eliminate_static_scc();
	sg.print_scc(dir);

	std::vector<std::pair<unsigned int, unsigned int>> edges;
	iread_insert_edge(dir + get_inc_f(insert_file, insert_type), edges, sg.n2y);
	int start_idx = 0, end_idx = 0;

	T dynamic_graph(sg.cg, sg.n2v, sg.pspaces, sg.ds, sg.n2y, sg.is_in_same_year, sg.is_in_unorder, sg.min_year, sg.max_year, edges, sg.np);
	size_t med = Memory::toKB(Memory::getCurrentRSS());

	if (percentage < 0) {
		for (int i = 0; i < inc_percentage_vec.size(); i++) {
			end_idx = (int)(edges.size() * inc_percentage_vec[i] / DELTAG_SIZE);
			dynamic_graph.add_edges(start_idx, end_idx);
			start_idx = end_idx;
		}
	}
	else if (percentage == 0) {
		std::vector<std::pair<unsigned int, unsigned int>> single_cand;
		int step = (edges.size() - 1) / EVAL_SINGLE_CNT;
		for (int i = 0; i < EVAL_SINGLE_CNT; i++) {
			single_cand.push_back(edges[step * i]);
		}
		
		dynamic_graph.inc_edges = single_cand;
		dynamic_graph.add_edges(0, EVAL_SINGLE_CNT);
	}
	else {
		end_idx = (int)(edges.size() * percentage / DELTAG_SIZE);
		dynamic_graph.add_edges(start_idx, end_idx);
	}

	size_t meinc = Memory::toKB(Memory::getCurrentRSS()) - med;
	myres.gen_dynamic_scc(dynamic_graph.ds.scc_nodes);
	myres.dynamic_info(edges.size(), B);
	myres.graph_memory = meg;
	myres.static_memory = mescc;
	myres.dynamic_memory = meinc;
	myres.inc_fn = get_inc_f(insert_file, insert_type);
}



template <typename T>
void do_dynamic_baseline(T & sg, const std::string& dir, const string & pyear, const string& pref,
	const string &insert_type, string & insert_file, const double percentage, const int B) {
	std::vector<std::pair<unsigned int, unsigned int>> edges;
	sg.read_citation_graph();
	sg.gen_citation_graph();
	iread_insert_edge(dir + get_inc_f(insert_file, insert_type), edges, sg.n2y);
	size_t meg = Memory::toKB(Memory::getCurrentRSS());
	int start_idx = 0, end_idx = 0;

	if (percentage < 0) {
		for (int i = 0; i < inc_percentage_vec.size(); i++) {
			end_idx = (int)(edges.size() * inc_percentage_vec[i] / DELTAG_SIZE);
			sg.add_inc_edges(edges, start_idx, end_idx);
			sg.topo_sort();
			sg.topo_sort_pearce2();
			record_cur_per_info(myres);
			start_idx = end_idx;
			
		}
	}
	else if (percentage == 0) {
		std::vector<std::pair<unsigned int, unsigned int>> single_cand;
		int step = (edges.size() - 1) / EVAL_SINGLE_CNT;
		for (int i = 0; i < EVAL_SINGLE_CNT; i++) {
			single_cand.push_back(edges[step * i]);
		}
		
		sg.add_inc_edges(single_cand, 0, EVAL_SINGLE_CNT);
		sg.topo_sort();
		sg.topo_sort_pearce2();
		record_cur_per_info(myres);

	}
	else {
		end_idx = (int)(edges.size() * percentage / DELTAG_SIZE);
		sg.add_inc_edges(edges, start_idx, end_idx);
		sg.topo_sort();
		sg.topo_sort_pearce2();
		record_cur_per_info(myres);
	}


	size_t meinc = Memory::toKB(Memory::getCurrentRSS()) -meg;
	sg.print_scc(dir);
	string dnname = "Dyn_Baseline-";
	myres.set_dynamic_name(dnname);
	myres.gen_dynamic_scc(sg.topo_scc);
	myres.dynamic_info(edges.size(), B);
	myres.graph_memory = meg;
	myres.dynamic_memory = meinc;
	myres.inc_fn = get_inc_f(insert_file, insert_type);
}


// bat_TC_ahrsz+
template <typename T, typename S>
void do_batch_tc_ahrsz(S & sg, const std::string& dir, const string & pyear, const string& pref,
	const string &insert_type, string & insert_file, const double percentage, const int B) {
	sg.read_citation_graph();
	sg.gen_citation_graph();
	//std::cout << "gen_citation_graph. Physical Process Memory Cost: " << Memory::toKB(Memory::getCurrentRSS()) << " MB" << std::endl;
	size_t meg = Memory::toKB(Memory::getCurrentRSS());
	sg.topo_sort();
	size_t mescc = Memory::toKB(Memory::getCurrentRSS()) - meg;
	sg.eliminate_static_scc();
	sg.print_scc(dir);
	
	std::vector<std::pair<unsigned int, unsigned int>> edges;
	iread_insert_edge(dir + get_inc_f(insert_file, insert_type), edges, sg.n2y);
	int start_idx = 0, end_idx = 0;
	//std::cout << "iread_insert_edge. Physical Process Memory Cost: " << Memory::toKB(Memory::getCurrentRSS()) << " MB" << std::endl;

	// to debug the results. 
	//std::random_device rd;
	//std::shuffle(edges.begin(), edges.end(), rd);
	//std::vector<std::pair<unsigned int, unsigned int>>::iterator it;
	//long cur = clock();
	//std::ofstream ofs(dir + "/log/rand-" + to_string(cur) + ".txt");
	//for (it = edges.begin(); it != edges.end(); it++) {
	//	ofs << (*it).first << "\t" << (*it).second << std::endl;
	//}
	//ofs.close();

	T dynamic_graph(sg.cg, sg.n2v, sg.pspaces, sg.ds, sg.n2y, sg.is_in_same_year, sg.is_in_unorder, sg.min_year, sg.max_year, edges, sg.np, B);
	size_t med = Memory::toKB(Memory::getCurrentRSS());

	if (percentage < 0) {
		for (int i = 0; i < inc_percentage_vec.size(); i++) {
			end_idx = (int)(edges.size() * inc_percentage_vec[i] / DELTAG_SIZE);
			dynamic_graph.add_edges(start_idx, end_idx);
			start_idx = end_idx;
		}
	}
	else if (percentage == 0) {
		std::vector<std::pair<unsigned int, unsigned int>> single_cand;
		int step = (edges.size() - 1) / EVAL_SINGLE_CNT;
		for (int i = 0; i < EVAL_SINGLE_CNT; i++) {
			single_cand.push_back(edges[step * i]);
		}
		dynamic_graph.bs = 1;
		dynamic_graph.inc_edges = single_cand;
		dynamic_graph.add_edges(0, EVAL_SINGLE_CNT);
	}
	else {
		end_idx = (int)(edges.size() * percentage / DELTAG_SIZE);
		dynamic_graph.add_edges(start_idx, end_idx);
	}

	size_t meinc = Memory::toKB(Memory::getCurrentRSS()) - med;

	myres.gen_dynamic_scc(dynamic_graph.ds.scc_nodes);
	myres.dynamic_info(edges.size(), B);
	myres.graph_memory = meg;
	myres.static_memory = mescc;
	myres.dynamic_memory = meinc;
	myres.inc_fn = get_inc_f(insert_file,  insert_type);
	//print_scc_infile(dir, dynamic_graph.ds.scc_nodes, dynamic_graph.self_name, edges.size());

}

// bat_pk2
template <typename T, typename S>
void do_batch_pk2(S & sg, const std::string& dir, const string & pyear, const string& pref,
	const string &insert_type, string & insert_file, const double percentage, const int B) {
	sg.read_citation_graph();
	sg.gen_citation_graph();
	size_t meg = Memory::toKB(Memory::getCurrentRSS());
	sg.topo_sort();
	size_t mescc = Memory::toKB(Memory::getCurrentRSS()) - meg;
	sg.eliminate_static_scc();
	sg.print_scc(dir);

	std::vector<std::pair<unsigned int, unsigned int>> edges;
	iread_insert_edge(dir + get_inc_f(insert_file, insert_type), edges, sg.n2y);
	int start_idx = 0, end_idx = 0;
	T dynamic_graph(sg.cg, sg.ds, sg.n2v, edges, sg.np, B);
	size_t med = Memory::toKB(Memory::getCurrentRSS());
	
	if (percentage < 0) {
		for (int i = 0; i < inc_percentage_vec.size(); i++) {
			end_idx = (int)(edges.size() * inc_percentage_vec[i] / DELTAG_SIZE);
			dynamic_graph.add_edges(start_idx, end_idx);
			start_idx = end_idx;
		}
	}
	else if (percentage == 0) {
		std::vector<std::pair<unsigned int, unsigned int>> single_cand;
		int step = (edges.size() - 1) / EVAL_SINGLE_CNT; 
		for (int i = 0; i < EVAL_SINGLE_CNT;i ++) {
			single_cand.push_back(edges[step * i]);
		}
		dynamic_graph.bs = 1; 
		dynamic_graph.inc_edges = single_cand;
		dynamic_graph.add_edges(0, EVAL_SINGLE_CNT);

	}
	else {
		end_idx = (int)(edges.size() * percentage / DELTAG_SIZE);
		dynamic_graph.add_edges(start_idx, end_idx);
	}
 

	size_t meinc = Memory::toKB(Memory::getCurrentRSS()) - med;

	myres.gen_dynamic_scc(dynamic_graph.ds.scc_nodes);
	myres.dynamic_info(edges.size(), B);
	myres.graph_memory = meg;
	myres.static_memory = mescc;
	myres.dynamic_memory = meinc;
	myres.inc_fn = get_inc_f(insert_file, insert_type);
}

//bat_ahrsz
template <typename T, typename S >
void do_batch_ahrsz(S & sg, const std::string& dir, const string & pyear, const string& pref,
	const string &insert_type, string & insert_file, const double percentage, const int B) {						// sg: stands for static graph.
	sg.read_citation_graph();
	sg.gen_citation_graph();
	size_t meg = Memory::toKB(Memory::getCurrentRSS());
	sg.topo_sort();
	size_t mescc = Memory::toKB(Memory::getCurrentRSS()) - meg;
	sg.eliminate_static_scc();
	sg.print_scc(dir);
	
	std::vector<std::pair<unsigned int, unsigned int>> edges;
	iread_insert_edge(dir + get_inc_f(insert_file, insert_type), edges, sg.n2y);

	int start_idx = 0, end_idx = 0;
	T dynamic_graph(sg.cg, sg.n2v, sg.p_space, sg.ds, edges, sg.np, B);
	size_t med = Memory::toKB(Memory::getCurrentRSS());
	if (percentage < 0) {
		for (int i = 0; i < inc_percentage_vec.size(); i++) {
			end_idx = (int)(edges.size() * inc_percentage_vec[i] / DELTAG_SIZE);
			dynamic_graph.add_edges(start_idx, end_idx);
			start_idx = end_idx;
		}
	}
	else if (percentage == 0) {
		std::vector<std::pair<unsigned int, unsigned int>> single_cand;
		int step = (edges.size() - 1) / EVAL_SINGLE_CNT;
		for (int i = 0; i < EVAL_SINGLE_CNT; i++) {
			single_cand.push_back(edges[step * i]);
		}
		dynamic_graph.bs = 1;
		dynamic_graph.inc_edges = single_cand;
		dynamic_graph.add_edges(0, EVAL_SINGLE_CNT);
	}
	else {
		end_idx = (int)(edges.size() * percentage / DELTAG_SIZE);
		dynamic_graph.add_edges(start_idx, end_idx);

	}
	size_t meinc = Memory::toKB(Memory::getCurrentRSS()) - med;

	myres.gen_dynamic_scc(dynamic_graph.ds.scc_nodes);
	myres.dynamic_info(edges.size(), B);
	myres.graph_memory = meg;
	myres.static_memory = mescc;
	myres.dynamic_memory = meinc;
	myres.inc_fn = get_inc_f(insert_file, insert_type);
}


void do_cycle_detection(Config &config) {
	cond_scc_in = std::vector<bool>(config.n, false);
	cond_scc_out = std::vector<bool>(config.n, false);
	int np = config.n;
	string dir = config.folder;
	string pyear = NODE_NAME;
	string pref = EDGE_NAME;
	string insert_file = "";

	if (config.node_name.length() != 0) {
		pyear = config.node_name;
	}
	if (config.edge_name.length() != 0) {
		pref = config.edge_name;
	}
	if (config.insert_file.length() != 0) {
		insert_file = config.insert_file;
	}

	if (config.alg_type == TYPE_STATIC) {
		if (config.alg_name == TARJAN) {
			static_scc_tarjan static_graph(np, dir + pyear, dir + pref);
			do_static_cd<static_scc_tarjan>(static_graph, dir);
		}
		else if (config.alg_name == GABOW) {
			static_scc_gabow static_graph(np, dir + pyear, dir + pref);
			do_static_cd<static_scc_gabow>(static_graph, dir);
		}
		else if (config.alg_name == KOSARAJU) {
			static_scc_kosaraju static_graph(np, dir + pyear, dir + pref);
			do_static_cd<static_scc_kosaraju>(static_graph, dir);
		}
		else if (config.alg_name == PEARCE) {
			//static_scc_pearce static_graph(np, dir + pyear);
			//do_static_cd<static_scc_pearce, subgraph1>(static_graph, np, dir + pyear, dir + pref, dir);
		}
		else if (config.alg_name == PEARCE2) {
			static_scc_pearce2 static_graph(np, dir + pyear, dir + pref);
			do_static_cd<static_scc_pearce2>(static_graph, dir);
		}
		else if (config.alg_name == TC_TARJAN) {
			static_scc_tarjan_tc static_graph(np, dir + pyear, dir + pref);
			do_static_cd<static_scc_tarjan_tc>(static_graph, dir);
		}
		else if (config.alg_name == TC_GABOW) {
			static_scc_gabow_tc static_graph(np, dir + pyear, dir + pref);
			do_static_cd<static_scc_gabow_tc>(static_graph, dir);
		}
		else if (config.alg_name == TC_KOSARAJU) {
			static_scc_kosaraju_tc static_graph(np, dir + pyear, dir + pref);
			do_static_cd<static_scc_kosaraju_tc>(static_graph, dir);
		}
		else if (config.alg_name == TC_PEARCE2) {
			static_scc_pearce2_tc static_tc(np, dir + pyear, dir + pref);
			do_static_cd<static_scc_pearce2_tc>(static_tc, dir);
			record_cur_per_info(myres);
		}
		myres.set_file_info(dir, pyear, pref );
		myres.to_file();
	}
	else if (config.alg_type == TYPE_STATIC_DEF) {
		//remove this type. 
	}
	else if (config.alg_type == TYPE_INC) {
		if (config.alg_name == MNR) {
			static_scc_pearce2 static_graph(np, dir + pyear, dir + pref);
			do_dynamic_cd<dn_mnr_scc, static_scc_pearce2>(static_graph, dir, pyear, pref, config.insert_type, insert_file, config.inc_percentage, config.bat_size);
		}
		else if (config.alg_name == PK) { // PK, i.e., IGC
			static_scc_pearce2 static_graph(np, dir + pyear, dir + pref);
			do_dynamic_cd<dn_igc_scc, static_scc_pearce2>(static_graph, dir, pyear, pref, config.insert_type, insert_file, config.inc_percentage, config.bat_size);
		}
		else if (config.alg_name == AHRSZ) {
			static_scc_pearce2_ol static_graph(np, dir + pyear, dir + pref);
			do_dynamic_cd2<dn_ahrsz1_scc, static_scc_pearce2_ol >(static_graph, dir, pyear, pref, config.insert_type, insert_file, config.inc_percentage, config.bat_size);
		}
		else if (config.alg_name == HKMST) {
			static_scc_pearce2_ol static_graph(np, dir + pyear, dir + pref);
			do_dynamic_cd2<dn_hkmst_scc, static_scc_pearce2_ol >(static_graph, dir, pyear, pref, config.insert_type, insert_file, config.inc_percentage, config.bat_size);
		}
		else if (config.alg_name == TC_MNR) {
			static_scc_pearce2_tc static_graph(np, dir + pyear, dir + pref);
			do_dynamic_cd3<dn_tc_mnr_scc, static_scc_pearce2_tc>(static_graph, dir, pyear, pref, config.insert_type, insert_file, config.inc_percentage, config.bat_size);
		}
		else if (config.alg_name == TC_PK) {
			static_scc_pearce2_tc static_graph(np, dir + pyear, dir + pref);
			do_dynamic_cd3<dn_tc_igc_scc, static_scc_pearce2_tc>(static_graph, dir, pyear, pref, config.insert_type, insert_file, config.inc_percentage, config.bat_size);
		}
		else if (config.alg_name == TC_AHRSZ) {
			static_scc_pearce2_tc_ol static_graph(np, dir + pyear, dir + pref);
			do_dynamic_cd4<sin_dscc, static_scc_pearce2_tc_ol >(static_graph, dir, pyear, pref, config.insert_type, insert_file, config.inc_percentage, config.bat_size);
		}
		else if (config.alg_name == TC_HKMST) {
			static_scc_pearce2_tc_ol static_graph(np, dir + pyear, dir + pref);
			do_dynamic_cd4<dn_tc_hkmst_scc, static_scc_pearce2_tc_ol >(static_graph, dir, pyear, pref, config.insert_type, insert_file, config.inc_percentage, config.bat_size);
		}
		else if (config.alg_name == DYN_PEARCE2 || config.alg_name == DYN_TC_PEARCE2) {
			static_scc_pearce2_tc static_graph(np, dir + pyear, dir + pref);
			do_dynamic_baseline<static_scc_pearce2_tc>(static_graph, dir, pyear, pref, config.insert_type, insert_file, config.inc_percentage, config.bat_size);
		}

		myres.set_file_info(dir, pyear, pref );
		myres.to_file();
	}
	else if (config.alg_type == TYPE_BAT) {
		if (config.alg_name == BAT_PK2) {
			//const std::string BAT_PK2 = "bat_pk2";
			static_scc_pearce2 static_graph(np, dir + pyear, dir + pref);
			do_batch_pk2<bat_pk2_scc, static_scc_pearce2>(static_graph, dir, pyear, pref, config.insert_type, insert_file, config.inc_percentage, config.bat_size);
		}
		else if (config.alg_name == BAT_AHRSZ) {
			//const std::string BAT_AHRSZ = "bat_ahrsz";
			static_scc_pearce2_ol static_graph(np, dir + pyear, dir + pref);
			do_batch_ahrsz<bat_ahrsz_scc, static_scc_pearce2_ol >(static_graph, dir, pyear, pref, config.insert_type, insert_file, config.inc_percentage, config.bat_size);
		}
		else if (config.alg_name == BAT_AHRSZP) {
			//const std::string BAT_AHRSZ = "bat_ahrsz_p";
			static_scc_pearce2_ol static_graph(np, dir + pyear, dir + pref);
			do_batch_ahrsz<bat_ahrszp_scc, static_scc_pearce2_ol >(static_graph, dir, pyear, pref, config.insert_type, insert_file, config.inc_percentage, config.bat_size);
		}
		else if (config.alg_name == BAT_CD4CG) {
			//const std::string BAT_CD4CG = "bat_cd4cg";
			static_scc_pearce2_tc_ol static_graph(np, dir + pyear, dir + pref);
			do_batch_tc_ahrsz<bat_dscc, static_scc_pearce2_tc_ol >(static_graph, dir, pyear, pref, config.insert_type, insert_file, config.inc_percentage, config.bat_size);
		}
		else if (config.alg_name == BAT_CD4CG2) {
			static_scc_pearce2_tc_ol static_graph(np, dir + pyear, dir + pref);
			do_batch_tc_ahrsz<bat_dscc_test, static_scc_pearce2_tc_ol >(static_graph, dir, pyear, pref, config.insert_type, insert_file, config.inc_percentage, config.bat_size);
		}
		myres.set_file_info(dir, pyear, pref );
		myres.to_file();
	}
	else {
	}
}


int main(int argc, char **argv) {
	// increase stack size 
#ifdef linux
	inc_stack_size(NEW_STACK_SIZE);
#endif // linux


	Config config;
	try {
		config = parseParams(argc, argv);
		config.check();
		config.display();
		myres.simplify_data = config.simple_data;
		myres.scale_factor = config.scale_factor; 
	}
	catch (const std::exception &ex) {
		cerr << ex.what() << '\n';
		return ERROR_IN_COMMAND_LINE;
	}

	do_cycle_detection(config);


	return SUCCESS;
}

