#ifndef UTILS_HPP
#define UTILS_HPP


#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <assert.h>
#include <algorithm> // std::find
#include <chrono>
#include <fstream>
#include <tuple>
#ifdef linux
#include <sys/resource.h>
#endif


using namespace std;

#ifndef TIMES_PER_SEC
#define TIMES_PER_SEC (1.0e9)
#endif

// static scc algorithms
const std::string TARJAN = "tarjan";
const std::string GABOW = "gabow";
const std::string KOSARAJU = "kosaraju";
const std::string PEARCE = "pearce";
const std::string PEARCE2 = "pearce2";
const std::string TC_TARJAN = "tc_tarjan";
const std::string TC_GABOW = "tc_gabow";
const std::string TC_KOSARAJU = "tc_kosaraju";
const std::string TC_PEARCE2 = "tc_pearce2";

// the baseline of the incremental alg. 
const std::string DYN_PEARCE2 = "dyn_pearce2";
const std::string DYN_TC_PEARCE2 = "dyn_tc_pearce2";
const std::string DYN_BASELINE= "dyn_baseline";


// inc scc algorithms
const std::string AHRSZ = "ahrsz";
const std::string HKMST = "hkmst";
const std::string PK = "pk";
const std::string MNR = "mnr";
const std::string TC_AHRSZ = "tc_ahrsz";
const std::string TC_HKMST = "tc_hkmst";
const std::string TC_PK = "tc_pk";
const std::string TC_MNR = "tc_mnr";

// bat-inc scc algorithms
const std::string BAT_PK2 = "bat_pk2";
const std::string BAT_AHRSZ = "bat_ahrsz";
const std::string BAT_AHRSZP = "bat_ahrsz_p";
const std::string BAT_CD4CG = "bat_cd4cg";
const std::string BAT_CD4CG2 = "bat_cd4cg2";

const std::string TYPE_STATIC = "static";
const std::string TYPE_STATIC_DEF = "static_def"; // run all of the static algorithms.
const std::string TYPE_INC = "inc";
const std::string TYPE_BAT = "batch";

const std::string SEQ_INSERT = "seq";
const std::string RAND_INSERT = "rand";
const std::string COMBINE_INSERT = "comb"; 
#define DELTAG_SIZE 30				// size of delta G 

struct Config {
	string folder;
	string insert_file;
	string alg_type;
	string alg_name;
	string insert_type;
	string node_name;
	string edge_name;
	int n;
	int bat_size;
	double inc_percentage;
	bool simple_data; 
	double scale_factor;


	void display() {
		std::cout << "====================Configurations==================" << std::endl;
		std::cout << "data folder:                    " << folder << '\n';
		std::cout << "algorithm type:                 " << alg_type << '\n';
		std::cout << "algorithm name:                 " << alg_name << '\n';
		std::cout << "type of inserted edge:          " << insert_type << '\n';
		std::cout << "#node of the graph:             " << n << '\n';
		std::cout << "batch size for bat-inc:         " << bat_size << '\n';
		std::cout << "percentage of edge increment:   " << inc_percentage << "\n";
		std::cout << "node file name of graph:        " << node_name << '\n';
		std::cout << "edge file name of graph:        " << edge_name << '\n';
		std::cout << "simplify the datasets:		  " << simple_data << std::endl;
		std::cout << "====================Configurations==================" << std::endl;
	}

	void check() {
		vector<string> alg_type_list = { TYPE_STATIC, TYPE_STATIC_DEF, TYPE_INC, TYPE_BAT };
		vector<string> static_algs = { TARJAN, GABOW, KOSARAJU, PEARCE, PEARCE2, TC_TARJAN, TC_GABOW, TC_KOSARAJU,TC_PEARCE2 };
		vector<string> dynamic_algs = { AHRSZ, HKMST, PK, MNR, TC_AHRSZ, TC_HKMST,TC_PK, TC_MNR, DYN_PEARCE2, DYN_TC_PEARCE2 };
		vector<string> bat_algs = { BAT_PK2 , BAT_AHRSZ,BAT_AHRSZP, BAT_CD4CG,BAT_CD4CG2 };
		
		auto f = std::find(alg_type_list.begin(), alg_type_list.end(), alg_type);
		assert(f != alg_type_list.end());
		if (*f == TYPE_STATIC) {
			auto static_f = std::find(static_algs.begin(), static_algs.end(), alg_name);
			assert(static_f != static_algs.end());
		}
		else if (*f == TYPE_INC) {
			auto  inc_f = std::find(dynamic_algs.begin(), dynamic_algs.end(), alg_name);
			assert(inc_f != dynamic_algs.end());
		}
		else if (*f == TYPE_BAT) {
			auto  bat_f = std::find(bat_algs.begin(), bat_algs.end(), alg_name);
			assert(bat_f != bat_algs.end());
		}
		else {
			assert(0 == 1 && "alg_type error!");
			throw runtime_error("alg type error!.");
		}
	}
};

class Timer {
public:
	static std::vector<double> timeUsed;
	static  std::vector<string> timeUsedDesc;
	int id;
	std::chrono::steady_clock::time_point startTime;
	bool showOnDestroy;

	Timer(int id, string desc = "", bool showOnDestroy = false) {
		this->id = id;
		while ((int)timeUsed.size() <= id) {
			timeUsed.push_back(0);
			timeUsedDesc.push_back("");
		}
		timeUsedDesc[id] = desc;
		startTime = std::chrono::steady_clock::now();
		this->showOnDestroy = showOnDestroy;
	}

	static double used(int id) {
		return timeUsed[id] / TIMES_PER_SEC;
	}

	~Timer() {
		auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - startTime).count();
		if (showOnDestroy) {
			std::cout << "time spend on " << timeUsedDesc[id] << ":" << duration / TIMES_PER_SEC << "s" << std::endl;
		}
		timeUsed[id] += duration;
	}

	static  void show(bool debug = false) {
		cout << "##### Timer #####" << endl;
		for (int i = 0; i < (int)timeUsed.size(); i++) {
			if (timeUsed[i] > 0) {
				//sprintf(str, "%.6lf", timeUsed[i] / TIMES_PER_SEC);
				std::cout << setprecision(6) << timeUsed[i] / TIMES_PER_SEC;
				// sprintf(t, "%4d %s %s", i, s.c_str(), timeUsedDesc[i].c_str());
				std::cout << setw(4) << i << timeUsedDesc[i].c_str() << std::endl;
			}
		}
	}

	static  void reset(int id) {
		if (id >= 0 && id < timeUsed.size()) {
			timeUsed[id] = 0;
		}
	}

	static  void clearAll() {
		timeUsed.clear();
		timeUsedDesc.clear();
	}
};


class my_timer {
private:
	std::chrono::system_clock::time_point start;

public:
	my_timer() {
		start = std::chrono::system_clock::now();
	}

public:
	string to_str(void) {
		std::time_t start_t = std::chrono::system_clock::to_time_t(start);
		return to_string(start_t);
	}

	double elapsed(void) {
		std::chrono::system_clock::time_point end;
		auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - start).count();
		return duration / TIMES_PER_SEC;
	}
};


class my_result {

public:
	// general info 
	string			dir;
	string			node_fn;
	string			edge_fn;
	string			inc_fn;
	unsigned int	same_year_cite;
	unsigned int	unorder_cite;
	unsigned int	normal_cite;
	size_t			graph_memory;
	size_t			static_memory;
	size_t			dynamic_memory;
	bool			simplify_data;
	double			scale_factor;

	// static algorithm results
	string			static_name;
	double			scc_time;
	unsigned int	paper_cnt;
	unsigned int	edge_cnt;
	int				max_year;
	int				min_year;
	unsigned int	access_nodes;
	unsigned int	access_edges;
	double			eliminate_scc_time;
	
	// dynamic algorithm results
	string			dynamic_name;
	double			dn_eval_time;
	double			dn_eval_time_AAN;
	unsigned int	inc_edges;
	int				batch;
	double			search_time;
	double			reassignment_time;
	unsigned int	aff_region; 
	unsigned int	invalid_edges;
	unsigned int	num_visited;					// evaluated for AHRSZ. number of visited, i.e., the number of nodes pushed into Min/Max heap.
	double			add_edge_time;
	double			condensation_time;
	double			tc_comp_time;					// the cost of our Time-Coupled Componnent 
	double			test_time1;					// test time
	double			test_time2;
	double			test_time3;
	double			test_time4;
	int				invalid_all_cnt; 
	int				invalid_part_cnt;
	int				invalid_local_cnt;


	std::vector<double> dn_eval_time_vec;
	std::vector<unsigned int> inc_edges_vec;
	std::vector<double> search_time_vec;
	std::vector<double> reassignment_time_vec;
	std::vector<unsigned int> aff_region_vec;
	std::vector<unsigned int> invalid_edges_vec;
	std::vector<double> scc_time_vec;

	bool is_shown;
	std::vector<std::pair<int, int >> res_scc;																		// static scc & scc size 
	std::vector<std::pair<int, int >> dynamic_res_scc;																// dynamic scc & scc size
	std::vector<std::tuple<bool, unsigned int, unsigned int, std::vector<unsigned int> >> invedge_affreg;			//invalid edges and affection region. bool is true: Gm, bool is false: Gsi 
	
public:
	my_result(bool _is_shown = false) {
		this->is_shown = _is_shown;
		// general info 
		dir					= "";
		node_fn				= "";
		edge_fn				= "";
		inc_fn				= "";
		same_year_cite		= 0;
		unorder_cite		= 0;
		normal_cite			= 0;
		graph_memory		= 0;
		static_memory		= 0;
		dynamic_memory		= 0;
		simplify_data		= false; 
		
		// static algorithm results
		static_name			= "";
		paper_cnt			= 0;
		scc_time			= 0;
		edge_cnt			= 0;
		max_year			= 0;
		min_year			= 0;
		access_nodes		= 0;
		access_edges		= 0;
		eliminate_scc_time	= 0;
	
		// dynamic algorithm results
		dynamic_name		= "";
		dn_eval_time		= 0;
		dn_eval_time_AAN	= 0; 
		inc_edges			= 0;
		batch				= 0;
		search_time			= 0;
		reassignment_time	= 0;
		aff_region			= 0;
		invalid_edges		= 0;
		num_visited			= 0;
		add_edge_time		= 0;
		condensation_time	= 0;
		tc_comp_time		= 0;
		
		test_time1 = 0;
		test_time2 = 0;
		test_time3 = 0;
		test_time4 = 0;
		invalid_all_cnt = 0;
		invalid_local_cnt = 0; 
		invalid_part_cnt = 0;

	}

	void to_file() {
		my_timer mt;
		
		std::stringstream stream;
		stream << std::fixed << std::setprecision(1) << scale_factor;
		std::string ssf = stream.str();

		std::ofstream ofs(dir + static_name + dynamic_name +
			"inc-" + to_string(inc_edges) + "-bat-" + to_string(batch) + "-r-" + ssf + "-" + mt.to_str() + ".txt");

		ofs << setiosflags(ios::fixed) << setprecision(6);
		ofs << "****************general info****************" << std::endl;
		ofs << "dir:                     " << dir				<< std::endl;
		ofs << "node file name:          " << node_fn			<< std::endl;
		ofs << "edge file name:          " << edge_fn			<< std::endl;
		ofs << "inc edge file name:      " << inc_fn			<< std::endl;
		ofs << "unorder year citations:  " << unorder_cite		<< std::endl;
		ofs << "same year citations:     " << same_year_cite	<< std::endl;
		ofs << "normal year citations:   " << normal_cite		<< std::endl;
		ofs << "space of graph (MB):     " << graph_memory		<< std::endl;
		ofs << "space of static(MB):     " << static_memory		<< std::endl;
		ofs << "space of inc   (MB):     " << dynamic_memory	<< std::endl;
		ofs << "simplify the datasets:   " << simplify_data		<< std::endl;


		ofs << std::endl;

		ofs << "***************static-alg info***************" << std::endl;
		ofs << "static alg name:         " << static_name		<< std::endl;
		ofs << "static SCC time:         " << to_str_double(scc_time_vec)    << std::endl;
		ofs << "node counts:             " << paper_cnt			<< std::endl;
		ofs << "edge counts:             " << edge_cnt			<< std::endl;
		ofs << "max year:                " << max_year			<< std::endl;
		ofs << "min year:                " << min_year			<< std::endl;
		ofs << "node access_ratio:       " << setprecision(6)	<< ((access_nodes * 1.0) / (paper_cnt * 1.0)) << std::endl;
		ofs << "access nodes:            " << access_nodes		<< std::endl;
		ofs << "edge access_ratio:       " << setprecision(6)	<< ((access_edges * 1.0) / (edge_cnt * 1.0)) << std::endl;
		ofs << "access edges:            " << access_edges		<< std::endl;
		ofs << "eliminate scc time:      " << eliminate_scc_time<< std::endl;
		ofs << std::endl;

		ofs << "***************dynamic-alg info***************" << std::endl;
		ofs << "dynamic alg name:        " << dynamic_name		<< std::endl;
		ofs << "dynamic alg time:        " << to_str_double(dn_eval_time_vec)  << std::endl;
		ofs << "dynamic alg time(AAN):   " << dn_eval_time_AAN	<< std::endl;
		ofs << "insert edges:            " << to_str_int(inc_edges_vec)			<< std::endl;
		ofs << "batch size:              " << batch				<< std::endl;
		ofs << "search time:             " << to_str_double(search_time_vec)		<< std::endl;
		ofs << "reassignment time:       " << to_str_double(reassignment_time_vec) << std::endl;
		ofs << "affect region:           " << to_str_int(aff_region_vec)		<< std::endl;
		ofs << "invalid edges:           " << to_str_int(invalid_edges_vec)		<< std::endl;
		ofs << "number of visited:       " << num_visited		<< std::endl;
		ofs << "condense SCC time:       " << condensation_time << std::endl;
		ofs << "all invalid edge cnt:    " << invalid_all_cnt << std::endl;
		ofs << "invalid partitioning cnt:" << invalid_part_cnt<< std::endl;
		ofs << "invalid local order cnt: " << invalid_local_cnt<< std::endl;

		ofs << std::endl;

		ofs << "************STATIC-SCC size & counts************" << std::endl;
		for (int i = 0; i < res_scc.size(); i++) {
			ofs << res_scc[i].first << "\t" << res_scc[i].second << std::endl;
		}

		ofs << std::endl;
		ofs << "************DYNAMIC-SCC size & counts************" << std::endl;
		for (int i = 0; i < dynamic_res_scc.size(); i++) {
			ofs << dynamic_res_scc[i].first << "\t" << dynamic_res_scc[i].second << std::endl;
		}

		ofs.close();

		/*
		std::ofstream ofs2(dir + static_name + dynamic_name + "-debug.txt");
		for (int i = 0; i < invedge_affreg.size(); i++) {
			if (std::get<0>(invedge_affreg[i]) == true) {
				ofs2 << "Gm " << std::get<1>(invedge_affreg[i]) << " -> " << std::get<2>(invedge_affreg[i]) << std::endl;
			}
			else {
				ofs2 << "Gsi" << std::get<1>(invedge_affreg[i]) << " -> " << std::get<2>(invedge_affreg[i]) << std::endl;
			}
			auto affreg = std::get<3>(invedge_affreg[i]);
			for (int j = 0; j < affreg.size(); j++) {
				ofs2 << affreg[j] << ", ";
			}
			ofs2 << std::endl << std::endl;
		}
		ofs2.close();
		*/
	}


public:
	void set_file_info(const string dir, const string node_fn, const string edge_fn ) {
		this->dir		= dir;
		this->node_fn	= node_fn;
		this->edge_fn	= edge_fn;
	}

	void set_read_data_info(const unsigned int pc, const unsigned int ec, const int _max_year, const int _min_year) {
		paper_cnt = pc; 
		edge_cnt = ec; 
		max_year = _max_year; 
		min_year = _min_year;
	}
	void set_cites(const int unorder_cites, const int same_year_cites, const int normal_cites) {
		this->unorder_cite = unorder_cites;
		this->same_year_cite = same_year_cites;
		this->normal_cite = normal_cites;
	}
	void set_static_name(const string static_name) {
		this->static_name = static_name;
	}
	void set_eliminate_scc_time(const double eli_time) {
		this->eliminate_scc_time = eliminate_scc_time;
	}

	void set_dynamic_name(const string dynamic_name) {
		this->dynamic_name = dynamic_name;
	}
	void dynamic_info(const int inc_edges, const int bat_size) {
		/*this->inc_edges = inc_edges;*/
		this->batch = bat_size;
	}
	void set_scc_info(const double scc_time) {
		this->scc_time = scc_time;
	}

	void set_search_time(const double search_time) {
		this->search_time = search_time;
	}
	double get_search_time() {
		return this->search_time;
	}	void set_reassign_time(const double reassign_time) {
		this->reassignment_time = reassign_time;
	}
	double get_reassign_time() {
		return this->reassignment_time;
	}
	void set_aff_region(const int aff_region) {
		this->aff_region = aff_region;
	}
	int get_aff_region() {
		return this->aff_region;
	}
	void set_invalid_edges(const int invalid_edges) {
		this->invalid_edges = invalid_edges;
	}
	int get_invalid_edges() {
		return this->invalid_edges;
	}

	void gen_sccs(std::vector<std::vector<unsigned int>> &scc_nodes) {
		//assert(!scc_nodes.empty() && "The scc beign printed is EMPTY.");
		int np = scc_nodes.size();
		std::vector<std::vector<unsigned int>>::iterator it;
		for (it = scc_nodes.begin(); it != scc_nodes.end(); it++) {
			if (it->size() > np) {
				np = it->size();
			}
		}
		std::vector<int>cc_size_cnt(np + 1, 0);
		for (it = scc_nodes.begin(); it != scc_nodes.end(); it++) {
			if ((*it).size() > 0) {
				cc_size_cnt[(*it).size()]++;
			}
		}
		for (int i = 1; i < np + 1; i++) {
			if (cc_size_cnt[i] > 0) {
				res_scc.push_back(std::make_pair(i, cc_size_cnt[i]));
			}
		}
	}

	void gen_dynamic_scc(std::vector<std::vector<unsigned int>> &scc_nodes) {
		//assert(!scc_nodes.empty() && "The scc beign printed is EMPTY.");
		int np = scc_nodes.size();
		std::vector<std::vector<unsigned int>>::iterator it;
		for (it = scc_nodes.begin(); it != scc_nodes.end(); it++) {
			if (it->size() > np) {
				np = it->size();
			}
		}
		std::vector<int>cc_size_cnt(np + 1, 0);
		for (it = scc_nodes.begin(); it != scc_nodes.end(); it++) {
			if ((*it).size() > 0) {
				cc_size_cnt[(*it).size()]++;
			}
		}
		for (int i = 1; i < np + 1; i++) {
			if (cc_size_cnt[i] > 0) {
				dynamic_res_scc.push_back(std::make_pair(i, cc_size_cnt[i]));
			}
		}
	}
	string to_str_int(std::vector<unsigned int> &res) {
		string s = "";
		for (int i = 0; i < res.size();i ++) {
			s += to_string(res[i]);
			if ( i != res.size() - 1 ) {
				s += " : ";
			}
		}
		return s; 
	}

	string to_str_double(std::vector<double> &res) {
		
		string s = "";
		for (int i = 0; i < res.size(); i++) {
			std::stringstream stream;
			stream << std::fixed << std::setprecision(6) << res[i];
			s += stream.str();
			if (i != res.size() - 1) {
				s += " : ";
			}
		}
		return s;
	}




};

inline void record_cur_per_info(my_result & myres) {
	myres.dn_eval_time_vec.push_back(myres.dn_eval_time);
	myres.inc_edges_vec.push_back(myres.inc_edges);
	myres.search_time_vec.push_back(myres.search_time);
	myres.reassignment_time_vec.push_back(myres.reassignment_time);
	myres.aff_region_vec.push_back(myres.aff_region);
	myres.invalid_edges_vec.push_back(myres.invalid_edges);
	myres.scc_time_vec.push_back(myres.scc_time);

}

extern my_result myres;
extern std::vector<bool> cond_scc_in;									// condense scc tmp variables. 
extern std::vector<bool> cond_scc_out;									

#ifdef linux
inline void inc_stack_size(const long ss) {
	const rlim_t kStackSize = ss * 1024L * 1024L;   // min stack size = 64 Mb
	struct rlimit rl;
	int result;
	result = getrlimit(RLIMIT_STACK, &rl);
	if (result == 0) {
		if (rl.rlim_cur < kStackSize) {
			rl.rlim_cur = kStackSize;
			result = setrlimit(RLIMIT_STACK, &rl);
			if (result != 0) {
				std::cout << "setrlimit returned result: " << result / 1024 / 1024 << " MB." << std::endl;
			}
		}
	}
}
#endif



#endif // !UTILS_HPP
