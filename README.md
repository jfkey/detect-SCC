# For ICDE Reviewers:
When linking the reference [34](https://github.com/jfkey/detect-SCC/blob/main/full-version.pdf) to open the full version, you might encounter the issue `unable to render code block`. This is likely caused by the large code block associated with the full version PDF, making it impossible for GitHub to render the content properly. 
You have two options to access the full version: (i) directly download the entire GitHub repository or (ii) send an email to liujunfeng@buaa.edu.cn


# SCC Detection
This repository implements strongly connected component(SCC) detection for schoarly data, including static algorithm (staDSCC), single incremental algorithm (sinDSCC) and batch incremental algorithms (batDSCC). 
Note that, the static tested algorithms Pearce[1], Tarjan[2], Gabow[3], Kosaraju[4], and incremental tested algorithms AHRSZ[5], HKMST[6], incSCC+[7], MNR[8] and PK2 [9] are also included. 

## Structure
* `dataset` Dataset is placed under this folder. Due to the large data collection size, only `AAN` is in the repository.
  * `AAN` collects the articles published at ACL from 1965 to 2011.
    * `indpp_refine.txt` The entire graph of AAN, and the data format is  `x \t y`, which means article `x` cites article `y`. Note that `x` and `y` are int values and start with `0`.
    * `indpyear.txt` The published time of articles, and the data format is  `x \t i`, which means article `x` is published at year `i`.
    * `indpp_startC.txt` The graph of `G` in incremental algorithms. Here the size of 
    * `inc_c30.txt` The increment updates to the graph `G`, i.e., `\Delta G`, here `|\Delta G| = 0.3 |G|`, and `|G| + |\Delta G|` is actually the entire datasets.   
* `header` header files of all test the algorithms, and we list some of them
  * `bat_AHRSZ_P_scc.h` The header file of batch version of `AHRSZP` algorithms, the SCC detection and maintenance use `batDSCC`.  
  * `bat_AHRSZ_scc.h` The header file of batch version of `AHRSZ` algorithms
  * `bat_DSCC.h` The header file of our batch incremental algorithm `batDSCC`
  * `dynamic_HKMST_scc.h` The header file of incremental algorithm `HKMST`. 
  * `static_scc_Gabow.h` The header file of static algorithm `Gabow`
  * `static_scc_Gabow_tc.h` The header file of static algorithm `Gabow_TC`, i.e., applying the strategy of `staDSCC`on `Gabow`
  * `static_scc_Pearce2.h`  The header file of static algorithm `Pearce`
  * `static_scc_Pearce2_tc.h` The header file of our static algorithm `staDSCC`
  * ...
* `source` source files of the test algorithms, and this corresponds to the header file. We also list some of them.
  * `bat_AHRSZ_P_scc.cpp` The source file of batch version of `AHRSZP` algorithms, the SCC detection and maintenance use `batDSCC`.  
  * `bat_AHRSZ_scc.cpp` The source file of batch version of `AHRSZ` algorithms
  * `bat_DSCC.cpp` The source file of our batch incremental algorithm `batDSCC`
  * `dynamic_HKMST_scc.cpp` The source file of incremental algorithm `HKMST`. 
  * `cd4cg.cpp` The main file. 
  * ...
* `full-version.pdf`  The full version of the paper of ICDE24 submission `168` "Incremental Detection of Strongly Connected Components for Scholarly Data". 
* `README.md` readme file. 
 

## Requirements

- Boost C++ Libraries (Version 1.71.0)
- Visual Studio 2017 (MSVC compiler)

# Usage of SCC detection
  
- Running for static algorithm `staDSCC`

```
-f  %dir%  -t "static" -b 0 -a "tc_pearce2" -n %nn%  -c 0  -p "comb"  -g "indpp_startC.txt"  -s true
```
- 
  - `-f` the dataset folder, and the datasets shall be organized same as `AAN`.
  - `-t` the type of algorithm to execute, is one of `static`, `inc`, and `batch`  
  - `-b` the batch size of batch incremental algorithms 
  - `-a` the algorithm name to execute, e.g., `Tarjan`, `Pearce`
  - `-n` the number of nodes of the citation graph 
  - `-c` The percentage of the increments in the original citation graph, i.e., `0.5`,`1`,`5`,`10`,`15`,`20`,`25`,`30` for tests. 
  - `-p` the insert type, default use `comb` to run the algorithm. 
  - `-g` the file name of the increments
  - `-s` it is utilized to test the performance on simplify datasets. set `true` for test the algorithms.

Note that, 
`%dir%` refer to the dataset folder, for example, it could be replaced by  `"D:\xxxx\detect-SCC\dataset\AAN\"`
`%nn%` refer to the node number of the citation graph, for example it equals `18041`. 
 
  
- Running single incremental `sinDSCC`

```-f  %dir%  -t "inc" -b 1 -a "tc_ahrsz" -n %nn%  -c 0  -p "comb"  -g "indpp_startC.txt"  -s true```

The parameters are the same as those described above. One can easily replace `tc_ahrsz` with `hkmst` to run `HKMST`.

- Running batch incremental `batDSCC`

```-f  %dir%  -t "batch" -b %batsize% -a "bat_cd4cg" -n %nn% -c  -1  -p "comb"  -g "indpp_startC.txt"  -s true```

The parameters are the same as those described above.


## Datasets
* [AAN](https://aclanthology.org/W09-3607/) collects the articles published at ACL from 1965 to 2011.
* [DBLP](https://www.aminer.cn/citation)    collects the publications at the DBLP Bibliography from 1936 to 2016.
* [ACM](https://www.aminer.cn/citation) collects in Computer Science from 1936 to 2016.
* [MAG](https://www.microsoft.com/en-us/research/project/microsoft-academic-graph/) gathers different types of publications e.g., books, conferences, journals and patents of various disciplines from 1800 to 2021.


# References
- [1] David J Pearce. 2016. A space-efficient algorithm for finding strongly connected
components. Inform. Process. Lett. 116, 1 (2016), 47–52.
- [2] Robert Tarjan. 1972. Depth-first search and linear graph algorithms. SIAM journal
on computing 1, 2 (1972), 146–160.
- [3] Harold N. Gabow. 2000. Path-based depth-first search for strong and biconnected
components. Inform. Process. Lett. 74, 3 (2000), 107–114.
- [4] Micha Sharir. 1981. A strong-connectivity algorithm and its applications in data
flow analysis. Computers & Mathematics with Applications 7, 1 (1981), 67–72.
- [5] Bowen Alpern, Roger Hoover, Barry K Rosen, Peter F Sweeney, and F Kenneth
Zadeck. 1990. Incremental evaluation of computational circuits. In Proceedings of
the first annual ACM-SIAM symposium on Discrete algorithms. 32–42.
- [6] Bernhard Haeupler, Telikepalli Kavitha, Rogers Mathew, Siddhartha Sen, and
Robert E Tarjan. 2012. Incremental cycle detection, topological ordering, and
strong component maintenance. ACM Transactions on Algorithms (TALG) 8, 1
(2012), 1–33.
- [7] Wenfei Fan, Chunming Hu, and Chao Tian. 2017. Incremental graph computations:
Doable and undoable. In Proceedings of the 2017 ACM International Conference
on Management of Data. 155–169.
- [8] Alberto Marchetti-Spaccamela, Umberto Nanni, and Hans Rohnert. 1996. Maintaining
a topological order under edge insertions. Inform. Process. Lett. 59, 1
(1996), 53–58.
- [9] David J Pearce and Paul HJ Kelly. 2010. A batch algorithm for maintaining a
topological order. In Proceedings of the Thirty-Third Australasian Conferenc on
Computer Science-Volume 102. 79–88.
